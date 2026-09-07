# ================================================================
#  GOLD PRICE FORECASTING  —  best ARIMA  vs  best ETS (trend only)
#  Daily data, modelled on the LOG of price (lambda = 0)
#
#  Pipeline (standard Box-Jenkins / ETS comparison):
#    1. Load + clean                       (Part 1)
#    2. Descriptive statistics + plot      (Part 2)
#    3. Decomposition: trend / seasonality (Part 3)  -> justifies non-seasonal model
#    4. Train / test split                 (Part 4)
#    5. Stationarity: ADF + KPSS -> d      (Part 5)
#    6. ACF / PACF plots                   (Part 6)
#    7. Best ARIMA (auto + manual grid)    (Part 7)
#    8. Best ETS with trend (Holt linear)  (Part 8)
#    9. Model ADEQUACY comparison          (Part 9)
#   10. FORECAST ERROR comparison          (Part 10)
#   11. Plots                              (Part 11)
#   12-18. ARIMA-GARCH extension           (added below)
# ================================================================

# ---- Part 0: packages -------------------------------------------------
pkgs <- c("forecast", "tseries", "FinTS")
to_install <- pkgs[!(pkgs %in% rownames(installed.packages()))]
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Part 1: load and clean data --------------------------------------
gold <- read.csv("GOLDDAILYPP.csv", stringsAsFactors = FALSE)
gold <- gold[, 1:2]
colnames(gold) <- c("Date", "Price")
gold$Date  <- as.Date(gold$Date, format = "%d-%b-%y")
gold$Price <- as.numeric(gsub(",", "", gold$Price))
gold <- gold[!is.na(gold$Date) & !is.na(gold$Price), ]
gold <- gold[order(gold$Date), ]
cat("Rows:", nrow(gold), "| From", format(min(gold$Date)), "to", format(max(gold$Date)), "\n\n")

# ---- Part 2: descriptive statistics + plot ----------------------------
plot(gold$Date, gold$Price, type = "l", col = "steelblue",
     main = "Gold Price - full series", xlab = "Date", ylab = "Price (INR/10g)")

cat("=== DESCRIPTIVE STATISTICS ===\n")
cat("N      :", nrow(gold), "\n")
cat("Mean   :", round(mean(gold$Price), 2), "\n")
cat("Median :", round(median(gold$Price), 2), "\n")
cat("SD     :", round(sd(gold$Price), 2), "\n")
cat("Min    :", min(gold$Price), "\n")
cat("Max    :", max(gold$Price), "\n")
cat("Range  :", max(gold$Price) - min(gold$Price), "\n")
cat("CV     :", round(sd(gold$Price) / mean(gold$Price) * 100, 2), "%\n\n")

# ---- Part 3: decomposition (trend / seasonality strength) -------------
# STL splits the series into trend + seasonal + remainder.
# Strength is 0 (none) to 1 (strong). For daily gold we expect a strong
# trend and almost no seasonality, which justifies a NON-seasonal model.
cat("=== DECOMPOSITION / SEASONALITY ===\n")
gold_ts_full <- ts(gold$Price, frequency = 252)        # ~252 trading days/year
stl_full <- stl(gold_ts_full, s.window = "periodic", robust = TRUE)
plot(stl_full, main = "STL Decomposition - Gold Price")

rem <- stl_full$time.series[, "remainder"]
sea <- stl_full$time.series[, "seasonal"]
tre <- stl_full$time.series[, "trend"]
trend_strength    <- max(0, 1 - var(rem) / var(tre + rem))
seasonal_strength <- max(0, 1 - var(rem) / var(sea + rem))
cat("Trend strength    :", round(trend_strength, 3), "\n")
cat("Seasonal strength :", round(seasonal_strength, 3), "\n")
cat("Decision:",
    ifelse(seasonal_strength > 0.40,
           "seasonality matters -> SARIMA",
           "seasonality weak -> non-seasonal ARIMA / ETS is appropriate"), "\n\n")

# ---- Part 4: train / test split ---------------------------------------
# TEST = the last 35 trading days. (To use an 80/20 split instead,
# replace the next line with:  test_size <- nrow(gold) - floor(0.80*nrow(gold)) )
test_size <- 35

n     <- nrow(gold)
train <- gold[1:(n - test_size), ]
test  <- gold[(n - test_size + 1):n, ]
train_ts <- ts(train$Price, frequency = 1)
actual   <- test$Price
h        <- length(actual)
cat("Train:", nrow(train), "obs  |  Test:", nrow(test), "obs  |  horizon h =", h, "\n\n")

# ---- Part 5: stationarity (ADF + KPSS) -> decide d --------------------
# Two tests, OPPOSITE logic, used together as a cross-check:
#   ADF  -> small p (< 0.05) means STATIONARY
#   KPSS -> small p (< 0.05) means NON-stationary
# A price is usually non-stationary; after one difference it becomes
# stationary, so d = 1.
cat("=== STATIONARITY (decides d) ===\n")
adf_lvl  <- suppressWarnings(adf.test(train_ts))
kpss_lvl <- suppressWarnings(kpss.test(train_ts))
cat("ADF  on price       p =", round(adf_lvl$p.value, 4),
    "->", ifelse(adf_lvl$p.value < 0.05, "stationary", "NON-stationary"), "\n")
cat("KPSS on price       p =", round(kpss_lvl$p.value, 4),
    "->", ifelse(kpss_lvl$p.value < 0.05, "NON-stationary", "stationary"), "\n")

adf_d1 <- suppressWarnings(adf.test(diff(train_ts)))
cat("ADF  on differenced p =", round(adf_d1$p.value, 4),
    "->", ifelse(adf_d1$p.value < 0.05, "stationary", "still NON-stationary"), "\n")

d_val <- if (adf_d1$p.value < 0.05) 1 else 2
cat("=> d =", d_val, "\n\n")

# ---- Part 6: ACF / PACF of the differenced series ---------------------
# These help justify the candidate ARIMA orders (p from PACF, q from ACF).
par(mfrow = c(1, 2))
acf (diff(train_ts), lag.max = 30, main = "ACF (differenced)",  col = "steelblue", lwd = 2)
pacf(diff(train_ts), lag.max = 30, main = "PACF (differenced)", col = "darkred",   lwd = 2)
par(mfrow = c(1, 1))

# ================================================================
# Part 7: BEST ARIMA
#   (a) auto.arima searches automatically
#   (b) a manual grid over p,q = 0..3, ranked by AICc (lower = better)
#   Both use lambda = 0 -> model the LOG of price.
#   A good sign: the two methods should agree on the same model.
# ================================================================
cat("=== BEST ARIMA ===\n")
fit_arima <- auto.arima(train_ts, d = d_val, seasonal = FALSE,
                        stepwise = FALSE, approximation = FALSE,
                        ic = "aicc", lambda = 0)
ord <- arimaorder(fit_arima)
cat(sprintf("auto.arima picked: ARIMA(%d,%d,%d)\n\n", ord[1], ord[2], ord[3]))
print(fit_arima)
cat("Manual grid (ranked by AICc, lowest = best):\n")
grid <- list()
for (p in 0:3) for (q in 0:3) {
  f <- tryCatch(Arima(train_ts, order = c(p, d_val, q), method = "ML", lambda = 0),
                error = function(e) NULL)
  if (!is.null(f))
    grid[[paste0(p, q)]] <- data.frame(
      Model = sprintf("ARIMA(%d,%d,%d)", p, d_val, q),
      AICc  = round(f$aicc, 2))
}
grid_tab <- do.call(rbind, grid)
grid_tab <- grid_tab[order(grid_tab$AICc), ]
print(grid_tab, row.names = FALSE)
cat(sprintf("\nManual grid best: %s  (matches auto.arima if same)\n\n", grid_tab$Model[1]))

fc_arima   <- forecast(fit_arima, h = h, level = 95)
pred_arima <- as.numeric(fc_arima$mean)

# ================================================================
# Part 8: BEST ETS WITH TREND (no seasonality)
#   Holt's linear method = ETS(A,A,N): a level + a trend, no season.
#   This is the exponential-smoothing model that accounts for trend
#   but not seasonality - exactly the comparison you want against ARIMA.
#   NOTE: fitted with lambda = 0, so this is ETS(A,A,N) on the LOG scale
#   (i.e. multiplicative trend on the original price scale). The label
#   ETS(A,A,N) refers to the model in log space.
# ================================================================
cat("=== BEST ETS (trend, no seasonality) ===\n")
fit_ets <- holt(train_ts, h = h, damped = FALSE, lambda = 0)
cat("Model: Holt's linear method = ETS(A,A,N) on log scale\n")
cat("AICc :", round(fit_ets$model$aicc, 1), "\n\n")
pred_ets <- as.numeric(fit_ets$mean)
summary(pred_ets)
# ---- Part 9/10: forecast accuracy on the test set ----------------------

cat("=== ETS forecast accuracy ===\n")
print(accuracy(pred_ets, actual))

cat("\n=== ARIMA forecast accuracy ===\n")
print(accuracy(pred_arima, actual))
# ================================================================
# Part 9: MODEL ADEQUACY COMPARISON  (residual tests)
#   A good model's leftover errors (residuals) should be:
#     Ljung-Box  : no leftover pattern   -> want p > 0.05 (PASS)
#     ARCH-LM    : steady swing size     -> want p > 0.05 (PASS)
#     Jarque-Bera: roughly bell-shaped   -> want p > 0.05 (PASS)
#   (Ljung-Box is the key Box-Jenkins adequacy test.)
# ================================================================
adeq <- function(fit, npar) {
  r  <- as.numeric(residuals(fit)); r <- r[!is.na(r)]
  lb <- Box.test(r, lag = 20, type = "Ljung-Box", fitdf = npar)$p.value
  ar <- ArchTest(r, lags = 12)$p.value
  jb <- jarque.bera.test(r)$p.value
  verdict <- function(p) ifelse(p > 0.05, "PASS", "FAIL")
  c(LjungBox   = sprintf("%.4f (%s)", lb, verdict(lb)),
    ARCH       = sprintf("%.4f (%s)", ar, verdict(ar)),
    JarqueBera = sprintf("%.4f (%s)", jb, verdict(jb)))
}
adeq_tab <- rbind(
  ARIMA = adeq(fit_arima,     npar = ord[1] + ord[3]),
  ETS   = adeq(fit_ets$model, npar = length(fit_ets$model$par))
)
cat("=== ADEQUACY (residual tests; PASS = p > 0.05) ===\n")
print(as.data.frame(adeq_tab))
cat("\n")

# ================================================================
# Part 10: FORECAST ERROR COMPARISON  (ARIMA vs ETS)
#   Lower RMSE / MAE / MAPE = more accurate on the test set.
# ================================================================
err <- function(pred, name) {
  e <- actual - pred
  data.frame(Model = name,
             RMSE = round(sqrt(mean(e^2)), 2),
             MAE  = round(mean(abs(e)), 2),
             MAPE_pct = round(mean(abs(e / actual)) * 100, 3))
}
acc <- rbind(err(pred_arima, sprintf("ARIMA(%d,%d,%d)", ord[1], ord[2], ord[3])),
             err(pred_ets,  "ETS(A,A,N) Holt"))
cat("=== FORECAST ERROR (lower = better) ===\n")
print(acc, row.names = FALSE)
cat("\nMore accurate model:", acc$Model[which.min(acc$RMSE)], "\n\n")

write.csv(acc,      "gold_error_comparison.csv", row.names = FALSE)
write.csv(grid_tab, "gold_arima_grid.csv",       row.names = FALSE)
# ---- Forecast vs Actual table (day by day) ----------------------------
compare <- data.frame(
  Date       = test$Date,
  Actual     = round(actual, 2),
  ARIMA      = round(pred_arima, 2),
  ETS        = round(pred_ets, 2),
  Err_ARIMA  = round(actual - pred_arima, 2),
  Err_ETS    = round(actual - pred_ets, 2),
  Lo95_ARIMA = round(as.numeric(fc_arima$lower), 2),
  Hi95_ARIMA = round(as.numeric(fc_arima$upper), 2)
)
# Did the actual price fall inside ARIMA's 95% prediction interval?
compare$In_95CI <- compare$Actual >= compare$Lo95_ARIMA &
  compare$Actual <= compare$Hi95_ARIMA

cat("=== FORECAST vs ACTUAL (all 35 test days) ===\n")
print(compare, row.names = FALSE)

cat("\nARIMA coverage (% of actuals inside 95% interval):",
    round(mean(compare$In_95CI) * 100, 1), "%\n\n")

write.csv(compare, "gold_forecast_vs_actual.csv", row.names = FALSE)

# ================================================================
# Part 11: PLOTS
# ================================================================
# (a) Actual vs the two forecasts over the test period
plot(test$Date, actual, type = "l", lwd = 2, col = "black",
     ylim = range(c(actual, pred_arima, pred_ets)),
     main = "Actual vs Forecast (test period)",
     xlab = "Date", ylab = "Price (INR/10g)")
lines(test$Date, pred_arima, col = "#e63946", lwd = 2, lty = 2)
lines(test$Date, pred_ets,   col = "#457b9d", lwd = 2, lty = 3)
legend("topleft", bty = "n", lwd = 2,
       legend = c("Actual", "ARIMA", "ETS (Holt)"),
       col    = c("black", "#e63946", "#457b9d"),
       lty    = c(1, 2, 3))

# (b) Error bar chart
barplot(acc$RMSE, names.arg = acc$Model, col = c("#e63946", "#457b9d"),
        border = "white", ylab = "RMSE", main = "Forecast error (RMSE)")

cat("=== DONE ===\n")
write.csv(gold[, c("Date","Price")], "gold_clean.csv", row.names = FALSE)


# ================================================================
# ================================================================
#  ARIMA-GARCH EXTENSION
#  ----------------------------------------------------------------
#  Reason: in Part 9 the ARIMA residuals were adequate in the MEAN
#  (Ljung-Box PASS) but FAILED ARCH-LM and Jarque-Bera. Those two
#  failures are the same feature of daily financial data — volatility
#  clustering, which also produces fat tails. The fix is to model the
#  variance explicitly: re-estimate an ARMA(p,q) mean JOINTLY with a
#  GARCH(1,1) variance equation, with the errors allowed to follow a
#  Student-t distribution (heavy tails) instead of the normal. The
#  (p, q) used for the mean are the orders identified in Part 7; the
#  integration order d = 1 is handled by working on log-RETURNS below.
#
#  HONEST FRAMING (state this in the report): rugarch estimates the
#  mean and variance parameters TOGETHER, so the ARMA coefficients here
#  are NOT identical to the stand-alone Part 7 ARIMA fit — only the
#  ORDERS (p, q) are carried over. The model is therefore an
#  ARMA(p,q)-GARCH(1,1) on log-returns, equivalent to ARIMA(p,1,q)-GARCH.
#
#  FITTING SAMPLE (deliberate, not convenience): this GARCH is fitted on the
#  TRAINING returns only -- diff(log(train$Price)). That is BECAUSE Parts 16-17
#  validate the resulting 95% interval against the held-out TEST days (the
#  coverage table). Fitting on the full sample would let the model see the very
#  days it is being scored on (leakage) and inflate the coverage, so train-only
#  is the methodologically correct choice here. This differs from the crude-oil
#  and exchange-rate chapters, where GARCH is purely DESCRIPTIVE (not scored
#  against held-out data) and is therefore fitted on the FULL sample, following
#  the reference papers. The rule is: fit on TRAIN when you validate the output
#  out-of-sample; fit on FULL when the model is descriptive.
# ================================================================
# ================================================================

# ---- Part 12: load rugarch + build the return series -----------------
if (!requireNamespace("rugarch", quietly = TRUE)) install.packages("rugarch")
library(rugarch)

# GARCH is fitted on RETURNS, not on the price level. diff(log(price))
# is the daily log-return; this IS the differenced (d = 1) log series,
# so the integration "I" of the ARIMA is handled here. *100 expresses
# the return in percent (same convention as the crude oil chapter) and
# keeps the optimiser numerically stable.
ret <- diff(log(train$Price)) * 100

cat("\n=== PART 12: RETURN SERIES FOR GARCH ===\n")
cat("Length  :", length(ret), " (one obs lost to differencing)\n")
cat("Mean    :", round(mean(ret), 4), "%\n")
cat("SD      :", round(sd(ret),   4), "%\n")
sk <- mean((ret - mean(ret))^3) / sd(ret)^3
ku <- mean((ret - mean(ret))^4) / sd(ret)^4
cat("Skewness:", round(sk, 4), "\n")
cat("Kurtosis:", round(ku, 4), " (>3 = heavy tails, GARCH-t appropriate)\n\n")

# ================================================================
# Part 13: FIT THE ARMA(p,q)-GARCH(1,1) MODEL  [= ARIMA(p,1,q)-GARCH]
#   Mean equation     : the ARMA orders (p, q) from Part 7
#   Variance equation : GARCH(1,1)  sigma^2_t = omega
#                                              + alpha * e^2_{t-1}
#                                              + beta  * sigma^2_{t-1}
#   Distributions compared: Normal vs Student-t (chosen by AIC),
#   exactly as in the crude oil and exchange rate chapters.
# ================================================================
cat("=== PART 13: FIT ARIMA(", ord[1], ",", ord[2], ",", ord[3],
    ")-GARCH(1,1) ===\n", sep = "")

# NOTE: 'armaOrder' below is rugarch's own argument name for the (p, q)
# mean part; we feed it the p and q of the ARIMA identified in Part 7.
# rugarch RE-ESTIMATES these mean coefficients jointly with the GARCH
# variance, so they will differ numerically from the Part 7 stand-alone
# ARIMA (only the orders are shared). The full model is ARIMA(p,1,q)-GARCH(1,1).
spec_norm <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(ord[1], ord[3]), include.mean = TRUE),
  distribution.model = "norm")

spec_t <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(ord[1], ord[3]), include.mean = TRUE),
  distribution.model = "std")

cat("Fitting ARIMA-GARCH with Normal distribution...\n")
fit_garch_norm <- ugarchfit(spec_norm, data = ret, solver = "hybrid")
cat("Fitting ARIMA-GARCH with Student-t distribution...\n\n")
fit_garch_t    <- ugarchfit(spec_t,    data = ret, solver = "hybrid")

# Compare the two distributions by AIC / BIC (lower = better)
ic_norm <- infocriteria(fit_garch_norm)
ic_t    <- infocriteria(fit_garch_t)
ic_tab  <- data.frame(
  Criterion = rownames(ic_norm),
  Normal    = round(as.numeric(ic_norm), 4),
  Student_t = round(as.numeric(ic_t),    4))
cat("=== DISTRIBUTION COMPARISON (Normal vs Student-t) ===\n")
print(ic_tab, row.names = FALSE)


# Select the better distribution by AIC (row 1)
if (ic_t[1] < ic_norm[1]) {
  best_garch <- fit_garch_t;    best_dist <- "Student-t"
  cat("\nSelected: Student-t (lower AIC) — expected for daily financial data\n\n")
} else {
  best_garch <- fit_garch_norm; best_dist <- "Normal"
  cat("\nSelected: Normal (lower AIC)\n\n")
}

cat("=== FITTED ARIMA(", ord[1], ",", ord[2], ",", ord[3],
    ")-GARCH(1,1) MODEL [", best_dist, "] ===\n", sep = "")
show(best_garch)

# ================================================================
# Part 14: PARAMETER ESTIMATION AND INTERPRETATION
# ================================================================
cat("\n=== PART 14: PARAMETER INTERPRETATION ===\n")
cf      <- coef(best_garch)
alpha1  <- cf["alpha1"]
beta1   <- cf["beta1"]
persist <- alpha1 + beta1
summary(best_garch)
print(best_garch)
cat("omega (long-run baseline variance):", round(cf["omega"], 6), "\n")
cat("alpha (reaction to last shock)    :", round(alpha1, 4), "\n")
cat("beta  (carry-over of past variance):", round(beta1, 4), "\n")
cat("alpha + beta (volatility persistence):", round(persist, 4), "\n")
if (persist > 0.9) {
  cat("  -> HIGH persistence: a volatility spike fades only slowly\n")
} else if (persist > 0.7) {
  cat("  -> MODERATE persistence: volatility settles gradually\n")
} else {
  cat("  -> LOW persistence: volatility calms quickly\n")
}
if ("shape" %in% names(cf)) {
  cat("shape (Student-t degrees of freedom):", round(cf["shape"], 3), "\n")
  cat("  -> smaller value = fatter tails; ~30+ behaves like the normal\n")
}
# Long-run (unconditional) average volatility implied by the model
if (persist < 1) {
  lr_var <- cf["omega"] / (1 - persist)
  cat("Implied long-run average volatility:", round(sqrt(lr_var), 4), "% per day\n")
}
garch_fc <- ugarchforecast(best_garch, n.ahead = h)
# Extract in-sample conditional standard deviation from fitted GARCH
actual_vol <- sigma(best_garch)
print(best_garch)






# Forecast volatility for test period
garch_forecast <- ugarchforecast(best_garch, n.ahead = h)

# Extract conditional standard deviation (volatility)
volatility_forecast <- sigma(garch_forecast)

# Plot
plot(volatility_forecast, type = "l", col = "steelblue", lwd = 2,
     main = "GARCH(1,1) Forecasted Volatility - Gold Price",
     xlab = "Days", ylab = "Conditional Standard Deviation")# ================================================================
# Part 15: RE-TEST THE STANDARDISED RESIDUALS
#   This directly answers the diagnostic failures from Part 9.
#   A good ARIMA-GARCH leaves standardised residuals that are:
#     - uncorrelated            (Ljung-Box on z)      -> mean still OK
#     - free of ARCH effects    (Ljung-Box on z^2,    -> variance now OK
#                                ARCH-LM on z)
#     - consistent with the chosen distribution (Student-t)
# ================================================================
cat("\n=== PART 15: STANDARDISED-RESIDUAL DIAGNOSTICS ===\n")
z <- as.numeric(residuals(best_garch, standardize = TRUE))

lb_z   <- Box.test(z,   lag = 20, type = "Ljung-Box")
lb_z2  <- Box.test(z^2, lag = 20, type = "Ljung-Box")
arch_z <- ArchTest(z, lags = 12)
vv <- function(p) ifelse(p > 0.05, "PASS", "FAIL")
lb_z
lb_z2
arch_z

cat(sprintf("Ljung-Box on z      p = %.4f (%s)  -> mean still adequate\n",
            lb_z$p.value,  vv(lb_z$p.value)))
cat(sprintf("Ljung-Box on z^2    p = %.4f (%s)  -> volatility now captured\n",
            lb_z2$p.value, vv(lb_z2$p.value)))
cat(sprintf("ARCH-LM(12) on z    p = %.4f (%s)  -> compare to plain ARIMA (was FAIL)\n",
            arch_z$p.value, vv(arch_z$p.value)))
skz <- mean((z - mean(z))^3) / sd(z)^3
kuz <- mean((z - mean(z))^4) / sd(z)^4 - 3
cat(sprintf("Skewness of z       = %.3f | Excess kurtosis = %.3f\n", skz, kuz))

cat("\nGoodness-of-fit to the chosen distribution (adjusted Pearson chi-square;\n")
cat("H0: standardised residuals follow the FITTED distribution, so large p = good fit):\n")
print(gof(best_garch, c(20, 30, 40, 50)))
cat("Once a Student-t is assumed, the relevant question is fit to the t,\n")
cat("not normality. (Jarque-Bera/Shapiro over-reject at large n like n~1497.)\n")

# Diagnostic plots for the standardised residuals
par(mfrow = c(2, 2))
plot(z, type = "l", col = "steelblue", main = "Standardised Residuals",
     ylab = "z", xlab = "Time"); abline(h = c(-2, 0, 2), col = c("red","black","red"), lty = 2)
# Histogram overlaid with the FITTED density: standardised Student-t (unit
# variance) when a t was selected, else the normal. This matches the model's
# distributional assumption rather than defaulting to the normal bell curve.
hist(z, breaks = 40, col = "lightblue", border = "white", probability = TRUE,
     main = "Histogram of z", xlab = "z")
if ("shape" %in% names(cf)) {
  curve(ddist("std", x, mu = 0, sigma = 1, shape = cf["shape"]),
        add = TRUE, col = "red", lwd = 2)               # fitted standardised-t density
} else {
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
}
acf(z,   main = "ACF of z",   lag.max = 24, col = "steelblue")
acf(z^2, main = "ACF of z^2", lag.max = 24, col = "darkred")
par(mfrow = c(1, 1))
qqnorm(z, main = "QQ Plot of Standardised Residuals (normal reference)", pch = 16,
       cex = 0.5, col = "steelblue"); qqline(z, col = "red", lwd = 2)

# Conditional volatility — what the ARIMA-GARCH adds over plain ARIMA
cond_vol <- sigma(best_garch)
plot(cond_vol, type = "l", col = "darkred", lwd = 1.3,
     main = "Conditional Volatility of GOLD",xlab="Time",
     ylab = "Volatility (% per day)")
abline(h = mean(cond_vol), col = "blue", lty = 2)

legend("topleft", c("Daily volatility", "Average"),
       col = c("darkred", "blue"), lty = c(1, 2), bty = "n", cex = 0.8)
###########################################################################
# Conditional volatility
cond_vol <- sigma(best_garch)

# Actual dates from your dataset
dates <- as.Date(train$Date)

length(cond_vol)
length(dates)
head(dates)
tail(dates)



# ================================================================
# Part 16: ARIMA-GARCH FORECAST WITH A VOLATILITY-ADAPTIVE 95% BAND
#   Point forecast: cumulate the forecast log-returns onto the last
#   price (the ARIMA mean). Interval: cumulative forecast variance,
#   using the quantiles of the fitted standardised residuals so the
#   band honours the heavy tails AND widens/narrows with volatility.
# ================================================================
cat("\n=== PART 16: ARIMA-GARCH FORECAST ===\n")
fc_g  <- ugarchforecast(best_garch, n.ahead = h)
mu_g  <- as.numeric(fitted(fc_g)) / 100      # forecast log-returns
sig_g <- as.numeric(sigma(fc_g))  / 100      # forecast daily volatility

logP0  <- log(tail(train$Price, 1))
cum_mu <- cumsum(mu_g)
cum_sd <- sqrt(cumsum(sig_g^2))
q_lo <- as.numeric(quantile(z, 0.025))       # as.numeric: drop the "2.5%" label
q_hi <- as.numeric(quantile(z, 0.975))       # as.numeric: drop the "97.5%" label

pred_garch <- exp(logP0 + cum_mu)
loG <- exp(logP0 + cum_mu + q_lo * cum_sd)
hiG <- exp(logP0 + cum_mu + q_hi * cum_sd)

garch_tab <- data.frame(
  Date    = test$Date,
  Actual  = round(actual, 2),
  Forecast= round(pred_garch, 2),
  Lo95    = round(loG, 2),
  Hi95    = round(hiG, 2),
  In_95CI = actual >= loG & actual <= hiG,
  CondVol = round(sig_g * 100, 4))           # forecast daily volatility (%)
print(garch_tab, row.names = FALSE)
write.csv(garch_tab, "gold_garch_forecast_vs_actual.csv", row.names = FALSE)

# ================================================================
# Part 17: INTERVAL-QUALITY COMPARISON (plain ARIMA vs ARIMA-GARCH)
#   GARCH is NOT a price competitor — the point forecast comes from the
#   ARIMA mean. So GARCH is judged on the INTERVAL, not on RMSE/MAE.
#   A good interval is NARROW while still covering ~95% of actuals.
# ================================================================
cat("\n=== PART 17: INTERVAL QUALITY (ARIMA vs ARIMA-GARCH) ===\n")
int_metric <- function(lo, hi, name) {
  data.frame(Model = name,
             MeanIntWidth = round(mean(hi - lo), 0),
             Coverage_pct = round(mean(actual >= lo & actual <= hi) * 100, 1))
}
int_tab <- rbind(
  int_metric(as.numeric(fc_arima$lower), as.numeric(fc_arima$upper),
             sprintf("ARIMA(%d,%d,%d) Normal band", ord[1], ord[2], ord[3])),
  int_metric(loG, hiG,
             sprintf("ARIMA(%d,%d,%d)-GARCH(1,1) %s", ord[1], ord[2], ord[3], best_dist)))
print(int_tab, row.names = FALSE)
write.csv(int_tab, "gold_interval_quality.csv", row.names = FALSE)
cat("\nInterpretation: the ARIMA mean already gives the price line; the\n")
cat("ARIMA-GARCH contribution is (1) clean residual diagnostics (Part 15)\n")
cat("and (2) a calibrated, volatility-adaptive interval — typically\n")
cat("narrower than the over-wide plain-ARIMA band while keeping coverage\n")
cat("near 95%.\n")

# ================================================================
# Part 18: ARIMA-GARCH FORECAST PLOT
# ================================================================
plot(test$Date, actual, type = "l", lwd = 2, col = "black",
     ylim = range(c(actual, loG, hiG, pred_arima)),
     main = paste0("Gold: ARIMA(", ord[1], ",", ord[2], ",", ord[3],
                   ")-GARCH(1,1) forecast with adaptive 95% band"),
     xlab = "Date", ylab = "Price (INR/10g)")
polygon(c(test$Date, rev(test$Date)), c(hiG, rev(loG)),
        col = rgb(0.7, 0.85, 1, 0.4), border = NA)
lines(test$Date, pred_garch, col = "#e63946", lwd = 2, lty = 2)
lines(test$Date, pred_arima, col = "#457b9d", lwd = 1.5, lty = 3)
legend("topleft", bty = "n", lwd = 2,
       legend = c("Actual", "ARIMA-GARCH forecast", "ARIMA forecast", "GARCH 95% band"),
       col = c("black", "#e63946", "#457b9d", rgb(0.7, 0.85, 1, 0.8)),
       lty = c(1, 2, 3, NA), pch = c(NA, NA, NA, 15))

cat("\n=== ARIMA-GARCH EXTENSION COMPLETE ===\n")












library(ggplot2)
library(dplyr)

df <- data.frame(
  Date = test$Date,
  Actual = actual,
  ARIMA_GARCH = pred_garch,
  ARIMA = pred_arima,
  lo = loG,
  hi = hiG
)

ord_lab <- paste0("(", ord[1], ",", ord[2], ",", ord[3], ")")

p <- ggplot(df, aes(x = Date)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = "GARCH 95% band"), alpha = 0.25) +
  geom_line(aes(y = Actual, color = "Actual"), linewidth = 0.9) +
  geom_line(aes(y = ARIMA_GARCH, color = "ARIMA-GARCH forecast"),
            linewidth = 0.7, linetype = "dashed") +
  geom_line(aes(y = ARIMA, color = "ARIMA forecast"),
            linewidth = 0.6, linetype = "dotted") +
  scale_color_manual(
    name = NULL,
    values = c("Actual" = "black",
               "ARIMA-GARCH forecast" = "#e63946",
               "ARIMA forecast" = "#1d3557")
  ) +
  scale_fill_manual(name = NULL, values = c("GARCH 95% band" = "#a8c8f0")) +
  labs(
    title = paste0("Gold: ARIMA", ord_lab, "-GARCH(1,1) Forecast"),
    subtitle = "with Adaptive 95% Confidence Band",
    x = "Date", y = "Price (INR/10g)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    legend.position = "top",
    legend.justification = "left",
    legend.box.spacing = unit(0.2, "lines"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "plain")
  ) +
  guides(fill = guide_legend(order = 2), color = guide_legend(order = 1))

print(p)
ggsave("gold_arima_garch_forecast.png", p, width = 8, height = 5.5, dpi = 300)



