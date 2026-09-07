# ================================================================
#  USD/INR EXCHANGE RATE — ARIMA + GARCH ANALYSIS
#  Scope : model & forecast the RATE (ARIMA) and its VOLATILITY (GARCH)
#  Data  : daily, ~2020 to 2026
#
#  WHY BOTH MODELS (literature-backed):
#    - ARIMA forecasts the LEVEL of the rate.
#    - Daily financial series usually show volatility clustering, which
#      ARIMA assumes away. The ARCH-LM test checks for it.
#    - If ARCH is present (p < 0.05), GARCH models the time-varying
#      volatility. This ARMA/ARIMA + GARCH design follows the crude-oil
#      paper (ARMA-EGARCH) and gold studies (Ping 2013; Ye 2014) in the
#      project library. If ARCH is absent, ARIMA alone is sufficient.
#
#  SEQUENCE:
#    1 packages | 2 load+clean | 3 describe | 4 log returns
#    5 STL (seasonality check) | 6 train/test split
#    7 stationarity ADF+KPSS+PP -> d | 8 ACF/PACF
#    9 ARIMA on train + diagnostics | 10 forecast + accuracy
#   11 ETS on train (comparison model) | 12 ARIMA vs ETS vs benchmarks
#   13 ARCH test (decides if GARCH is needed)
#   14 ARIMA-GARCH on returns + diagnostics + interpretation
#   15 volatility forecast | 16 final 30-day rate forecast | 17 save
# ================================================================

# Gold is modelled on the log (its variance grows with the level).
# USD/INR spans a small range, so the log is OPTIONAL here.
USE_LOG <- FALSE     # FALSE = model the rate directly (recommended for FX)

# ---- 1. PACKAGES -----------------------------------------------
packages <- c("forecast", "tseries", "FinTS", "rugarch", "moments")
to_install <- packages[!packages %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
invisible(lapply(packages, library, character.only = TRUE))

# ---- 2. LOAD + CLEAN -------------------------------------------
df <- read.csv("USDINR.csv", stringsAsFactors = FALSE)
df <- df[, 1:2]
colnames(df) <- c("Date", "Rate")

parse_dates <- function(x) {                 # data mixes mm-dd-yyyy AND mm/dd/yyyy
  x <- trimws(as.character(x))
  x <- gsub("-", "/", x)                     # unify separators: "-" -> "/"
  d <- as.Date(x, format = "%m/%d/%Y")       # both formats are MONTH-day-year
  na <- is.na(d)                             # 2-digit-year fallback, just in case
  if (any(na)) d[na] <- as.Date(x[na], format = "%m/%d/%y")
  d
}
df$Date <- parse_dates(df$Date)
df$Rate <- as.numeric(gsub("[^0-9.]", "", as.character(df$Rate)))
df <- df[!is.na(df$Date) & !is.na(df$Rate), ]
df <- df[!duplicated(df$Date), ]
df <- df[order(df$Date), ]
cat("Data range :", format(min(df$Date)), "to", format(max(df$Date)), "| Obs:", nrow(df), "\n\n")

# ---- 3. DESCRIBE THE WHOLE SERIES ------------------------------
plot(df$Date, df$Rate, type = "l", col = "#1f77b4",
     main = "USD/INR Daily Exchange Rate", xlab = "Date", ylab = "Rate")
cat("=== DESCRIPTIVE STATISTICS ===\n")
cat("Mean :", round(mean(df$Rate),4), " Median:", round(median(df$Rate),4),
    " SD:", round(sd(df$Rate),4), "\n")
cat("Min  :", round(min(df$Rate),4), " Max   :", round(max(df$Rate),4),
    " CV:", round(sd(df$Rate)/mean(df$Rate)*100,2), "%\n")
cat("Skewness:", round(skewness(df$Rate),4), " Kurtosis:", round(kurtosis(df$Rate),4), "\n\n")

# ---- 4. LOG RETURNS (for GARCH) --------------------------------
df$LogReturn <- c(NA, diff(log(df$Rate)) * 100)   # percentage log returns

# ---- 5. STL DECOMPOSITION (seasonality check) ------------------
cat("=== DECOMPOSITION / SEASONALITY ===\n")
stl_fit <- stl(ts(df$Rate, frequency = 5), s.window = "periodic", robust = TRUE)
plot(stl_fit, main = " Decomposition - USD/INR")
rem <- stl_fit$time.series[, "remainder"]; sea <- stl_fit$time.series[, "seasonal"]; tre <- stl_fit$time.series[, "trend"]
cat("Trend strength    :", round(max(0, 1 - var(rem)/var(tre + rem)), 3), "\n")
cat("Seasonal strength :", round(max(0, 1 - var(rem)/var(sea + rem)), 3),
    "-> weak seasonality, use non-seasonal models\n\n")

# ---- 6. TRAIN / TEST SPLIT -------------------------------------
train <- df[df$Date <= as.Date("2025-12-31"), ]
test  <- df[df$Date >  as.Date("2025-12-31"), ]
actual <- test$Rate
h <- nrow(test)
cat("=== TRAIN / TEST SPLIT ===\n")
cat("Train:", nrow(train), "| Test:", nrow(test), "| h =", h, "\n\n")

LAMBDA       <- if (USE_LOG) 0 else NULL
model_series <- if (USE_LOG) log(train$Rate) else train$Rate
train_ts     <- ts(train$Rate, frequency = 1)    # frequency = 1 -> non-seasonal

# ---- 7. STATIONARITY (ADF + KPSS + PP) -> d --------------------
cat("=== STATIONARITY (train) ===\n")
adf_l <- suppressWarnings(adf.test(model_series))
kpss_l<- suppressWarnings(kpss.test(model_series))
pp_l  <- suppressWarnings(pp.test(model_series))
cat(sprintf("Level      ADF=%.4f(%s) KPSS=%.4f(%s) PP=%.4f(%s)\n",
            adf_l$p.value, ifelse(adf_l$p.value<0.05,"stat","NON"),
            kpss_l$p.value,ifelse(kpss_l$p.value<0.05,"NON","stat"),
            pp_l$p.value,  ifelse(pp_l$p.value<0.05,"stat","NON")))
diff_series <- diff(model_series)
adf_d <- suppressWarnings(adf.test(diff_series))
cat(sprintf("Difference ADF=%.4f(%s)\n", adf_d$p.value,
            ifelse(adf_d$p.value<0.05,"stat","NON")))
d_val <- if (adf_d$p.value < 0.05) 1 else 2
cat("=> d =", d_val, "\n\n")

# ---- 8. ACF / PACF ---------------------------------------------
par(mfrow = c(1, 2))
ggtsdisplay(
  diff(train_ts),
  lag.max =12 ,
  main = "ACF and PACF of Differenced "
)
acf (diff_series, lag.max = 40, main = "ACF (differenced)",  col = "#1f77b4", lwd = 2)
pacf(diff_series, lag.max = 40, main = "PACF (differenced)", col = "#d62728", lwd = 2)
par(mfrow = c(1, 1))

# ================================================================
# 9. ARIMA ON TRAIN
# ================================================================
cat("=== ARIMA (train) ===\n")
arima_fit <- auto.arima(train_ts, d = d_val, seasonal = FALSE,
                        stepwise = FALSE, approximation = FALSE,
                        ic = "aicc", lambda = LAMBDA)
print(arima_fit)
ord <- arimaorder(arima_fit)
cat(sprintf("auto.arima selected: ARIMA(%d,%d,%d)\n\n", ord[1], ord[2], ord[3]))

# Manual grid: fit EVERY ARIMA(p,1,q) for p,q in 0:3 and rank by AICc.
# This documents that many models were compared, not just the winner.
cat("Manual grid (all ARIMA(p,1,q), ranked by AICc; lowest = best):\n")
grid <- list()
for (p in 0:3) for (q in 0:3) {
  f <- tryCatch(Arima(train_ts, order = c(p, d_val, q), method = "ML", lambda = LAMBDA),
                error = function(e) NULL)
  if (!is.null(f))
    grid[[paste0(p, q)]] <- data.frame(Model = sprintf("ARIMA(%d,%d,%d)", p, d_val, q),
                                       AICc = round(f$aicc, 2))
}
grid_tab <- do.call(rbind, grid)
grid_tab <- grid_tab[order(grid_tab$AICc), ]
print(grid_tab, row.names = FALSE)
cat(sprintf("Grid best: %s  (should match auto.arima above)\n\n", grid_tab$Model[1]))

lb_a <- Box.test(residuals(arima_fit), lag = 20, type = "Ljung-Box", fitdf = ord[1] + ord[3])
cat(sprintf("Ljung-Box p = %.4f -> %s\n\n", lb_a$p.value,
            ifelse(lb_a$p.value > 0.05, "PASS (adequate)", "FAIL")))
shapiro.test(residuals(arima_fit))
print(lb_a)
ArchTest(residuals(arima_fit))
checkresiduals(arima_fit)
# ---- 10. ARIMA forecast + accuracy -----------------------------
arima_fc   <- forecast(arima_fit, h = h, level = 95)
arima_pred <- as.numeric(arima_fc$mean)
cmp_arima <- data.frame(
  Date = test$Date, Actual = round(actual, 4), Forecast = round(arima_pred, 4),
  Error = round(actual - arima_pred, 4),
  Pct_Err = round((actual - arima_pred)/actual*100, 4),
  Lo95 = round(as.numeric(arima_fc$lower[,1]),4),
  Hi95 = round(as.numeric(arima_fc$upper[,1]),4))
cmp_arima$In_95CI <- cmp_arima$Actual >= cmp_arima$Lo95 & cmp_arima$Actual <= cmp_arima$Hi95
cat("=== ARIMA: ACTUAL vs FORECAST (head) ===\n")
print(head(cmp_arima, 10), row.names = FALSE)
mae_a <- mean(abs(cmp_arima$Error)); rmse_a <- sqrt(mean(cmp_arima$Error^2)); mape_a <- mean(abs(cmp_arima$Pct_Err))
cat(sprintf("\nARIMA  MAE=%.4f RMSE=%.4f MAPE=%.4f%% | coverage=%.1f%%\n\n",
            mae_a, rmse_a, mape_a, mean(cmp_arima$In_95CI)*100))


# ================================================================
# 11. ETS ON TRAIN (comparison model)
# ================================================================
cat("=== ETS (train) ===\n")
ets_fit  <- ets(train_ts, lambda = LAMBDA)
ets_pred <- as.numeric(forecast(ets_fit, h = h)$mean)
cat("ETS picked:", ets_fit$method, "\n")
mae_e <- mean(abs(actual-ets_pred)); rmse_e <- sqrt(mean((actual-ets_pred)^2))
mape_e <- mean(abs((actual-ets_pred)/actual))*100
cat(sprintf("ETS    MAE=%.4f RMSE=%.4f MAPE=%.4f%%\n\n", mae_e, rmse_e, mape_e))

# ================================================================
# 12. ARIMA vs ETS vs BENCHMARKS
# ================================================================
last_val <- tail(train$Rate, 1); drift <- mean(diff(train$Rate))
naive_pred <- rep(last_val, h); rw_pred <- last_val + drift*(1:h)
acc_one <- function(p) c(MAE=mean(abs(actual-p)), RMSE=sqrt(mean((actual-p)^2)),
                         MAPE=mean(abs((actual-p)/actual))*100)
comparison <- round(rbind(ARIMA=acc_one(arima_pred), ETS=acc_one(ets_pred),
                          Naive=acc_one(naive_pred), RandomWalk=acc_one(rw_pred)), 4)
comparison <- comparison[order(comparison[,"RMSE"]), ]
cat("=== MODEL COMPARISON (ranked by RMSE) ===\n")
print(comparison); cat("\nBest:", rownames(comparison)[1], "\n\n")

plot(test$Date, actual, type = "l", lwd = 2, col = "black",
     ylim = range(c(actual, arima_pred, ets_pred)),
     main = "USD/INR: Actual vs Forecast", xlab = "Date", ylab = "Rate")
lines(test$Date, arima_pred, col = "#1f77b4", lwd = 2, lty = 2)
lines(test$Date, ets_pred,   col = "red", lwd = 2, lty = 3)
legend("topleft", bty="n", lwd=2, legend=c("Actual","ARIMA","ETS"),
       col=c("black","#1f77b4","red"), lty=c(1,2,3))

# ================================================================
# 13. ARCH TEST  -> decides whether GARCH is needed
# ================================================================
cat("=== ARCH-LM TEST (on ARIMA residuals) ===\n")
arch_test <- ArchTest(residuals(arima_fit), lags = 12)
garch_needed <- arch_test$p.value < 0.05
cat(sprintf("ARCH-LM p = %.4f -> %s\n\n", arch_test$p.value,
            ifelse(garch_needed, "ARCH PRESENT: fit GARCH (below)",
                   "NO ARCH: ARIMA alone is sufficient, GARCH optional")))

# ================================================================
# 14. ARIMA-GARCH(1,1) ON RETURNS  (only meaningful if ARCH present)
# ----------------------------------------------------------------
# NOTE ON SCOPE (state this in the report): this GARCH is a DESCRIPTIVE
# volatility model fitted on the FULL-sample log-returns. Its job is to
# characterise volatility clustering over the whole history and to produce
# the forward 30-day volatility forecast in Section 15 -- it is NOT evaluated
# out-of-sample the way the ARIMA mean model is in Sections 10-12. So fitting
# it on all the data is deliberate, not a train/test "leakage" error: there is
# no test-set accuracy being claimed for the GARCH, only a description of risk
# and a forward forecast. The ARCH-LM test in Section 13 (on TRAIN residuals)
# is what justifies fitting it at all.
# ================================================================
cat("=== ARIMA-GARCH(1,1) ===\n")
returns <- train$LogReturn[!is.na(train$LogReturn)]
garch_spec <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(ord[1], ord[3]), include.mean = TRUE),
  distribution.model = "std")                    # Student-t for fat tails
garch_fit <- ugarchfit(spec = garch_spec, data = returns, solver = "hybrid")
mean(returns)
cat("\n=== RETURN SERIES (for volatility modelling) ===\n")
cat("Mean :", round(mean(returns), 3), "%   SD :",
    round(sd(returns), 3), "%\n")
cat("Skewness:", round(moments::skewness(returns), 3),
    " Kurtosis:", round(moments::kurtosis(returns), 3), "\n")
cat("INTERPRETATION: Kurtosis far above 3 = HEAVY TAILS (extreme months happen\n")
cat("more often than a normal bell curve predicts) -> use a Student-t below.\n")

cf <- coef(garch_fit); omega <- cf["omega"]; alpha <- cf["alpha1"]; beta <- cf["beta1"]
persistence <- alpha + beta
cat(sprintf("alpha=%.4f beta=%.4f persistence=%.4f\n", alpha, beta, persistence))
cat("Volatility persistence:",
    if (persistence > 0.95) "very high (shocks linger for a long time)."
    else if (persistence > 0.85) "high (shocks fade over weeks)."
    else "moderate (volatility settles quickly).", "\n")
cat("Long-run annualised volatility:",
    round(sqrt(omega/(1 - alpha - beta)) * sqrt(252), 2), "%\n")     # x sqrt(252), corrected

std_resid <- as.numeric(residuals(garch_fit, standardize = TRUE))
z <- as.numeric(residuals(garch_fit, standardize = TRUE))

lb_z   <- Box.test(z,   lag = 20, type = "Ljung-Box")
lb_z2  <- Box.test(z^2, lag = 20, type = "Ljung-Box")
arch_z <- ArchTest(z, lags = 12)
lb_z
lb_z2
arch_z
print(garch_fit)
par(mfrow = c(2, 2))
plot(z, type = "l", col = "steelblue", main = "Standardised Residuals",
     ylab = "z", xlab = "Time"); abline(h = c(-2, 0, 2), col = c("red","black","red"), lty = 2)
# Histogram overlaid with the F
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


# ================================================================
# 15. GARCH 30-DAY VOLATILITY FORECAST
# ================================================================
gfc <- ugarchforecast(garch_fit, n.ahead = 30)
vol_fore <- data.frame(
  Day = 1:30, Date = seq(max(df$Date)+1, by="day", length.out=30),
  DailyVol_pct = round(as.numeric(sigma(gfc)), 4),
  AnnualisedVol = round(as.numeric(sigma(gfc)) * sqrt(252), 4))
cat("=== 30-DAY VOLATILITY FORECAST (weekly) ===\n")
print(vol_fore[seq(7, 28, by = 7), c("Date","DailyVol_pct","AnnualisedVol")], row.names = FALSE)
cat("\n")

# ================================================================
# 16. FINAL ARIMA -> 30-DAY RATE FORECAST (refit on full data)
# ================================================================
cat("=== 30-DAY RATE FORECAST ===\n")
arima_final <- auto.arima(ts(df$Rate, frequency = 1), d = d_val, seasonal = FALSE,
                          stepwise = FALSE, approximation = FALSE, ic = "aicc", lambda = LAMBDA)
ffc <- forecast(arima_final, h = 30, level = 95)
fut_df <- data.frame(
  Day = 1:30, Date = seq(max(df$Date)+1, by="day", length.out=30),
  Forecast = round(as.numeric(ffc$mean),4),
  Lo95 = round(as.numeric(ffc$lower[,1]),4), Hi95 = round(as.numeric(ffc$upper[,1]),4))
print(fut_df[seq(7, 28, by = 7), ], row.names = FALSE)


# ================================================================
# 17. SAVE
# ================================================================
write.csv(cmp_arima,  "usdinr_arima_test_comparison.csv", row.names = FALSE)
write.csv(comparison, "usdinr_model_accuracy.csv")
write.csv(vol_fore,   "usdinr_garch_volatility_forecast.csv", row.names = FALSE)
write.csv(fut_df,     "usdinr_30day_forecast.csv", row.names = FALSE)
cat("\nSaved 4 CSV files.\n=== DONE ===\n")
write.csv(df[, c("Date","Rate")], "usdinr_clean.csv", row.names = FALSE)

