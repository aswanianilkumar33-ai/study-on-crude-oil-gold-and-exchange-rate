# ============================================================
#  COMPLETE ANALYSIS: CRUDE OIL, GOLD & USD/INR EXCHANGE RATE
#  Monthly Data | Financial Years 2000-01 to 2025-26
#  All packages explicitly namespaced to avoid conflicts
# ============================================================


# ============================================================
# STEP 0: INSTALL & LOAD PACKAGES
# ============================================================

pkgs <- c("dplyr","tidyr","stringr","readxl","ggplot2",
          "tseries","urca","vars","lmtest","sandwich","car")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(ggplot2)
library(tseries)   # adf.test
library(urca)      # ca.jo (Johansen)
library(vars)      # VAR, causality, irf, fevd, vec2var
library(lmtest)    # grangertest, bptest, dwtest, coeftest
library(sandwich)  # vcovHC (robust SE)
library(car)       # vif


# ============================================================
# STEP 1: LOAD & RESHAPE DATA
# ============================================================
# Excel file with 3 sheets. Each sheet has 13 columns:
#   Col 1 : Year  (e.g. "2000-01", "2001-02", ...)
#   Col 2-13: April May June July August September
#             October November December January February March
#
# Financial year logic:
#   April–December → calendar year = first 4 digits of Year
#   January–March  → calendar year = first 4 digits + 1

# *** CHANGE FILE PATH AND SHEET NAMES TO MATCH YOUR FILE ***
crude_raw <- read.csv("crudeoilmonthly.csv", stringsAsFactors = FALSE)
gold_raw <- read.csv("goldmonthly.csv",  stringsAsFactors = FALSE)
fx_raw <- read.csv("usdinrmonthly.csv", stringsAsFactors = FALSE)
# Month names exactly as they appear in your Excel headers
month_order  <- c("April","May","June","July","August","September",
                  "October","November","December","January","February","March")

# Map each month name to its calendar month number
month_to_num <- c(April=4L, May=5L, June=6L, July=7L,
                  August=8L, September=9L, October=10L,
                  November=11L, December=12L,
                  January=1L, February=2L, March=3L)
gold_raw <- read.csv("goldmonthly.csv", stringsAsFactors = FALSE)
colnames(gold_raw) <- trimws(colnames(gold_raw))
print(colnames(gold_raw))
print(head(gold_raw))
reshape_to_monthly <- function(raw_df, value_name) {
  colnames(raw_df) <- trimws(colnames(raw_df))
  colnames(raw_df)[1] <- "Year"
  
  long <- raw_df %>%
    dplyr::mutate(Year = as.character(Year)) %>%
    dplyr::filter(!is.na(Year), trimws(Year) != "") %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(month_order),
      names_to  = "month_name",
      values_to = "val"
    ) %>%
    dplyr::mutate(
      # Strip commas before parsing as numeric — fixes "4,460.00" -> 4460.00
      val       = suppressWarnings(as.numeric(gsub(",", "", val))),
      base_year = as.integer(stringr::str_extract(Year, "\\d{4}")),
      month_num = month_to_num[month_name],
      cal_year  = dplyr::if_else(month_num >= 4L, base_year, base_year + 1L),
      date      = as.Date(paste(cal_year,
                                formatC(month_num, width=2, flag="0"),
                                "01", sep = "-"))
    ) %>%
    dplyr::filter(!is.na(date), !is.na(val)) %>%
    dplyr::select(date, val) %>%
    dplyr::rename(!!value_name := val) %>%
    dplyr::arrange(date)
  
  return(long)
}

crude_df <- reshape_to_monthly(crude_raw, "crude")
gold_df  <- reshape_to_monthly(gold_raw,  "gold")
fx_df    <- reshape_to_monthly(fx_raw,    "fx")

# Diagnostic
cat("Crude:", nrow(crude_df), "rows |",
    as.character(min(crude_df$date)), "to", as.character(max(crude_df$date)), "\n")
cat("Gold: ", nrow(gold_df),  "rows |",
    as.character(min(gold_df$date)),  "to", as.character(max(gold_df$date)),  "\n")
cat("FX:   ", nrow(fx_df),    "rows |",
    as.character(min(fx_df$date)),    "to", as.character(max(fx_df$date)),    "\n")

# Merge all three on date
df <- crude_df %>%
  dplyr::inner_join(gold_df, by = "date") %>%
  dplyr::inner_join(fx_df,   by = "date") %>%
  dplyr::arrange(date)

cat("\n Merged rows:", nrow(df), "\n")
cat("Date range:", as.character(min(df$date)),
    "to", as.character(max(df$date)), "\n")
print(head(df))


# ============================================================
# STEP 2: CREATE TIME SERIES & LOG TRANSFORM
# ============================================================
# Log transform stabilises variance and allows elasticity interpretation.
# A 1% change in crude → β1 % change in FX rate.

start_yr  <- as.integer(format(min(df$date), "%Y"))
start_mon <- as.integer(format(min(df$date), "%m"))

crude_ts <- ts(df$crude, start = c(start_yr, start_mon), frequency = 12)
gold_ts  <- ts(df$gold,  start = c(start_yr, start_mon), frequency = 12)
fx_ts    <- ts(df$fx,    start = c(start_yr, start_mon), frequency = 12)

lcrude <- log(crude_ts)
lgold  <- log(gold_ts)
lfx    <- log(fx_ts)


# ============================================================
# STEP 3: DESCRIPTIVE STATISTICS
# ============================================================
cat("\n========== DESCRIPTIVE STATISTICS ==========\n")
cat("Gives basic properties of each variable before any tests.\n\n")

for (nm in c("crude","gold","fx")) {
  x <- df[[nm]]
  cat(sprintf("%-6s | Mean: %8.2f | Median: %8.2f | SD: %7.2f | Min: %7.2f | Max: %9.2f\n",
              nm, mean(x,na.rm=T), median(x,na.rm=T),
              sd(x,na.rm=T), min(x,na.rm=T), max(x,na.rm=T)))
}


# ============================================================
# STEP 4: PLOT TIME SERIES
# ============================================================
# Visual check for trends, breaks, and volatility clusters.
# Upward trend = likely non-stationary (important for Step 5).

df_long <- df %>%
  tidyr::pivot_longer(cols = c(crude, gold, fx),
                      names_to  = "variable",
                      values_to = "value") %>%
  dplyr::mutate(variable = dplyr::recode(variable,
                  crude = "Crude Oil (USD/bbl)",
                  gold  = "Gold (INR/10g or USD/oz)",
                  fx    = "USD/INR Exchange Rate"))

print(
  ggplot(df_long, aes(x = date, y = value, colour = variable)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~variable, scales = "free_y", ncol = 1) +
    labs(title = "Monthly Crude Oil, Gold & USD/INR (2000–2026)",
         x = "Date", y = "Value") +
    theme_minimal() +
    theme(legend.position = "none")
)


# ============================================================
# STEP 5: CORRELATION ANALYSIS
# ============================================================
# Quick measure of linear association.
# Interpretation: +1 = perfect positive, -1 = perfect negative, 0 = none.
# NOTE: correlation does not imply causation.

cat("\n========== PEARSON CORRELATION (Raw Levels) ==========\n")
print(round(cor(df[, c("crude","gold","fx")], use = "complete.obs"), 4))

cat("\n--- Correlation (Log Levels) ---\n")
log_df <- data.frame(lcrude = as.numeric(lcrude),
                     lgold  = as.numeric(lgold),
                     lfx    = as.numeric(lfx))
print(round(cor(log_df, use = "complete.obs"), 4))


# ============================================================
# STEP 6: UNIT ROOT TESTS (Stationarity)
# ============================================================
# MUST DO before regression or VAR.
# Non-stationary series give spurious (fake) regression results.
#
# ADF Test:
#   H0: Series has a unit root (non-stationary)
#   p < 0.05 → stationary
#   p > 0.05 → non-stationary (has unit root)
#
# I(1) = non-stationary at levels, stationary after 1st difference.
# Most price/FX series are I(1).

cat("\n========== ADF UNIT ROOT TESTS ==========\n")

adf_summary <- function(x, label) {
  cat("\n---", label, "---\n")
  t_level <- tseries::adf.test(x, alternative = "stationary")
  t_diff  <- tseries::adf.test(diff(x), alternative = "stationary")
  cat("Level    p-value:", round(t_level$p.value, 4),
      "→", ifelse(t_level$p.value < 0.05, "STATIONARY I(0)", "NON-STATIONARY"), "\n")
  cat("1st Diff p-value:", round(t_diff$p.value, 4),
      "→", ifelse(t_diff$p.value < 0.05, "STATIONARY I(1)", "STILL NON-STATIONARY"), "\n")
}

adf_summary(lcrude, "Log Crude Oil")
adf_summary(lgold,  "Log Gold Price")
adf_summary(lfx,    "Log USD/INR")

cat("\nINTERPRETATION: If all three are I(1),")
cat(" proceed to Johansen cointegration test.\n")


# ============================================================
# STEP 7: JOHANSEN COINTEGRATION TEST (Long-Run Relationship)
# ============================================================
# Tests whether I(1) variables share a long-run equilibrium.
# Even if they drift short-term, cointegration = tied together long-run.
#
# Trace test & Max-Eigen test:
#   H0 (r=0): No cointegrating vectors
#   Reject if test statistic > critical value at 5%
#
# If cointegrated → use VECM (Step 11)
# If NOT cointegrated → use VAR in first differences (Step 8)

cat("\n========== JOHANSEN COINTEGRATION TEST ==========\n")

log_mat <- cbind(lcrude, lgold, lfx)

lag_sel <- vars::VARselect(log_mat, lag.max = 12, type = "const")
cat("Lag selection (AIC / BIC / HQ / FPE):\n")
print(lag_sel$selection)
opt_lag <- max(as.integer(lag_sel$selection["AIC(n)"]), 2L)
cat("Using lag:", opt_lag, "\n")

joh_trace <- urca::ca.jo(log_mat, type = "trace", ecdet = "const", K = opt_lag)
joh_eigen <- urca::ca.jo(log_mat, type = "eigen", ecdet = "const", K = opt_lag)

cat("\n--- Trace Test ---\n")
print(summary(joh_trace))

cat("\n--- Max-Eigenvalue Test ---\n")
print(summary(joh_eigen))

cat("\nINTERPRETATION:")
cat("\n  Test stat > 10pct critical value → reject H0 → cointegration exists.")
cat("\n  If r>=1 confirmed → long-run relationship exists → use VECM.\n")


# ============================================================
# STEP 8: VAR MODEL (Vector Autoregression)
# ============================================================
# Models each variable as a function of its own and others' past values.
# No variable is pre-assigned as dependent — all are endogenous.
# Run on FIRST DIFFERENCES of log series (growth rates).
#
# Used as basis for:
#   - Granger causality (Step 9)
#   - Impulse Response Functions (Step 10)
#   - Variance Decomposition (Step 10)

cat("\n========== VAR MODEL (First Differences) ==========\n")

dlcrude <- diff(lcrude)
dlgold  <- diff(lgold)
dlfx    <- diff(lfx)

diff_mat <- cbind(dlcrude, dlgold, dlfx)
colnames(diff_mat) <- c("dlcrude", "dlgold", "dlfx")

lag_sel2 <- vars::VARselect(diff_mat, lag.max = 12, type = "const")
cat("Optimal lags:\n"); print(lag_sel2$selection)
best_lag <- max(as.integer(lag_sel2$selection["AIC(n)"]), 1L)
cat("Using lag:", best_lag, "\n")

var_model <- vars::VAR(diff_mat, p = best_lag, type = "const")
print(summary(var_model))

# Stability: all eigenvalue roots must be < 1
r <- vars::roots(var_model)
cat("\nVAR roots:", round(r, 4), "\n")
cat(ifelse(all(r < 1), "VAR is STABLE\n", "WARNING: VAR is UNSTABLE\n"))


# ============================================================
# STEP 9: GRANGER CAUSALITY TESTS
# ============================================================
# Tests whether past values of X improve prediction of Y.
# This is predictive causality, not true economic causation.
#
# H0: X does NOT Granger-cause Y
# p < 0.05 → X DOES Granger-cause Y
#
# We test all 6 directions among the three variables.

cat("\n========== GRANGER CAUSALITY TESTS ==========\n")
cat("p < 0.05 = causality exists\n\n")

diff_df <- as.data.frame(diff_mat)

# Using VAR-based causality (multivariate — controls for all variables)
cat("--- Crude causes [Gold, FX]? ---\n")
print(vars::causality(var_model, cause = "dlcrude")$Granger)

cat("\n--- Gold causes [Crude, FX]? ---\n")
print(vars::causality(var_model, cause = "dlgold")$Granger)

cat("\n--- FX causes [Crude, Gold]? ---\n")
print(vars::causality(var_model, cause = "dlfx")$Granger)

# Bivariate Granger tests (pairwise, more specific)
cat("\n--- Bivariate: Crude → FX ---\n")
print(lmtest::grangertest(dlfx ~ dlcrude, order = best_lag, data = diff_df))

cat("\n--- Bivariate: FX → Crude ---\n")
print(lmtest::grangertest(dlcrude ~ dlfx, order = best_lag, data = diff_df))

cat("\n--- Bivariate: Gold → FX ---\n")
print(lmtest::grangertest(dlfx ~ dlgold, order = best_lag, data = diff_df))

cat("\n--- Bivariate: FX → Gold ---\n")
print(lmtest::grangertest(dlgold ~ dlfx, order = best_lag, data = diff_df))

cat("\n--- Bivariate: Crude → Gold ---\n")
print(lmtest::grangertest(dlgold ~ dlcrude, order = best_lag, data = diff_df))

cat("\n--- Bivariate: Gold → Crude ---\n")
print(lmtest::grangertest(dlcrude ~ dlgold, order = best_lag, data = diff_df))


# ============================================================
# STEP 10: IMPULSE RESPONSE & VARIANCE DECOMPOSITION
# ============================================================
# IRF: How does a shock in one variable affect another over time?
#   x-axis = months ahead; y-axis = response size
#   If confidence band includes zero → effect not significant
#
# FEVD: What % of a variable's forecast variance is due to each variable?
#   Helps identify dominant drivers in the system.

cat("\n========== IMPULSE RESPONSE FUNCTIONS ==========\n")

irf1 <- vars::irf(var_model, impulse="dlcrude", response="dlfx",
                  n.ahead=24, boot=TRUE, ci=0.95)
plot(irf1, main="IRF: USD/INR response to Crude Oil shock")

irf2 <- vars::irf(var_model, impulse="dlgold", response="dlfx",
                  n.ahead=24, boot=TRUE, ci=0.95)
plot(irf2, main="IRF: USD/INR response to Gold Price shock")

irf3 <- vars::irf(var_model, impulse="dlfx", response="dlcrude",
                  n.ahead=24, boot=TRUE, ci=0.95)
plot(irf3, main="IRF: Crude Oil response to USD/INR shock")

cat("\n========== FORECAST ERROR VARIANCE DECOMPOSITION ==========\n")
cat("Shows % of forecast variance in each variable explained by others.\n")
fevd_res <- vars::fevd(var_model, n.ahead = 24)
print(fevd_res)
plot(fevd_res, main = "Variance Decomposition")


# ============================================================
# STEP 11: OLS REGRESSION (USD/INR as Dependent Variable)
# ============================================================
# Static long-run regression.
# Model: log(FX) = b0 + b1*log(Crude) + b2*log(Gold) + error
#
# Coefficients are ELASTICITIES (log-log model):
#   b1 = % change in FX for 1% rise in crude oil
#   b2 = % change in FX for 1% rise in gold
#   R² = % of FX variation explained by crude + gold
#
# Valid if series are cointegrated (Engle-Granger sense).

cat("\n========== OLS REGRESSION ==========\n")
cat("Dependent: Log(USD/INR) | Independent: Log(Crude), Log(Gold)\n\n")

ols <- lm(lfx ~ lcrude + lgold)
print(summary(ols))

cat("\n--- Robust Standard Errors (HC1) ---\n")
cat("Corrects SE for heteroscedasticity.\n")
print(lmtest::coeftest(ols, vcov = sandwich::vcovHC(ols, type = "HC1")))

cat("\n--- VIF (Multicollinearity Check) ---\n")
cat("VIF > 10 = serious multicollinearity between crude and gold.\n")
print(car::vif(ols))

cat("\n--- Breusch-Pagan Test (Heteroscedasticity) ---\n")
cat("H0: Constant variance. p < 0.05 → heteroscedasticity present.\n")
print(lmtest::bptest(ols))

cat("\n--- Durbin-Watson Test (Autocorrelation) ---\n")
cat("DW near 2 = no autocorrelation. DW < 1.5 or > 2.5 = problem.\n")
print(lmtest::dwtest(ols))

par(mfrow = c(2,2))
plot(ols, main = "OLS Diagnostics")
par(mfrow = c(1,1))


# ============================================================
# STEP 12: ENGLE-GRANGER COINTEGRATION (Pairwise)
# ============================================================
# Simpler alternative to Johansen for 2-variable cointegration.
# Method: Run OLS, then test if residuals are stationary.
# Stationary residuals → variables are cointegrated (long-run link).
#
# H0: Residuals are non-stationary (no cointegration)
# p < 0.05 → cointegration exists

cat("\n========== ENGLE-GRANGER COINTEGRATION (Pairwise) ==========\n")

cat("\n--- FX ~ Crude Oil ---\n")
eg1  <- lm(lfx ~ lcrude)
cat("ADF on residuals (p < 0.05 = cointegrated):\n")
print(tseries::adf.test(residuals(eg1)))

cat("\n--- FX ~ Gold ---\n")
eg2  <- lm(lfx ~ lgold)
cat("ADF on residuals (p < 0.05 = cointegrated):\n")
print(tseries::adf.test(residuals(eg2)))

cat("\n--- Crude ~ Gold ---\n")
eg3  <- lm(lcrude ~ lgold)
cat("ADF on residuals (p < 0.05 = cointegrated):\n")
print(tseries::adf.test(residuals(eg3)))


# ============================================================
# STEP 13: VECTOR ERROR CORRECTION MODEL (VECM)
# ============================================================
# USE ONLY IF Johansen test (Step 7) confirmed cointegration.
#
# VECM adds the long-run equilibrium error (ECT) into the VAR.
# ECT (error correction term):
#   Must be NEGATIVE and SIGNIFICANT for valid long-run correction.
#   e.g., ECT = -0.12 → 12% of last period's disequilibrium
#   is corrected each month.
#
# Short-run coefficients: immediate effects of changes
# Long-run (cointegrating vector): structural equilibrium

cat("\n========== VECM (run only if cointegration confirmed) ==========\n")

# r = number of cointegrating vectors from Johansen test
# Change r = 1 if Johansen showed exactly 1 cointegrating vector
vecm      <- urca::cajorls(joh_trace, r = 1)
print(vecm)

vecm_var  <- vars::vec2var(joh_trace, r = 1)
print(summary(vecm_var))

cat("\nINTERPRETATION:")
cat("\n  ect1 coefficient should be negative and significant.")
cat("\n  Negative ect1 confirms long-run adjustment back to equilibrium.")
cat("\n  Magnitude = speed of adjustment per month.\n")


# ============================================================
# STEP 14: SUMMARY
# ============================================================
cat("\n")
cat("================================================================\n")
cat("                    RESULTS DECISION GUIDE\n")
cat("================================================================\n")
cat("
Step 6  Unit Root     → All I(1)? → proceed to cointegration
Step 7  Johansen      → Cointegrated? → use VECM (Step 13)
                      → Not cointegrated? → use VAR in diffs (Step 8)
Step 8  VAR           → Foundation for Steps 9 & 10
Step 9  Granger       → Which variable predicts which?
Step 10 IRF / FEVD    → How long do shocks last? Who drives whom?
Step 11 OLS           → Elasticity estimates (if cointegrated)
Step 12 Engle-Granger → Pairwise long-run link check
Step 13 VECM          → Speed of long-run adjustment

EXPECTED FINDING (based on literature):
  Dependent variable   = USD/INR Exchange Rate
  Independent variables= Crude Oil Price, Gold Price
  Both crude oil and gold have significant long-run effects on Rupee.
  Granger causality likely runs: Crude → FX and Gold → FX.
")
cat("================================================================\n")

