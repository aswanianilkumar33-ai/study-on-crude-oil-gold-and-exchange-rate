##Aim / Objectives

The broad aim of this study is to model and forecast the prices of crude oil, gold, and the USD/INR exchange rate for the Indian economy, and to examine how these three series are related to one another over time. The specific objectives are:

**To identify the most suitable time series model for forecasting the monthly price of crude oil, evaluated through out-of-sample forecast accuracy.

**To determine the best-performing forecasting model for the daily price of gold by comparing ARIMA and exponential smoothing (ETS) approaches.

**To identify the most appropriate model for forecasting the daily USD/INR exchange rate, assessed against naive and random-walk benchmarks.

**To examine the long-run and short-run relationships among crude oil, gold, and the exchange rate, and the direction of influence between them.

**To compare the forecasting performance of classical time series models with machine learning models, and to assess whether the differences in accuracy are meaningful.

##Theoretical Framework of Topic

The study draws on the following core statistical and econometric concepts:

**Stationarity and Unit Root Testing: A time series is stationary if its statistical properties do not change over time. The Augmented Dickey-Fuller (ADF) and KPSS tests were used to check whether each series needed differencing before modelling.

**Box-Jenkins ARIMA/SARIMA Framework: A widely used approach for modelling and forecasting a series based on its own past values and past forecast errors, identified using the autocorrelation and partial autocorrelation functions and information criteria (AIC/BIC).

**Exponential Smoothing (ETS): An alternative forecasting approach that gives more weight to recent observations, used here as a benchmark against ARIMA.

**GARCH Modelling: Since financial series often show periods of high and low volatility clustered together, a GARCH(1,1) model with a Student-t distribution was used to describe this time-varying risk.

**Cointegration and the Vector Error Correction Model (VECM): Johansen's cointegration test was used to check whether crude oil, gold, and the exchange rate share a stable long-run relationship, and a VECM was used to study how each series adjusts back toward this equilibrium.

**Granger Causality and Impulse Response Analysis: Used to study the direction of influence between the three series and how a shock to one series affects the others over time.

**Machine Learning Models: Random Forest, XGBoost, and Long Short-Term Memory (LSTM) neural networks were used as modern alternatives to classical models, to test whether they could forecast gold and the exchange rate more accurately.

##Data Sources

The study uses four datasets, each drawn from an authoritative Indian government or exchange database:

**Crude oil price (monthly): Indian Basket crude oil price, in USD per barrel, from the Petroleum Planning and Analysis Cell (PPAC), Ministry of Petroleum and Natural Gas — April 2000 to May 2026 (314 observations).

**Gold price (daily): From the National Stock Exchange of India (NSE), in INR per 10 grams — January 2020 to May 2026 (1,537 observations).

**Gold price (monthly, for multivariate analysis): From the Reserve Bank of India (RBI) database, in INR per 10 grams — April 2000 to May 2026.

**USD/INR exchange rate (daily and monthly): From the historical data database on Investing.com — daily data from January 2020, and monthly data from April 2000, both to 2026.

##Software Used

**All classical time series modelling — ARIMA, SARIMA, GARCH, the Johansen cointegration test, the VECM, Granger causality, and impulse response analysis — was carried out in R, using the forecast, tseries, urca, vars, and rugarch packages. The machine learning models (Random Forest, XGBoost, and LSTM) were implemented in Python using scikit-learn, XGBoost, and TensorFlow/Keras.

