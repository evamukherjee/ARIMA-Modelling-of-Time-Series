rm(list=ls(all=TRUE))
library(ggplot2) 
library(zoo)
library(readr) 
library(lubridate)
library(forecast)
library(tsoutliers) # Load the tsoutliers package
library(tseries)
library(astsa) 
library(portes)
library(dplyr)

# Load and prepare data
valeurs_mensuelles <- read_delim("valeurs_mensuelles.csv", delim = ";")
View(valeurs_mensuelles) 
base <- valeurs_mensuelles[-c(1:3),c(1:2)] 
colnames(base) <- c("Period","Index_CVS_CJO") 
View(base)
base$Index_CVS_CJO <- as.numeric(base$Index_CVS_CJO)
base$Period <- as.Date(paste0(base$Period, "-01"))

qplot(y=base$Index_CVS_CJO,  x  =as.yearmon(base[[1]]),  geom  =  c("point","smooth","line"),  xlab  =  "Period",  ylab  =  "Index of production of cocoa, chocolate and confectionary products") 
boxplot(base$Index_CVS_CJO ~ month(base$Period),
        col = "lightblue", # Color of the boxes
        pch = 20,          # Plotting character for outliers
        cex = 0.5,         # Size of outliers
        main = "Monthly Boxplots of Production Index", # Title
        ylab = "Index (Base 100 in 2021)",            # Y-axis label
        xlab = "Month") 

# Convert to time series object
start_year <- year(base$Period[1])
start_month <- month(base$Period[1])
ts_data <- ts(base$Index_CVS_CJO, start = c(start_year, start_month), frequency = 12)

# Decompose time series
ts_data_dec <- decompose(ts_data, type = "multiplicative")
autoplot(ts_data_dec, main = "Decomposition of the CVS-CJO Series", xlab = "Year")

# Question 2

# Stationarity tests
adf_test_result <- adf.test(base$Index_CVS_CJO)
print(adf_test_result)
kpss_test_result <- kpss.test(base$Index_CVS_CJO)
print(kpss_test_result)

# Removing deterministic trend based on ADF and KPSS test results
time_index <- time(ts_data)
trend_model <- lm(ts_data ~ time_index)
detrended_numeric <- residuals(trend_model) # Store as numeric first

# Convert detrended_numeric back to a time series object
detrended_ts_data <- ts(detrended_numeric, start = start(ts_data), frequency = frequency(ts_data))

adf.test(detrended_ts_data)  # Expect p < 0.05 (stationary)
kpss.test(detrended_ts_data, null = "Level")  # Expect p > 0.05 (stationary)
kpss.test(detrended_ts_data, null = "Trend")

#Question 3

par(mfrow=c(2,1))

plot(y=base$Index_CVS_CJO, x =as.yearmon(base[[1]]),type = "line", main="Original Series", xlab = "", ylab = "Index of production")

plot(y=detrended_ts_data, x =as.yearmon(base[[1]]),type = "line", main="Stationary Series", xlab = "", ylab = "Index(linear detrended) of production")

par(mfrow=c(1,1))

##PART II

# Visualize ACF and PACF of the transformed series to help identify ARMA orders (useful for Part 2)
# Visualize ACF and PACF of the transformed (detrended) series to help identify ARMA orders (useful for Part 2)
acf(detrended_ts_data, main = "ACF of Detrended Series", lag.max = 36)
pacf(detrended_ts_data, main = "PACF of Detrended Series", lag.max = 36)

# ACF and PACF remind us of an AR(3) and an MA(2). We will test all ARMA models such that
#p<=3 and q<=2
model1=sarima(base$Index_CVS_CJO, 3, 0, 0)
model2=sarima(base$Index_CVS_CJO, 3, 0, 1)
model3=sarima(base$Index_CVS_CJO, 3, 0, 2)
model4=sarima(base$Index_CVS_CJO, 2, 0, 0)
model5=sarima(base$Index_CVS_CJO, 2, 0, 1)
model6=sarima(base$Index_CVS_CJO, 2, 0, 2)
model7=sarima(base$Index_CVS_CJO, 1, 0, 0)
model8=sarima(base$Index_CVS_CJO, 1, 0, 1)
model9=sarima(base$Index_CVS_CJO, 1, 0, 2)
model10=sarima(base$Index_CVS_CJO, 0, 0, 1)
model11=sarima(base$Index_CVS_CJO, 0, 0, 2)

# Estimated coefficients 

model1$ttable
model2$ttable
model3$ttable
model4$ttable
model5$ttable
model6$ttable
model7$ttable
model8$ttable
model9$ttable
model10$ttable
model11$ttable


# Ljung-Box test (non-autocorrelation of residuals) 
#firgure out how many lags to do
LjungBox(model1$fit)
LjungBox(model2$fit)
LjungBox(model3$fit)
LjungBox(model4$fit)
LjungBox(model5$fit)
LjungBox(model6$fit)
LjungBox(model7$fit)
LjungBox(model8$fit)
LjungBox(model9$fit)
LjungBox(model10$fit)
LjungBox(model11$fit)



# Tests joints et visualisation des résidus
checkresiduals(model1$fit)
checkresiduals(model2$fit)
checkresiduals(model3$fit)
checkresiduals(model4$fit)
checkresiduals(model5$fit)
checkresiduals(model6$fit)
checkresiduals(model7$fit)
checkresiduals(model8$fit)
checkresiduals(model9$fit)
checkresiduals(model10$fit)
checkresiduals(model11$fit)

# Selecting the best model based on the information criteria 
# AIC
aic <- AIC(model1$fit, model2$fit, model3$fit, model4$fit, model5$fit,
           model6$fit, model7$fit, model8$fit, model9$fit, model10$fit, model11$fit)
aic
which.min(aic$AIC)

#BIC
bic <- BIC(model1$fit, model2$fit, model3$fit, model4$fit, model5$fit,
           model6$fit, model7$fit, model8$fit, model9$fit, model10$fit, model11$fit)
bic
which.min(bic$BIC)

#Model 8 which is ARIMA(1,0,1) shows the best results for AIC and BIC criterion

#Question 5

# Fit selected ARIMA model 
arima_finalmodel <- arima(detrended_ts_data, order = c(1, 0, 1)) 

final_model <- arima_finalmodel

# To find final model coefficients
summary(final_model)

#PART III

#Question 6 and 7
library(ellipse)

# Step 1: Forecast 2 steps ahead (using detrended series and d=0)
final_model_clean <- Arima(detrended_ts_data, order = c(0, 0, 1))

# Now forecast 2 steps ahead safely
Pred <- predict(final_model_clean, n.ahead = 2, se.fit = TRUE)

# Step 2: Compute variance-covariance matrix of (X_{T+1}, X_{T+2})
sigma2 <- final_model$sigma2
theta1 <- ifelse("ma1" %in% names(final_model$coef), final_model$coef["ma1"], 0)

# Variance of 1-step and 2-step ahead forecasts
sigma_g1 <- sqrt(sigma2)
sigma_g2 <- sqrt(sigma2 * (1 + theta1^2))

# Covariance between X_{T+1} and X_{T+2}
rho <- sigma2 * theta1
Sigma <- matrix(c(sigma_g1^2, rho, rho, sigma_g2^2), nrow = 2)
print("Variance-Covariance Matrix:")
print(Sigma)


# Step 3: Confidence ellipse
ell <- ellipse(Sigma, centre = Pred$pred, level = 0.95, npoints = 1000)

# Step 5: Plotting
plot(ell,xlab = "Forecast for T+1",
     ylab = "Forecast for T+2",
     main = "95% Confidence Ellipse for 2-Step Forecasts",type = "l")
points(x = Pred$pred[1], y = Pred$pred[2], pch = 19, col = "red", cex = 1.5)


