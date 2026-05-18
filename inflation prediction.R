rm(list = ls())

install.packages("readxl")
install.packages("dplyr")
install.packages("openxlsx")
install.packages("forecast")
install.packages("zoo")
install.packages("vars")
install.packages("urca")


# Load the library
library(readxl)
library(dplyr)
library(openxlsx)
library(forecast)

library(zoo)
library(vars)
library(urca)

mydata <- read.csv("C:/Users/miles/OneDrive/Documents/Data projects/inflation_data.csv", header=TRUE, sep= "," )


############NEW Method 
t <- length(date)
n <- 4


IF <- mydata$inflation
UN <- mydata$unemployment
CA <- mydata$cash.rate
FU <- mydata$fuel

# Remove commas and convert blank strings to NA
FU <- gsub(",", "", FU)
FU[FU == ""] <- NA

# Convert to numeric
FU <- as.numeric(FU)

str(FU)


x <- cbind(IF, UN, CA, FU)

x <- na.omit(x)


vecm_ts <- ts(x, start = c(1993, 4), frequency = 4)


lag_select <- VARselect(vecm_ts, lag.max = 20, type = "const")
lag_select$selection

K_lags <- 10

length(K_lags)

i=0
for (i in 1:length(K_lags)){
johansen_test <- ca.jo(
  vecm_ts,
  type = "trace",
  ecdet = "const",
  K = K_lags[i],
  spec = "transitory"
)

print(summary(johansen_test))

}

rank <- 1
vecm_var <- vec2var(johansen_test, r = rank)
vecm_forecast <- predict(vecm_var, n.ahead = 12, ci = 0.95)

# Pull CPI forecast results
cpi_fc <- vecm_forecast$fcst$IF

forecast_mean <- cpi_fc[, "fcst"]
forecast_lower <- cpi_fc[, "lower"]
forecast_upper <- cpi_fc[, "upper"]

# Historical CPI
actual_cpi <- vecm_ts[, "IF"]

# Create x-axis
actual_time <- time(actual_cpi)
forecast_time <- seq(
  from = max(actual_time) + 1/frequency(vecm_ts),
  by = 1/frequency(vecm_ts),
  length.out = length(forecast_mean)
)

# Plot actual CPI
plot(
  actual_cpi,
  type = "l",
  col = "black",
  lwd = 1.5,
  xlim = c(2010, max(forecast_time)),
  ylim = range(c(actual_cpi, forecast_lower, forecast_upper), na.rm = TRUE),
  main = "Forecast for CPI - VECM, Autolags: 10",
  xlab = "Years",
  ylab = "CPI"
)

# Vertical forecast line
abline(v = max(actual_time), lty = 2, col = "grey")

# Add forecast mean
lines(forecast_time, forecast_mean, col = "blue", lwd = 1.5, lty = 2)

# Add confidence intervals
lines(forecast_time, forecast_upper, col = "red", lwd = 1.2, lty = 3)
lines(forecast_time, forecast_lower, col = "red", lwd = 1.2, lty = 3)


##################################################################
# Creating an impulse response funtion 
##################################################################

fuel_ts <- ts(FU ,start = c(1993, 3), frequency = 4)

plot(
  fuel_ts,
  type = "l",
  xlim= c(1990,2025),
  main = "Fuel Prices Over Time",
  xlab = "Time",
  ylab = "Fuel Price Index"
)



x <- cbind(IF, UN, CA, FU)
x <- na.omit(x)

x_diff <- diff(x)


var_model <- VAR(x_diff, p = 2, type = "const")

irf_fuel <- irf(
  var_model,
  impulse = "FU",
  response = "IF",
  n.ahead = 12,
  boot = TRUE,
  ci = 0.95
)

plot(irf_fuel)


log_FU <- log(FU)

######################## creating a percentage shock in fuel 

d_log_FU <- diff(log_FU)
d_IF <- diff(IF)
d_UN <- diff(UN)
d_CA <- diff(CA)

dx <- cbind(d_log_FU, d_IF, d_UN, d_CA)
dx <- na.omit(dx)
var_model <- VAR(dx, p = 2, type = "const" )
log(1996 / 1838)
log(780/645)

shock_size <- 0.1900436/ sd(d_log_FU, na.rm = TRUE)


irf_fuel <- irf(
  var_model,
  impulse = "d_log_FU",
  response = "d_IF",
  n.ahead = 12,
  boot = TRUE,
  ci = 0.95
)

irf_fuel$irf$d_log_FU <- irf_fuel$irf$d_log_FU * shock_size
irf_fuel$Lower$d_log_FU <- irf_fuel$Lower$d_log_FU * shock_size
irf_fuel$Upper$d_log_FU <- irf_fuel$Upper$d_log_FU * shock_size

plot(irf_fuel,
     main = "19.0% Shock IRF",
     xlab = "Number of Quaters",
     ylab = "Percentage Point Change")

#Data starts from 1993  Q4 when the Rba annouced a 2-3% target on inflation 

# Generate ACF and PACF to see if we can identify trends to see if the lags in the data are signifcant 
i=1
for (i in 1:10)
{ 
  acf(mydata[1 + i], main= colnames(mydata[1+i]))
  pacf(mydata[1 + i], main= colnames(mydata[1+i]))
}

# create variables 
IF <- mydata$inflation[-t]
UN <- mydata$unemployment[-t]
CA <- mydata$cash.rate[-t]

t <- length(date)
n <- 3

x <- cbind(IF, UN, CA)

# 
VAR_est <- list()
ic_var <- matrix(nrow = 20, ncol = 3)
colnames(ic_var) <- c("p", "aic", "bic")
for (p in 1:20)
{
  VAR_est[[p]] <- VAR(x, p)
  ic_var[p,] <- c(p, AIC(VAR_est[[p]]),
                  BIC(VAR_est[[p]]))
}
# ordering best AIC and BIC models 
ic_aic_var <- ic_var[order(ic_var[,2]),]
ic_bic_var <- ic_var[order(ic_var[,3]),]


# we will use the LM test by setting type = "BG", but other tests are just
# as valid.

adq_set_var <- as.matrix(ic_var[c(2, 9, 10), ])
adq_idx_var <- c(2, 9, 10)
nmods <- length(adq_idx_var)
for (i in 1:nmods)
{
  p <- adq_idx_var[i]
  print(paste0("Checking VAR(", p, ")"))
  print(serial.test(VAR_est[[p]], lags.bg = 1,
                    type = "BG"))
}


# Test Eigenvalue 
nmods <- length(adq_idx_var)
for (i in 1:nmods)
{
  p <- adq_idx_var[i]
  print(paste0("VAR(", p,
               "): Maximum absolute eigenvalue is ",
               max(vars::roots(VAR_est[[p]]))))
}



# VECM models 
MO <- c(2,9,10)
VECM_est <- list()
ic_vecm <- matrix(nrow = 4* (1 + n), ncol = 4)
colnames(ic_vecm) <- c("p", "r", "aic", "bic")
i <- 0
for (p in MO)
{
  for (r in 0:n)
  {
    i <- i + 1
    if (r == n)
    {
      VECM_est[[i]] <- VAR(x, p)
    }
    else if (r == 0)
    {
      VECM_est[[i]] <- VAR(diff(x), p - 1)
    }
    else
    {
      VECM_est[[i]] <- vec2var(ca.jo(x, K = p), r)
    }
    ic_vecm[i,] <- c(p, r, AIC(VECM_est[[i]]),
                     BIC(VECM_est[[i]]))
  }
}
ic_aic_vecm <- ic_vecm[order(ic_vecm[,3]),][1:7,]
ic_bic_vecm <- ic_vecm[order(ic_vecm[,4]),][1:7,]
ic_int_vecm <- intersect(as.data.frame(ic_aic_vecm),
                         as.data.frame(ic_bic_vecm))


adq_set_vecm <- as.matrix(arrange(as.data.frame(
  ic_int_vecm), p, r))
adq_idx_vecm <- match(data.frame(t(adq_set_vecm[, 1:2])),
                      data.frame(t(ic_vecm[, 1:2])))



nmods <- length(adq_idx_vecm)
for (i in 1:nmods)
{
  p <- adq_set_vecm[i, 1]
  r <- adq_set_vecm[i, 2]
  print(paste0("Checking VECM(", p, ", ", r, ")"))
  print(serial.test(VECM_est[[adq_idx_vecm[i]]],
                    lags.bg = 1,
                    type = "BG"))
}



hrz = 12
Vecm_fcst <- list()
xlim <- c(length(t) - 3 * hrz,
          length(t) + hrz)
ylim <- c(mydata$inflation[xlim[1]],
          max(mydata$inflation) + 0.2)
ylim <- range(mydata$inflation, na.rm = TRUE)
ylim <- c(ylim[1], ylim[2] + 0.2)
for (i in 1:nmods)
{
  p <- adq_idx_vecm[i]
  Vecm_fcst[[i]] <- predict(VECM_est[[p]],
                           n.ahead = hrz)
  plot(Vecm_fcst[[i]], names = "lrgdp",
       xlim = xlim, ylim = ylim,
       main = "Forecast for Log Real GDP",
       xlab = "Horizon",
       ylab = "RRP")
}
