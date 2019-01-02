
library("SX.weather"); library(VineCopula); library(copula); library(sn)

dat <- apply(offset.forecast(ecmwf)[1:2,,,15,-1], 1, c)

plot(dat, pch = 20, col = adjustcolor("violetred4", alpha = 0.3))
hist(dat[,1], breaks = "fd", col = adjustcolor("violetred4", alpha = 0.3), border = "violetred4", prob = T)
lines((-100:150)/10, dnorm((-100:150)/10, mean(dat[,1]), sd(dat[,1])), col = "red3", lwd = 2)

hist(dat[,2], breaks = "fd", col = adjustcolor("violetred4", alpha = 0.3), border = "violetred4", prob = T)
lines((-100:150)/10, dnorm((-100:150)/10, mean(dat[,2]), sd(dat[,2])), col = "red3", lwd = 2)
lines((-100:150)/10, dsn((-100:150)/10, xi = mean(dat[,2]), omega = sd(dat[,2]), alpha = -skewness(dat[,2])), col = "blue", lwd = 2)
lines((-100:150)/10, dcauchy((-100:150)/10, mean(dat[,2]), sd(dat[,2])), col = "blue", lwd = 2)

cor(dat, method = "kendall")    # 0.54 for year 1, 0.60 for all years

# fit copula over 1 year, 1 leadtime only
temp.copula <- BiCopSelect(pobs(dat)[,1], pobs(dat)[,2], familyset = NA)  

# over year one:  t, rho = 0.74, df 16
# over all years: t, rho = 0.79, df 14

plot((-50:50)/100, dt((-50:50)/100, 13), type = "l")
lines((-50:50)/100, dt((-50:50)/100, 16), col = "red")

mvn <- rCopula(2000, normalCopula(param = 0.8, dim = 2))
st <- rCopula(2000, tCopula(param = 0.8, dim = 2, df = 14))
plot(mvn, pch = 20, col = adjustcolor("steelblue", alpha = 0.4), xlab = "u", ylab = "v", 
     main = expression(paste("Normal copula, ", rho,  " = 0.8")))
plot(st, pch = 20, col = adjustcolor("steelblue", alpha = 0.4), xlab = "u", ylab = "v", 
     main = expression(paste("Student-t copula, ", rho,  " = 0.8, df = 14")))

