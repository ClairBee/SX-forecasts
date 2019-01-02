
# comparison of modelling methods over SRFT data

library("SX.weather"); library("ensembleBMA"); library("ensembleMOS")
setwd("~/Documents/PhD/Miniprojects/03-Sichun-paper/Code")

# 8 different ensembles
# 52 unique dates
# 969 stations
# 48-hr lag (so, this is a 2-day leadtime forecast)

# select only complete cases: 130 stations
data(srft)
cc <- ddply(srft, .(station), summarise, n = length(date))
dat <- srft[srft$station %in% cc$station[cc$n == 52],]

dat[,1:9] <- dat[,1:9] - 273.15

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Fit Bayesian model                                                                            ####

temps <- array(dim = c(9, 52, 130),
               dimnames = list("model" = colnames(dat)[1:9],
                               "date" = sapply(unique(dat$date), substr,1,8),
                               "station" = unique(dat$station)))

invisible(sapply(levels(dat$date), function(dt) {
    invisible(sapply(colnames(dat)[1:9], function(mod) {
        invisible(sapply(dimnames(temps)$station, function(stn) {
            temps[mod, substr(dt,1,8), stn] <<- dat[dat$date == dt & dat$station == stn, mod]
        }))
    }))
}))

# data provides Y_bar. n_i and C_i are missing and need to be estimated somehow.
# Sigma can be calculated.
# estimate eta and lambda from historic data.

# fit model only over first 5 stations in list for speed (also for covariance estimation...)
stn <- 1:5

bayes.fit <- abind(sapply(22:52, function(dt) {
    
    tr <- sweep(temps[-9,,], 2:3, temps[9,,], "-")[, dt-(21:2), stn]
    eta <- apply(tr, 3, mean)
    Lambda <- cov(apply(tr, 2:3, mean))
    
    Y.bar <- temps[, dt, stn]
    C.i <- aaply(temps[,,stn], 1, cov) * 0      # dummy C in absence of individual forecasts
    Sigma <- cov(Y.bar)

    D.i <- aaply(sweep(C.i, 2:3, Sigma, "+"), 1, solve)
        
    S <- Lambda + solve(apply(D.i, 2:3, sum))
    tau <- (S %*% ((solve(diag(ncol(Y.bar)) + (apply(D.i, 2:3, sum) %*% 
                                         Lambda))) %*% apply(aaply(abind(Y.bar, D.i, along = 2), 
                                                                   1, function(arr) arr[-1, ] %*% arr[1, ]), 2, sum))) - 
        eta
    
    abind(tau = tau, S, along = 2)
}, simplify = F), along = 0)

dimnames(bayes.fit)[[1]] <- dimnames(temps)$date[22:52]

# predictions
bayes.out <- bayes.fit[,,"tau"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# Fit BMA & MOS models                                                                          ####

# extract data, truncate to only 5 station locations
temp.data <- ensembleData(forecasts = dat[, 1:8],
                          dates = levels(dat$date)[dat$date], 
                          observations = dat$obs,
                          latitude = dat$lat,
                          longitude = dat$lon,
                          station = levels(dat$station)[dat$station],
                          forecastHour = 48,
                          initializationTime = "00")
temp.data <- temp.data[temp.data$station %in% dimnames(temps)$station[1:5],]

bma.fit <- ensembleBMA(temp.data, model = "normal", trainingDays = 20)
mos.fit <- ensembleMOS(temp.data, model = "normal", trainingDays = 20)

bma.out <- quantileForecast(bma.fit, temp.data)
mos.out <- array(quantileForecast(mos.fit, temp.data), dim = c(5,31))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Visual comparison of bias-correction terms                                                    ####

# eta & Lambda (for bias-correction comparison)
eta.fit <- abind(sapply(22:52, function(dt) {
    
    tr <- sweep(temps[-9,,], 2:3, temps[9,,], "-")[, dt-(21:2), stn]
    apply(tr, c(1,3), mean)
}, simplify = F), along = 0)

Lambda.fit <- abind(sapply(22:52, function(dt) {
    
    tr <- sweep(temps[-9,,], 2:3, temps[9,,], "-")[, dt-(21:2), stn]
    diag(cov(t(apply(tr, 1:2, mean))))
    
}, simplify = F), along = 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# fitted errors by station & by model
pdf("../Plots/SRFT-error-consistency.pdf", height = 4); {
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(eta.fit)[[3]], function(varb) {
        matplot(eta.fit[,,varb], type = "l", lty = 1, xlab = "", ylab = "", main = varb, ylim = c(-3,4))
    }))
    mtext("Training errors by station", outer = T)
    
    par(mfrow = c(3,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(eta.fit)[[2]], function(varb) {
        matplot(eta.fit[,varb,], type = "l", lty = 1, xlab = "", ylab = "", main = varb, ylim = c(-3,4))
    }))
    mtext("Training errors by ensemble", outer = T)
}; dev.off()
# per-variable errors are more consistent: all ensembles make similar mistakes to one another
# (although not consistent across this  training set)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# bias offset coefficients (currently not aligned due to possible non-contiguous dates?)
{
    matplot(t(bma.fit$biasCoefs[1,,]), type = "l", lty = 1, ylab = "", main = "Bias offset coefficients",
            ylim = range(bma.fit$biasCoefs[1,,], mos.fit$a, 0))
    matplot(t(mos.fit$a), type = "l", add = T, lwd = 2)
    matplot(-apply(eta.fit, 1, mean), type = "l", add = T, lwd = 3, col = "steelblue")
    
    # are there holes in the date record?
    
    legend("bottomright", bty = "n", lwd = 2, col = c("black", "steelblue"),
           c("MOS overall offset", expression(paste(eta, " offset"))))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# bias slope coefficients
matplot(t(bma.fit$biasCoefs[2,,]), type = "l", ylab = "", main = "BMA slope coefficients")

matplot(Lambda.fit, type = "l")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# Compare model output


matplot(t(bma.fit$biasCoefs[2,,]), type = "l", ylab = "", main = "BMA slope coefficients")

matplot(t(bma.fit$weights), type = "l", ylab = "", main = "BMA weights")
matplot(t(mos.fit$B[,]), type = "l", ylab = "", main = "MOS weights")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Verification                                                                                  ####

# plots of predicted vs observed temperatures
pred <- abind("bayes" = bayes.out,
              "bma" = t(array(bma.out, dim = c(5,31))),
              "mos" = t(mos.out),
              "obs" = temps[9,22:52,stn], along = 0)
err <- sweep(pred[1:3,,], 2:3, pred[4,,], "-")
rmse <- sqrt(apply(err^2, c(1,3), mean))


pdf("../Plots/SRFT-models.pdf", height = 4); {
    
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(pred)[[3]], function(varb) {
        matplot(t(pred[,,varb]), type = "l", lty = c(1,2,2,1), lwd = c(1,1,1,2), 
                col = c("red3", "steelblue", "green3", "black"), xlab = "", ylab = "Temp", main = varb)
    })) 
    mtext("Predicted & observed temps for each model", outer = T)
    plot.new()
    legend("center", bty = "n", lty = c(1,2,2,1), lwd = c(1,1,1,2), 
           col = c("red3", "steelblue", "green3", "black"), dimnames(pred)[[1]])
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    invisible(sapply(dimnames(pred)[[3]], function(varb) {
        matplot(t(err[,,varb]), type = "l", lty = c(1,2,2,1), lwd = c(1,1,1,2), 
                col = c("red3", "steelblue", "green3", "black"), xlab = "", ylab = "Temp", main = varb)
        abline(h = 0, col = "orange")
    })) 
    plot.new()
    legend("center", bty = "n", lty = c(1,2,2), lwd = c(1,1,1), 
           col = c("red3", "steelblue", "green3"), dimnames(pred)[[1]][1:3])
    mtext("Prediction errors for each model", outer = T)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    par(mfrow = c(1,1))
    matplot(t(rmse), ylim = range(0, rmse), pch = c(20,1,1), col = c("red3", "steelblue", "green3"),
            main = "RMSE for each model")
    legend("topleft", bty = "n", pch = c(20,1,1), 
           col = c("red3", "steelblue", "green3"), dimnames(pred)[[1]][1:3])
}; dev.off()


# CRPS function - has mysteriously stopped working
{
    bayes.crps <- aaply(abind("o" = temps[9,22:52,stn], bayes.fit, along = 3), 1,
                        function(arr) es.crps(o = arr[,"o"], mu = arr[,"tau"], sig = arr[,3:7]) )
    
    bma.crps2 <- array(crps(bma.fit, temp.data)[,"BMA"], dim = c(31,5))
    mos.crps <- crps(mos.fit, temp.data)
    
    bma.crps1 <- crps(bma.fit, temp.data, dates = "2004012300")
}


# verification rank histograms
bma.pit <- array(pit(bma.fit, temp.data), dim = c(5, 31))
#mos.pit <- pit(mos.fit, temp.data)

bayes.pit <- aaply(abind("o" = temps[9,22:52,stn], bayes.fit, along = 3), 1, function(arr) {
    pnorm(arr[,"o"], mean = arr[,"tau"], sd = sqrt(diag(arr[,3:7])))
})

pdf("../Plots/SRFT-PIT-hist.pdf", height = 4); {
    hist(bma.pit, breaks = c(0:20)/20, prob = T, border = "steelblue", col = adjustcolor("steelblue", alpha = 0.3),
         ylim = c(0,3), main = "PIT for BMA (blue) & Bayes (red)", xlab = "", ylab = "")
    hist(bayes.pit, breaks = c(0:20)/20, prob = T, add = T, border = "red3", col = adjustcolor("red3", alpha = 0.3))
    abline(h = 1, col = "red")
}; dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Consider a single day: first available date, 2004-01-24                                       ####

# look at biavariate case first (easier to visualise)

plot(t(dat[dat$date == "2004012400" & dat$station %in% c("KMYL ", "KMWH "), 1:8]), pch = 20, 
     col = "green3", main = "Ensemble forecasts", xlab = "KMYL", ylab = "KMWH")
points(t(dat[dat$date == "2004012400" & dat$station %in% c("KMYL ", "KMWH "), "observation"]), pch = 4, lwd = 2)

ncol(dat[dat$date %in% (which(levels(dat$date) == "2004012400") - (22:3)) & dat$station %in% c("KMYL ", "KMWH "),]
plot()