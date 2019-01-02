
# replicate output of BMA & MOS models

library("SX.weather"); library("ensembleBMA"); library("ensembleMOS")
setwd("~/Documents/PhD/Miniprojects/03-Sichun-paper/Code")

lt <- 5; d <- 32; y <- 5 
timestamps <- gsub("-", "", as.Date(load.data("../Data/ECMWF_europe_starttime.rda")[16:105,1:7]/ 24, "1800-01-01"))
dt <- timestamps[90*(y-1) + d]

exch <- c(rep(1, 50), rep(2, 20), rep(3, 23))
names(exch) <- c(paste0("ec", 1:50), paste0("nc", 1:20), paste0("mo", 1:23))

par(pch = 20, mar = c(2,2,3,1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Create necessary data sets                                                                    ####

fc <- cbind(apply(apply(offset.forecast(ecmwf)[1:2,,,lt,-1], c(1,4), c), 3, c),
            apply(apply(offset.forecast(ncep)[1:2,,,lt,-1], c(1,4), c), 3, c),
            apply(apply(offset.forecast(ukmo)[1:2,,,lt,-1], c(1,4), c), 3, c))
colnames(fc) <- names(exch)

o <- c(apply(obs[1:2,,], 1, c))

temp.dat <- ensembleData(forecasts = fc,
                         dates = rep(timestamps, 2),
                         observations = o,
                         forecastHour = 0,
                         station = rep(c("N", "S"), each = 630),
                         latitude = rep(c(56, 52), each = 630),
                         longitude = rep(c(-3, 0), each = 630),
                         initializationTime = "00",
                         exchangeable = exch)

temp.n <- temp.dat[temp.dat$station == "N",]
temp.s <- temp.dat[temp.dat$station == "S",]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# extract training data for the target date
tr <- trainingData(temp.dat, trainingDays = 25, date = dt)
tr.n <- trainingData(temp.n, trainingDays = 25, date = dt)
tr.s <- trainingData(temp.s, trainingDays = 25, date = dt)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# datasets using ensemble mean values only
tm.dat <- ensembleData(forecasts = cbind("ec" = apply(fc[,1:50], 1, mean),
                                         "nc" = apply(fc[,51:70], 1, mean),
                                         "mo" = apply(fc[,71:93], 1, mean)),
                       dates = rep(timestamps, 2),
                       observations = o,
                       forecastHour = 0,
                       station = rep(c("N", "S"), each = 630),
                       latitude = rep(c(56, 52), each = 630),
                       longitude = rep(c(-3, 0), each = 630),
                       initializationTime = "00")

tm.n <- tm.dat[tm.dat$station == "N",]
tm.s <- tm.dat[tm.dat$station == "S",]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Fit models using packages                                                                     ####

# models fitted using individual ensemble members
bma.fit <- ensembleBMA(temp.dat, model = "normal", trainingDays = 25, dates = dt)
bma.n <- ensembleBMA(temp.n, model = "normal", trainingDays = 25, dates = dt)
bma.s <- ensembleBMA(temp.s, model = "normal", trainingDays = 25, dates = dt)

mos.fit <- ensembleMOS(temp.dat, model = "normal", trainingDays = 25, dates = dt)
mos.n <- ensembleMOS(temp.n, model = "normal", trainingDays = 25, dates = dt)
mos.s <- ensembleMOS(temp.s, model = "normal", trainingDays = 25, dates = dt)


bma.mean.fit <- ensembleBMA(tm.dat, model = "normal", trainingDays = 25, dates = dt)
mos.mean.fit <- ensembleMOS(tm.dat, model = "normal", trainingDays = 25, dates = dt)

bma.mean.n <- ensembleBMA(tm.n, model = "normal", trainingDays = 25, dates = dt)
bma.mean.s <- ensembleBMA(tm.s, model = "normal", trainingDays = 25, dates = dt)

mos.mean.n <- ensembleMOS(tm.n, model = "normal", trainingDays = 25, dates = dt)
mos.mean.s <- ensembleMOS(tm.s, model = "normal", trainingDays = 25, dates = dt)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# compare to results using ensemble means only
mos.mean.fit$a; mos.fit$a
mos.mean.fit$B; ddply(data.frame(model = sapply(rownames(mos.fit$B), substr, 1, 2), b = mos.fit$B)
                      , .(model), summarise, sum(X20120101))
mos.mean.fit$c; mos.fit$c
mos.mean.fit$d; mos.fit$d

bma.mean.fit$biasCoefs; bma.fit$biasCoefs[, c(1,51,71),]
bma.mean.fit$sd; bma.fit$sd
bma.mean.fit$weights; bma.fit$weights[c(1,51,71)] * c(50, 20, 23)

# MOS fits same model regardless of inclusion of individual ensemble members. 
# BMA produces quite different fit with & without individual members.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Bias correction replicated                                                                    ####

# ensemble BMA starts with linear regression over each ensemble's training data
# each ensemble member, time & location is treated as drawn from same PDF (one per ensemble source)
bma.coefs <- lapply(list("ec" = cbind(c(obs[1:2,d-(24:0),y]), 1, c(offset.forecast(ecmwf)[1:2,d-(24:0),y,lt,-1])),
                         "nc" = cbind(c(obs[1:2,d-(24:0),y]), 1, c(offset.forecast(ncep)[1:2,d-(24:0),y,lt,-1])),
                         "mo" = cbind(c(obs[1:2,d-(24:0),y]), 1, c(offset.forecast(ukmo)[1:2,d-(24:0),y,lt,-1]))),
                    function(dat) {
                        Y.n <- as.matrix(dat[,1])
                        X.n <- as.matrix(dat[,2:3])
                        
                        solve(t(X.n) %*% X.n) %*% t(X.n) %*% Y.n
                    })

# ensemble MOS does essentially the same, but over means & in single model

    Y.n <- as.matrix(c(obs[1:2,d-(24:0),y]))
    X.n <- apply(abind("C" = array(1, dim = c(2,25)),
                       "ec" = apply(offset.forecast(ecmwf)[1:2,d-(24:0),y,lt,-1], 1:2, mean),
                       "nc" = apply(offset.forecast(ncep)[1:2,d-(24:0),y,lt,-1], 1:2, mean),
                       "mo" = apply(offset.forecast(ukmo)[1:2,d-(24:0),y,lt,-1], 1:2, mean), along = 3),
                 3, c)
    solve(t(X.n) %*% X.n) %*% t(X.n) %*% Y.n
}
