
# check whether offsetting is correct - is method reasonable for online operation?

library("SX.weather")

# fit.mimic doesn't actually use the observation at all. Can therefore use it as prior.
# And rewrite the function to be a bit less stupid.

####################################################################################################

# SANITY CHECK                                                                                  ####

# does this whole 'recursive' approach actually make sense, in terms of leadtimes & chronology?
# step-by-step forecast of a given date at decreasing leadtimes, instead

y <- 5; d.t <- 70; lt <- 15; d.fc <- d.t-lt+1
# Forecast for y5, d50; begin at lt 14 (forecast made on d36)

# offset the observations, instead of the forecasts
obs.os <- abind(array(NA, dim = c(5,15,7)), obs, along = 2)


# analogue search
{
    cand <- abind("y.o" = obs.os[,2:105,],
                  "y.ecmwf" = apply(ecmwf[,2:105,,lt,-1], 1:3, mean),
                  "y.ncep" = apply(ncep[,2:105,,lt,-1], 1:3, mean),
                  "y.ukmo" = apply(ukmo[,2:105,,lt,-1], 1:3, mean),
                  "c.ecmwf" = apply(ecmwf[,1:104,,lt,-1], 1:3, mean),
                  "c.ncep" = apply(ncep[,1:104,,lt,-1], 1:3, mean),
                  "c.ukmo" = apply(ukmo[,1:104,,lt,-1], 1:3, mean),
                  along = 0)
    cand <- abind(array(NA, dim = dim(cand[,,1,])), cand, along = 3)
    cand.sd <- apply(cand, 1:2, sd, na.rm = T)
    
    targ <- cand[,,d.fc, y]
    
    dist <- sweep(sqrt(sweep(cand, 1:2, targ, "-")^2), 1:2, cand.sd, "/")
    m.dist <- apply(dist, 3:4, mean)
    
    an <- which(m.dist <= sort(m.dist)[26] & m.dist > 0, arr.ind = T)
}

# create training set
{
    # ensemble mean forecasts
    em <- abind("ecmwf" = apply(ecmwf[,,,lt,-1], 1:3, mean),
                "ncep" = apply(ncep[,,,lt,-1], 1:3, mean),
                "ukmo" = apply(ukmo[,,,lt,-1], 1:3, mean), along = 0)
    
    tr.obs <- apply(obs.os, 1, "[", an)
    tr.err <- apply(sweep(apply(em, 1:2, "[", an), c(1,3), apply(obs.os, 1, "[", an), "-"), c(1,3), mean)
}

# data available on day of forecast:
#   - current observation
#   - current ensemble forecasts
#   - training data for day of forecast
#   - prior: alpha = latest observation
#            Gamma = covariance of latest observations

d.14 <- cbind("o" = obs.os[1:2, d.fc, y],
              "ec" = ecmwf[1:2, d.fc, y, lt,-1],
              "nc" = ncep[1:2, d.fc, y, lt,-1], 
              "mo" = ukmo[1:2, d.fc, y, lt,-1],
              "tr" = t(tr.err[,1:2]))

uninf <- fit.mimic(fc = list("ecmwf" = ecmwf[, d.fc, y, lt,-1],
                             "ncep" = ncep[, d.fc, y, lt,-1],
                             "ukmo" = ukmo[, d.fc, y, lt,-1]),
                   tr = t(tr.err))

persist <- fit.mimic(fc = list("ecmwf" = ecmwf[, d.fc, y, lt,-1],
                            "ncep" = ncep[, d.fc, y, lt,-1],
                            "ukmo" = ukmo[, d.fc, y, lt,-1]),
                  tr = t(tr.err),
                  prior = list(alpha = obs.os[, d.fc, 1],
                               Gamma = cov(t(obs.os[,d.fc -(24:0),y]))))

analogue <- fit.mimic(fc = list("ecmwf" = ecmwf[, d.fc, y, lt,-1],
                               "ncep" = ncep[, d.fc, y, lt,-1],
                               "ukmo" = ukmo[, d.fc, y, lt,-1]),
                     tr = t(tr.err),
                     prior = list(alpha = apply(tr.obs, 2, mean),
                                  Gamma = cov(tr.obs)))

ens.cols <- c("steelblue", "coral3", "green3")
ens.cols2 <- c(rep(ens.cols[1], 50), rep(ens.cols[2], 20), rep(ens.cols[3], 23))

plot(t(d.14[,2:94]), pch = 1, 
     col = sapply(ens.cols2, adjustcolor, alpha = 0.5),
     main = paste0("Raw forecasts at leadtime ", lt-1))
points(em[, 1:2, d.fc, y], pch = 20, col = sapply(ens.cols, adjustcolor, alpha = 0.7))
points(t(obs.os[1:2, d.t, y]), pch = 4)
points(t(uninf[1:2,"Tau"]), col = "black")
points(t(persist[1:2,"Tau"]), col = "red")
points(t(analogue[1:2,"Tau"]), col = "blue")


