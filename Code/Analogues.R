
# BEFORE CONTINUING...
    # CODE ENERGY SCORE (USING MC APPROXIMATION)
    # CODE BOX ORDINAL TRANSFORM?
    # CODE MULTIVARIATE PIT
    # CODE MST
    # CODE DETERMINANT SHARPNESS
    # CODE DEVIATION FROM UNIFORMITY
    
    # ... AND APPLY ALL OF THE ABOVE TO THE BASE MODELS ALREADY FITTED
    # (maybe print some BW plots to avoid having to carry laptop?)

# THEN...
    # look at distribution of verifying observations vs target for closest analogues
    # consider sequential filtering: threshold at maximum distance over prior n days?

library("SX.weather"); library("CB.Misc")

require(fpc)            # Bhattcharyya distance 
require(distrEx)        # Hellinger distance, total variation etc

# IDENTIFYING ANALOGUES TO THE CURRENT FORECAST

# get analogues (= training set) per ensemble mean
# calculate ens. mean error for each training instance
# Eta & Lambda are the mean & covariance matrix over this training set.

#--------------------------------------------------------------------------------------------------

# ideally, search for analogues that are similar in terms of current forecast & preceding weather
# (may have to reserve preceding weather for prior on alpha/Gamma?)

# look for analogues to each model's forecast, then average, or find analogue to average?
# try individual analogues first (mainly to understand similarities & code)

# how close were the verifying observations on those days to the current verifying obs?

####################################################################################################

# MAHALANOBIS DISTANCE FROM ALL ENS. MEANS TO CURRENT FORECAST DIST                             ####

fc.all <- abind("ecmwf" = apply(offset.forecast(ecmwf)[1:2,,,,-1], 1:4, mean),
                "ncep" = apply(offset.forecast(ncep)[1:2,,,,-1], 1:4, mean),
                "ukmo" = apply(offset.forecast(ukmo)[1:2,,,,-1], 1:4, mean),
                along = 0)

# consider analogues to day 5, year 1 only - longest leadtime
# Mahalanobis distance to single nearest forecast, per model

# model current ensemble forecasts as distribution, find distance to all previous points
{
    Y.bar <- abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
                   "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
                   "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
                   along = 0)[,1:2,,,]
    
    C <- abind("ecmwf" = square.mat(apply(aperm(offset.forecast(ecmwf)[,,,,-1], c(5,1:4)), 3:5, cov)),
               "ncep" = square.mat(apply(aperm(offset.forecast(ncep)[,,,,-1], c(5,1:4)), 3:5, cov)),
               "ukmo" = square.mat(apply(aperm(offset.forecast(ukmo)[,,,,-1], c(5,1:4)), 3:5, cov)), 
               along = 0)[,1:2,1:2,,,]
}

md <- abind("ecmwf" = mahalanobis(x = apply(apply(offset.forecast(ecmwf)[1:2,,,,-1], 1:4, mean), 1, rbind),
                                  center = apply(offset.forecast(ecmwf)[1:2,5,1,14,-1], 1, mean),
                                  cov = cov(t(offset.forecast(ecmwf)[1:2,5,1,14,-1]))),
            "ncep" = mahalanobis(x = apply(apply(offset.forecast(ncep)[1:2,,,,-1], 1:4, mean), 1, rbind),
                                 center = apply(offset.forecast(ncep)[1:2,5,1,14,-1], 1, mean),
                                 cov = cov(t(offset.forecast(ncep)[1:2,5,1,14,-1]))),
            "ukmo" = mahalanobis(x = apply(apply(offset.forecast(ukmo)[1:2,,,,-1], 1:4, mean), 1, rbind),
                                   center = apply(offset.forecast(ukmo)[1:2,5,1,14,-1], 1, mean),
                                   cov = cov(t(offset.forecast(ukmo)[1:2,5,1,14,-1]))),
            along = 0)

plot(apply(apply(offset.forecast(ecmwf)[1:2,,,,-1], 1:4, mean), 1, rbind)[md["ecmwf",] < 1,],
     col = adjustcolor("steelblue4", alpha = 0.5), pch = 20, 
     main = "Mahalanobis distance of forecasts")
points(apply(apply(offset.forecast(ncep)[1:2,,,,-1], 1:4, mean), 1, rbind)[md["ncep",] < 1,],
       col = adjustcolor("darkolivegreen4", alpha = 0.3), pch = 20)
points(apply(apply(offset.forecast(ukmo)[1:2,,,,-1], 1:4, mean), 1, rbind)[md["ukmo",] < 1,],
       col = adjustcolor("goldenrod3", alpha = 0.1), pch = 20)
points(t(apply(offset.forecast(ecmwf)[1:2,5,1,14,-1], 1, mean)), col = "darkred", pch = 20)

sum(apply(md, 2, max) < 1)

plot(apply(apply(offset.forecast(ecmwf)[1:2,,,,-1], 1:4, mean), 1, rbind)[apply(md, 2, max) < 1,],
     col = adjustcolor("steelblue4", alpha = 0.5), pch = 20, 
     main = "Mahalanobis distance of forecasts - 2")
points(apply(apply(offset.forecast(ncep)[1:2,,,,-1], 1:4, mean), 1, rbind)[apply(md, 2, max) < 1,],
       col = adjustcolor("indianred2", alpha = 0.3), pch = 20)
points(apply(apply(offset.forecast(ukmo)[1:2,,,,-1], 1:4, mean), 1, rbind)[apply(md, 2, max) < 1,],
       col = adjustcolor("goldenrod3", alpha = 0.2), pch = 20)
points(t(apply(offset.forecast(ecmwf)[1:2,5,1,14,-1], 1, mean)), col = "darkred", pch = 20)

####################################################################################################

# HELLINGER / BHATTCHARYYA DISTANCES                                                            ####

require(fpc)

models <- abind("Y.bar" = apply(abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
                                      "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
                                      "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
                                      along = 0)[,1:2,,,], c(1:2, 5), c),
                "C" = apply(abind("ecmwf" = square.mat(apply(aperm(offset.forecast(ecmwf)[,,,,-1], c(5,1:4)), 3:5, cov)),
                                  "ncep" = square.mat(apply(aperm(offset.forecast(ncep)[,,,,-1], c(5,1:4)), 3:5, cov)),
                                  "ukmo" = square.mat(apply(aperm(offset.forecast(ukmo)[,,,,-1], c(5,1:4)), 3:5, cov)), 
                                  along = 0)[,1:2,1:2,,,], c(1:3, 6), c),
                along = 3)

current.fc <- list("mu" = Y.bar["ecmwf", ,1,1,14],
                   "sigma" = C["ecmwf",,,1,1,14])

# around 1s to compare 630 potential analogues at all 15 leadtimes
bd <- apply(models[,"ecmwf",,,], c(1,4), function(arr) {
        bhattacharyya.dist(mu1 = current.fc$mu,
                           mu2 = arr[1,],
                           Sigma1 = current.fc$sigma,
                           Sigma2 = arr[-1,])
})

# find 50 smallest distances & plot verifying observations
