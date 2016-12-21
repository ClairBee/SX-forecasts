

# regarding the problem of calculating the determinant over increasing numbers of variables...
# (which is required to calculate the Bhattacharyya distance to identify analogues)
# IF we're looking for the closest match to all three ensemble means, then covariance may not be vital.

# each ensemble consists of a set of random perturbations, propagated over a certain leadtime.
# Therefore if we limit search to same lead time (justified for physical reasons),
# each ensemble should have a reasonably similar covariance matrix (up to a certain point).
# Therefore using the covariance assumption doesn't necessarily add anything.

# Basically: if we use the positions of the ensemble means + the observation as point estimates,
# is that more informative than the spreads of the individual perturbations?

# may allow more efficient scaling up to a spatial field?

####################################################################################################

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
boxplot(bd)

# plot determinants?
dt <- apply(models[,,-1,,], c(1:2,5), det)

boxplot(dt[,"ecmwf",], border = "green3", col = adjustcolor("gold", alpha = 0.2))
boxplot(dt[,"ncep",], add = T, border = "red", col = adjustcolor("coral", alpha = 0.2))
boxplot(dt[,"ukmo",], add = T, border = "blue", col = adjustcolor("skyblue", alpha = 0.2))

####################################################################################################

# SIMPLE METRIC - MAX DISTANCE OF Q OF INTEREST                                                 ####

# search for analogues over 8pts: obs + 3 forecast means for two variables
# use Euclidean distance for now (all same variable, so no normalisation needed)

target <- rbind("o" = obs[1:2,1,1], 
                "ecmwf" = apply(offset.forecast(ecmwf)[1:2,1,1,1,-1], 1, mean),
                "ncep" = apply(offset.forecast(ncep)[1:2,1,1,1,-1], 1, mean),
                "ukmo" = apply(offset.forecast(ukmo)[1:2,1,1,1,-1], 1, mean))
search <- abind("o" = obs[1:2,,], 
                "ecmwf" = apply(offset.forecast(ecmwf)[1:2,,,1,-1], 1:3, mean),
                "ncep" = apply(offset.forecast(ncep)[1:2,,,1,-1], 1:3, mean),
                "ukmo" = apply(offset.forecast(ukmo)[1:2,,,1,-1], 1:3, mean),
                along = 0)

# get Euclidean norm of each obs/ensemble
t.dist <- apply(sweep(search, 1:2, target, "-"), c(1,3,4), function(v) sqrt(sum(v^2)))
t.max <- apply(t.dist, 2:3, max)
t.mean <- apply(t.dist, 2:3, mean)

analogues.max <- which(t.max <= sort(t.max)[26], arr.ind = T)[2:26,]
analogues.mean <- which(t.mean <= sort(t.mean)[26], arr.ind = T)[2:26,]
# 8 differences.

# ideally remove all dates after target in the current year.

plot(target, pch = 20, col = c("black", "red", "green3", "blue"), xlim = c(2,7), ylim = c(5.5,10),
     main = "Analogues selected by max distance from target")
points(search["ecmwf", 1, , ][analogues.max], search["ecmwf", 2, , ][analogues.max], col = adjustcolor("red", alpha = 0.4), pch = 4)
points(search["ncep", 1, , ][analogues.max], search["ncep", 2, , ][analogues.max], col = adjustcolor("green3", alpha = 0.4), pch = 4)
points(search["ukmo", 1, , ][analogues.max], search["ukmo", 2, , ][analogues.max], col = adjustcolor("blue", alpha = 0.4), pch = 4)
points(search["o", 1, , ][analogues.max], search["o", 2, , ][analogues.max], col = adjustcolor("black", alpha = 0.4), pch = 3)

plot(target, pch = 20, col = c("black", "red", "green3", "blue"), xlim = c(2,7), ylim = c(5.5,10),
     main = "Analogues selected by mean distance from target")
points(search["ecmwf", 1, , ][analogues.mean], search["ecmwf", 2, , ][analogues.mean], col = adjustcolor("red", alpha = 0.4), pch = 4)
points(search["ncep", 1, , ][analogues.mean], search["ncep", 2, , ][analogues.mean], col = adjustcolor("green3", alpha = 0.4), pch = 4)
points(search["ukmo", 1, , ][analogues.mean], search["ukmo", 2, , ][analogues.mean], col = adjustcolor("blue", alpha = 0.4), pch = 4)
points(search["o", 1, , ][analogues.mean], search["o", 2, , ][analogues.mean], col = adjustcolor("black", alpha = 0.4), pch = 3)

# now use these analogues to fit a model...
# use mean distance as initial decider
training.set <- apply(search, 1:2, "[", analogues.mean)
tr.error <- sweep(training.set[,-1,], c(1,3), training.set[,1,], "-")

Eta <- apply(tr.error, 2:3, mean)

Lam <- square.mat(apply(tr.error, 2:3, cov))

Y.bar <- abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
               "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
               "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
               along = 0)
C <- abind("ecmwf" = square.mat(apply(aperm(offset.forecast(ecmwf)[,,,,-1], c(5,1:4)), 3:5, cov)),
           "ncep" = square.mat(apply(aperm(offset.forecast(ncep)[,,,,-1], c(5,1:4)), 3:5, cov)),
           "ukmo" = square.mat(apply(aperm(offset.forecast(ukmo)[,,,,-1], c(5,1:4)), 3:5, cov)), 
           along = 0)

Sig <- square.mat(apply(abind(lapply(list(ecmwf, ncep, ukmo),
                                     function (model) {
                                         apply(offset.forecast(model)[,,,,-1], 1:4, mean)
                                     }), along = 0), 3:5, cov))

D <- sweep(sweep(C, 1, table(attr(superensemble(), "m"))[2:4], "/"), 2:6, Sig, "+")

# truncate everything
Y.bar <- Y.bar[,1:2,1,1,1]
C <- C[,1:2,1:2,1,1,1]
Sig <- Sig[1:2,1:2,1,1,1]
D <- D[,1:2,1:2,1,1,1]

# now fit model
S <- solve(apply(array.solve(D, 1), 1:2, sum)) + Lam
S <- array.solve(apply(array.solve(D, c(1, 4:5)), c(1:2, 4:5), sum), c(3:4)) + Lam

####################################################################################################

# CHANGEPOINT IN NCEP MODEL                                                                     ####

# NCEP model was upgraded to correct for GFS bias in March 2009 (end of year 3)
# this may tally with a change in the annual forecast error at this time...

err <- apply(forecast.errors(ncep)[,,,,-1], 1:4, mean)

matplot(err[1,,,1], type = "l", lty = 1, col = c(rep("black", 3), rep("red", 4)))
abline(0,0,col = "blue")

plot(c(err[1,,,1]), type = "l")
abline(0,0,col = "blue")
abline(v = c(0:7)*90, lty = 2, col = "cyan3")
abline(v = 90*3, col = "red")

library(changepoint)
{
    cp <- apply(err[,,,1], 1, function(o) cpt.mean(c(o)))
    
    cp[[1]]@cpts; cp[[2]]@cpts
    
    par(mfrow = c(2,1), mar = c(2,2,3,1))
    plot(cp[[1]], main = "temp.n", xlab = "", ylab = "")
    abline(v = 90*3, col = "orange", lwd = 2, lty = 2)
    plot(cp[[2]], main = "temp.s", xlab = "", ylab = "")
    abline(v = 90*3, col = "orange", lwd = 2, lty = 2)
    
    # December 2010 was unusually cold (down to -10);
    # this may explain the apparent extension of the warm bias.
    
    m.err <- cbind("pre" = apply(err[,,1:3,1], 1, mean),
                   "post" = apply(err[,,4:7,1], 1, mean))
    
    err.var <- cbind("pre" = apply(err[,,1:3,1], 1, function(v) var(c(v))),
                     "post" = apply(err[,,4:7,1], 1, function(v) var(c(v))))
    
    plot(cp[[3]], main = "pc1", xlab = "", ylab = "")
    lines(x = c(0,270, NA, 270, 630, NA), y = rep(m.err["pc1", ], each = 3), col = "blue")
    
    plot(cp[[4]], main = "pc2", xlab = "", ylab = "")
    lines(x = c(0,270, NA, 270, 630, NA), y = rep(m.err["pc2", ], each = 3), col = "blue")
    
    plot(cp[[5]], main = "pc3", xlab = "", ylab = "")
    lines(x = c(0,270, NA, 270, 630, NA), y = rep(m.err["pc3", ], each = 3), col = "blue")
    
    tt <- apply(err[,,,1], 1, function(dat) t.test(c(dat[,1:3]), c(dat[,4:7])))
    lapply(lapply(tt, "[[", "p.value"), round, 5)
}

# okay, we have a changepoint. Is the effect noticeable in the fitted models?

fitted <- readRDS("./Models/lambda-hist.rds")
plot(fitted$tau[1, 1,], type = "l")

lines(c(obs[1,,])[26:630], col = "red")

err <- fitted$tau[1, 1,] - c(obs[1,,])[26:630]
plot(fitted$tau[1, 1,] - c(obs[1,,])[26:630], type = "l")
abline(v = 245, col = "red")

t.test(err[1:245], err[246:605])
# no evidence of any difference in mean error before/after changepoint.

# What about effect on weights in BMA/MOS fitted models? 
m.bma <- readRDS("./Models/ensBMA-temps.rds")[[1]]
m.mos <- readRDS("./Models/ensMOS-temps.rds")[[1]]
ens.data <- readRDS("./Models/Ens-data.rds")[[1]]

# certainly seems to appear in the BMA bias coefficients, not obvious in the weights
{
    matplot(t(m.bma$biasCoefs[1,,]), type = "l", lty = c(2,1,2), col = c("skyblue", "black", "coral"),
            xlab = "", ylab = "", main = "BMA bias correction offset")
    abline(v = (0:7)*90 - 25, col = "cyan3", lty = 2)
    abline(v = 245, col = "red")
    abline(h = 0, col = "blue")
    legend("bottomright", col = c("skyblue", "black", "coral", NA, "red"), lty = c(2,1,2,1,1), bty = "n",
           legend = c("ECMWF", "NCEP", "UKMO", NA, "Changepoint"))
    
    matplot(t(m.bma$biasCoefs[2,,]), type = "l", lty = c(2,1,2), col = c("skyblue", "black", "coral"),
            xlab = "", ylab = "", main = "BMA bias correction slope")
    abline(v = (0:7)*90 - 25, col = "cyan3", lty = 2)
    abline(v = 245, col = "red")
    abline(h = 1.05, col = "blue")
    legend("bottomright", col = c("skyblue", "black", "coral", NA, "red"), lty = c(2,1,2,1,1), bty = "n",
           legend = c("ECMWF", "NCEP", "UKMO", NA, "Changepoint"))
    
    matplot(t(m.bma$weights), type = "l", lty = c(2,1,2), col = c("skyblue", "black", "coral"),
            xlab = "", ylab = "", main = "BMA model weights")
}

# also not so obvious in MOS.
{
    plot(c(m.mos$a), type = "l")
    abline(v = (0:7)*90 - 25, col = "cyan3", lty = 2)
    abline(v = 245, col = "red")
    abline(h = 0, col = "blue")
    
    matplot(t(m.mos$B), type = "l", lty = c(2,1,2), col = c("skyblue", "black", "coral"),
            xlab = "", ylab = "", main = "MOS per-model coefficients")
    abline(v = (0:7)*90 - 25, col = "cyan3", lty = 2)
    abline(v = 245, col = "red")
    abline(h = 0, col = "blue")
    
    qf <- quantileForecast(m.mos, ens.data)
    err <- cbind("temp.n" = qf[1:606,] - ens.data$observations[25:630],
                 "temp.s" = qf[(1:606)+606,] - ens.data$observations[(25:630)+30])
    plot(err[,1], type = "l")
    abline(v = 295, col = "red")
    plot(err[,2], type = "l")
    
}

# should also check for changepoints in any of the other models!

err <- apply(forecast.errors(ukmo)[,,,,-1], 1:4, mean)
cp <- apply(err[,,,1], 1, function(o) cpt.mean(c(o)))
cp[[1]]@cpts; cp[[2]]@cpts; cp[[3]]@cpts; cp[[4]]@cpts; cp[[5]]@cpts

cp[[1]]@param.est$mean; cp[[2]]@param.est$mean; cp[[3]]@param.est$mean; cp[[4]]@param.est$mean; cp[[5]]@param.est$mean

par(mfrow = c(5,1), mar = c(2,2,1,1))
lapply(cp, plot)


# ECMWF & UKO:
# no evidence at shorter leadtimes. May be a bias shift after ~300days in temp at longest leads.
# However shift seems to be halfway through y4, so probably artificial.
# also, moved around (at lt 8, temp.n is 430 and temp.s 163, suggesting this is not a true CP)
# shift appears at similar point in both models, and mean diff is v small, suggesting artefact

sum(cp[[1]]@param.est$mean * c(1,-1)); sum(cp[[2]]@param.est$mean * c(1,-1))

# NCEP bias shift is 1.75 & 1 at LT 0; 2.9 & 2.1 at LT 15.
# ECMWF is 0 & 0; 2 & 1.2
# UKMO is 0 & 0; 1.9 & 1.2

# seems that shift only occurs in NCEP.

####################################################################################################

# MAHALANOBIS DISTANCE VS EUCLIDEAN NORM                                                        ####

# consider single forecast; treat as 8d multivariate

lt <- "14"; d <- 90; y <- 7

dat <- abind("o" = obs[1:2,,],
             "ecmwf" = )

fc.mu <- cbind(obs[1:2,d,y], 
                  apply(offset.forecast(ecmwf)[1:2,d,y,lt,-1], 1, mean),
           apply(offset.forecast(ncep)[1:2,d,y,lt,-1], 1, mean),
           apply(offset.forecast(ukmo)[1:2,d,y,lt,-1], 1, mean))

fc.sig <- cov(fc.mu)

lt <- "14"
Y.bar <- abind("ecmwf" = apply(offset.forecast(ecmwf)[1:2,,,lt,-1], 1:3, mean),
                "ncep" = apply(offset.forecast(ncep)[1:2,,,lt,-1], 1:3, mean),
                "ukmo" = apply(offset.forecast(ukmo)[1:2,,,lt,-1], 1:3, mean),
                along = 0)

# model current ensemble forecasts as distribution, find distance to all previous points
{
    Y.bar <- abind("ecmwf" = apply(offset.forecast(ecmwf)[1:2,,,lt,-1], 1:3, mean),
                   "ncep" = apply(offset.forecast(ncep)[1:2,,,lt,-1], 1:3, mean),
                   "ukmo" = apply(offset.forecast(ukmo)[1:2,,,lt,-1], 1:3, mean),
                   along = 0)

    sig <- square.mat(apply(Y.bar, 3:4, cov))
    
    C <- abind("ecmwf" = square.mat(apply(aperm(offset.forecast(ecmwf)[,,,lt,-1], c(4,1:3)), 3:4, cov)),
               "ncep" = square.mat(apply(aperm(offset.forecast(ncep)[,,,lt,-1], c(4,1:3)), 3:4, cov)),
               "ukmo" = square.mat(apply(aperm(offset.forecast(ukmo)[,,,lt,-1], c(4,1:3)), 3:4, cov)), 
               along = 0)
}

md <- abind("ecmwf" = mahalanobis(x = apply(apply(offset.forecast(ecmwf)[1:2,,,lt,-1], 1:3, mean), 1, rbind),
                                  center = apply(offset.forecast(ecmwf)[1:2,5,1,lt,-1], 1, mean),
                                  cov = cov(t(offset.forecast(ecmwf)[1:2,5,1,lt,-1]))),
            "ncep" = mahalanobis(x = apply(apply(offset.forecast(ncep)[1:2,,,lt,-1], 1:3, mean), 1, rbind),
                                 center = apply(offset.forecast(ncep)[1:2,5,1,lt,-1], 1, mean),
                                 cov = cov(t(offset.forecast(ncep)[1:2,5,1,lt,-1]))),
            "ukmo" = mahalanobis(x = apply(apply(offset.forecast(ukmo)[1:2,,,lt,-1], 1:3, mean), 1, rbind),
                                 center = apply(offset.forecast(ukmo)[1:2,5,1,lt,-1], 1, mean),
                                 cov = cov(t(offset.forecast(ukmo)[1:2,5,1,lt,-1]))),
            along = 0)

####################################################################################################

# PRINCIPAL COMPONENTS IN ANALOGUE SELECTION                                                    ####

# start with checking analogues to a single forecast.
require("SX.weather")
lt <- "3"; d <- 2; y <- 1

# target forecast: single day. Use yesterday's observation as part of vector
t <- rbind("o" = obs[,d-1,y],
           "ecmwf" = apply(offset.forecast(ecmwf)[,d,y,lt,-1], 1, mean),
           "ncep" = apply(offset.forecast(ncep)[,d,y,lt,-1], 1, mean),
           "ukmo" = apply(offset.forecast(ukmo)[,d,y,lt,-1], 1, mean))

# candidate analogues (offset observations by 1 day)
cand <- abind("o" = obs[,1:89,],
              "ecmwf" = apply(offset.forecast(ecmwf)[,2:90,,lt,-1], 1:3, mean),
              "ncep" = apply(offset.forecast(ncep)[,2:90,,lt,-1], 1:3, mean),
              "ukmo" = apply(offset.forecast(ukmo)[,2:90,,lt,-1], 1:3, mean),
              along = 0)

# standard deviation of each variable, for normalisation
cand.sd <- apply(cand, 1:2, sd)

# calculate normalised point distances for all possible analogues
dist <- sweep(sqrt(sweep(cand, 1:2, t, "-")^2), 1:2, cand.sd, "/")

t.dist.2 <- apply(dist[,1:2,,], 3:4, mean)
t.dist.5 <- apply(dist, 3:4, mean)

o.2 <- order(t.dist.2)
o.5 <- order(t.dist.5)

plot(c(t.dist.2), c(t.dist.5)[order(t.dist.2)])


####################################################################################################

# DB CLUSTERING OF FORECASTS                                                                    ####

db.plot <- function(d, y, lt, model, eps = 2, MinPts = 2) {
    
    require(plyr)
    yy <- formatC(7:14, width = 2, flag = "0")[y]
    
    fc <- t(offset.forecast(model)[1:2,d,y,lt,-1])
    db <- dbscan(fc, eps = eps, MinPts = MinPts)
    
    rng <- apply(model[1:2,,,,-1], 1, range)
    
    ccols <- c("grey", "red", "blue", "gold", "green3")
    
    plot(fc, col = adjustcolor(ccols[db$cluster+1], alpha = 0.5), pch = 20, 
         xlim = rng[,1], ylim = rng[,2], main = "")
    legend("topleft", legend = paste0(yy, "d", d, ", LT ", lt-1), bty = "n")
    points(t(apply(fc, 2, mean)), pch = 20, lwd = 2)
    
    # include verifying observation
        points(t(obs[1:2, d, y]), pch = 4, lwd = 2, cex = 1)
    
    # also add individual cluster means (IF multiple clusters identified)
    if(max(db$cluster) > 1) {
        
        df <- ddply(data.frame(fc, cl = db$cluster), .(cl),
                    summarise, tn.mean = mean(temp.n), ts.mean = mean(temp.s))
        points(df[df$cl > 0,2:3], bg = ccols[df$cl[df$cl > 0]+1], pch = 21, lwd = 1)
    }
}

db <- apply(offset.forecast(ukmo)[1:2,,,,-1], 2:4, function(fc) dbscan(t(fc), eps = 3, MinPts = 2)$cluster)

db.max <- apply(db, 2:4, max)
db.clusters <- which(apply(db, 2:4, max) > 1, arr.ind = T)
db.clusters <- db.clusters[order(db.clusters[,3], db.clusters[,2], db.clusters[,1]),]

db.plot(10, 5, 11, ukmo)
db.plot(20, 1, 14, ukmo)
db.plot(7, 4, 8, ukmo)

pdf("./Documents/PhD/Miniprojects/03-Sichun-paper/Plots/Ensemble-clustering-UKMO.pdf", height = 15, width = 15); {
    par(mfrow = c(7, 5), mar = c(2,2,1,1), oma = c(0,0,2,0))
    
    for (p in 0:2) {
        invisible(apply(db.clusters[(1:35) + (p*35),], 1, function(dt) db.plot(dt[1], dt[2], dt[3], ukmo)))
    }
    
    mtext("Clusters in ensemble forecasts: ukmo", outer = T)
}; dev.off()


# numbers by leadtime
table(db.clusters[,3])

#        2   3   4   5   6   7   8   9  10  11  12  13  14  
# ECMWF      4  10  20  25  29  40  47  44  51  66  76          412
# NCEP       1   7  22  34  53  73 103 111 112 121 128 139      904
# UKMO   1   7  10  21  34  63  93 111 146 137 157 175 169      1124

# interesting. ECMWF seems much more confident in its predictions.
# how does this stack up against forecast error?

ec <- apply(offset.forecast(ecmwf)[1:2,,,,-1], 4, function(ens) {
    apply(abind("o" = obs[1:2,,], "fc" = ens, along = 4), 2:3, 
          function(fc) {es.crps(o = fc[,1], efc = fc[,-1])})})

nc <- apply(offset.forecast(ncep)[1:2,,,,-1], 4, function(ens) {
    apply(abind("o" = obs[1:2,,], "fc" = ens, along = 4), 2:3, 
          function(fc) {es.crps(o = fc[,1], efc = fc[,-1])})})
    
mo <- apply(offset.forecast(ukmo)[1:2,,,,-1], 4, function(ens) {
    apply(abind("o" = obs[1:2,,], "fc" = ens, along = 4), 2:3, 
          function(fc) {es.crps(o = fc[,1], efc = fc[,-1])})})

matplot(cbind(apply(ec, 2, mean), apply(nc, 2, mean), apply(mo, 2, mean)),
        type = "l", lty = 1, col = c("red", "blue", "green3"))

matplot(ec, type = "l", lty = 1, col = adjustcolor("steelblue", alpha = 0.4), ylim = c(0,12))
lines(ec[,1], col = "black")
lines(ec[,15], col = "darkred")

matplot(nc, type = "l", lty = 1, add = T, col = adjustcolor("green3", alpha = 0.4), ylim = c(0,12))

matplot(mo, type = "l", lty = 1, add = T, col = adjustcolor("red3", alpha = 0.4), ylim = c(0,12))

matplot(cbind(apply(ec, 1, mean), apply(nc, 1, mean), apply(mo, 1, mean)), type = "l", lty = 1,
        col = c("red", "blue", "green3"), ylab = "")

fbma <- readRDS("./Documents/PhD/Miniprojects/03-Sichun-paper/Models/ensBMA-temps.rds")

qq <- sapply(fbma, "[[", "weights")
sum(qq < 0)

# no negative weights observed.

####################################################################################################

# BMA COUNTER-EXAMPLE                                                                           ####

# look for an example where the observation is definitely outside the forecast region
hull.plot <- function(d, y, lt) {
    
    full.rng <- apply(abind(array(rep(obs[1:2,,],15), dim = c(2,90,7,15)),
                      offset.forecast(ecmwf)[1:2,,,,-1],
                      offset.forecast(ncep)[1:2,,,,-1],
                      offset.forecast(ukmo)[1:2,,,,-1], along = 5), 1, range)
        
    lt <- toString(lt)
    yy <- formatC(7:14, width = 2, flag = "0")[y]
    c.cols = c("black", sapply(c("coral3", "green3", "steelblue"), adjustcolor, alpha = 0.4))
    
    o <- obs[1:2,d,y]
    ec <- t(offset.forecast(ecmwf)[1:2,d,y,lt,-1])
    nc <- t(offset.forecast(ncep)[1:2,d,y,lt,-1])
    mo <- t(offset.forecast(ukmo)[1:2,d,y,lt,-1])
    
    ind <- c(1, rep(2, 50), rep(3, 20), rep(4, 23))
    
    plot(abind(o, ec, nc, mo, along = 1), col = c.cols[ind], pch = c(4, rep(20, 93)), lwd = 2,
         xlim = full.rng[,1], ylim = full.rng[,2],
         main = paste0("Ensemble overlap - 20", yy, "-", d, "; lt ", lt))
    
    ch.ec <- chull(ec); ch.nc <- chull(nc); ch.mo <- chull(mo)
    ch <- chull(abind(o, ec, nc, mo, along = 1))

    lines(ec[c(ch.ec, ch.ec[1]),], col = c.cols[2])
    lines(nc[c(ch.nc, ch.nc[1]),], col = c.cols[3])
    lines(mo[c(ch.mo, ch.mo[1]),], col = c.cols[4])
    
    lines(abind(o, ec, nc, mo, along = 1)[c(ch, ch[1]),])
    
    legend("topleft", pch = c(rep(20, 3), NA, 4), col = c(c.cols[2:4], NA, "black"), bty = "n",
           legend = c("ECMWF", "NCEP", "UKMO", NA, "Obs"))
    
    win <- owin(full.rng[,1], full.rng[,2])
    
    # check for crossings
    ec.psp <- as.psp(cbind(ec[ch.ec,], ec[c(ch.ec[-1], ch.ec[1]),]), window = win)
    nc.psp <- as.psp(cbind(nc[ch.nc,], nc[c(ch.nc[-1], ch.nc[1]),]), window = win)
    mo.psp <- as.psp(cbind(mo[ch.mo,], mo[c(ch.mo[-1], ch.mo[1]),]), window = win)
    
    ch.psp <- as.psp(cbind(abind(o, ec, nc, mo, along = 1)[ch,],
                           abind(o, ec, nc, mo, along = 1)[c(ch[-1], ch[1]),]),
                     window = win)
    
    # line from upper corner of window to observation
    obs.psp <- as.psp(list(x0 = full.rng[1,1], y0 = full.rng[2,2], x1 = o[1], y1 = o[2]), window = win)
    
    # check crossings
    cr <- sapply(list("super" = crossing.psp(ch.psp, obs.psp),
               "ecmwf" = crossing.psp(ec.psp, obs.psp),
               "ncep" = crossing.psp(nc.psp, obs.psp),
               "ukmo" = crossing.psp(mo.psp, obs.psp)),
               function(px) length(px$x) %% 2 == 1)
    
    print(cr)
}

hull.plot(d, y, lt)
invisible(sapply(14:0, function(ll) hull.plot(1,2,ll)))

ch.ecmwf <- apply(offset.forecast(ecmwf)[1:2,,,,-1], 2:4, function(fc) {
        ch <- chull(t(fc))
        t(fc[,c(ch, ch[1])])
    })

library(spatstat)
full.rng <- apply(abind(array(rep(obs[1:2,,],15), dim = c(2,90,7,15)),
                        offset.forecast(ecmwf)[1:2,,,,-1],
                        offset.forecast(ncep)[1:2,,,,-1],
                        offset.forecast(ukmo)[1:2,,,,-1], along = 5), 1, range)



# and how frequently within each of the others?
cr.ecmwf <- apply(offset.forecast(ecmwf)[1:2,,,,-1], 4, function(fc) {
    dat <- abind("o" = obs[1:2,,], fc, along = 4)
    apply(dat, 2:3, function(ens) {
        
        ch <- chull(t(ens[,-1]))
        ch.o <- chull(t(ens))
        
        if(length(ch) == length(ch.o)) {all((ch+1) == ch.o)} else {F}
        
        #plot(t(ens), col = adjustcolor("steelblue", alpha = 0.5), pch = 20)
        #points(t(ens[,1]), pch = 4, lwd = 2)
        
        #lines(t(ens[,c(ch, ch[1])+1]), col = "steelblue")
        #lines(t(ens[,c(ch.o, ch.o[1])]), col = "red")
    })
})
cr.ncep <- apply(offset.forecast(ncep)[1:2,,,,-1], 4, function(fc) {
    dat <- abind("o" = obs[1:2,,], fc, along = 4)
    apply(dat, 2:3, function(ens) {
        
        ch <- chull(t(ens[,-1]))
        ch.o <- chull(t(ens))
        
        if(length(ch) == length(ch.o)) {all((ch+1) == ch.o)} else {F}
    })
})
cr.ukmo <- apply(offset.forecast(ukmo)[1:2,,,,-1], 4, function(fc) {
    dat <- abind("o" = obs[1:2,,], fc, along = 4)
    apply(dat, 2:3, function(ens) {
        
        ch <- chull(t(ens[,-1]))
        ch.o <- chull(t(ens))
        
        if(length(ch) == length(ch.o)) {all((ch+1) == ch.o)} else {F}
    })
})

cr <- rbind("ecmwf" = c(cr.ecmwf), "ncep" = c(cr.ncep), "ukmo" = c(cr.ukmo))

# ECMWF 4700  (49.7%)
# NCEP  3763  (39.8%)
# UKMO  4716  (49.9%)

sum(colSums(cr) == 0)       
# not in any ensemble:  3015  (32%)
# a single ensemble:    2117  (22%)      ECMWF 804      NCEP  482       UKMO 831
# exactly two:          1892  (20%)     !ECMWF 422     !NCEP 1037      !UKMO 433    
# all three:            2426  (26%)

{
    cr.ens.mean <- apply(abind("ecmwf" = apply(offset.forecast(ecmwf)[1:2,,,,-1], 1:4, mean),
                               "ncep" = apply(offset.forecast(ncep)[1:2,,,,-1], 1:4, mean),
                               "ukmo" = apply(offset.forecast(ukmo)[1:2,,,,-1], 1:4, mean),
                               rev.along = 0),
                         4, function(fc) {
                             dat <- abind("o" = obs[1:2,,], fc, along = 4)
                             apply(dat, 2:3, function(ens) {
                                 
                                 ch <- chull(t(ens[,-1]))
                                 ch.o <- chull(t(ens))
                                 
                                 if(length(ch) == length(ch.o)) {all((ch+1) == ch.o)} else {F}
                             })
                         })
}   # within the convex hull of the ensemble means:   392  (4%)
{
    cr.ens <- apply(abind("ecmwf" = offset.forecast(ecmwf)[1:2,,,,-1],
                      "ncep" = offset.forecast(ncep)[1:2,,,,-1], 
                      "ukmo" = offset.forecast(ukmo)[1:2,,,,-1],
                      along = 5),
                4, function(fc) {
                    dat <- abind("o" = obs[1:2,,], fc, along = 4)
                    apply(dat, 2:3, function(ens) {
                        
                        ch <- chull(t(ens[,-1]))
                        ch.o <- chull(t(ens))
                        
                        if(length(ch) == length(ch.o)) {all((ch+1) == ch.o)} else {F}
                    })
                })}   # within the convex hull of the superensemble:   7020 (74%)

# after NCEP bias correction? ()
cr.post <- rbind("ecmwf" = c(cr.ecmwf[271:630,]), "ncep" = c(cr.ncep[271:630,]), "ukmo" = c(cr.ukmo[271:630,]))

sum(colSums(cr.post) == 3) 
# not in any ensemble:  2031  (38%)
# a single ensemble:     999  (19%)      ECMWF 369      NCEP 224       UKMO 406
# exactly two:           922  (17%)     !ECMWF 220     !NCEP 396      !UKMO 206    
# all three:            1448  (26%)