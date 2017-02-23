
# all models trained over analogue set (rather than historic data)

library("SX.weather")

l.cols <- c("black", "blue", "red3"); l.wd <- c(2,1,1)

####################################################################################################

# CREATE TRAINING SET                                                                           ####

# get array of 25 analogues per day/year/leadtime (~75s)
an.indx <- abind(lapply(1:15, function(lt) {
    
    # search space - all candidates at that leadtime
    cand <- abind("y.o" = obs[,1:89,],
                  "y.ecmwf" = apply(offset.forecast(ecmwf)[,1:89,,lt,-1], 1:3, mean),
                  "y.ncep" = apply(offset.forecast(ncep)[,1:89,,lt,-1], 1:3, mean),
                  "y.ukmo" = apply(offset.forecast(ukmo)[,1:89,,lt,-1], 1:3, mean),
                  "c.ecmwf" = apply(offset.forecast(ecmwf)[,2:90,,lt,-1], 1:3, mean),
                  "c.ncep" = apply(offset.forecast(ncep)[,2:90,,lt,-1], 1:3, mean),
                  "c.ukmo" = apply(offset.forecast(ukmo)[,2:90,,lt,-1], 1:3, mean),
                  along = 0)
    cand <- abind(array(NA, dim = dim(cand[,,1,])), cand, along = 3)

    # indices of 25 closest analogues
    an <- aaply(cand, 3:4, function(targ) {
        
        if (sum(!is.na(targ)) == 0) {
            # if target is padding (as when running over array), return an NA array of expected size
            return(array(NA, dim = c(25, 2)))
        }
        
        # standard deviation for normalisation
        cand.sd <- apply(cand, 1:2, sd, na.rm = T)
        
        # distance from each point on target to equivalent in candidate space
        dist <- sweep(sqrt(sweep(cand, 1:2, targ, "-")^2), 1:2, cand.sd, "/")
        
        # mean distance across all variables
        m.dist <- apply(dist, 3:4, mean)
        
        # identify n analogues
        an <- which(m.dist <= sort(m.dist)[26], arr.ind = T)
        
        # find target vector and remove from analogue list
        t <- which(m.dist == 0, arr.ind = T)
        an <- an[!(an[,1] == t[1] & an[,2] == t[2]),]
    })
}), along = 2.5)
   

# extract ensemble means & verifying observations for each analogue
em <- abind("ecmwf" = apply(offset.forecast(ecmwf)[,,,,-1], 1:4, mean),
            "ncep" = apply(offset.forecast(ncep)[,,,,-1], 1:4, mean),
            "ukmo" = apply(offset.forecast(ukmo)[,,,,-1], 1:4, mean),
            along = 0)

tr.fc <- abind(sapply(1:15, function(lt) {
    aaply(an.indx, 1:2, function(a) {
        apply(em[,,,,lt], 1:2, "[", a[lt,,])
    })}, simplify = F), along = 2.5)

tr.obs <- aaply(an.indx, 1:3, function(a) {
    apply(obs, 1, "[", a)
})

####################################################################################################

# FIND BAYESIAN POSTERIOR                                                                       ####

# training data: forecast errors
tr.error <- aperm(apply(sweep(tr.fc, c(1:4, 6), tr.obs, "-"), c(1:4, 6), mean), c(5, 1:4))

# fit mimic at all days/years/leadtimes (~60s)
mimic <- abind(sapply(1:15, function(lt) {
    src <- abind("o" = obs,
                 "ecmwf" = offset.forecast(ecmwf)[,,,lt,-1],
                 "ncep" = offset.forecast(ncep)[,,,lt,-1],
                 "ukmo" = offset.forecast(ukmo)[,,,lt,-1],
                 "tr" = tr.error[,,,lt,],
                 along = 4)
    
    suppressWarnings(aaply(src, 2:3, function(arr) {
        fit.mimic(o = arr[,1],
                  fc = list("ecmwf" = arr[,2:51], "ncep" = arr[,52:71], "ukmo" = arr[,72:94]),
                  tr = arr[,95:119])
    }))
}, simplify = F), along = 2.5)

saveRDS(mimic, "../Models/mimic-an25.rds")

####################################################################################################

# FIT ENSEMBLE BMA & ENSEMBLE MOS                                                               ####

require(ensembleMOS); require(ensembleBMA)

ts <- load.data("../Data/ECMWF_europe_starttime.rda")[16:105,1:7]/24

# combine daily forecast/obs with analogues for easier processing
all.fc <- apply(abind(aperm(tr.fc, c(5:6,1:3,4)), em, along = 6)[,1:2,,,,], c(1,3:5),
                function(vv) c(t(vv)))

d <- 5; y <- 5; lt <- 1;

# prep arrays to hold all results
bma.fitted <- array(dim = c(90,7,15,11), 
                    dimnames = list(NULL, NULL, NULL,
                                    c("ec1", "ec2", "nc1", "nc2", "mo1", "mo2",
                                      "sd", "w.ec", "w.nc", "w.mo", "ll")))
mos.fitted <- array(dim = c(90,7,15,6),
                    dimnames = list(NULL, NULL, NULL,
                                    c("a", "b1", "b2", "b3", "c", "d")))

q.fc <- array(dim = c(90,7,15,2,2,2),
              dimnames = list(NULL, NULL, NULL, c("bma", "mos"), c("temp.n", "temp.s"), c("fc", "crps")))

# ran over all dates & leadtimes: ~55m (constructing data & fitting models)
# currently unable to calculate CRPS automatically. Code is ready if I fix this
system.time({
invisible(sapply(1:15, function(lt) {
    invisible(sapply(2:90, function(d) {
        invisible(sapply(1:7, function(y) {
            
            # create ensemble object to run model
            eo <- ensembleData(forecasts = all.fc[,,d,y,lt],
                               observations = c(rbind(tr.obs[d,y,lt,,1:2], obs[1:2,d,y])), 
                               dates = rep(gsub("-", "", as.Date(ts[d,y]-(25:0), "1800-01-01")), 2),
                               forecastHour = 00,
                               station = rep(c("N", "S"), each = 26),
                               initializationTime = "00")
            
            tr <- trainingData(eo, 25, date = max(ensembleValidDates(eo)))
#            targ <- ensembleData(forecasts = t(em[,1:2,d,y,lt]),
#                                 observations = obs[1:2,d,y],
#                                 dates = max(ensembleValidDates(eo)),
#                                 forecastHour = 00,
#                                 initializationTime = "00")
                
            bma <- fitBMA(tr, model = "normal")
            mos <- fitMOS(tr, model = "normal")
            
            bma.fitted[d,y,lt,] <<- unlist(bma)[-11]
            mos.fitted[d,y,lt,] <<- unlist(mos)[1:6]
            
            q.fc[d,y,lt,"bma",,"fc"] <<- quantileForecast(bma, eo)[c(26, 52)]
            q.fc[d,y,lt,"mos",,"fc"] <<- quantileForecast(mos, eo)[c(26, 52)]
            
#            q.fc[d,y,lt,"mos",,"crps"] <<- crps(mos, targ)
        }))
    }))
}))
})

saveRDS(bma.fitted, "../Models/an25-bma-fitted.rds")
saveRDS(mos.fitted, "../Models/an25-mos-fitted.rds")
saveRDS(q.fc, "../Models/an25-bma-mos-fc.rds")


####################################################################################################

# CREATE MOS & BMA SIMULATIONS FOR ASSESSMENT                                                   ####

# PIT and CRPS functions are strangely unreliable. 
# Use simulations instead to get verification ranks & energy score
bma.fitted <- readRDS("../Models/an25-bma-fitted.rds")
mos.fitted <- readRDS("../Models/an25-mos-fitted.rds")

require(mixtools)

fcm <- abind("ecmwf" = apply(offset.forecast(ecmwf)[,,,,-1], 1:4, mean),
             "ncep" = apply(offset.forecast(ncep)[,,,,-1], 1:4, mean),
             "ukmo" = apply(offset.forecast(ukmo)[,,,,-1], 1:4, mean),
             along = 0)

bma.mu <- abind("ecmwf" = sweep(sweep(fcm[1,1:2,,,], 2:4, bma.fitted[,,,"ec2"], "*"),
                                2:4, bma.fitted[,,,"ec1"], "+"),
                "ncep" = sweep(sweep(fcm[2,1:2,,,], 2:4, bma.fitted[,,,"nc2"], "*"),
                               2:4, bma.fitted[,,,"nc1"], "+"),
                "ukmo" = sweep(sweep(fcm[3,1:2,,,], 2:4, bma.fitted[,,,"mo2"], "*"),
                               2:4, bma.fitted[,,,"mo1"], "+"),
                along = 0)
bma.pred <- apply(sweep(bma.mu, c(1, 3:5), aperm(bma.fitted[,,,8:10], c(4,1:3)), "*"), 2:5, sum)

mos.mu <- abind("ecmwf" = sweep(fcm[1,1:2,,,], 2:4, mos.fitted[,,,"b1"], "*"),
                "ncep" = sweep(fcm[2,1:2,,,], 2:4, mos.fitted[,,,"b2"], "*"),
                "ukmo" = sweep(fcm[3,1:2,,,], 2:4, mos.fitted[,,,"b3"], "*"),
                along = 0)
mos.pred <- sweep(apply(mos.mu, 2:5, sum), 2:4, mos.fitted[,,,"a"], "+")

# c + ds^2, S^2 is ensemble variance
mos.var <- sweep(sweep(apply(fcm[,1:2,,,], 2:5, var), 
                       2:4, mos.fitted[,,,"d"], "*"),
                 2:4,  mos.fitted[,,,"c"], "+")

# not sure why this is coming out differently. Rounding perhaps?
plot(bma.pred, aperm(q.fc[,,,"bma",], c(4,1:3)), pch = ".")
abline(0,1,col = "red")
bma.diff <- bma.pred - aperm(q.fc[,,,"bma",], c(4,1:3))
quantile(abs(bma.diff), 0.95, na.rm = T)
# 95th quantile of errors is < 0.03. Hopefully accurate enough for current purposes.

plot(mos.pred, aperm(q.fc[,,,"mos",], c(4,1:3)), pch = ".")
abline(0,1,col = "red")
max(abs(mos.pred - aperm(q.fc[,,,"mos",], c(4,1:3))), na.rm = T)
# Max difference is 3.6e-15. Close enough
# what about variance? checked agaisnt quantileForecast function, this calculation is correct.

quantileForecast(mos, eo, c(0.05, 0.5, 0.95))[c(26, 52),]
qnorm(c(0.05, 0.5, 0.95), mean = mos.pred[1,5,5,1], sd = sqrt(mos.var[1,5,5,1]))
qnorm(c(0.05, 0.5, 0.95), mean = mos.pred[2,5,5,1], sd = sqrt(mos.var[2,5,5,1]))


sim <- array(dim = c(2,2,90,7,15,1000), 
             dimnames = list(c("bma", "mos"), c("temp.n", "temp.s"), NULL, NULL, NULL, NULL))

# over all leadtimes & years: 10s for both models.
set.seed("24747")

invisible(sapply(1:15, function(lt) {
    invisible(sapply(1:7, function(y) {
        invisible(sapply(2:90, function(d) {
            
            sim["bma","temp.n",d,y,lt,] <<- rnormmix(1000, 
                                                     lambda = bma.fitted[d,y,lt,8:10], 
                                                     mu = bma.mu[,"temp.n",d,y,lt], 
                                                     sigma = bma.fitted[d,y,lt,7])
            sim["bma","temp.s",d,y,lt,] <<- rnormmix(1000, 
                                                     lambda = bma.fitted[d,y,lt,8:10],
                                                     mu = bma.mu[,"temp.s",d,y,lt], 
                                                     sigma = bma.fitted[d,y,lt,7])
            
            sim["mos","temp.n",d,y,lt,] <<- rnorm(1000, 
                                                  mean = mos.pred["temp.n",d,y,lt],
                                                  sd = sqrt(mos.var["temp.n",d,y,lt]))
            sim["mos","temp.s",d,y,lt,] <<- rnorm(1000, 
                                                  mean = mos.pred["temp.s",d,y,lt],
                                                  sd = sqrt(mos.var["temp.s",d,y,lt]))
        }))
    }))
}))

saveRDS(sim, "../Models/an25-bma-mos-sim1000.rds")

# use same source data to get predicted probability of freezing on each day
mos.fprob <- apply(abind("mean" = mos.pred, "var" = mos.var, along = 0), 2:5,
                   function(mos) pnorm(0, mos["mean"], sqrt(mos["var"])))
                   

mos.fprob <- pnorm(0, mean = mos.pred["temp.n",d,y,lt], sd = sqrt(mos.var["temp.n",d,y,lt]))

bma.lambda <- aperm(bma.fitted[,,,8:10], c(4,1:3))
bma.sigma <- bma.fitted[,,,7]

fp.bma <- array(dim = dim(bma.pred), dimnames = dimnames(bma.pred))

f0.bma <- apply(sweep(aaply(bma.mu, 1:2, 
                            function(arr) pnorm(0, mean = arr, sd = bma.fitted[,,,7])),
                      c(1,3:5), aperm(bma.fitted[,,,8:10], c(4,1:3)), "*"),
                2:5, sum)
                
    
####################################################################################################

# COMPARE OUTPUT                                                                                ####

q.fc <- readRDS("../Models/an25-bma-mos-fc.rds")
mimic <- readRDS("../Models/an25-mimic.rds")
sim <- readRDS("../Models/an25-bma-mos-sim1000.rds")

#===============================================================================================

# --> multivariate CRPS - energy score                                                         ####
{
    v.dat <- aaply(mimic, 3, function(arr) abind("o" = aperm(obs, c(2:3, 1)), arr, along = 4))
    es <- apply(v.dat[,-1,,,], 1:3, function(arr) es.crps(o = arr[,1], mu = arr[,2], sig = arr[,-(1:2)]))
    boxplot(apply(es, 1, cbind))
}

# energy score from simulations
es <- function(x, o, k = length(x)) {
    norm.1 <- mean(apply(t(x) - o, 2, function(v) sqrt(sum(v^2))))
    norm.2 <- sum(sapply(x[1:k - 1] - x[2:k], function(v) sqrt(sum(v^2))))/(2 * (k - 1))
    norm.1 - norm.2
}

es.auto <- aaply(sim[,,-1,,,], c(1, 5), function(arr) {
    apply(abind(obs[1:2,-1,], arr), 1:3, function(v) es(x = v[-1], o = v[1]))
})

es.mim <- aaply(mim.sim, 3, function(arr) {
    apply(abind(obs[1:2,-1,], aperm(arr, c(4,1:3))), 1:3, function(v) es(x = v[-1], o = v[1]))
})

es.hist <- aaply(hist.sim, 1, function(arr) {
    apply(abind(o.605, aperm(arr, c(3,1:2))), 1:2, function(v) es(x = v[-1], o = v[1]))
})

crps.es <- abind("mimic" = apply(es.mim, 1:2, mean),
            apply(es.auto, 1:3, mean),
            "mimic.hist" = apply(es.hist, 1:2, mean), along = 1)
dimnames(crps.es)[[3]] <- dimnames(obs)[[1]][1:2]


#===============================================================================================

# --> univariate RMSE for all models                                                        ####
{
    err <- abind("mimic" = sweep(aperm(mimic[,,,1:2,1], c(4, 1:3)), 1:3, obs[1:2,,], "-"),
                 sweep(aperm(q.fc, c(4:5,1:3)), 2:4, obs[1:2,,], "-"), 
                 along = 1)
    
    rmse.lt <- sqrt(apply(err^2, c(1,2,5), mean, na.rm = T))
    
    # add RMSE of original models
    rmse.lt <- abind(rmse.lt,
                     "mimic.org" = t(sqrt(apply(sweep(readRDS("../Models/Mimic-hist.rds")[,,1:2,1], 2:3, 
                                                      apply(obs[1:2,,], 1, c)[-(1:25),], "-")^2,
                                                c(1,3), mean))),
                     "bma.org" = sqrt(apply(sweep(readRDS("../Models/ensBMA-fc.rds")[,"0.5",,], 1:2, 
                                                  t(apply(obs[1:2,,], 1, c)[-(1:24),]), "-")^2,
                                            c(1,3), mean)),
                     "bma.mos" = sqrt(apply(sweep(readRDS("../Models/ensMOS-fc.rds")[,"0.5",,], 1:2, 
                                                  t(apply(obs[1:2,,], 1, c)[-(1:24),]), "-")^2,
                                            c(1,3), mean)),
                     along = 1)
}

# plots
pdf("../Plots/Analogue-perf/RMSE.pdf", height = 4, width = 7); {
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0)); 
    
    matplot(0:14, t(rmse.lt[,"temp.n",]), type = "l", lwd = rep(c(2,1,1), 2),
            col = rep(c("black", "blue", "red3"), 2), lty = rep(c(1, 2), each = 3),
            xlab = "", ylab = "", main = "temp.n", ylim = range(0, rmse.lt))
    legend("bottomright", col = c("black", "blue", "red3"), lty = 1,
           legend = dimnames(rmse.lt)[[1]][1:3], bty = "n")
    legend("topleft", lty = c(1,2), bty = "n", legend = c("Analogues", "Historic"))
    abline(v = 7, col = "cyan3", lty = 3)
    
    matplot(0:14, t(rmse.lt[,"temp.s",]), type = "l", lwd = rep(c(2,1,1), 2),
            col = rep(c("black", "blue", "red3"), 2), lty = rep(c(1, 2), each = 3),
            xlab = "", ylab = "", main = "temp.s", ylim = range(0, rmse.lt))
    legend("bottomright", col = c("black", "blue", "red3"), lty = 1,
           legend = dimnames(rmse.lt)[[1]][1:3], bty = "n")
    legend("topleft", lty = c(1,2), bty = "n", legend = c("Analogues", "Historic"))
    abline(v = 7, col = "cyan3", lty = 3)
    
    mtext("RMSE: analogue vs historic training set", outer = T)
    
    par(mfrow = c(1,1))
}; dev.off()

# boxplots
{
    bp.err <- apply(abs(err), c(1:2, 5), c)
    bp.m <- apply(bp.err, 2:4, median, na.rm = T)
    bp.cols <- adjustcolor(c("coral", "green3", "skyblue"), alpha = 0.4)
    m.cols <- c("red", "green3", "blue")
        
    par(mfrow = c(2,3))
    
    invisible(sapply(dimnames(bp.err)[[3]], function(varb) {
        invisible(sapply(1:3, function(model) {
            
            boxplot(bp.err[,model,varb,], main = dimnames(bp.err)[[2]][model], names = 0:14, 
                    col = bp.cols[model], ylim = range(0, bp.err, na.rm = T))
            
            matplot(t(bp.m[,varb,]), type = "l", add = T, lwd = 2, lty = 1, col = m.cols)
            
            abline(v = 7.5, col = "cyan3", lty = 3)
            legend("topleft", legend = varb, bty = "n")
        }))
    }))
    mtext("Mean asolute error", outer = T)
}

#===============================================================================================

# --> verification rank histograms (based on simulations where necessary)                   ####

vr <- aaply(sim, c(1, 5), function(ens) {
    apply(abind(obs[1:2,,], ens), 1:3, function(v) which(sort(v) == v[1]))})

# marginal PIT
pit <- aaply(mimic[,,,1:2,1:3], 3, function(mim) {
    abind("temp.n" = pnorm(obs["temp.n",,], mim[,,"temp.n", "Tau"], sqrt(mim[,,"temp.n", "temp.n"])),
          "temp.s" = pnorm(obs["temp.s",,], mim[,,"temp.s", "Tau"], sqrt(mim[,,"temp.s", "temp.s"])),
          along = 0)
})

# multivariate Box DOT is not comparable - centre-outward ordering instead
pdf("../Plots/Analogue-perf/VR-hists.pdf", height = 4, width = 9); {
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(1:15, function(lt) {
        invisible(sapply(dimnames(vr)[[3]], function(varb) {
            
            hist(vr["bma",lt,varb,,] /1001, prob = T, xlab = "", ylab = "", ylim = c(0,1.9),
                 main = varb, col = adjustcolor("skyblue", alpha = 0.4), 
                 border = adjustcolor("skyblue", alpha = 0.7))
            
            hist(vr["mos",lt,varb,,]/1001, prob = T, add = T, col = adjustcolor("gold", alpha = 0.4), 
                 border = adjustcolor("gold", alpha = 0.7))
            
            hist(pit[lt,varb,,], prob = T, add = T)
           # hist(vr.mim[lt, varb,,]/1001, prob = T, add = T, border = "red", lty = 2)
            
            abline(h = 1, col = "red")
        }))
        mtext(paste0("Verification rank histograms (blue = BMA, gold = MOS, black = mimic); LT ", lt-1),
              outer = T)
    }))
}; dev.off()

#===============================================================================================

# --> deviation from uniformity                                                             ####

# simulate from mimic & use to obtain verification ranks
require(mvtnorm)
mim.sim <- aaply(mimic[-1,,,1:2,1:3], 1:3, function(mim) rmvnorm(1000, mean = mim[,"Tau"], sigma = mim[,-1]))
vr.mim <- aaply(aperm(mim.sim, c(5,1:4)), 4, function(ens) {
    apply(abind(obs[1:2,-1,], ens), 1:3, function(v) which(sort(v) == v[1]))})
dimnames(vr.mim)[[2]] <- dimnames(obs)[[1]][1:2]

ud <- abind("mimic" = apply(vr.mim, 1:2, u.dev, l = 1001),
            apply(vr, 1:3, u.dev, l = 1001), along = 1)

pdf("../Plots/Analogue-perf/deviation from uniformity.pdf", height = 4, width = 9); {
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(ud)[[3]], function(varb) {
        matplot(0:14, t(ud[,,varb]), type = "l", ylim = range(ud, 0, na.rm = T),
                xlab = "", ylab = "", main = varb, col = l.cols, lwd = l.wd)
        legend("bottomright", bty = "n", col = l.cols, lty = 1, lwd = l.wd, legend = dimnames(ud)[[1]])
    }))
    mtext("Deviation from uniformity", outer = T)
}; dev.off()


#===============================================================================================

# --> skewness of verification ranks                                                        ####

require(moments)
sk <- abind("mimic" = apply(vr.mim, 1:2, function(arr) skewness(c(arr))),
            apply(vr, 1:3, function(arr) skewness(c(arr))), along = 1)

matplot(t(sk[,,"temp.n"]), type = "l")
abline(h = 0, col = "cyan3")

pit.skew <- apply(pit[,,-1,], 1:2, function(arr) skewness(c(arr)))
lines(pit.skew[,"temp.n"], col = "magenta3")

#===============================================================================================

# --> Brier score for probability of freezing                                               ####
{
    fp.mim <- aaply(mimic[,,,1:2,1:3], 3, function(mim) {
        abind("temp.n" = pnorm(0, mim[,,"temp.n", "Tau"], sqrt(mim[,,"temp.n", "temp.n"])),
              "temp.s" = pnorm(0, mim[,,"temp.s", "Tau"], sqrt(mim[,,"temp.s", "temp.s"])),
              along = 0)
    })
    
    fp.mos <- apply(abind("mean" = mos.pred, "var" = mos.var, along = 0), 2:5,
                    function(mos) pnorm(0, mos["mean"], sqrt(mos["var"])))
    
    fp.bma <- apply(sweep(aaply(bma.mu, 1:2, 
                                function(arr) pnorm(0, mean = arr, sd = bma.fitted[,,,7])),
                          c(1,3:5), aperm(bma.fitted[,,,8:10], c(4,1:3)), "*"),
                    2:5, sum)
    
    f.prob <- abind("mimic" = fp.mim,
                    "bma" = aperm(fp.bma, c(4,1:3)),
                    "mos" = aperm(fp.mos, c(4,1:3)),
                    along = 0)
    
    brier <- apply(sweep(f.prob, 3:5, obs[1:2,,] <= 0, "-")^2, 1:3, mean, na.rm = T)
}


pdf("../Plots/Analogue-perf/Brier-score.pdf", height = 4, width = 9); {
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(dimnames(brier)[[3]], function(varb) {
        matplot(0:14, t(brier[,,varb]), type = "l", col = c("black", "blue", "red3"), lwd = c(2,1,1),
                main = varb, xlab = "", ylab = "", ylim = range(0, brier))
        legend("bottomright", legend = dimnames(brier)[[1]], col = c("black", "blue", "red3"), 
               lwd = c(2,1,1), bty = "n")
    }))
    mtext("Brier score across all leadtimes - probability of freezing", outer = T)
}; dev.off()

####################################################################################################

# --> Brier scores for probabilities of other temps                                         ####

# (-5, -2.5, 0, 2.5, 5, 7.5, 10)
get.brier <- function(temp = 0) {
    fp.mim <- aaply(mimic[,,,1:2,1:3], 3, function(mim) {
        abind("temp.n" = pnorm(temp, mim[,,"temp.n", "Tau"], sqrt(mim[,,"temp.n", "temp.n"])),
              "temp.s" = pnorm(temp, mim[,,"temp.s", "Tau"], sqrt(mim[,,"temp.s", "temp.s"])),
              along = 0)
    })
    
    fp.mos <- apply(abind("mean" = mos.pred, "var" = mos.var, along = 0), 2:5,
                    function(mos) pnorm(temp, mos["mean"], sqrt(mos["var"])))
    
    fp.bma <- apply(sweep(aaply(bma.mu, 1:2, 
                                function(arr) pnorm(temp, mean = arr, sd = bma.fitted[,,,7])),
                          c(1,3:5), aperm(bma.fitted[,,,8:10], c(4,1:3)), "*"),
                    2:5, sum)
    

    
    
    f.prob <- abind("mimic" = fp.mim,
                    "bma" = aperm(fp.bma, c(4,1:3)),
                    "mos" = aperm(fp.mos, c(4,1:3)),
                    along = 0)
    
    fp.hist <- apply(sweep(aaply(mimic.hist[,,1:2,1:3], 1, function(mim) {
        abind("temp.n" = pnorm(temp, mim[,"temp.n", "Tau"], sqrt(mim[,"temp.n", "temp.n"])),
              "temp.s" = pnorm(temp, mim[,"temp.s", "Tau"], sqrt(mim[,"temp.s", "temp.s"])),
              along = 0)}), 2:3, o.605 <= temp, "-")^2, 1:2, mean, na.rm = T)
    
    abind(apply(sweep(f.prob, 3:5, obs[1:2,,] <= temp, "-")^2, 1:3, mean, na.rm = T),
          "mimic.hist" = fp.hist, along = 1)
}

bb <- get.brier()

brier.scores <- abind(invisible(sapply(c(-2:4) * 2.5, function(t) {
    get.brier(temp = t)
}, simplify = F)), along = 0)

dimnames(brier.scores)[[1]] <- c(-2:4) * 2.5

pdf("../Plots/Analogue-perf/Brier-scores-all-temps.pdf", height = 4, width = 9); {
    
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(dimnames(brier.scores)[[1]], function(tmp) {
        invisible(sapply(dimnames(brier.scores)[[4]], function(varb) {
        matplot(0:14, t(brier.scores[tmp,,,varb]), type = "l", col = l.cols, lwd = l.wd, lty = l.ty,
                main = varb, xlab = "", ylab = "", ylim = range(0, brier.scores))
        abline(h = 0, col = "cyan3")
        legend("topleft", col = l.cols, lwd = l.wd, lty = l.ty,
               legend = dimnames(brier.scores)[[2]], bty = "n")
        }))
        mtext(paste0("Brier scores when predicting temp of ", tmp), outer = T)
    }))
    
}; dev.off()

#===============================================================================================

####################################################################################################

# SINGLE COMBINED PLOT OF ALL METRICS                                                           ####

# calculate scores for historically-trained mimic, for comparison
{
    mimic.hist <- readRDS("../Models/Mimic-hist.rds")
    o.605 <- t(apply(obs[1:2,,], 1, c)[26:630,])

    # RMSE
    {
        err <- abind("mimic" = sweep(aperm(mimic[,,,1:2,1], c(4, 1:3)), 1:3, obs[1:2,,], "-"),
                     sweep(aperm(q.fc, c(4:5,1:3)), 2:4, obs[1:2,,], "-"), 
                     along = 1)
        
        rmse <- abind(sqrt(apply(err^2, c(1,2,5), mean, na.rm = T)), 
                      "mimic.hist" = sqrt(apply(sweep(aperm(mimic.hist[,,1:2,1], c(3:1)),
                                                      1:2, o.605, "-")^2, c(1,3), mean)),
                      along = 1)
    }

    # Brier score (freezing)
    {
        brier.h <- apply(sweep(aaply(mimic.hist[,,1:2,1:3], 1, function(mim) {
            abind("temp.n" = pnorm(0, mim[,"temp.n", "Tau"], sqrt(mim[,"temp.n", "temp.n"])),
                  "temp.s" = pnorm(0, mim[,"temp.s", "Tau"], sqrt(mim[,"temp.s", "temp.s"])),
                  along = 0)}), 2:3, o.605 <= 0, "-")^2, 1:2, mean, na.rm = T)
        brier <- abind(brier, "mimic.hist" = brier.h, along = 1)
    }

    # deviation from uniformity
    {
        hist.sim <- aaply(mimic.hist[,,1:2,1:3], 1:2, function(mim) rmvnorm(1000, mean = mim[,"Tau"], sigma = mim[,-1]))
        
        vr.mim <- aaply(aperm(mim.sim, c(5,1:4)), 4, function(ens) {
            apply(abind(obs[1:2,-1,], ens), 1:3, function(v) which(sort(v) == v[1]))})
        
        vr.hist <- aaply(aperm(hist.sim, c(4,2,1,3)), 3, function(ens) {
            apply(abind(o.605, ens), 1:2, function(v) which(sort(v) == v[1]))})
        dimnames(vr.hist)[[2]] <- dimnames(obs)[[1]][1:2]
        
        ud.hist <- apply(vr.hist, 1:2, u.dev, l = 1001)
        ud <- abind(ud, "mimic.hist" = ud.hist, along = 1)
    }
    
    # skewness of verification ranks
    {
        sk <- abind(sk, "mimic.hist" = apply(vr.hist, 1:2, skewness), along = 1)
    }
}

l.cols <- c("black", "blue", "red3", "grey"); l.wd <- c(2,1,1,1); l.ty <- c(1,1,1,2)

pdf("../Plots/Analogue-perf/Mimic-performance.pdf", height = 4, width = 9); {
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    # Energy score (univariate -> CRPS)
    {
        invisible(sapply(dimnames(crps.es)[[3]], function(varb) {
            matplot(0:14, t(crps.es[,,varb]), type = "l", col = l.cols, lwd = l.wd, lty = l.ty,
                    main = varb, xlab = "", ylab = "", ylim = range(0, crps.es))
            abline(h = 0, col = "cyan3")
            legend("bottomright", col = l.cols, lwd = l.wd, lty = l.ty,
                   legend = dimnames(brier)[[1]], bty = "n")
        }))
        mtext("CRPS (energy score) from simulations", outer = T)
    }
    
    # RMSE
    {
        invisible(sapply(dimnames(brier)[[3]], function(varb) {
            matplot(0:14, t(rmse[,varb,]), type = "l", col = l.cols, lwd = l.wd, lty = l.ty,
                    main = varb, xlab = "", ylab = "", ylim = range(0, rmse))
            abline(h = 0, col = "cyan3")
            legend("bottomright", col = l.cols, lwd = l.wd, lty = l.ty,
                   legend = dimnames(brier)[[1]], bty = "n")
        }))
        mtext("RMSE", outer = T)
    }
    
    # Brier score
    {
        invisible(sapply(dimnames(brier)[[3]], function(varb) {
            matplot(0:14, t(brier[,,varb]), type = "l", col = l.cols, lwd = l.wd, lty = l.ty,
                    main = varb, xlab = "", ylab = "", ylim = range(0, brier))
            abline(h = 0, col = "cyan3")
            legend("bottomright", col = l.cols, lwd = l.wd, lty = l.ty,
                   legend = dimnames(brier)[[1]], bty = "n")
        }))
        mtext("Brier score - probability of freezing", outer = T)
    }

    # skewness of verification ranks
    {
        invisible(sapply(dimnames(sk)[[3]], function(varb) {
            matplot(0:14, t(sk[,,varb]), type = "l", xlab = "", ylab = "", main = varb,
                    col = l.cols, lwd = l.wd, lty = l.ty, ylim = range(sk, 0))
            abline(h = 0, col = "cyan3")
            legend("bottomright", col = l.cols, lwd = l.wd, lty = l.ty,
                   legend = dimnames(sk)[[1]], bty = "n")
        }))
        mtext("Skewness of verification ranks", outer = T)
    }
    
    # deviation from uniformity
    {
        invisible(sapply(dimnames(ud)[[3]], function(varb) {
            matplot(0:14, t(ud[,,varb]), type = "l", ylim = range(ud, 0, na.rm = T),
                    xlab = "", ylab = "", main = varb, col = l.cols, lwd = l.wd, lty = l.ty)
            legend("bottomright", bty = "n", col = l.cols, lwd = l.wd, lty = l.ty, legend = dimnames(ud)[[1]])
            abline(h = 0, col = "cyan3")
        }))
        mtext("Deviation from uniformity", outer = T)
    }
    
}; dev.off()

# CRPS, Brier-other?


####################################################################################################

# SINGLE-ENSEMBLE BAYESIAN POSTERIORS                                                           ####

# check plateau point for each ensemble - where does model error take over?

ens.training.set <- function(model) {
    
    # get array of 25 analogues per day/year/leadtime (~75s)
    an.indx <- abind(lapply(1:15, function(lt) {
        
        # search space - all candidates at that leadtime
        cand <- abind("y.o" = obs[,1:89,],
                      "y.fc" = apply(offset.forecast(model)[,1:89,,lt,-1], 1:3, mean),
                      "c.fc" = apply(offset.forecast(model)[,2:90,,lt,-1], 1:3, mean),
                      along = 0)
        cand <- abind(array(NA, dim = dim(cand[,,1,])), cand, along = 3)
        
        # indices of 25 closest analogues
        an <- aaply(cand, 3:4, function(targ) {
            
            if (sum(!is.na(targ)) == 0) {
                # if target is padding (as when running over array), return an NA array of expected size
                return(array(NA, dim = c(25, 2)))
            }
            
            # standard deviation for normalisation
            cand.sd <- apply(cand, 1:2, sd, na.rm = T)
            
            # distance from each point on target to equivalent in candidate space
            dist <- sweep(sqrt(sweep(cand, 1:2, targ, "-")^2), 1:2, cand.sd, "/")
            
            # mean distance across all variables
            m.dist <- apply(dist, 3:4, mean)
            
            # identify n analogues
            an <- which(m.dist <= sort(m.dist)[26], arr.ind = T)
            
            # find target vector and remove from analogue list
            t <- which(m.dist == 0, arr.ind = T)
            an <- an[!(an[,1] == t[1] & an[,2] == t[2]),]
        })
    }), along = 2.5)
    
    # extract analogue forecast errors
    em <- apply(forecast.errors(ecmwf)[,,,,-1], 1:4, mean)
    abind(sapply(1:15, function(lt) {
        aaply(an.indx[,,lt,,], 1:2, function(a) {
            apply(em[,,,lt], 1, "[", a)
        })}, simplify = F), along = 2.5)
}

tr.ec <- ens.training.set(ecmwf)
tr.nc <- ens.training.set(ncep)
tr.mo <- ens.training.set(ukmo)

o <- obs[1:2,5,5]
ef <- offset.forecast(ecmwf)[1:2,5,5,1,-1]
tr <- t(tr.ec[5,5,1,,1:2])

# use analogues to find eta & lambda
ens.mimic <- function(o, ef, tr) {
    
    Y.bar <- apply(ef, 1, mean)
    C <- cov(t(ef))
    
    Eta <- apply(tr, 1, mean)
    Lambda <- cov(t(tr))
    
    S <- Lambda + C
    Tau <- S %*% (solve(diag(2) + solve(C) %*% Lambda) %*% (solve(C) %*% Y.bar)) - Eta
    
    abind(Tau = Tau, S, along = 2)
}

# run over each model set
ec.mimic <- abind(invisible(sapply(1:15, function(lt) {
    suppressWarnings(aaply(abind(obs[1:2,,],
                                 aperm(tr.ec[,,lt,,1:2], c(4,1:3)),
                                 offset.forecast(ecmwf)[1:2,,,lt,-1], 
                                 along = 4), 2:3, function(arr) {
                                     ens.mimic(o = arr[,1],
                                               ef = arr[,-(1:26)],
                                               tr = arr[,2:26])}))
}, simplify = F)), along = 2.5)

nc.mimic <- abind(invisible(sapply(1:15, function(lt) {
    suppressWarnings(aaply(abind(obs[1:2,,],
                                 aperm(tr.nc[,,lt,,1:2], c(4,1:3)),
                                 offset.forecast(ncep)[1:2,,,lt,-1], 
                                 along = 4), 2:3, function(arr) {
                                     ens.mimic(o = arr[,1],
                                               ef = arr[,-(1:26)],
                                               tr = arr[,2:26])}))
}, simplify = F)), along = 2.5)

mo.mimic <- abind(invisible(sapply(1:15, function(lt) {
    suppressWarnings(aaply(abind(obs[1:2,,],
                                 aperm(tr.mo[,,lt,,1:2], c(4,1:3)),
                                 offset.forecast(ukmo)[1:2,,,lt,-1], 
                                 along = 4), 2:3, function(arr) {
                                     ens.mimic(o = arr[,1],
                                               ef = arr[,-(1:26)],
                                               tr = arr[,2:26])}))
}, simplify = F)), along = 2.5)


# now calculate & compare RMSE (for starters)

rmse <- abind("mimic" = sqrt(apply(sweep(mimic[,,,1:2,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ecmwf" = sqrt(apply(sweep(ec.mimic[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ncep" = sqrt(apply(sweep(nc.mimic[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ukmo" = sqrt(apply(sweep(mo.mimic[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              along = 0)

l.cols <- c("black", "blue", "red3", "green3"); l.wd <- c(2,1,1,1); l.ty <- 1

pdf("../Plots/Analogue-perf/Individual-ensemble-mimic-RMSE.pdf", height = 4, width = 9); {
   
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(dimnames(rmse)[[3]], function(varb) {
        matplot(0:14, t(rmse[,,varb]), type = "l", col = l.cols, lwd = l.wd, lty = l.ty,
                main = varb, xlab = "", ylab = "", ylim = range(0, rmse))
        abline(h = 0, col = "cyan3")
        abline(v = 7, col = "gold")
        legend("bottomright", col = l.cols, lwd = l.wd, lty = l.ty,
               legend = dimnames(rmse)[[1]], bty = "n")
    }))
    mtext("RMSE - MME vs individual ensembles", outer = T)
    
}; dev.off()



####################################################################################################

# SINGLE-ENSEMBLE POSTERIORS TRAINED ON MME ANALOGUES                                           ####

# create training set as per first section above
tr.errors <- sweep(tr.fc, c(1:4, 6), tr.obs, "-")

ens.mimic <- function(o, ef, tr) {
    
    Y.bar <- apply(ef, 1, mean)
    C <- cov(t(ef))
    
    Eta <- apply(tr, 1, mean)
    Lambda <- cov(t(tr))
    
    S <- Lambda + C
    Tau <- S %*% (solve(diag(2) + solve(C) %*% Lambda) %*% (solve(C) %*% Y.bar))
    
    abind(Tau = Tau, S, along = 2)
}

ec.mimic.mma <- abind(invisible(sapply(1:15, function(lt) {
    suppressWarnings(aaply(abind(obs[1:2,,],
                                 aperm(tr.errors[,,lt,,"ecmwf",1:2], c(4,1:3)),
                                 offset.forecast(ecmwf)[1:2,,,lt,-1], 
                                 along = 4), 2:3, function(arr) {
                                     ens.mimic(o = arr[,1],
                                               ef = arr[,-(1:26)],
                                               tr = arr[,2:26])}))
}, simplify = F)), along = 2.5)

nc.mimic.mma <- abind(invisible(sapply(1:15, function(lt) {
    suppressWarnings(aaply(abind(obs[1:2,,],
                                 aperm(tr.errors[,,lt,,"ncep",1:2], c(4,1:3)),
                                 offset.forecast(ncep)[1:2,,,lt,-1], 
                                 along = 4), 2:3, function(arr) {
                                     ens.mimic(o = arr[,1],
                                               ef = arr[,-(1:26)],
                                               tr = arr[,2:26])}))
}, simplify = F)), along = 2.5)

mo.mimic.mma <- abind(invisible(sapply(1:15, function(lt) {
    suppressWarnings(aaply(abind(obs[1:2,,],
                                 aperm(tr.errors[,,lt,,"ukmo",1:2], c(4,1:3)),
                                 offset.forecast(ukmo)[1:2,,,lt,-1], 
                                 along = 4), 2:3, function(arr) {
                                     ens.mimic(o = arr[,1],
                                               ef = arr[,-(1:26)],
                                               tr = arr[,2:26])}))
}, simplify = F)), along = 2.5)


rmse <- abind("mimic" = sqrt(apply(sweep(mimic[,,,1:2,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ecmwf" = sqrt(apply(sweep(ec.mimic[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ncep" = sqrt(apply(sweep(nc.mimic[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ukmo" = sqrt(apply(sweep(mo.mimic[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ecmwf.mma" = sqrt(apply(sweep(ec.mimic.mma[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ncep.mma" = sqrt(apply(sweep(nc.mimic.mma[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ukmo.mma" = sqrt(apply(sweep(mo.mimic.mma[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              along = 0)

l.cols <- c("black", rep(c("blue", "red3", "green3"), 2)); l.wd <- c(2,rep(1,6)); l.ty <- c(rep(1,4), rep(2, 3))

invisible(sapply(dimnames(rmse)[[3]], function(varb) {
    matplot(0:14, t(rmse[,,varb]), type = "l", col = l.cols, lwd = l.wd, lty = l.ty,
            main = varb, xlab = "", ylab = "", ylim = range(0, rmse))
    abline(h = 0, col = "cyan3")
    abline(v = 7, col = "gold")
    legend("bottomright", col = l.cols, lwd = l.wd, lty = l.ty,
           legend = dimnames(rmse)[[1]], bty = "n")
}))
mtext("RMSE - MME vs individual ensembles", outer = T)

####################################################################################################

# CHECK ANALOGUE SELECTION FOR SINGLE-ENSEMBLE FIT VS ALL MODELS                                ####

lt <- 1

cand.mma <- abind("y.o" = obs[,1:89,],
              "y.ecmwf" = apply(offset.forecast(ecmwf)[,1:89,,lt,-1], 1:3, mean),
              "y.ncep" = apply(offset.forecast(ncep)[,1:89,,lt,-1], 1:3, mean),
              "y.ukmo" = apply(offset.forecast(ukmo)[,1:89,,lt,-1], 1:3, mean),
              "c.ecmwf" = apply(offset.forecast(ecmwf)[,2:90,,lt,-1], 1:3, mean),
              "c.ncep" = apply(offset.forecast(ncep)[,2:90,,lt,-1], 1:3, mean),
              "c.ukmo" = apply(offset.forecast(ukmo)[,2:90,,lt,-1], 1:3, mean),
              along = 0)
cand.mma <- abind(array(NA, dim = dim(cand.mma[,,1,])), cand.mma, along = 3)
cand.sd.mma <- apply(cand.mma, 1:2, sd, na.rm = T)
cand.sd.ec <- apply(cand.ec, 1:2, sd, na.rm = T)

cand.ec <- abind("y.o" = obs[,1:89,],
                  "y.ecmwf" = apply(offset.forecast(ecmwf)[,1:89,,lt,-1], 1:3, mean),
                  "c.ecmwf" = apply(offset.forecast(ecmwf)[,2:90,,lt,-1], 1:3, mean),
                  along = 0)
cand.ec <- abind(array(NA, dim = dim(cand.ec[,,1,])), cand.ec, along = 3)

d <- 5; y <- 5
targ.mma <- cand.mma[,,d,y]
targ.ec <- cand.ec[,,d,y]

# calculate & normalise distances
dist.mma <- sweep(sqrt(sweep(cand.mma, 1:2, targ.mma, "-")^2), 1:2, cand.sd.mma, "/")
dist.ec <- sweep(sqrt(sweep(cand.ec, 1:2, targ.ec, "-")^2), 1:2, cand.sd.ec, "/")

m.dist.mma <- apply(dist.mma, 3:4, mean)
m.dist.ec <- apply(dist.ec, 3:4, mean)

an.mma <- which(m.dist.mma <= sort(m.dist.mma)[26], arr.ind = T)
# find target vector and remove from analogue list
t <- which(m.dist.mma == 0, arr.ind = T)
an.mma <- an.mma[!(an.mma[,1] == t[1] & an.mma[,2] == t[2]),]

an.ec <- which(m.dist.ec <= sort(m.dist.ec)[26], arr.ind = T)
t <- which(m.dist.ec == 0, arr.ind = T)
an.ec <- an.ec[!(an.ec[,1] == t[1] & an.ec[,2] == t[2]),]

# use analogue indices to identify training set

fc.errors <- abind("ecmwf" = apply(offset.forecast(ecmwf)[,,,,-1], 1:4, mean),
                   "ncep" = apply(offset.forecast(ncep)[,,,,-1], 1:4, mean),
                   "ukmo" = apply(offset.forecast(ukmo)[,,,,-1], 1:4, mean),
                   along = 0)

tr.mma <- suppressWarnings(aaply(an.mma, 1, function(a) apply(fc.errors[,,,,1], 1:2, "[", t(a))))
tr.ec <- suppressWarnings(aaply(an.ec, 1, function(a) apply(fc.errors[,,,,1], 1:2, "[", t(a))))

mimic.mma <- fit.mimic(o = obs[,5,5],
                       fc = list("ecmwf" = offset.forecast(ecmwf)[,5,5,1,-1],
                                 "ncep" = offset.forecast(ncep)[,5,5,1,-1],
                                 "ukmo" = offset.forecast(ukmo)[,5,5,1,-1]),
                       tr = t(apply(tr.mma, c(1,3), mean)))

mimic.ec <- ens.mimic(o = obs[1:2,5,5],
                      ef = offset.forecast(ecmwf)[1:2,5,5,1,-1],
                      tr = t(tr.ec[,"ecmwf",1:2]))

mimic.ec
mimic.mma

points(mimic.ec[1:2,1], col = "red", pch = 4)
mimic.mma[1:2, 1]
obs[1:2,5,5]

####################################################################################################

# MME VS SINGLE-ENSEMBLE MIMICS                                                                 ####

# find analogue indices based on mma & individual ensembles
analogue.indices <- function(lt = 1, n = 25) {
    
    cand <- abind("y.o" = obs[,1:89,],
                  "y.ecmwf" = apply(offset.forecast(ecmwf)[,1:89,,lt,-1], 1:3, mean),
                  "c.ecmwf" = apply(offset.forecast(ecmwf)[,2:90,,lt,-1], 1:3, mean),
                  "y.ncep" = apply(offset.forecast(ncep)[,1:89,,lt,-1], 1:3, mean),
                  "c.ncep" = apply(offset.forecast(ncep)[,2:90,,lt,-1], 1:3, mean),
                  "y.ukmo" = apply(offset.forecast(ukmo)[,1:89,,lt,-1], 1:3, mean),
                  "c.ukmo" = apply(offset.forecast(ukmo)[,2:90,,lt,-1], 1:3, mean),
                  along = 0)
    cand <- abind(array(NA, dim = dim(cand[,,1,])), cand, along = 3)
    
    # standard deviation for normalisation
    cand.sd <- apply(cand, 1:2, sd, na.rm = T)
    
    # normalised per-model and per-variable distances between target & other candidates
    dist <- aaply(cand, 3:4, function(targ) {
        sweep(sqrt(sweep(cand, 1:2, targ, "-")^2), 1:2, cand.sd, "/")
    })
    
    # analogues for multi-model mimic
    {
        # calculate mean distances over each subset
        m.dist.mma <- apply(dist, c(1:2, 5:6), mean)
        
        # identify closest analogues
        an.mma <- suppressWarnings(aaply(m.dist.mma[-1,,,], 1:2, 
                                         function(md) which(md <= sort(md)[n+1] & md > 0, arr.ind = T)))
        an.mma <- abind(array(NA, dim = dim(an.mma[1,,,])), an.mma, along = 1)
    }

    # analogues for ecmwf
    {
        # calculate mean distances over each subset
        m.dist.ec <- apply(dist[,,1:3,,,], c(1:2, 5:6), mean)
        
        # identify closest analogues
        an.ec <- suppressWarnings(aaply(m.dist.ec[-1,,,], 1:2, 
                                        function(md) which(md <= sort(md)[26] & md > 0, arr.ind = T)))
        an.ec <- abind(array(NA, dim = dim(an.ec[1,,,])), an.ec, along = 1)
    }
    
    # analogues for ncep
    {
        # calculate mean distances over each subset
        m.dist.nc <- apply(dist[,,c(1, 4:5),,,], c(1:2, 5:6), mean)
        
        # identify closest analogues
        an.nc <- suppressWarnings(aaply(m.dist.nc[-1,,,], 1:2, 
                                        function(md) which(md <= sort(md)[26] & md > 0, arr.ind = T)))
        an.nc <- abind(array(NA, dim = dim(an.nc[1,,,])), an.nc, along = 1)
    }
    
    # analogues for ukmo
    {
        # calculate mean distances over each subset
        m.dist.mo <- apply(dist[,,c(1, 6:7),,,], c(1:2, 5:6), mean)
        
        # identify closest analogues
        an.mo <- suppressWarnings(aaply(m.dist.mo[-1,,,], 1:2, 
                                        function(md) which(md <= sort(md)[26] & md > 0, arr.ind = T)))
        an.mo <- abind(array(NA, dim = dim(an.mo[1,,,])), an.mo, along = 1)
    }
    
    abind("mimic" = an.mma,
          "ecmwf" = an.ec,
          "ncep" = an.nc,
          "ukmo" = an.mo,
          along = 0)
}

# find analogues at all leadtimes (~5 minutes) 
an.indx <- abind(invisible(sapply(1:15, analogue.indices, simplify = F)), along = 3.5)

saveRDS(an.indx, "../Models/Analogue-indices.rds")

# create training sets
{
    fc.errors <- abind("ecmwf" = apply(forecast.errors(ecmwf)[,,,,-1], 1:4, mean),
                       "ncep" = apply(forecast.errors(ncep)[,,,,-1], 1:4, mean),
                       "ukmo" = apply(forecast.errors(ukmo)[,,,,-1], 1:4, mean),
                       along = 0)
    
    tr.mma <- abind(sapply(1:15, function(lt) {
        suppressWarnings(aaply(an.indx["mimic",,,lt,,], 1:2, function(a) {
            apply(fc.errors[,,,,lt], 1:2, "[", a)
        }))
    }, simplify = F), along = 2.5)
    
    tr.ec <- abind(sapply(1:15, function(lt) {
        suppressWarnings(aaply(an.indx["ecmwf",,,lt,,], 1:2, function(a) {
            apply(fc.errors["ecmwf",,,,lt], 1, "[", a)
        }))
    }, simplify = F), along = 2.5)
    
    tr.nc <- abind(sapply(1:15, function(lt) {
        suppressWarnings(aaply(an.indx["ncep",,,lt,,], 1:2, function(a) {
            apply(fc.errors["ncep",,,,lt], 1, "[", a)
        }))
    }, simplify = F), along = 2.5)
    
    tr.mo <- abind(sapply(1:15, function(lt) {
        suppressWarnings(aaply(an.indx["ukmo",,,lt,,], 1:2, function(a) {
            apply(fc.errors["ukmo",,,,lt], 1, "[", a)
        }))
    }, simplify = F), along = 2.5)
}

# now, to fit the mimics
{
    ens.mimic <- function(o, ef, tr) {
        
        Y.bar <- apply(ef, 1, mean)
        C <- cov(t(ef))
        
        Eta <- apply(tr, 1, mean)
        Lambda <- cov(t(tr))
        
        S <- Lambda + C
        Tau <- S %*% (solve(diag(2) + solve(C) %*% Lambda) %*% (solve(C) %*% Y.bar)) - Eta
        
        abind(Tau = Tau, S, along = 2)
    }
    
    ec.mimic <- abind(invisible(sapply(1:15, function(lt) {
        suppressWarnings(aaply(abind(obs[1:2,,],
                                     aperm(tr.ec[,,lt,,1:2], c(4,1:3)),
                                     offset.forecast(ecmwf)[1:2,,,lt,-1], 
                                     along = 4), 2:3, function(arr) {
                                         ens.mimic(o = arr[,1],
                                                   ef = arr[,-(1:26)],
                                                   tr = arr[,2:26])}))
    }, simplify = F)), along = 2.5)
    
    nc.mimic <- abind(invisible(sapply(1:15, function(lt) {
        suppressWarnings(aaply(abind(obs[1:2,,],
                                     aperm(tr.nc[,,lt,,1:2], c(4,1:3)),
                                     offset.forecast(ncep)[1:2,,,lt,-1], 
                                     along = 4), 2:3, function(arr) {
                                         ens.mimic(o = arr[,1],
                                                   ef = arr[,-(1:26)],
                                                   tr = arr[,2:26])}))
    }, simplify = F)), along = 2.5)
    
    mo.mimic <- abind(invisible(sapply(1:15, function(lt) {
        suppressWarnings(aaply(abind(obs[1:2,,],
                                     aperm(tr.mo[,,lt,,1:2], c(4,1:3)),
                                     offset.forecast(ukmo)[1:2,,,lt,-1], 
                                     along = 4), 2:3, function(arr) {
                                         ens.mimic(o = arr[,1],
                                                   ef = arr[,-(1:26)],
                                                   tr = arr[,2:26])}))
    }, simplify = F)), along = 2.5)
    
    mma.mimic <- abind(sapply(1:15, function(lt) {
        src <- abind("o" = obs,
                     "ecmwf" = offset.forecast(ecmwf)[,,,lt,-1],
                     "ncep" = offset.forecast(ncep)[,,,lt,-1],
                     "ukmo" = offset.forecast(ukmo)[,,,lt,-1],
                     "tr" = aperm(apply(tr.mma[,,lt,,,], c(1:3, 5), mean), c(4,1:3)),
                     along = 4)
        
        suppressWarnings(aaply(src, 2:3, function(arr) {
            fit.mimic(o = arr[,1],
                      fc = list("ecmwf" = arr[,2:51], "ncep" = arr[,52:71], "ukmo" = arr[,72:94]),
                      tr = arr[,95:119])
        }))
    }, simplify = F), along = 2.5)
    
    ec.mimic.mma <- abind(invisible(sapply(1:15, function(lt) {
        suppressWarnings(aaply(abind(obs[1:2,,],
                                     aperm(tr.mma[,,lt,,"ecmwf",1:2], c(4,1:3)),
                                     offset.forecast(ecmwf)[1:2,,,lt,-1], 
                                     along = 4), 2:3, function(arr) {
                                         ens.mimic(o = arr[,1],
                                                   ef = arr[,-(1:26)],
                                                   tr = arr[,2:26])}))
    }, simplify = F)), along = 2.5)
    
    nc.mimic.mma <- abind(invisible(sapply(1:15, function(lt) {
        suppressWarnings(aaply(abind(obs[1:2,,],
                                     aperm(tr.mma[,,lt,,"ncep",1:2], c(4,1:3)),
                                     offset.forecast(ncep)[1:2,,,lt,-1], 
                                     along = 4), 2:3, function(arr) {
                                         ens.mimic(o = arr[,1],
                                                   ef = arr[,-(1:26)],
                                                   tr = arr[,2:26])}))
    }, simplify = F)), along = 2.5)
    
    mo.mimic.mma <- abind(invisible(sapply(1:15, function(lt) {
        suppressWarnings(aaply(abind(obs[1:2,,],
                                     aperm(tr.mma[,,lt,,"ukmo",1:2], c(4,1:3)),
                                     offset.forecast(ukmo)[1:2,,,lt,-1], 
                                     along = 4), 2:3, function(arr) {
                                         ens.mimic(o = arr[,1],
                                                   ef = arr[,-(1:26)],
                                                   tr = arr[,2:26])}))
    }, simplify = F)), along = 2.5)
}

# calculate RMSE for each (simplest metric, may add to this later)

mimic.hist <- readRDS("../Models/Mimic-hist.rds")

mh.err <- sweep(mimic.hist[,,1:2,1], 2:3, apply(obs, 1, c)[26:630,1:2], "-")

rmse <- abind("mimic" = sqrt(apply(sweep(mma.mimic[,,,1:2,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ecmwf.mma" = sqrt(apply(sweep(ec.mimic.mma[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ncep.mma" = sqrt(apply(sweep(nc.mimic.mma[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ukmo.mma" = sqrt(apply(sweep(mo.mimic.mma[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ecmwf" = sqrt(apply(sweep(ec.mimic[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ncep" = sqrt(apply(sweep(nc.mimic[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "ukmo" = sqrt(apply(sweep(mo.mimic[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 3:4, mean, na.rm = T)),
              "mimic.hist" = sqrt(apply(mh.err^2, c(1, 3), mean, na.rm = T)),
              along = 0)


yr <- 1:3
rmse.1 <- abind("mimic" = sqrt(apply(sweep(mma.mimic[,yr,,1:2,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
              "ecmwf.mma" = sqrt(apply(sweep(ec.mimic.mma[,yr,,,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
              "ncep.mma" = sqrt(apply(sweep(nc.mimic.mma[,yr,,,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
              "ukmo.mma" = sqrt(apply(sweep(mo.mimic.mma[,yr,,,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
              "ecmwf" = sqrt(apply(sweep(ec.mimic[,yr,,,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
              "ncep" = sqrt(apply(sweep(nc.mimic[,yr,,,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
              "ukmo" = sqrt(apply(sweep(mo.mimic[,yr,,,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
              "mimic.hist" = sqrt(apply(mh.err[,1:245,]^2, c(1, 3), mean, na.rm = T)),
              along = 0)

yr <- 4:7
rmse.2 <- abind("mimic" = sqrt(apply(sweep(mma.mimic[,yr,,1:2,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
                "ecmwf.mma" = sqrt(apply(sweep(ec.mimic.mma[,yr,,,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
                "ncep.mma" = sqrt(apply(sweep(nc.mimic.mma[,yr,,,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
                "ukmo.mma" = sqrt(apply(sweep(mo.mimic.mma[,yr,,,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
                "ecmwf" = sqrt(apply(sweep(ec.mimic[,yr,,,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
                "ncep" = sqrt(apply(sweep(nc.mimic[,yr,,,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
                "ukmo" = sqrt(apply(sweep(mo.mimic[,yr,,,1], c(4,1,2), obs[1:2,,yr], "-")^2, 3:4, mean, na.rm = T)),
                "mimic.hist" = sqrt(apply(mh.err[,246:605,]^2, c(1, 3), mean, na.rm = T)),
                along = 0)

mh <- sweep(mimic.hist[,hr,1:2,1], 2:3, apply(obs,1,c)[hr,1:2], "-")^2
mh <- sqrt(apply(sweep(mimic.hist[,hr,1:2,1], 2:3,
                       apply(obs,1,c)[hr,1:2], "-")^2, c(1, 3), mean, na.rm = T))

l.cols <- c("black", rep(c("red3", "blue", "green3"), 2), "black")
l.ty <- c(rep(1, 4), rep(2, 3), 3); l.wd <- c(2, rep(1,6), 2) 

pdf("../Plots/Analogue-perf/Individual-ensemble-mimic-RMSE.pdf", height = 4, width = 9); {
    
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(dimnames(rmse)[[3]], function(varb) {
        matplot(0:14, t(rmse[,,varb]), type = "l", col = l.cols, lty = l.ty, lwd = l.wd, 
                xlab = "", ylab = "", ylim = range(0,rmse), main = varb)
        legend("bottomright", col = l.cols, lty = l.ty, lwd = l.wd, bty = "n",
               legend = dimnames(rmse)[[1]], cex = 0.7)
        abline(v = 7, col = "cyan3", lty = 2)
    }))
    mtext("RMSE for invividual ensembles vs MME (all years)", outer = T)
    
    invisible(sapply(dimnames(rmse.1)[[3]], function(varb) {
        matplot(0:14, t(rmse.1[,,varb]), type = "l", col = l.cols, lty = l.ty, lwd = l.wd, 
                xlab = "", ylab = "", ylim = range(0, rmse.1), main = varb)
        legend("bottomright", col = l.cols, lty = l.ty, lwd = l.wd, bty = "n",
               legend = dimnames(rmse.1)[[1]], cex = 0.7)
        abline(v = 7, col = "cyan3", lty = 2)
    }))
    mtext("RMSE for invividual ensembles vs MME (Dec 2007 - Feb 2010)", outer = T)
    
    
    invisible(sapply(dimnames(rmse.2)[[3]], function(varb) {
        matplot(0:14, t(rmse.2[,,varb]), type = "l", col = l.cols, lty = l.ty, lwd = l.wd, 
                xlab = "", ylab = "", ylim = range(0, rmse.2), main = varb)
        legend("bottomright", col = l.cols, lty = l.ty, lwd = l.wd, bty = "n",
               legend = dimnames(rmse.2)[[1]], cex = 0.7)
        abline(v = 7, col = "cyan3", lty = 2)
    }))
    mtext("RMSE for invividual ensembles vs MME (Dec 2010 - Feb 2014)", outer = T)
    
}; dev.off()

ncep.err <- sqrt(apply(sweep(nc.mimic.mma[,,,,1], c(4,1,2), obs[1:2,,], "-")^2, 2:4, mean, na.rm = T))

matplot(0:14, t(ncep.err[,,1]), type = "l", main = "temp.n", col = c(rep("red3", 3), rep("black", 4)))
abline(h = 0, col = "cyan3")

matplot(0:14, t(ncep.err[,,2]), type = "l", main = "temp.s", col = c(rep("red3", 3), rep("black", 4)))
