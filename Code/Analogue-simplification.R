
library("SX.weather"); library("CB.Misc")

# NEXT STEP: CHECK MARGINAL METRICS & COMPARE TO BMA/MOS RESULTS
# NEED A BOT/PIT EQUIVALENT TO TEST OF DEVIATION FROM UNIFORMITY

# REFACTOR FUNCTIONS (need single-day function for all elements, which can be applied to array)

# need to set up a generic function to create candidate array for any number of forecasts & obs

# could we include control forecasts when calculating Sigma?
# (would double size of estimate of ensemble spread)

####################################################################################################

# FUNCTION TO IDENTIFY ANALOGUES                                                                ####

lt <- "0"

find.analogues <- function(targ, cand, n = 25) {
    
    if (sum(!is.na(targ)) == 0) {
        # if target is padding (as when running over array), return an NA array of expected size
        return(array(NA, dim = c(n, 2)))
    }
    
    # standard deviation for normalisation
    cand.sd <- apply(cand, 1:2, sd, na.rm = T)
    
    # distance from each point on target to equivalent in candidate space
    dist <- sweep(sqrt(sweep(cand, 1:2, targ, "-")^2), 1:2, cand.sd, "/")
    
    # mean distance across all variables
    m.dist <- apply(dist, 3:4, mean)
    
    # identify n analogues
    an <- which(m.dist <= sort(m.dist)[n+1], arr.ind = T)
    
    # find target vector and remove from analogue list
    t <- which(m.dist == 0, arr.ind = T)
    an <- an[!(an[,1] == t[1] & an[,2] == t[2]),]
    
    # return indices of analogues identified
    return(an)
}

# candidate analogues Based on two days' forecast, one day's observations
{
    cand <- abind("y.o" = obs[,1:89,],
                  "y.ecmwf" = apply(offset.forecast(ecmwf)[,1:89,,lt,-1], 1:3, mean),
                  "y.ncep" = apply(offset.forecast(ncep)[,1:89,,lt,-1], 1:3, mean),
                  "y.ukmo" = apply(offset.forecast(ukmo)[,1:89,,lt,-1], 1:3, mean),
                  "c.ecmwf" = apply(offset.forecast(ecmwf)[,2:90,,lt,-1], 1:3, mean),
                  "c.ncep" = apply(offset.forecast(ncep)[,2:90,,lt,-1], 1:3, mean),
                  "c.ukmo" = apply(offset.forecast(ukmo)[,2:90,,lt,-1], 1:3, mean),
                  along = 0)
    cand <- abind(array(NA, dim = dim(cand[,,1,])), cand, along = 3)
}


####################################################################################################

# FIT MODEL WITH ANALOGUE TRAINING SET FOR SINGLE DAY                                           ####

# inputs: obs, ensemble fc, training set for a single day
d <- 5; y <- 5; lt <- "0"

o <- obs[,d,y]
fc <- list("ecmwf" = offset.forecast(ecmwf)[,d,y,lt,-1],
           "ncep" = offset.forecast(ncep)[,d,y,lt,-1], 
           "ukmo" = offset.forecast(ukmo)[,d,y,lt,-1])

tr <- apply(apply(abind("ecmwf" = apply(forecast.errors(ecmwf)[,,,lt,-1], 1:3, mean),
                        "ncep" = apply(forecast.errors(ncep)[,,,lt,-1], 1:3, mean),
                        "ukmo" = apply(forecast.errors(ukmo)[,,,lt,-1], 1:3, mean), along = 0),
                  2:4, mean), 1, "[", find.analogues(cand[,,d,y], cand))

# calculate components of mimic
Y.bar <- abind(lapply(fc, apply, 1, mean), along = 0)
C <- abind(lapply(lapply(fc, t), cov), along = 0)

D.i <- aaply(sweep(sweep(C, 1, sapply(fc, dim)[2,], "/"), 2:3, cov(Y.bar), "+"), 1, solve)

Eta <- apply(tr, 2, mean)
Lambda <- cov(tr)

S <- Lambda + solve(apply(D.i, 2:3, sum))
Tau <- (S %*% 
            ((solve(diag(5) + (apply(D.i, 2:3, sum) %*% Lambda))) %*%
                 apply(aaply(abind(Y.bar, D.i, along = 2), 1, 
                             function(arr) arr[-1,] %*% arr[1,]), 2, sum))) - Eta

plot(rbind(o, t(Tau), apply(Y.bar, 2, mean), Y.bar),
     pch = c(rep(20, 3), rep(4, 3)), col = c("black", "red", "grey", adjustcolor(c("steelblue", "coral3", "green3"))))

####################################################################################################

# RUN AS FUNCTION OVER ARRAY OF INPUT DATA                                                      ####

fc <- list("ecmwf" = offset.forecast(ecmwf)[,,,lt,-1],
           "ncep" = offset.forecast(ncep)[,,,lt,-1], 
           "ukmo" = offset.forecast(ukmo)[,,,lt,-1])

an <- aaply(cand, 3:4, find.analogues, cand)        # around 5s to identify all analogues, single LT

mean.errors <- apply(abind("ecmwf" = apply(forecast.errors(ecmwf)[,,,lt,-1], 1:3, mean),
                           "ncep" = apply(forecast.errors(ncep)[,,,lt,-1], 1:3, mean),
                           "ukmo" = apply(forecast.errors(ukmo)[,,,lt,-1], 1:3, mean), along = 0),
                     2:4, mean)
tr <- aaply(an, 1:2, function(a) aaply(mean.errors, 1, "[", a))


fit.mimic <- function(o, fc, tr, arr.out = T) {
    
        Y.bar <- abind(lapply(fc, apply, 1, mean), along = 0)
        C <- abind(lapply(lapply(fc, t), cov), along = 0)
        
        D.i <- aaply(sweep(sweep(C, 1, sapply(fc, dim)[2,], "/"), 2:3, cov(Y.bar), "+"), 1, solve)
        
        Eta <- apply(t(tr), 2, mean)
        Lambda <- cov(t(tr))
        
        S <- Lambda + solve(apply(D.i, 2:3, sum))
        Tau <- (S %*% 
                    ((solve(diag(5) + (apply(D.i, 2:3, sum) %*% Lambda))) %*%
                         apply(aaply(abind(Y.bar, D.i, along = 2), 1, 
                                     function(arr) arr[-1,] %*% arr[1,]), 2, sum))) - Eta
        
    if (arr.out) {
        # return results in array form (more useful for array operations)
        return(abind("Tau" = Tau, S, along = 2))
    } else {
        # otherwise, return output as list
        return(list(Tau = Tau, S = S))
    }
}

# adjust this to include forecasts & training data at all leadtimes
all.dat <- abind(obs, fc[[1]], fc[[2]], fc[[3]], aperm(tr, c(3,1:2,4)), along = 4)


# ~4s to run over all dates for single leadtime
system.time({
zz <- suppressWarnings(aaply(all.dat, 2:3, function(arr) {
    fit.mimic(o = arr[,1],
              fc = list("ecmwf" = arr[,2:51], "ncep" = arr[,52:71], "ukmo" = arr[,72:94]),
              tr = arr[,95:119])
}))})

plot(zz[,1,1,1], obs[1,,1], pch = 20)

ec.mean <- apply(offset.forecast(ecmwf)[,,,,-1], 1:3, mean)
points(ec.mean[1,,1], obs[1,,1], pch = 20, col = adjustcolor("steelblue", alpha = 0.4))
points(apply(offset.forecast(ncep)[,,,,-1], 1:3, mean)[1,,1], 
       obs[1,,1], pch = 20, col = adjustcolor("coral3", alpha = 0.4))
points(apply(offset.forecast(ukmo)[,,,,-1], 1:3, mean)[1,,1], 
       obs[1,,1], pch = 20, col = adjustcolor("green3", alpha = 0.4))
points(apply(abind(apply(offset.forecast(ecmwf)[,,,,-1], 1:3, mean)[1,,1],
                   apply(offset.forecast(ncep)[,,,,-1], 1:3, mean)[1,,1],
                   apply(offset.forecast(ukmo)[,,,,-1], 1:3, mean)[1,,1],
                   along = 0), 2, mean),
       obs[1,,1], pch = 4, col = "red", lwd = 2)

# verification

v.dat <- abind("o" = aperm(obs, c(2:3,1)), zz, along = 4)

boxed <- apply(v.dat, 1:2, function(arr) box.dot(arr[,"o"], arr[,"Tau"], arr[,-(1:2)]))

hist(boxed, breaks = "fd", prob = T)
abline(h = 1, col = "red", lty = 2)

es <- apply(v.dat[-1,,,], 1:2, function(arr) es.crps(o = arr[,"o"], mu = arr[,"Tau"], sig = arr[,-(1:2)]))
mean(es)    

hist.fitted <- readRDS("../Models/lambda-hist.rds")

dens.mv.crps <- apply(abind("o" = array(t(apply(obs[1:2,,], 1, rbind))[,26:630], dim = c(2,1, 605)),
                            "mu" = hist.fitted$tau[1:2,lt,,drop = F],
                            "sig" = hist.fitted$s[1:2,1:2,lt,,drop = F],
                            along = 1),
                      4, function(arr) {
                          es.crps(o = arr["o",,], mu = arr["mu",,], sig = arr[-(1:2),,], k = 1000)
                      })

####################################################################################################

# AUTOMATE FITTING OVER TRAINING DATA                                                           ####

# currently fixed at 2 days' forecasts, 1 day's observations. Generalise later
get.candidates <- function(lt) {

    cand <- abind("y.o" = obs[,1:89,],
                  "y.ecmwf" = apply(offset.forecast(ecmwf)[,1:89,,lt,-1], 1:3, mean),
                  "y.ncep" = apply(offset.forecast(ncep)[,1:89,,lt,-1], 1:3, mean),
                  "y.ukmo" = apply(offset.forecast(ukmo)[,1:89,,lt,-1], 1:3, mean),
                  "c.ecmwf" = apply(offset.forecast(ecmwf)[,2:90,,lt,-1], 1:3, mean),
                  "c.ncep" = apply(offset.forecast(ncep)[,2:90,,lt,-1], 1:3, mean),
                  "c.ukmo" = apply(offset.forecast(ukmo)[,2:90,,lt,-1], 1:3, mean),
                  along = 0)
    cand <- abind(array(NA, dim = dim(cand[,,1,])), cand, along = 3)
    cand
}

find.analogues <- function(targ, cand, n = 25) {
    
    if (sum(!is.na(targ)) == 0) {
        # if target is padding (as when running over array), return an NA array of expected size
        return(array(NA, dim = c(n, 2)))
    }
    
    # standard deviation for normalisation
    cand.sd <- apply(cand, 1:2, sd, na.rm = T)
    
    # distance from each point on target to equivalent in candidate space
    dist <- sweep(sqrt(sweep(cand, 1:2, targ, "-")^2), 1:2, cand.sd, "/")
    
    # mean distance across all variables
    m.dist <- apply(dist, 3:4, mean)
    
    # identify n analogues
    an <- which(m.dist <= sort(m.dist)[n+1], arr.ind = T)
    
    # find target vector and remove from analogue list
    t <- which(m.dist == 0, arr.ind = T)
    an <- an[!(an[,1] == t[1] & an[,2] == t[2]),]
    
    # return indices of analogues identified
    return(an)
}

prep.data <- function(lt) {
    
    fc <- list("ecmwf" = offset.forecast(ecmwf)[,,,lt,-1],
               "ncep" = offset.forecast(ncep)[,,,lt,-1], 
               "ukmo" = offset.forecast(ukmo)[,,,lt,-1])
    
    cand <- get.candidates(lt)
    
    an <- aaply(cand, 3:4, find.analogues, cand)        # around 5s to identify all analogues, single LT
    
    mean.errors <- apply(abind("ecmwf" = apply(forecast.errors(ecmwf)[,,,lt,-1], 1:3, mean),
                               "ncep" = apply(forecast.errors(ncep)[,,,lt,-1], 1:3, mean),
                               "ukmo" = apply(forecast.errors(ukmo)[,,,lt,-1], 1:3, mean), along = 0),
                         2:4, mean)
    tr <- aaply(an, 1:2, function(a) aaply(mean.errors, 1, "[", a))
    
    abind(obs, fc[[1]], fc[[2]], fc[[3]], aperm(tr, c(3,1:2,4)), along = 4)
}

fit.mimic <- function(o, fc, tr, arr.out = T) {
    
    Y.bar <- abind(lapply(fc, apply, 1, mean), along = 0)
    C <- abind(lapply(lapply(fc, t), cov), along = 0)
    
    D.i <- aaply(sweep(sweep(C, 1, sapply(fc, dim)[2,], "/"), 2:3, cov(Y.bar), "+"), 1, solve)
    
    Eta <- apply(t(tr), 2, mean)
    Lambda <- cov(t(tr))
    
    S <- Lambda + solve(apply(D.i, 2:3, sum))
    Tau <- (S %*% 
                ((solve(diag(5) + (apply(D.i, 2:3, sum) %*% Lambda))) %*%
                     apply(aaply(abind(Y.bar, D.i, along = 2), 1, 
                                 function(arr) arr[-1,] %*% arr[1,]), 2, sum))) - Eta
    
    if (arr.out) {
        # return results in array form (more useful for array operations)
        return(abind("Tau" = Tau, S, along = 2))
    } else {
        # otherwise, return output as list
        return(list(Tau = Tau, S = S))
    }
}

# ~ 2m to prep data over all leadtimes & dates
zz <- abind(lapply(1:15, prep.data), rev.along = 3.5)

# ~4s to run over all dates for single leadtime; ~1m for all leadtimes, all dates
system.time({
    ff <- suppressWarnings(aaply(zz, 2:4, function(arr) {
        fit.mimic(o = arr[,1],
                  fc = list("ecmwf" = arr[,2:51], "ncep" = arr[,52:71], "ukmo" = arr[,72:94]),
                  tr = arr[,95:119])
    }))})

saveRDS(ff, "../Models/fitted-an.rds")
# so: less than 4 minutes to fit over all leadtimes, including finding training analogues

####################################################################################################

# VERIFICATION                                                                                  ####

an.fitted <- readRDS("../Models/fitted-an.rds")

# multivariate CRPS

v.dat <- aaply(an.fitted, 1, function(arr) abind("o" = aperm(obs, c(2:3, 1)), arr, along = 4))

# ~1m to run over all leadtimes, all dates
es <- apply(v.dat[,-1,,,], 1:3, function(arr) es.crps(o = arr[,1], mu = arr[,2], sig = arr[,-(1:2)]))

plot(0:14, apply(es, 1, mean), type = "l", main = "Multivariate CRPS (energy score)",
     xlab = "", ylab = "")



hist.fitted <- readRDS("../Models/lambda-hist.rds")
f.hist <- abind(lapply(1:15, function(lt) {
    abind("o" = t(apply(obs, 1, cbind))[,26:630],
          "Tau" = hist.fitted$tau[,lt,],
          hist.fitted$s[,,lt,],
          along = 1)
}), along = 0)

es.hist <- apply(f.hist, c(1, 4), function(arr) es.crps(o = arr[1,], mu = arr[2,], sig = arr[-(1:2),]))

matplot(0:14, t(abind(apply(es, 1, mean), apply(es.hist, 1, mean), along = 0)), type = "l",
        col = c("red", "black"), lty = c(1,2), xlab = "", ylab = "",
        main = "Energy scores for posterior distribution")

legend("bottomright", bty = "n", lty = c(2,1), col = c("black", "red"), 
       legend = c("Immediate history", "2-day analogues"))

# can't get energy score for ensemble BMA (MOS should be possible but probably not worth coding)... compare univariate instead
bma.fitted <- readRDS("../Models/ensBMA-temps.rds")
mos.fitted <- readRDS("../Models/ensMOS-temps.rds")
ens.data <- readRDS("../Models/Ens-data.rds")

require(ensembleBMA)
crps.bma <- abind(sapply(1:15, function(lt) crps(bma.fitted[[lt]], ens.data[[lt]])[,"BMA"], simplify = F),
                  along = 0)

require(ensembleMOS)
crps.mos <- abind(sapply(1:15, function(lt) crps(mos.fitted[[lt]], ens.data[[lt]]), simplify = F),
                  along = 0)

crps.an <- apply(v.dat[,-1,,,], 1:3, 
                 function(arr) verification::crps(obs = arr[,"o"], 
                                                  pred = cbind(arr[,"Tau"],
                                                               sqrt(diag(arr[,-(1:2)]))))$crps)

dim(apply(crps.an, 1:2, mean)[1:2,])

crps.all <- abind("bma" = apply(array(crps.bma, dim = c(15, 606, 2)), c(1,3), mean),
                  "mos" = apply(array(crps.mos, dim = c(15, 606, 2)), c(1,3), mean),
                  "an" = t(apply(crps.an, 1:2, mean)[1:2,]),
                  along = 0)

pdf("../Plots/Analogue-fitting-first-attempt.pdf"); {
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    matplot(0:14, t(crps.all[,,1]), type = "l", ylim = range(crps.all, 0), lty = 1,
            col = c("red", "blue", "black"), lwd = c(1,1,2),
            xlab = "", ylab = "", main = "temp.n")
    legend("bottomright", bty = "n", lty = 1, lwd = c(1,1,2), col = c("red", "blue", "black"),
           legend = c("BMA", "MOS", "Analogues"))
    
    matplot(0:14, t(crps.all[,,2]), type = "l", ylim = range(crps.all, 0), lty = 1,
            col = c("red", "blue", "black"), lwd = c(1,1,2),
            xlab = "", ylab = "", main = "temp.s")
    legend("bottomright", bty = "n", lty = 1, lwd = c(1,1,2), col = c("red", "blue", "black"),
           legend = c("BMA", "MOS", "Analogues"))
    
    mtext("Marginal CRPS for analogue-fitted model vs BMA/MOS", outer = T)
}; dev.off()

####################################################################################################

# INVESTIGATION: S-MATRIX                                                                       ####

# boxplot of S over all models fitted
pdf("../Plots/Analogue-inv/Analogue-S.pdf"); {
    par(mfrow = c(5,1), mar = c(1,1,1,1), oma = c(0,0,2,0))
    invisible(sapply(1:15, function(lt) {
        invisible(sapply(dimnames(v.dat)[[4]], function(varb) {
            boxplot(apply(v.dat[lt,,,varb, -(1:2)], 3, cbind), ylim = c(-1,10))
            abline(h = c(-1,0,1), col = c("cyan3", "red", "cyan3"), lty = 2)
        }))
        mtext(paste0("Elements of S fitted with analogues - LT ", lt-1), outer = T)
    }))
}; dev.off()

# compare to S with historical training data
hist.fit <- readRDS("../Models/lambda-hist.rds")

pdf("../Plots/Analogue-inv/Hist-S.pdf"); {
    par(mfrow = c(5,1), mar = c(1,1,1,1), oma = c(0,0,2,0))
    
    invisible(sapply(1:15, function(lt) {
        invisible(sapply(dimnames(v.dat)[[4]], function(varb) {
            boxplot(apply(hist.fit$s[varb,,lt,], 1, cbind), ylim = c(-1,10))
            abline(h = c(-1,0,1), col = c("cyan3", "red", "cyan3"), lty = 2)
        }))
        mtext(paste0("Elements of S fitted with historical data - LT ", lt-1), outer = T)
    }))
}; dev.off()

pdf("../Plots/Analogue-inv/Compare-S.pdf"); {
    par(mfrow = c(5,1), mar = c(2,2,1,1), oma = c(0,0,2,0))
    
    invisible(sapply(1:15, function(lt) {
        invisible(sapply(dimnames(v.dat)[[4]], function(varb) {
            boxplot(apply(hist.fit$s[varb,,lt,], 1, cbind), ylim = c(-1,10), border = "orange")
            boxplot(apply(v.dat[lt,,,varb, -(1:2)], 3, cbind), add = T)
            abline(h = c(-1,0,1), col = c("cyan3", "red", "cyan3"), lty = 2)
        }))
        mtext(paste0("Elements of S fitted with historical data (gold) & analogues (black) - LT ", lt-1), outer = T)
    }))
}; dev.off()

# from LT 8, S matrix elements for temp.n & temp.s are slightly higher in the analogue-trained model
# than in the historically-trained model. Not sure why this should be, though...
# try comparing Lambda & Eta to investigate.

####################################################################################################

# INVESTIGATION: ETA & LAMBDA                                                                   ####

plot.path <- "../Plots/Analogue-inv/"
# improved performance after 7 days is presumably due largely to Eta & Lambda?
{
    # retrieve analogues over which Eta & Lambda were trained
    get.training.data <- function(lt) {
        cand <- get.candidates(lt)
        
        an <- aaply(cand, 3:4, find.analogues, cand)        # around 5s to identify all analogues, single LT
        
        mean.errors <- apply(abind("ecmwf" = apply(forecast.errors(ecmwf)[,,,lt,-1], 1:3, mean),
                                   "ncep" = apply(forecast.errors(ncep)[,,,lt,-1], 1:3, mean),
                                   "ukmo" = apply(forecast.errors(ukmo)[,,,lt,-1], 1:3, mean), along = 0),
                             2:4, mean)
        tr <- aaply(an, 1:2, function(a) aaply(mean.errors, 1, "[", a))
        tr
    }
    
    tr <- abind(sapply(formatC(0:14), get.training.data, simplify = F), along = 0)
    
    Eta <- apply(tr, 1:4, mean)
    Lambda <- aaply(tr, 1:3, function(ts) cov(t(ts)))
    
    pdf(paste0(plot.path, "Eta-fitted.pdf")); {
        par(mfrow = c(5,1), mar = c(2,2,1,1), oma = c(0,0,2,0))
        invisible(sapply(dimnames(Eta)[[4]], function(varb) {
            boxplot(apply(Eta[,,,varb], 1, rbind), ylim = c(-4,4))
            abline(h = c(0,1,-1), col = c("red", rep("cyan3", 2)), lty = 2)
            legend("topleft", legend = varb, bty = "n")
        }))
        mtext(expression(paste("Values of ", eta, " at all leadtimes")), outer = T)
    }; dev.off()

    pdf(paste0(plot.path, "Lambda-fitted.pdf")); {
        par(mfrow = c(5,1), mar = c(2,2,1,1), oma = c(0,0,2,0))
        invisible(sapply(1:15, function(lt) {
            invisible(sapply(dimnames(Lambda)[[4]], function(varb) {   
                boxplot(apply(Lambda[lt,,,varb,], 3, rbind), ylim = c(-1,10))
                abline(h = c(0,1,-1), col = c("red", rep("cyan3", 2)), lty = 2)
                legend("topright", legend = varb, bty = "n")
            }))
            mtext(bquote("Values of" ~ Lambda ~ "at leadtime" ~ .(lt)), outer = T)
        }))
    }; dev.off()
    
####################################################################################################

# INVESTIGATION: PATTERNS IN FORECAST ERRORS                                                    ####
    
fce <- abind("ecmwf" = apply(forecast.errors(ecmwf)[,,,,-1], 1:4, mean),
             "ncep" = apply(forecast.errors(ncep)[,,,,-1], 1:4, mean),
             "ukmo" = apply(forecast.errors(ukmo)[,,,,-1], 1:4, mean),
             along = 0)
    
yrng <- c(-5,5)   
pdf(paste0(plot.path, "FC-errors.pdf")); {
    par(mfrow = c(4,1), mar = c(2,2,1,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(fce)[[2]], function(varb) {
        invisible(sapply(dimnames(fce)[[1]], function(mm) {
            boxplot(apply(fce[mm,varb,,,], 3, rbind), ylim = yrng)
            abline(h = 0, col = "red", lty = 2)
            abline(h = mean(fce[mm, varb,,,]), col = "cyan3", lty = 2)
            legend("topleft", bty = "n", legend = mm)
        }))
        boxplot(apply(apply(fce[,varb,,,], 2:4, mean), 3, rbind), ylim = yrng)
        abline(h = 0, col = "red", lty = 2)
        abline(h = mean(fce[, varb,,,]), col = "cyan3", lty = 2)
        legend("topleft", bty = "n", legend = "Mean")
        
        mtext(paste0("Forecast errors for ", varb, " in each model, at each leadtime"), outer = T)
    }))
}; dev.off()

####################################################################################################

# ANALOGUES WITH & WITHOUT PRINCIPAL COMPONENTS                                                 ####

# currently fixed at 2 days' forecasts, 1 day's observations. Generalise later
get.candidates <- function(lt) {
    
    cand <- abind("y.o" = obs[,1:89,],
                  "y.ecmwf" = apply(offset.forecast(ecmwf)[,1:89,,lt,-1], 1:3, mean),
                  "y.ncep" = apply(offset.forecast(ncep)[,1:89,,lt,-1], 1:3, mean),
                  "y.ukmo" = apply(offset.forecast(ukmo)[,1:89,,lt,-1], 1:3, mean),
                  "c.ecmwf" = apply(offset.forecast(ecmwf)[,2:90,,lt,-1], 1:3, mean),
                  "c.ncep" = apply(offset.forecast(ncep)[,2:90,,lt,-1], 1:3, mean),
                  "c.ukmo" = apply(offset.forecast(ukmo)[,2:90,,lt,-1], 1:3, mean),
                  along = 0)
    cand <- abind(array(NA, dim = dim(cand[,,1,])), cand, along = 3)
    cand[,1:2,,]
}

find.analogues <- function(targ, cand, n = 25) {
    
    if (sum(!is.na(targ)) == 0) {
        # if target is padding (as when running over array), return an NA array of expected size
        return(array(NA, dim = c(n, 2)))
    }
    
    # standard deviation for normalisation
    cand.sd <- apply(cand, 1:2, sd, na.rm = T)
    
    # distance from each point on target to equivalent in candidate space
    dist <- sweep(sqrt(sweep(cand, 1:2, targ, "-")^2), 1:2, cand.sd, "/")
    
    # mean distance across all variables
    m.dist <- apply(dist, 3:4, mean)
    
    # identify n analogues
    an <- which(m.dist <= sort(m.dist)[n+1], arr.ind = T)
    
    # find target vector and remove from analogue list
    t <- which(m.dist == 0, arr.ind = T)
    an <- an[!(an[,1] == t[1] & an[,2] == t[2]),]
    
    # return indices of analogues identified
    return(an)
}

prep.data <- function(lt) {
    
    fc <- list("ecmwf" = offset.forecast(ecmwf)[,,,lt,-1],
               "ncep" = offset.forecast(ncep)[,,,lt,-1], 
               "ukmo" = offset.forecast(ukmo)[,,,lt,-1])
    
    cand <- get.candidates(lt)
    
    an <- aaply(cand, 3:4, find.analogues, cand)        # around 5s to identify all analogues, single LT
    
    mean.errors <- apply(abind("ecmwf" = apply(forecast.errors(ecmwf)[,,,lt,-1], 1:3, mean),
                               "ncep" = apply(forecast.errors(ncep)[,,,lt,-1], 1:3, mean),
                               "ukmo" = apply(forecast.errors(ukmo)[,,,lt,-1], 1:3, mean), along = 0),
                         2:4, mean)
    tr <- aaply(an, 1:2, function(a) aaply(mean.errors, 1, "[", a))
    
    abind(obs, fc[[1]], fc[[2]], fc[[3]], aperm(tr, c(3,1:2,4)), along = 4)
}

fit.mimic <- function(o, fc, tr, arr.out = T) {
    
    Y.bar <- abind(lapply(fc, apply, 1, mean), along = 0)
    C <- abind(lapply(lapply(fc, t), cov), along = 0)
    
    D.i <- aaply(sweep(sweep(C, 1, sapply(fc, dim)[2,], "/"), 2:3, cov(Y.bar), "+"), 1, solve)
    
    Eta <- apply(t(tr), 2, mean)
    Lambda <- cov(t(tr))
    
    S <- Lambda + solve(apply(D.i, 2:3, sum))
    Tau <- (S %*% 
                ((solve(diag(5) + (apply(D.i, 2:3, sum) %*% Lambda))) %*%
                     apply(aaply(abind(Y.bar, D.i, along = 2), 1, 
                                 function(arr) arr[-1,] %*% arr[1,]), 2, sum))) - Eta
    
    if (arr.out) {
        # return results in array form (more useful for array operations)
        return(abind("Tau" = Tau, S, along = 2))
    } else {
        # otherwise, return output as list
        return(list(Tau = Tau, S = S))
    }
}

# ~3m to run from scratch for all leadtimes, all dates
system.time({
    
    zz <- abind(lapply(1:15, prep.data), rev.along = 3.5)
    
    ff <- suppressWarnings(aaply(zz, 2:4, function(arr) {
        fit.mimic(o = arr[,1],
                  fc = list("ecmwf" = arr[,2:51], "ncep" = arr[,52:71], "ukmo" = arr[,72:94]),
                  tr = arr[,95:119])
    }))})

saveRDS(ff, "../Models/fitted-an-NO-PC.rds")

v.npc <- aaply(ff, 1, function(arr) abind("o" = aperm(obs, c(2:3, 1)), arr, along = 4))
es.npc <- apply(v.npc[,-1,,,], 1:3, function(arr) es.crps(o = arr[,1], mu = arr[,2], sig = arr[,-(1:2)]))

pdf(paste0(plot.path, "posterior-ES.pdf")); {
    matplot(0:14,
            t(abind(apply(es, 1, mean),
                    apply(es.hist, 1, mean),
                    apply(es.npc, 1, mean),
                    along = 0)), 
            type = "l", col = c("red", "black", "blue"), lty = c(1,2,1), xlab = "", ylab = "",
            main = "Energy scores for posterior distribution")
    
    legend("bottomright", bty = "n", lty = c(2,1,1), col = c("black", "red", "blue"), 
           legend = c("Immediate history", "2-day analogues with PC", "2-day analogues - no PC"),
           title = "Training set", title.adj = 0)
}; dev.off()

####################################################################################################

# ANALOGUES FOUND OVER 3 DAYS' FORECAST & OBS                                                   ####

get.3d.candidates <- function(lt) {
    
    models <- list("ecmwf" = ecmwf, "ncep" = ncep, "ukmo" = ukmo)
    
    dn <- c(paste("obs", 1:2, sep = "."), sapply(names(models), paste, 0:2, sep = "."))
    cand <- array(NA, dim = c(11, 5, 90, 7),
                  dimnames = append(dimnames(obs), list(dn), 0))
    
    cand["obs.1",,2:90,] <- obs[, 1:89,]
    cand["obs.2",,3:90,] <- obs[, 1:88,]
    
    # switch to function later if necessary
    offset.ec <- apply(offset.forecast(ecmwf)[,,,lt,-1], 1:3, mean)
        cand["ecmwf.0",,,] <- offset.ec
        cand["ecmwf.1",,2:90,] <- offset.ec[, 1:89,]
        cand["ecmwf.2",,3:90,] <- offset.ec[, 1:88,]
        
    offset.nc <- apply(offset.forecast(ncep)[,,,lt,-1], 1:3, mean)
        cand["ncep.0",,,] <- offset.nc
        cand["ncep.1",,2:90,] <- offset.nc[, 1:89,]
        cand["ncep.2",,3:90,] <- offset.nc[, 1:88,]
        
    offset.mo <- apply(offset.forecast(ukmo)[,,,lt,-1], 1:3, mean)
        cand["ukmo.0",,,] <- offset.mo
        cand["ukmo.1",,2:90,] <- offset.mo[, 1:89,]
        cand["ukmo.2",,3:90,] <- offset.mo[, 1:88,]
        
    cand
}

# NOTE CHANGE TO ERROR CHECKING IN FIRST LINE OF FUNCTION 
find.analogues <- function(targ, cand, n = 25) {
    
    if (sum(is.na(targ)) > 0) {
        # if target is padding (as when running over array), return an NA array of expected size
        return(array(NA, dim = c(n, 2)))
    }
    
    # standard deviation for normalisation
    cand.sd <- apply(cand, 1:2, sd, na.rm = T)
    
    # distance from each point on target to equivalent in candidate space
    dist <- sweep(sqrt(sweep(cand, 1:2, targ, "-")^2), 1:2, cand.sd, "/")
    
    # mean distance across all variables
    m.dist <- apply(dist, 3:4, mean)
    
    # identify n analogues
    an <- which(m.dist <= sort(m.dist)[n+1], arr.ind = T)
    
    # find target vector and remove from analogue list
    t <- which(m.dist == 0, arr.ind = T)
    an <- an[!(an[,1] == t[1] & an[,2] == t[2]),]
    
    # return indices of analogues identified
    return(an)
}

# prep data
{
    fc <- list("ecmwf" = offset.forecast(ecmwf)[,,,lt,-1],
               "ncep" = offset.forecast(ncep)[,,,lt,-1], 
               "ukmo" = offset.forecast(ukmo)[,,,lt,-1])
    
    cand <- get.3d.candidates(lt)
    
    an <- aaply(cand.3d, 3:4, find.analogues, cand.3d)        # around 5s to identify all analogues, single LT
    
    mean.errors <- apply(abind("ecmwf" = apply(forecast.errors(ecmwf)[,,,lt,-1], 1:3, mean),
                               "ncep" = apply(forecast.errors(ncep)[,,,lt,-1], 1:3, mean),
                               "ukmo" = apply(forecast.errors(ukmo)[,,,lt,-1], 1:3, mean), along = 0),
                         2:4, mean)
    tr <- aaply(an, 1:2, function(a) aaply(mean.errors, 1, "[", a))
    
dat <- abind(obs, fc[[1]], fc[[2]], fc[[3]], aperm(tr, c(3,1:2,4)), along = 4)
}

fit.mimic <- function(o, fc, tr, arr.out = T) {
    
    Y.bar <- abind(lapply(fc, apply, 1, mean), along = 0)
    C <- abind(lapply(lapply(fc, t), cov), along = 0)
    
    D.i <- aaply(sweep(sweep(C, 1, sapply(fc, dim)[2,], "/"), 2:3, cov(Y.bar), "+"), 1, solve)
    
    Eta <- apply(t(tr), 2, mean)
    Lambda <- cov(t(tr))
    
    S <- Lambda + solve(apply(D.i, 2:3, sum))
    Tau <- (S %*% 
                ((solve(diag(5) + (apply(D.i, 2:3, sum) %*% Lambda))) %*%
                     apply(aaply(abind(Y.bar, D.i, along = 2), 1, 
                                 function(arr) arr[-1,] %*% arr[1,]), 2, sum))) - Eta
    
    if (arr.out) {
        # return results in array form (more useful for array operations)
        return(abind("Tau" = Tau, S, along = 2))
    } else {
        # otherwise, return output as list
        return(list(Tau = Tau, S = S))
    }
}

cand.3d <- get.3d.candidates(5)
an <- find.analogues(cand.3d[,,5,5], cand.3d)
dat <- prep.data(5)

    zz <- abind(lapply(1:15, prep.data), rev.along = 3.5)
    
    ff.3d <- suppressWarnings(aaply(zz, 2:4, function(arr) {
        fit.mimic(o = arr[,1],
                  fc = list("ecmwf" = arr[,2:51], "ncep" = arr[,52:71], "ukmo" = arr[,72:94]),
                  tr = arr[,95:119])
    }))

saveRDS(ff.3d, "../Models/fitted-an-3d.rds")

# energy score for comparison...
v.3d <- aaply(ff.3d, 1, function(arr) abind("o" = aperm(obs, c(2:3, 1)), arr, along = 4))
es.3d <- apply(v.3d[,-(1:2),,,], 1:3, function(arr) es.crps(o = arr[,1], mu = arr[,2], sig = arr[,-(1:2)]))

pdf(paste0(plot.path, "posterior-ES.pdf")); {
    matplot(0:14,
            t(abind(apply(es, 1, mean),
                    apply(es.hist, 1, mean),
                    apply(es.npc, 1, mean),
                    apply(es.3d, 1, mean),
                    along = 0)), 
            type = "l", col = c("red", "black", "blue", "green3"), lty = c(1,2,1,1), xlab = "", ylab = "",
            main = "Energy scores for posterior distribution")
    
    legend("bottomright", bty = "n", lty = c(2,1,1, 1), col = c("black", "red", "blue", "green3"), 
           legend = c("Immediate history", "2-day analogues with PC", "2-day analogues - no PC", "3-day analogues with PC"),
           title = "Training set", title.adj = 0)
}; dev.off()

####################################################################################################

# RETRAINING HISTORICAL MODEL                                                                   ####

# refitting the historically-trained model should confirm the code -
# if result is the same, we have a success.

ens.mean.fc <- apply(apply(abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
                                 "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
                                 "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
                                 along = 0), 2:5, mean), c(1, 4), rbind)

ens.mean.error <- sweep(ens.mean.fc, 1:2, apply(obs, 1, rbind), "-")

# get offset error: 25 prior forecast errors for 5 vars, 5 leadtimes, 605 forecasts
hist.mean.error <- abind(invisible(lapply(1:(630-25), 
                                          function(i) ens.mean.error[i+(0:24),,])),
                         rev.along = 0)

# create single array of all data to fit models
zz <- abind(lapply(1:15, function(lt) {
    abind("o" = apply(obs, 1, cbind)[26:630,],
          apply(offset.forecast(ecmwf)[,,,lt,-1], c(1, 4), cbind)[26:630,,],
          apply(offset.forecast(ncep)[,,,lt,-1], c(1, 4), cbind)[26:630,,],
          apply(offset.forecast(ukmo)[,,,lt,-1], c(1, 4), cbind)[26:630,,],
          aperm(hist.mean.error[,,lt,], c(3,2,1)),
          along = 3)
}), along = 0)

# ~80s to fit all 605 timesteps at all leadtimes
ff <- suppressWarnings(aaply(zz, 1:2, function(arr) {
        fit.mimic(o = arr[,1],
                  fc = list("ecmwf" = arr[,2:51], "ncep" = arr[,52:71], "ukmo" = arr[,72:94]),
                  tr = arr[,95:119])
    }))
v.hist <- aaply(ff, 1, function(arr) abind("o" = apply(obs, 1, cbind)[26:630,], arr, along = 3))
es.hist <- apply(v.hist, 1:2, function(arr) es.crps(o = arr[,1], mu = arr[,2], sig = arr[,-(1:2)]))
    
saveRDS(ff, "../Models/Mimic-hist.rds")

# energy score for model fitted with 3-day analogues
{
    ff.3d <- readRDS("../Models/fitted-an-3d.rds")
    v.3d <- aaply(ff.3d, 1, function(arr) abind("o" = aperm(obs, c(2:3, 1)), arr, along = 4))
    es.3d <- apply(v.3d[,-(1:2),,,], 1:3, function(arr) es.crps(o = arr[,1], mu = arr[,2], sig = arr[,-(1:2)]))
}

# energy score for model fitted with 2-day analogues
{
    ff.2d <- readRDS("../Models/fitted-an.rds")
    v.2d <- aaply(ff.2d, 1, function(arr) abind("o" = aperm(obs, c(2:3, 1)), arr, along = 4))
    es.2d <- apply(v.2d[,-(1:2),,,], 1:3, function(arr) es.crps(o = arr[,1], mu = arr[,2], sig = arr[,-(1:2)]))
}

# reload & check energy score for historical model
{
    hist.org <- readRDS("../Models/lambda-hist.rds")
    f.hist.org <- abind(lapply(1:15, function(lt) {
        abind("o" = t(apply(obs, 1, cbind))[,26:630],
              "Tau" = hist.org$tau[,lt,],
              hist.org$s[,,lt,],
              along = 1)
    }), along = 0)
    
    es.hist.org <- apply(f.hist.org, c(1, 4), function(arr) es.crps(o = arr[1,], mu = arr[2,], sig = arr[-(1:2),]))
}

matplot(0:14,
        t(abind(apply(es.hist, 1, mean),
                apply(es.2d, 1, mean),
                apply(es.3d, 1, mean),
                apply(es.hist.org, 1, mean),
                along = 0)), 
        type = "l", col = c("black", "red", "blue", "orange"), lty = c(1,2,1,1), xlab = "", ylab = "",
        main = "Energy scores for posterior distribution")

legend("bottomright", bty = "n", lty = c(2,1,1,2), col = c("black", "red", "blue", "black"), 
       legend = c("Immediate history", "2-day analogues", "3-day analogues", "History (org model)"),
       title = "Training set", title.adj = 0)

####################################################################################################

# MULTIVARIATE VERIFICATION METRICS                                                             ####

# energy scores, rmse / spread, histograms, determinant sharpness etc

# import all models of interest & create verification arrays
{
    v.arr <- list("an" = apply(aaply(readRDS("../Models/fitted-an.rds"), 1,
                                     function(arr) abind("o" = aperm(obs, c(2:3, 1)), arr, along = 4))[,-1,,,],
                               c(1,4:5), cbind),
                  "npc" = apply(aaply(readRDS("../Models/fitted-an-NO-PC.rds"), 1,
                                      function(arr) abind("o" = aperm(obs, c(2:3, 1)), arr, along = 4))[,-1,,,],
                                c(1,4:5), cbind),
                  "3d" = apply(aaply(readRDS("../Models/fitted-an-3d.rds"), 1, 
                                     function(arr) abind("o" = aperm(obs, c(2:3, 1)), arr, along = 4))[,-(1:2),,,],
                               c(1,4:5), cbind),
                  "hist" = aperm(aaply(readRDS("../Models/Mimic-hist.rds"), 1,
                                       function(arr) abind("o" = apply(obs, 1, cbind)[26:630,], arr, along = 3)),
                                 c(2,1,3:4)))
}

# set up plotting functions & parameters
{
    plot.path <- "../Plots/Analogue-perf/"
    
    ll <- list("lty" = 1, "lwd" = c(2,1,1,1), "col" = c("black", "blue", "green3", "red"))
    
    res.plot <- function(res, main = "", ...) {
        
        res.arr <- abind(lapply(res, apply, 2, mean), rev.along = 0)
        matplot(0:14, res.arr, type = "l", xlab = "", ylab = "",
                lty = ll$lty, lwd = ll$lwd, col = ll$col, 
                main = main, ylim = range(0, res.arr), ...)
        
        legend("bottomright", bty = "n", lty = ll$lty, lwd = ll$lwd, col = ll$col, 
               legend = c("2-day analogues", "2-day, no PC", "3-day analogues", "Recent history"),
               title = "Training set")
    }
}

# energy score
{
    es <- lapply(v.arr, apply, 1:2, function(arr) es.crps(o = arr[,1], mu = arr[,2], sig = arr[,-(1:2)]))
    res.plot(es, main = "Energy score")
}

# determinant sharpness
{
    ds <- lapply(v.arr, apply, 1:2, function(arr) det.sharpness(arr[,-(1:2)]))
    res.plot(ds, main = "Determinant sharpness")
}

# RMSE
# is multivariate RMSE same as sqrt(vector norm)? check by rearranging.
{
    rmse <- lapply(v.arr, apply, 1:2, function(arr) sqrt(mean((arr[,2] - arr[,1])^2)))
    res.plot(rmse, main = "RMSE")
}

# deviation from uniformity
# NEED TO CALCULATE ANALYTICALLY
# currently calculated over synthetic ensemble of size 100 from each fitted model
{
    require(mvtnorm)
    ud <- abind("an" = apply(aaply(readRDS("../Models/fitted-an.rds")[,-1,,,], 1:3, 
                                   function(arr) rmvnorm(100, mean = arr[,1], sigma = arr[,-1])),
                             1, function(se) u.dev(verif.ranks(o = obs[,-1,], ens = aperm(se, c(4,1:3))), l = 101)),
                "npc" = apply(aaply(readRDS("../Models/fitted-an-NO-PC.rds")[,-1,,,], 1:3, 
                                    function(arr) rmvnorm(100, mean = arr[,1], sigma = arr[,-1])),
                              1, function(se) u.dev(verif.ranks(o = obs[,-1,], ens = aperm(se, c(4,1:3))), l = 101)),
                "3d" = apply(aaply(readRDS("../Models/fitted-an-3d.rds")[,-(1:2),,,], 1:3, 
                                   function(arr) rmvnorm(100, mean = arr[,1], sigma = arr[,-1])),
                             1, function(se) u.dev(verif.ranks(o = obs[,-(1:2),], ens = aperm(se, c(4,1:3))), l = 101)),
                "hist" = apply(aaply(readRDS("../Models/Mimic-hist.rds"), 1:2, 
                                     function(arr) rmvnorm(100, mean = arr[,1], sigma = arr[,-1])), 1,
                               function(se) u.dev(verif.ranks(o = apply(obs, 1, cbind)[26:630,],
                                                                        ens = aperm(se, c(1,3,2))))),
                along = 0)
}

pdf(paste0(plot.path, "mimic-evaluation.pdf")); {
    
    par(mfrow = c(2,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    res.plot(es, main = "Energy score")
    matplot(0:14, t(ud), type = "l", xlab = "", ylab = "",
            lty = ll$lty, lwd = ll$lwd, col = ll$col, 
            main = "Deviation from uniformity", ylim = range(0, ud))
    legend("topright", bty = "n", lty = ll$lty, lwd = ll$lwd, col = ll$col, 
           legend = c("2-day analogues", "2-day, no PC", "3-day analogues", "Recent history"),
           title = "Training set")
    
    res.plot(rmse, main = "RMSE")
    res.plot(ds, main = "Determinant sharpness")
    
    mtext("Performance of mimic variants", outer = T)
    
}; dev.off()

# can determinant sharpness be used as multivariate equivalent to spread?  
# doesn't look like it. Confirm by checking over marginals.

# ordinate transform histograms
bot <- lapply(v.arr, apply, 1:2, function(arr) box.dot(o = arr[,1], mu = arr[,2], s = arr[,-(1:2)]))


bot.max <- apply(abind(lapply(invisible(sapply(1:15, function(lt) {
    invisible(sapply(bot, function(or) {
        hist(or[,lt], breaks = c(0:10)/10, plot = F)$density
    }, simplify = F))}, simplify = F)), abind, along = 0), along = 0), 1, max)


pdf(paste0(plot.path, "mimic-BOT-hists.pdf")); {
    hist.cols <- adjustcolor(c("gold", "cyan3", "steelblue", "green3"), alpha = 0.4)
    
    par(mfrow = c(15, 4), mar = c(0.1,0.1,0.1,0.1), oma = c(0,0,2,0))
    
    invisible(sapply(1:15, function(lt) {
        invisible(sapply(1:4, function(or) {
            
            hist(bot[[or]][,lt], breaks = c(0:10)/10, prob = T, ylim = c(0, bot.max[lt]),
                 xlab = "", ylab = "", main = "", col = hist.cols[or], xaxt = "n", yaxt = "n")
            abline(h = 1, col = "red3", lty = 2)
            legend("top", legend = paste0(names(bot)[or], " - LT ", lt-1), bty = "n")
        }))
    }))
    mtext("Box ordinate density transform histograms for each mimic", outer = T)
}; dev.off()


####################################################################################################

# MARGINAL METRICS                                                                              ####

# (allows comparison to ensemble MOS, ensemble BMA)