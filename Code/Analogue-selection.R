
# analogue distances: how to select the 'best' analogues?

# final version must only look at PAST REFORECASTS. Not future models from same year.

# quantiles are essentially just choosing the closest n values - still essentially arbitrary
# but still probably only interested in values of eg. 4wks' data (say 28)

# test difference of weighting scheme for 2wks, 3wks, 4wks, 5wks.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

library("SX.weather")

# search space - all candidates, single day's system, search only one leadtime
cand <- abind(invisible(sapply(1:15, function(lt) {
    cand <- abind("y.o" = obs[,1:89,],
                  "c.ecmwf" = apply(offset.forecast(ecmwf)[,2:90,,lt,-1], 1:3, mean),
                  "c.ncep" = apply(offset.forecast(ncep)[,2:90,,lt,-1], 1:3, mean),
                  "c.ukmo" = apply(offset.forecast(ukmo)[,2:90,,lt,-1], 1:3, mean),
                  along = 0)
    abind(array(NA, dim = dim(cand[,,1,])), cand, along = 3)
}, simplify = F)), along = 5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# analogues identified across all 5 variables                                                   ####

# normalise search space
matplot(0:14, apply(cand[-1,,,,], c(5, 2), sd, na.rm = T), type = "l", 
        xlab = "Leadtime", ylab = "SD", ylim = range(0, cand.sd), main = "SD across all forecasts")
c.norm <- sweep(cand, c(2,5), apply(cand[-1,,,,], c(2, 5), sd, na.rm = T), "/")    

# get all normalised distances and weights

n.dist <- invisible(aaply(c.norm, 5, function(ssp) {
    invisible(aaply(ssp[,,-1,], 3:4, function(targ) {
        apply(abs(sweep(ssp, 1:2, targ, "-")), 3:4, mean, na.rm = T)
    }))
})) 

# get all analogue indices & distances (~20s for all leadtimes)
n.an <- aaply(n.dist, c(2:3,1), function(dd) {
    abind(sapply(1:6, function(wk) {
        d.lim <- sort(dd)[7 * wk + 1]
        an <- which(dd > sort(dd)[7 * (wk-1) + 1] & dd <= sort(dd)[7 * wk + 1], arr.ind = T)
        cbind(an, "d" = dd[an])
    }, simplify = F), along = 3)
})

# fit all errors based on simplified Bayesian model
model.errors <- function(fc, o) {
    
    # convert forecasts to errors
    fc <- lapply(fc, sweep, 1, o)
    Y.bar <- abind(lapply(fc, apply, 1, mean), along = 0)
    C <- abind(lapply(lapply(fc, t), cov), along = 0)
    D.i <- aaply(sweep(sweep(C, 1, sapply(fc, dim)[2, ], "/"), 
                       2:3, cov(Y.bar), "+"), 1, solve)
    
    S <- solve(apply(D.i, 2:3, sum))
    
    tau <- S %*% apply(aaply(abind(Y.bar, D.i, along = 2), 1, function(arr) arr[-1, ] %*% arr[1, ]), 2, sum)
    
    return(abind(Tau = tau, S, along = 2))
}

# should not take longer than running normal mimic model! ~75 seconds for all leadtimes
err <- abind(sapply(1:15, function(lt) {
    aaply(abind(offset.forecast(ecmwf)[,,,lt,-1],
          offset.forecast(ncep)[,,,lt,-1],
          offset.forecast(ukmo)[,,,lt,-1],
          obs, along = 4), 2:3, function(arr) model.errors(fc = list(arr[,1:50], arr[,51:70], arr[,71:93]),
                                                           o = arr[,94]))
    }, simplify = F), along = 2.5)

saveRDS(err, "../fitted-errors.rds")

# extract all analogue errors, with weights (~50s for all leadtimes, 6wks' analogues)
an.err <- abind(sapply(1:15, function(lt) {
    apply(n.an[,,lt,,,], c(1:3, 5), function(a) c(apply(err[,,lt,,1], 3, "[", t(a[1:2])), a[3]))
}, simplify = F), along = 3.5)

saveRDS(an.err, "../fitted-error-analogues.rds")

extract.errors <- function(an.err, n.weeks = 2, normalise.weights = T) {
    tr <- apply(an.err[,,,,,1:n.weeks, drop = F], 1:4, cbind)
    tr <- aaply(tr, c(3:5,1), function(er) c(er, "w" = er["d"]^-1))
    
    if (normalise.weights) {
        cw <- which(dimnames(tr)[[5]] == "w.d")
        tr[,,,,cw] <- aaply(tr[,,,,cw], c(1:3), function(er) er / sum(er))
    }
    
    names(dimnames(tr)) <- c("day", "year", "LT", "analogue", "values")
    # returns array with single layer of observations, weights etc
    return(tr)
}

tr.1w <- extract.errors(an.err, n.weeks = 1)       # ~  9s for 1wk
tr.2w <- extract.errors(an.err, n.weeks = 2)       # ~ 17s for 2wk
tr.3w <- extract.errors(an.err, n.weeks = 3)       # ~ 26s for 3wk
tr.4w <- extract.errors(an.err, n.weeks = 4)       # ~ 38s for 4wk
tr.5w <- extract.errors(an.err, n.weeks = 5)       # ~ 56s for 5wk
tr.6w <- extract.errors(an.err, n.weeks = 6)       # ~ 70s for 6wk

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# weighted eta & lambda                                                                         ####

# ~7-15s for all leadtimes, regardless of number of analogues used (with weights normalised)
wel.1wk <- aaply(tr.1w, 1:3, function(ts) abind(cov.wt(ts[,1:5], wt = ts[,7])[2:1], along = 1))
wel.2wk <- aaply(tr.2w, 1:3, function(ts) abind(cov.wt(ts[,1:5], wt = ts[,7])[2:1], along = 1))
wel.3wk <- aaply(tr.3w, 1:3, function(ts) abind(cov.wt(ts[,1:5], wt = ts[,7])[2:1], along = 1))
wel.4wk <- aaply(tr.4w, 1:3, function(ts) abind(cov.wt(ts[,1:5], wt = ts[,7])[2:1], along = 1))
wel.5wk <- aaply(tr.5w, 1:3, function(ts) abind(cov.wt(ts[,1:5], wt = ts[,7])[2:1], along = 1))

system.time({
    wel.6wk <- aaply(tr.6w, 1:3, function(ts) abind(cov.wt(ts[,1:5], wt = ts[,7])[2:1], along = 1))
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# and now, to fit the models using weighted eta and lambda                                      ####

fit.model <- function(fc, eta, Lambda) {
    
    d <- dim(fc[[1]])[[1]]
    Y.bar <- abind(lapply(fc, apply, 1, mean), along = 0)
    C <- abind(lapply(lapply(fc, t), cov), along = 0)
    D.i <- aaply(sweep(sweep(C, 1, sapply(fc, dim)[2, ], "/"), 
                       2:3, cov(Y.bar), "+"), 1, solve)
    
    S <- Lambda + solve(apply(D.i, 2:3, sum))
    Tau <- (S %*% ((solve(diag(d) + (apply(D.i, 2:3, sum) %*% 
                                         Lambda))) %*% apply(aaply(abind(Y.bar, D.i, along = 2), 
                                                                   1, function(arr) arr[-1, ] %*% arr[1, ]), 2, sum))) - 
        eta
    
    return(abind(Tau = Tau, S, along = 2))
}

# ~2m to fit all models, regardless of training set used
fit.1wk <- aaply(abind(offset.forecast(ecmwf)[,-1,,,-1],
                      offset.forecast(ncep)[,-1,,,-1],
                      offset.forecast(ukmo)[,-1,,,-1],
                      aperm(wel.1wk, c(5,1:3,4)), along = 5), 2:4, function(arr) {
                          fit.model(fc = list("ecmwf" = arr[,1:50],
                                              "ncep" = arr[,51:70],
                                              "ukmo" = arr[,71:93]),
                                    eta = arr[,94],
                                    Lambda = arr[,95:99])
                      })

fit.2wk <- aaply(abind(offset.forecast(ecmwf)[,-1,,,-1],
                       offset.forecast(ncep)[,-1,,,-1],
                       offset.forecast(ukmo)[,-1,,,-1],
                       aperm(wel.2wk, c(5,1:3,4)), along = 5), 2:4, function(arr) {
                           fit.model(fc = list("ecmwf" = arr[,1:50],
                                               "ncep" = arr[,51:70],
                                               "ukmo" = arr[,71:93]),
                                     eta = arr[,94],
                                     Lambda = arr[,95:99])
                       })

fit.3wk <- aaply(abind(offset.forecast(ecmwf)[,-1,,,-1],
                       offset.forecast(ncep)[,-1,,,-1],
                       offset.forecast(ukmo)[,-1,,,-1],
                       aperm(wel.3wk, c(5,1:3,4)), along = 5), 2:4, function(arr) {
                           fit.model(fc = list("ecmwf" = arr[,1:50],
                                               "ncep" = arr[,51:70],
                                               "ukmo" = arr[,71:93]),
                                     eta = arr[,94],
                                     Lambda = arr[,95:99])
                       })

fit.4wk <- aaply(abind(offset.forecast(ecmwf)[,-1,,,-1],
                       offset.forecast(ncep)[,-1,,,-1],
                       offset.forecast(ukmo)[,-1,,,-1],
                       aperm(wel.4wk, c(5,1:3,4)), along = 5), 2:4, function(arr) {
                           fit.model(fc = list("ecmwf" = arr[,1:50],
                                               "ncep" = arr[,51:70],
                                               "ukmo" = arr[,71:93]),
                                     eta = arr[,94],
                                     Lambda = arr[,95:99])
                       })

fit.5wk <- aaply(abind(offset.forecast(ecmwf)[,-1,,,-1],
                       offset.forecast(ncep)[,-1,,,-1],
                       offset.forecast(ukmo)[,-1,,,-1],
                       aperm(wel.5wk, c(5,1:3,4)), along = 5), 2:4, function(arr) {
                           fit.model(fc = list("ecmwf" = arr[,1:50],
                                               "ncep" = arr[,51:70],
                                               "ukmo" = arr[,71:93]),
                                     eta = arr[,94],
                                     Lambda = arr[,95:99])
                       })

fit.6wk <- aaply(abind(offset.forecast(ecmwf)[,-1,,,-1],
                       offset.forecast(ncep)[,-1,,,-1],
                       offset.forecast(ukmo)[,-1,,,-1],
                       aperm(wel.6wk, c(5,1:3,4)), along = 5), 2:4, function(arr) {
                           fit.model(fc = list("ecmwf" = arr[,1:50],
                                               "ncep" = arr[,51:70],
                                               "ukmo" = arr[,71:93]),
                                     eta = arr[,94],
                                     Lambda = arr[,95:99])
                       })

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# compare RMSE                                                                                  ####

org.mimic <- readRDS("../Models/an25-mimic.rds")
all.fitted <- abind("25, unweighted" = org.mimic[-1,,,,],
                    "7, weighted" = fit.1wk, "14, weighted" = fit.2wk, "21, weighted" = fit.3wk, 
                    "28, weighted" = fit.4wk, "35, weighted" = fit.5wk, "42, weighted" = fit.6wk, along = 0)

rmse <- sqrt(apply(sweep(all.fitted[,,,,,1], c(5,2:3), obs[,-1,], "-")^2, c(1, 4:5), mean))

pdf("../Plots/Training-set-size/rmse.pdf", height = 4); {
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    mp.cols <- c("black", "blue", "cyan3", "green2", "gold", "red", "magenta3")
    invisible(sapply(dimnames(rmse)[[3]][1:2], function(varb) {
        matplot(0:14, t(rmse[,,varb]), type = "l", main = varb, xlab = "", ylab = "", lty = 1, col = mp.cols,
                ylim = range(0, rmse))
        legend("bottomright", bty = "n", cex = 0.7, lty = 1, col = mp.cols, dimnames(rmse)[[1]])
    }))
    mtext(expression(paste("RMSE with different training sets for ", eta , " and ", Lambda)), outer = T)
}; dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# determinant sharpness: what do mean/max/min covariance matrices look like?                    ####

ds <- apply(apply(all.fitted[,,,,1:2,2:3], 1:4, det)^0.25, c(1,4), mean)

pdf("../Plots/Training-set-size/det-sharpness.pdf", height = 4); {
    par(mfrow = c(1,1), mar = c(2,2,3,1))
    
    matplot(0:14, t(ds), type = "l", col = mp.cols, lty = 1, xlab = "", ylab = "", 
            main = "Mean determinant sharpness of posterior", ylim = range(0, ds))
    legend("bottomright", bty = "n", cex = 0.7, lty = 1, col = mp.cols, dimnames(rmse)[[3]])
}; dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# shape of Lambda generated by each method                                                      ####
tr.org <- readRDS("../Models/tr-set-default.rds")
el.25 <- abind("tau" = aperm(apply(tr.org, 1:4, mean), c(2:4, 1)), aaply(tr.org, 2:4, function(tr) cov(t(tr))),
            along = 4)[-1,,,,]
            

all.tr <- abind("25, unweighted" = el.25, "7, weighted" = wel.1wk, "14, weighted" = wel.2wk, "21, weighted" = wel.3wk, 
                "28, weighted" = wel.4wk, "35, weighted" = wel.5wk, "42, weighted" = wel.6wk, along = 0)
    
mean.covmat <- apply(all.tr[,,,,-1,], c(1,4:6), mean)

pdf("../Plots/Training-set-size/covariance-matrices.pdf", height = 7); {
    par(mfrow = c(2,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(c(1,6,11,15), function(lt) {
        plot(0, type = "n", xlim = c(-8,8), ylim = c(-8,8), xlab = "", ylab = "",
             main = paste0("Leadtime ", lt-1))
        
        invisible(sapply(1:7, function(i) {
            lines(ellipse(mean.covmat[i,lt,1:2,1:2]), col = mp.cols[i])
        }))
        legend("bottomright", bty = "n", cex = 0.7, lty = 1, col = mp.cols, dimnames(mean.covmat)[[1]])
    }))
    mtext(expression(paste("Mean ", Lambda, " for temp.n and temp.s for each training set")), outer = T)
}; dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# compare verification rank histograms                                                          ####

# use PIT to compare marginal performance over both temps

pit <- aaply(all.fitted[,,,,1:2,1:3], c(1,4), function(mim)
    abind("temp.n" = pnorm(obs["temp.n",-1,], mim[,,"temp.n", "Tau"], sqrt(mim[,,"temp.n", "temp.n"])),
          "temp.s" = pnorm(obs["temp.s",-1,], mim[,,"temp.s", "Tau"], sqrt(mim[,,"temp.s", "temp.s"])),
          along = 0))

lt <- 1
pdf(paste0("../Plots/Training-set-size/PIT-hist-lt", lt-1, ".pdf"), height = 14); {
    par(mfrow = c(7, 2), mar = c(1,1,2,1), oma = c(0,0,2,0))
    invisible(sapply(1:7, function(i) {
        invisible(sapply(dimnames(pit)[[3]], function(varb) {
            hist(pit[i,lt,varb,,], breaks = "fd", prob = T, xaxt = "none", yaxt = "none", 
                 border = mp.cols[i], main = paste0(dimnames(pit)[[1]][i],", ", varb),
                 col = adjustcolor(mp.cols[i], alpha = 0.3))
            abline(h = 1, col = "red3")
        }))
    }))
    mtext(paste0("Marginal PIT histograms for each training set - leadtime ", lt-1), outer = T)
}; dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# discretise & get deviation from uniformity
# treat all as an ensemble of 49
vr <- apply(pit, 1:5, function(p) findInterval(p, (0:50)/50))
ud <- apply(vr, 1:3, u.dev, l = 50)
require(moments)
skew <- apply(vr, 1:3, function(r) skewness(c(r)))

pdf("../Plots/Training-set-size/PIT-nonuniformity.pdf", height = 4); {
par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))

invisible(sapply(dimnames(ud)[[3]], function(varb) {
    matplot(0:14, t(ud[,,varb]), type = "l", lty = 1, col = mp.cols, xlab = "", ylab = "", main = varb,
            ylim = range(0, ud))
    legend("bottomright", bty = "n", cex = 0.7, lty = 1, col = mp.cols, dimnames(ud)[[1]])
}))
mtext("Deviation from uniformity across candidate training sets", outer = T)

invisible(sapply(dimnames(skew)[[3]], function(varb) {
    matplot(0:14, t(skew[,,varb]), type = "l", lty = 1, col = mp.cols, xlab = "", ylab = "", main = varb,
            ylim = c(-0.2,0.2))
    abline(h = 0, col = "grey")
    legend("bottomright", bty = "n", cex = 0.7, lty = 1, col = mp.cols, dimnames(skew)[[1]])
}))
mtext("Skewness across candidate training sets", outer = T)
}; dev.off()



####################################################################################################
####################################################################################################
# unnormalised weights for 75 closest analogues to each
n <- 75; lt <- 15
qq <- aaply(n.dist[lt,,,,], 1:2, function(dd) {
    d.lim <- sort(dd)[n+1]
    an <- which(dd > 0 & dd <= d.lim, arr.ind = T)
    w <- dd[an]^-1
    
    ens <- abind("ec" = apply(aaply(offset.forecast(ecmwf)[,,,lt,-1], 1:3, mean, na.rm = T), 1, "[", an),
                 "nc" = apply(aaply(offset.forecast(ncep)[,,,lt,-1], 1:3, mean, na.rm = T), 1, "[", an),
                 "mo" = apply(aaply(offset.forecast(ukmo)[,,,lt,-1], 1:3, mean, na.rm = T), 1, "[", an),
                 along = 0)
    em <- apply(ens, 2:3, mean)
    
    cw <- cov.wt(em, wt = w)
    abind("tau" = cw$center, cw$cov)
})

setwd("~/Documents/PhD/Miniprojects/03-Sichun-paper/Code")
saveRDS(qq, "./weighted-covmat-tmp.rds")

# unnormalised weights for 14 closest analogues to each: ~ 12m (single leadtime)
# unnormalised weights for 21 closest analogues to each: ~ 12m (single leadtime)

n <- 21; lt <- 15
system.time({qq <- aaply(n.dist[lt,,,,], 1:2, function(dd) {
    d.lim <- sort(dd)[n+1]
    an <- which(dd > 0 & dd <= d.lim, arr.ind = T)
    w <- dd[an]^-1
    
    ens <- abind("ec" = apply(aaply(offset.forecast(ecmwf)[,,,lt,-1], 1:3, mean, na.rm = T), 1, "[", an),
                 "nc" = apply(aaply(offset.forecast(ncep)[,,,lt,-1], 1:3, mean, na.rm = T), 1, "[", an),
                 "mo" = apply(aaply(offset.forecast(ukmo)[,,,lt,-1], 1:3, mean, na.rm = T), 1, "[", an),
                 along = 0)
    em <- apply(ens, 2:3, mean)
    
    cw <- cov.wt(em, wt = w)
    abind("tau" = cw$center, cw$cov)
})})
saveRDS(qq, "./weighted-covmat-21.rds")

# GENERIC PROCESS FOR ANALOGUE SELECTION
#   - get all distances
#   - get all analogue indices (array of weekly batches: 1,2,3,4,5,6 weeks available) & weights

#   - get all Bayesian-inferred errors (simple model without lambda)
#   - find weighted mean & covariance for each batch

wc.14 <- readRDS("./weighted-covmat-14.rds")
wc.21 <- readRDS("./weighted-covmat-21.rds")
wc.75 <- readRDS("./weighted-covmat-75.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# fit model with each value of eta and lambda                                                   ####

fit.model <- function(fc, eta, Lambda) {
    
    d <- dim(fc[[1]])[[1]]
    Y.bar <- abind(lapply(fc, apply, 1, mean), along = 0)
    C <- abind(lapply(lapply(fc, t), cov), along = 0)
    D.i <- aaply(sweep(sweep(C, 1, sapply(fc, dim)[2, ], "/"), 
                       2:3, cov(Y.bar), "+"), 1, solve)
    
    S <- Lambda + solve(apply(D.i, 2:3, sum))
    Tau <- (S %*% ((solve(diag(d) + (apply(D.i, 2:3, sum) %*% 
                                         Lambda))) %*% apply(aaply(abind(Y.bar, D.i, along = 2), 
                                                                   1, function(arr) arr[-1, ] %*% arr[1, ]), 2, sum))) - 
        eta
    
    return(abind(Tau = Tau, S, along = 2))
}

fit.14 <- aaply(abind(offset.forecast(ecmwf)[,-1,,15,-1],
                      offset.forecast(ncep)[,-1,,15,-1],
                      offset.forecast(ukmo)[,-1,,15,-1],
                      aperm(wc.14, c(3,1:2,4)), along = 4), 2:3, function(arr) {
                          fit.model(fc = list("ecmwf" = arr[,1:50],
                                              "ncep" = arr[,51:70],
                                              "ukmo" = arr[,71:93]),
                                    eta = arr[,94],
                                    Lambda = arr[,95:99])
                      })

fit.21 <- aaply(abind(offset.forecast(ecmwf)[,-1,,15,-1],
                      offset.forecast(ncep)[,-1,,15,-1],
                      offset.forecast(ukmo)[,-1,,15,-1],
                      aperm(wc.21, c(3,1:2,4)), along = 4), 2:3, function(arr) {
                          fit.model(fc = list("ecmwf" = arr[,1:50],
                                              "ncep" = arr[,51:70],
                                              "ukmo" = arr[,71:93]),
                                    eta = arr[,94],
                                    Lambda = arr[,95:99])
                      })

fit.75 <- aaply(abind(offset.forecast(ecmwf)[,-1,,15,-1],
                      offset.forecast(ncep)[,-1,,15,-1],
                      offset.forecast(ukmo)[,-1,,15,-1],
                      aperm(wc.75, c(3,1:2,4)), along = 4), 2:3, function(arr) {
                          fit.model(fc = list("ecmwf" = arr[,1:50],
                                              "ncep" = arr[,51:70],
                                              "ukmo" = arr[,71:93]),
                                    eta = arr[,94],
                                    Lambda = arr[,95:99])
                      })

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Assess performance across different training sets                                             ####

rmse <- abind("14" = sqrt(apply(sweep(fit.14[,,,1], c(3,1:2), obs[,-1,], "-")^2, 3, mean)),
              "21" = sqrt(apply(sweep(fit.21[,,,1], c(3,1:2), obs[,-1,], "-")^2, 3, mean)),
              "75" = sqrt(apply(sweep(fit.75[,,,1], c(3,1:2), obs[,-1,], "-")^2, 3, mean)), along = 0)

matplot(t(rmse))

# this seems very wrong. Try plotting & comparing
tr.org <- readRDS("../Models/tr-set-default.rds")

d <- 5; y <- 5; lt <- 15; tr <- t(tr.org[1:2,d+1,y,lt,])
plot(tr, pch = 20, col = "grey", xlim = c(-10,10), ylim = c(-10,10))
points(t(apply(tr, 2, mean)), pch = 20)
lines(ellipse(cov(tr), centre = t(apply(tr, 2, mean))))

lines(ellipse(wc.14[d,y,1:2,2:3], centre = wc.14[d,y,1:2,1]), col = "red")
lines(ellipse(wc.21[d,y,1:2,2:3], centre = wc.21[d,y,1:2,1]), col = "blue")
lines(ellipse(wc.75[d,y,1:2,2:3], centre = wc.75[d,y,1:2,1]), col = "green3")

# suspect that this method requires normalised weights....
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# analogues identified over temperature only                                                    ####
# temperature distances only
t.dist <- invisible(sapply(1:15, function(lt) {
    aaply(cand[,1:2,,,lt], 3:4, function(targ) {
    abs(sweep(cand[,1:2,,,lt], 1:2, targ, "-"))
})
}, simplify = F))

###################################################################################################
# temporary working area
# normalised distances over all variables

mean.t.dist <- aaply(t.dist, c(1:2, 5:6), mean, na.rm = T)

q.t.dist <- apply(m.dist, 1:2, quantile, 0.1, na.rm = T)
td.small <- apply(m.dist, 1:2, function(d) sum(d < 0.75, na.rm = T))

plot(sort(td.small))
min(td.small[td.small > 0])

dist <- aaply(cand[,1:2,,], 3:4, function(targ) {
    
    # standard deviation for normalisation
    cand.sd <- apply(cand, 1:2, sd, na.rm = T)
    
    # normalised distance from each point on target to equivalent in candidate space
    sweep(sqrt(sweep(cand, 1:2, targ, "-")^2), 1:2, cand.sd, "/")
})

m.dist <- aaply(dist, c(1:2, 5:6), mean, na.rm = T)
max.dist <- aaply(dist, c(1:2, 5:6), max, na.rm = T)

d <- 10; y <- 5
plot(sort(m.dist[d,y,,]))
hist(m.dist[d,y,,], breaks = "fd")

plot(sort(max.dist[d,y,,]))
q.1 <- apply(m.dist, 1:2, quantile, 0.1, na.rm = T)

d.small <- apply(m.dist, 1:2, function(d) sum(d < 0.75, na.rm = T))

plot(sort(d.small))
min(d.small[d.small > 0])
sort(d.small)[1:15]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
