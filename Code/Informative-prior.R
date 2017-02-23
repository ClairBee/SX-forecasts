
library("SX.weather")

# train models over analogues & historic data for reference
# now use informative prior

# load mimic & get RMSE for comparison
mimic <- readRDS("../Models/an25-mimic.rds")
mimic.rmse <- sqrt(apply(sweep(mimic[,,,,1], c(4,1:2), obs, "-")^2, c(3:4), mean, na.rm = T))

####################################################################################################

# REVISED FUNCTION                                                                              ####

inf.mimic <- function(o, fc, tr, prior, arr.out = T) {
        
    d <- length(o)
    
        Y.bar <- abind(lapply(fc, apply, 1, mean), along = 0)
        C <- abind(lapply(lapply(fc, t), cov), along = 0)
        
        D.i <- aaply(sweep(sweep(C, 1, sapply(fc, dim)[2,], "/"), 2:3, cov(Y.bar), "+"), 1, solve)
        
        Eta <- apply(t(tr), 2, mean)
        Lambda <- cov(t(tr))
        
        if (missing(prior)) {
            S <- Lambda + solve(apply(D.i, 2:3, sum))
            Tau <- (S %*% 
                        ((solve(diag(d) + (apply(D.i, 2:3, sum) %*% Lambda))) %*%
                             apply(aaply(abind(Y.bar, D.i, along = 2), 1, 
                                         function(arr) arr[-1,] %*% arr[1,]), 2, sum))) - Eta
        } else {
            S <- solve(solve(prior$Gamma) + solve(Lambda + solve(apply(D.i, 2:3, sum))))
            Tau <- (S %*% 
                        ((solve(prior$Gamma) %*% (prior$alpha - Eta)) + 
                             (solve(diag(d) + (apply(D.i, 2:3, sum) %*% Lambda))) %*%
                             apply(aaply(abind(Y.bar, D.i, along = 2), 1, 
                                         function(arr) arr[-1,] %*% arr[1,]), 2, sum))) - Eta
        }
       
        if (arr.out) {
            # return results in array form (more useful for array operations)
            return(abind("Tau" = Tau, S, along = 2))
        } else {
            # otherwise, return output as list
            return(list(Tau = Tau, S = S))
        }
}

tr.error <- readRDS("../Models/tr-set-default.rds")

mimic[5,5,5,,]
#  Tau     temp.n     temp.s         pc1         pc2         pc3
#  temp.n  1.6001802 1.37023176 0.82216801 0.056409491 0.034205879 0.103422174
#  temp.s  4.4146433 0.82216801 1.20883957 0.046451564 0.060412545 0.109585777
#  pc1     0.8364430 0.05640949 0.04645156 0.039395010 0.002889492 0.006275338
#  pc2     0.7422758 0.03420588 0.06041255 0.002889492 0.028904276 0.007628357
#  pc3    -2.0256275 0.10342217 0.10958578 0.006275338 0.007628357 0.072004068

inf.mimic(o = obs[,5,5],
          fc = list("ecmwf" = offset.forecast(ecmwf)[,5,5,5,-1],
                    "ncep" = offset.forecast(ncep)[,5,5,5,-1],
                    "ukmo" = offset.forecast(ukmo)[,5,5,5,-1]),
          tr = tr.error[,5,5,5,])

# how to set prior on actual temperature?

#--------------------------------------------------------------------------------

persistence.prior <- list(alpha = obs[,4,5], Gamma = cov(apply(obs, 1, cbind)[364 -(50:1),]))
# use yesterday's observations & covariance of preceding 50 days


persistence.mimic <- inf.mimic(o = obs[,5,5],
                               fc = list("ecmwf" = offset.forecast(ecmwf)[,5,5,5,-1],
                                         "ncep" = offset.forecast(ncep)[,5,5,5,-1],
                                         "ukmo" = offset.forecast(ukmo)[,5,5,5,-1]),
                               tr = tr.error[,5,5,5,],
                               prior = persistence.prior)

#--------------------------------------------------------------------------------

# --> Persistence prior                                                             ####

o.seq <- apply(obs, 1, cbind)

p.data <- abind(invisible(sapply(1:15, function(lt) {
    abind("o" = o.seq[26:630,],
          "ec" = apply(offset.forecast(ecmwf)[,,,lt,-1], c(1,4), cbind)[26:630,,],
          "nc" = apply(offset.forecast(ncep)[,,,lt,-1], c(1,4), cbind)[26:630,,],
          "mo" = apply(offset.forecast(ukmo)[,,,lt,-1], c(1,4), cbind)[26:630,,],
          "tr" = apply(tr.error[,,,lt,], c(1,4), cbind)[26:630,,],
          "alpha" = o.seq[25:629,],
          "Gamma" = abind(invisible(sapply(26:630, 
                                           function (i) cov(o.seq[i - (25:1),]), simplify = F)), along = 0),
          along = 3)}, simplify = F)), along = 0)

persistence.mimic <- aaply(p.data, 1:2, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,2:51],
                        "ncep" = arr[,52:72],
                        "ukmo" = arr[,73:94]),
              tr = arr[,95:119],
              prior = list(alpha = arr[,120],
                           Gamma = arr[,121:125]))
})

#--------------------------------------------------------------------------------

# --> Analogue prior                                                                ####

tr.obs <- aaply(an.indx, 1:3, function(a) {
    apply(obs, 1, "[", a)
})

# alpha & Gamma taken from analogue observations (NOT errors)

a.data <- abind(invisible(sapply(1:15, function(lt) {
    abind("ec" = offset.forecast(ecmwf)[,,,lt,-1],
          "nc" = offset.forecast(ncep)[,,,lt,-1],
          "mo" = offset.forecast(ukmo)[,,,lt,-1], 
          "tr" = apply(tr.error[,,lt,,,], c(5, 1:3), mean), 
          "alpha" = apply(tr.obs[,,lt,,], c(4, 1:2), mean),
          "Gamma" = aperm(aaply(tr.obs[,,lt,,], c(1:2), cov), c(3, 1:2, 4)),
          along = 4)}, simplify = F)), along = 0)

analogue.mimic <- aaply(a.data, c(1,3:4), function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

#--------------------------------------------------------------------------------

# --> Offset persistence prior                                                      ####

o.seq <- apply(obs, 1, cbind)

os.p.data <- abind(invisible(sapply(1:15, function(lt) {
    abind("o" = o.seq[26:630,],
          "ec" = apply(offset.forecast(ecmwf)[,,,lt,-1], c(1,4), cbind)[26:630,,],
          "nc" = apply(offset.forecast(ncep)[,,,lt,-1], c(1,4), cbind)[26:630,,],
          "mo" = apply(offset.forecast(ukmo)[,,,lt,-1], c(1,4), cbind)[26:630,,],
          "tr" = apply(tr.error[,,,lt,], c(1,4), cbind)[26:630,,],
          "alpha" = o.seq[25:629-(lt+1),],
          "Gamma" = abind(invisible(sapply(26:630, 
                                           function (i) cov(o.seq[i - (25:1),]), simplify = F)), along = 0),
          along = 3)}, simplify = F)), along = 0)

os.persistence.mimic <- aaply(os.p.data, 1:2, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,2:51],
                        "ncep" = arr[,52:72],
                        "ukmo" = arr[,73:94]),
              tr = arr[,95:119],
              prior = list(alpha = arr[,120],
                           Gamma = arr[,121:125]))
})

#--------------------------------------------------------------------------------

# --> quick performance check: RMSE                                                 ####

persistence.rmse <- sqrt(apply(sweep(persistence.mimic[,,,"Tau"], 2:3, o.seq[26:630,], "-")^2, c(1,3), mean, na.rm = T))
analogue.rmse <- sqrt(apply(sweep(analogue.mimic[,,,,"Tau"], c(4,2:3), obs, "-")^2, c(1,4), mean, na.rm = T))
os.persistence.rmse <- sqrt(apply(sweep(os.persistence.mimic[,,,"Tau"], 2:3, o.seq[26:630,], "-")^2, c(1,3), mean, na.rm = T))


pdf("../Plots/Analogue-perf/Priors.pdf", height = 4, width = 9); {
    
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    plot(mimic.rmse[,"temp.n"], type = "l", ylim = range(mimic.rmse, 0, persistence.rmse), 
         xlab = "", ylab = "", main = "temp.n")
    lines(analogue.rmse[,"temp.n"], col = "blue")
    lines(persistence.rmse[,"temp.n"], col = "red")
    abline(v = 7, col = "cyan3", lty = 2)
    
    legend("bottomright", col = c("black", "red", "blue"), lty = 1, bty = "n",
           legend = c("Uninformative", "Persistence", "Analogues"))
    
    
    plot(mimic.rmse[,"temp.s"], type = "l", ylim = range(mimic.rmse, 0, persistence.rmse), 
         xlab = "", ylab = "", main = "temp.s")
    lines(analogue.rmse[,"temp.s"], col = "blue")
    lines(persistence.rmse[,"temp.s"], col = "red")
    abline(v = 7, col = "cyan3", lty = 2)
    
    legend("bottomright", col = c("black", "red", "blue"), lty = 1, bty = "n",
           legend = c("Uninformative", "Persistence", "Analogues"))
    
    mtext("RMSE of models with different priors", outer = T)
    
}; dev.off()

#--------------------------------------------------------------------------------

# --> Compare to BMA/MOS as well....                                                ####

q.fc <- readRDS("../Models/an25-bma-mos-fc.rds")

rmse <- abind("uninformative" = mimic.rmse[,1:2],
#            "persistence" = persistence.rmse[,1:2],
            "persistence" = os.persistence.rmse[,1:2],
            "analogue" = analogue.rmse[,1:2],
            aperm(sqrt(apply(sweep(q.fc, c(5, 1:2), obs[1:2,,], "-")^2, 3:5, mean, na.rm = T)), c(2,1,3)),
            along = 1)

m.col <- c("black","blue", "red", "green3", "purple"); m.lty <- c(1,1,1,2,2)

pdf("../Plots/Analogue-perf/Priors.pdf", height = 4, width = 9); {
    
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(dimnames(rmse)[[3]], function(varb) {
        matplot(t(rmse[,,varb]), type = "l", main = varb, col = m.col, lty = m.lty)
        legend("bottomright", bty = "n", legend = dimnames(rmse)[[1]], col = m.col, lty = m.lty, cex = 0.7)
    }))
    
    mtext("RMSE of models with different priors", outer = T)
    
}; dev.off()

####################################################################################################

# 'BUCKET' APPROACH                                                                             ####

# add least sophisticated methods for comparison, eg. mean of means

super.fc <- apply(abind("ecmwf" = apply(offset.forecast(ecmwf)[,,,,-1], 1:4, mean),
                        "ncep" = apply(offset.forecast(ncep)[,,,,-1], 1:4, mean),
                        "ukmo" = apply(offset.forecast(ukmo)[,,,,-1], 1:4, mean),
                        along = 0), 2:5, mean)

bucket.fc <- apply(offset.forecast(superensemble())[,,,,-(1:3)], 1:4, mean)

rmse <- abind("mimic" = mimic.rmse, 
              "super" = t(sqrt(apply(sweep(super.fc, 1:3, obs, "-")^2, c(1,4), mean))),
              "bucket" = t(sqrt(apply(sweep(bucket.fc, 1:3, obs, "-")^2, c(1,4), mean))),
              along = 0)

par(mfrow = c(2,1), mar = c(2,2,3,1), oma = c(0,0,2,0))
invisible(sapply(dimnames(rmse)[[3]][1:2], function(varb) {
    matplot(0:14, t(rmse[,,varb]), type = "l", xlab = "", ylab = "", main = varb,
            col = c("black", "blue", "green3"), lty = c(1,2,2), ylim = range(rmse))
    legend("bottomright", bty = "n", lty = c(1,2,2), col = c("black", "blue", "green3"),
           legend = dimnames(rmse)[[1]])
}))
mtext("RMSE", outer = T)

####################################################################################################

# COMPARE ERRORS BY BIAS                                                                        ####

err <- abind("mimic" = sweep(mimic[,,,1:2,"Tau"], c(4,1:2), obs[1:2,,], "-"),
                 sweep(q.fc, c(5, 1:2), obs[1:2,,], "-"),
             along = 4)

qq <- apply(err, 3:5, function(e) mean(e > 0, na.rm = T))

rmse.pv <- sqrt(apply((err * (err > 0))^2, 3:5, mean, na.rm = T))

pdf("../Plots/Analogue-perf/Error-bias.pdf"); {
    par(mar = c(2,2,3,1))
    
    matplot(0:14, rmse.pv[,,1] * qq[,,1], type = "l", col = c("black", "red", "blue"), lty = 1,
            xlab = "", ylab = "", main = "RMSE of positive errors, weighted by proportion of errors")
    matplot(0:14, rmse.pv[,,2] * qq[,,2], type = "l", col = c("black", "red", "blue"), lty = 2, add = T)
    
    legend("bottomright", col = rep(c("black", "red", "blue"), 2), lty = rep(c(1,2), each = 3), cex = .7,
           legend = paste(dimnames(rmse.pv)[[2]], rep(c("temp.n", "temp.s"), each = 3)), bty = "n")
}; dev.off()

#---------------------------------------------------------------------------------

# bias of each forecast vs bias of ensemble BMA/MOS

q.fc <- readRDS("../Models/an25-bma-mos-fc.rds")

err <- sweep(abind("ecmwf" = apply(offset.forecast(ecmwf)[1:2,,,,-1], 1:4, mean),
             "ncep" = apply(offset.forecast(ncep)[1:2,,,,-1], 1:4, mean),
             "ukmo" = apply(offset.forecast(ukmo)[1:2,,,,-1], 1:4, mean),
             "bma" = aperm(q.fc[,,,"bma",], c(4,1:3)),
             "mos" = aperm(q.fc[,,,"mos",], c(4,1:3)),
             "mimic" = aperm(mimic[,,,1:2,"Tau"], c(4,1:3)), 
             along = 0), 2:4, obs[1:2,,], "-")

m.err <- apply(err, 1:2, mean, na.rm = T)

plot(err["ecmwf","temp.n",,,], err["bma","temp.n",,,], pch = ".", xlab = "ecmwf", ylab = "bma",
     col = adjustcolor("steelblue", alpha = 0.5))
points(err["ncep","temp.n",,,], err["bma","temp.n",,,], pch = ".",
       col = adjustcolor("red3", alpha = 0.5))
points(err["ukmo","temp.n",,,], err["bma","temp.n",,,], pch = ".", 
       col = adjustcolor("green3", alpha = 0.5))
points(err["mimic","temp.n",,,], err["bma","temp.n",,,], pch = ".")

abline(h = m.err["bma","temp.n"], lty = 2)
abline(v = m.err[1:3,"temp.n"], lty = 2, col = c("steelblue", "red3", "green3"))
abline(v = m.err["mimic", "temp.n"], lty = 2)

####################################################################################################

# ANALOGUES AS SUPPLEMENTARY ENSEMBLE, RATHER THAN INFORMATIVE PRIOR                            ####

tr.fc <- abind(sapply(1:15, function(lt) {
    aaply(an.indx, 1:2, function(a) {
        apply(em[,,,,lt], 1:2, "[", a[lt,,])
    })}, simplify = F), along = 2.5)

an.ens.data <- abind("ec" = offset.forecast(ecmwf)[,,,,-1],
                "nc" = offset.forecast(ncep)[,,,,-1],
                "mo" = offset.forecast(ukmo)[,,,,-1], 
                "an" = apply(tr.fc, c(6, 1:4), mean),
                "tr" = apply(tr.error, c(6, 1:4), mean),
                rev.along = 1)

an.ens.mimic <- aaply(an.ens.data, c(2:4), function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,2:51],
                        "ncep" = arr[,52:72],
                        "ukmo" = arr[,73:94],
                        "an" = arr[,95:119]),
              tr = arr[,120:143])
})

rmse <- sqrt(apply(sweep(abind("No prior" = mimic,
              "Analogue prior" = aperm(analogue.mimic, c(2:3, 1, 4:5)),
              "Analogue ensemble" = an.ens.mimic,
              along = 0)[,,,,1:2,"Tau"], c(5,2:3), obs[1:2,,], "-")^2, c(5,4,1), mean, na.rm = T))

pdf("../Plots/Analogue-perf/Analogue-application.pdf", height = 4, width = 7); {
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(dimnames(rmse)[[1]], function(varb) {
        matplot(0:14, rmse[varb,,], type = "l", lty = c(1,2,2), col = c("black", "blue", "red"),
                main = varb, ylim = range(0,rmse))
        legend("bottomright", lty = c(1,2,2), col = c("black", "blue", "red"), bty = "n", cex = 0.7,
               legend = dimnames(rmse)[[3]])
    }))
    mtext("Use of analogues as synthetic ensemble vs informative prior", outer = T)
}; dev.off()

