
library("CB.Misc"); library("SX.weather")

# STILL TO DO...

#   - don't forget CompStat homework still to do (independence sampler)

#   - then look at superensemble & comparitive plots - once everything else is fully aligned
#   - rework the forecast error calculation to use offset-corrected forecast. Probably easier


# calculation used is for RMSE of ensemble mean - not mean RMSE of ensemble
# find some references on this - what is used in published papers?

####################################################################################################

# SUPERENSEMBLE PLOTS                                                                           ####

# do we need to distinguish between controls and perturbations?

super <- superensemble(list(ecmwf, ncep, ukmo))

s.err <- forecast.errors(super)
s.mean.error <- apply(s.err, c(1, 4, 5), mean)

ens.cols <- c("blue", "red", "green3")
pdf("./Plots/Superensemble mean error.pdf"); {
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(s.mean.error)[[1]][c(3:5, 1:2)], function(varb) {
        
        # plot all ensemble members
        matplot(s.mean.error[varb, ,], type = "l", lty = 1,
                col = c(rep(NA, 3), adjustcolor(ens.cols[attributes(super)$m], alpha = 0.1)),
                xlab = "", ylab = "", main = varb)
        
        # ensemble means
        invisible(sapply(1:3, function(n) lines(apply(s.mean.error[varb,,attributes(super)$m == n], 1, mean),
                                                col = ens.cols[n], lty = 2)))
        
        # controls
        matplot(s.mean.error[varb,,1:3], lwd = 1, col = ens.cols, add = T, type = "l", lty = 1)
        
        # overall ens. mean
        lines(apply(s.err[varb,,,,-(1:3)], 3, mean), lty = 2, lwd = 2)
        
        # overall control mean
        lines(apply(s.err[varb,,,,(1:3)], 3, mean), lty = 1, lwd = 2)
        
        abline(h = 0, col = "darkgrey", lty = 3)
    }))
    
    plot.new()
    legend("left", lty = c(1,1,1,1,NA,1,2), lwd = c(1,1,1,1,NA,2,2), col = c(ens.cols, rep("black", 4)), bty = "n",
           legend = c("ECMWF", "NCEP", "UKMO", "Superensemble", NA, "Control forecast", "Perturbations"))
    
    mtext("Superensemble - mean forecast error at each leadtime", outer = T)

}; dev.off()

#----------------------------------------------------------------------------------------

# plot of RMSE for each forecast
ctrl.mean.rmse <- sqrt(apply(s.err[,,,,1:3]^2, c(1,4), mean))
pert.mean.rmse <- sqrt(apply(apply(s.err[,,,,-(1:3)], 1:4, mean)^2, c(1,4), mean))

member.rmse <- forecast.rmse(super)
ens.means <- abind(invisible(sapply(1:3, 
                                    function(n) {
                                        sqrt(apply(apply(s.err[,,,,attributes(super)$m == n], 1:4, mean)^2, c(1,4), mean))
                                    }, simplify = F)),
                   along = 3)

pdf("./Plots/Superensemble RMSE.pdf"); {
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(dimnames(s.mean.error)[[1]][c(3:5, 1:2)], function(varb) {
        
        # plot all ensemble member rmses
        matplot(member.rmse[varb, ,], type = "l", lty = 1,
                col = c(rep(NA, 3), adjustcolor(ens.cols[attributes(super)$m], alpha = 0.1)),
                xlab = "", ylab = "", main = varb)
        
        # individual control rmse
        matplot(member.rmse[varb,,1:3], lwd = 1, col = ens.cols, add = T, type = "l", lty = 1)
        
        # individual ensemble rmse
        matplot(ens.means[varb,,], type = "l", add = T, lwd = 1, lty = 2, col = ens.cols)
        
        # overall ens. mean
        lines(pert.mean.rmse[varb,], lty = 2, lwd = 2)
        
        # overall control mean
        lines(ctrl.mean.rmse[varb,], lty = 1, lwd = 2)
        
        abline(h = 0, col = "darkgrey", lty = 3)
    }))
    
    plot.new()
    legend("left", lty = c(1,1,1,1,NA,1,2), lwd = c(1,1,1,1,NA,2,2), col = c(ens.cols, rep("black", 4)), bty = "n",
           legend = c("ECMWF", "NCEP", "UKMO", "Superensemble", NA, "Control forecast", "Perturbations"))
    
    mtext("Superensemble - RMSE at each leadtime", outer = T)
    
}; dev.off()

#----------------------------------------------------------------------------------------

# spread
ens.spread <- ensemble.spread(super[,,,,-(1:3)])

member.spread <- abind(invisible(sapply(1:3,
                                        function(n) ensemble.spread(super[,,,,attributes(super)$m == n]),
                                        simplify = F)),
                       along = 3)

pdf("./Plots/Superensemble RMSE vs spread.pdf"); {
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(dimnames(s.mean.error)[[1]][c(3:5, 1:2)], function(varb) {
        
        ens.cols <- c("goldenrod1", "chartreuse3", "cyan3")
        matplot(abind(ens.means, member.spread, pert.mean.rmse, ens.spread, along = 3)[varb,,],
                col = c(rep(ens.cols, 2), rep("black", 2)),
                lty = c(rep(1, 3), rep(2, 3), 1, 2),
                lwd = c(rep(1, 6), rep(2, 2)),
                type = "l", main = varb, xlab = "", ylab = "")
    }))
    #type 2 = spread, tye 1 = rmse
    
    plot.new()
    legend("left", lty = c(1,1,1,1,NA,1,2), lwd = c(1,1,1,1,NA,2,2), col = c(ens.cols, rep("black", 4)), bty = "n",
           legend = c("ECMWF", "NCEP", "UKMO", "Superensemble", NA, "RMSE of ensemble mean", "Ensemble spread"))
    
    mtext("Superensemble - RMSE vs spread at each leadtime", outer = T)
}; dev.off()
    
####################################################################################################

# COVARIANCES                                                                                   ####

ecmwf.cov <- se.covariances(offset.forecast(ecmwf))
ncep.cov <- se.covariances(offset.forecast(ncep))
ukmo.cov <- se.covariances(offset.forecast(ukmo))

c.sigma <- se.sigma()

c.lambda <- se.lambda()

c.eta <- apply(apply(forecast.errors(superensemble())[,,,,-(1:3)], 1:4, mean), c(1,2,4), mean)

y.bar <- abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
               "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
               "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
               along = 0)

d <- abind("ecmwf" = ecmwf.cov / (dim(ecmwf)[[5]] - 1) + c.sigma,
           "ncep" = ncep.cov / (dim(ncep)[[5]] - 1) + c.sigma,
           "ukmo" = ukmo.cov / (dim(ukmo)[[5]] - 1) + c.sigma,
           along = 0)

d.inv.sum <-  apply(array.solve(d, c(1, 4:6)), c(1:2, 4:6), sum)

c.precision <- array.solve(sweep(array.solve(d.inv.sum, c(3:5)),      # invert d.sum
                                 c(1:3, 5), c.lambda, "+"),           # add lambda
                           3:5)                                       # invert total

s <- array.solve(c.precision, 3:5)

c.tau <- se.tau(D = d, L = c.lambda, Y = y.bar, S = s, Eta = c.eta)

####################################################################################################

# RESULTS - METHOD 1                                                                            ####

# RMSE & spread, CRPS

spr <- sqrt(apply(apply(s, 3:5, diag), c(1,4), mean))
rmse <- se.rmse(c.tau)

c.crps <- array(crps(c(rep(obs, dim(c.tau)[[4]])),
               cbind(c(c.tau),
                     sqrt(apply(s, 3:5, diag))))$crps,
               dim = dim(c.tau), dimnames = dimnames(c.tau))

crps.mean <- apply(c.crps, c(1,4), mean)

pdf("./Plots/Superensemble spread & RMSE.pdf"); {
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(rmse)[[1]][c(3:5, 1:2)], function(varb) {
        matplot(cbind(rmse[varb,], spr[varb,]), ylim = range(0, rmse[varb,], spr[varb,]),
                type = "l", xlab = "", ylab = "", col = c("darkred", "darkblue"), main = varb)
    }))
    
    plot.new()
    legend("left", lty = c(1,2), lwd = 1, col = c("darkred", "darkblue"), bty = "n",
           legend = c("Superensemble RMSE", "Superensemble spread"))
    
    mtext("Superensemble - RMSE and spread at each leadtime", outer = T)
}; dev.off()

####################################################################################################

# COVARIANCE - METHOD 2                                                                         ####

# calculates covariance using previous 30 days, treating all observations as consecutive.
# Observations from Feb to December are not consecutive, so this seems flawed...
# (could conceivably work if run over each year separately, but would massively reduce savailale data)

####################################################################################################

# COVARIANCE - TEMP ONLY                                                                        ####

ecmwf.cov <- se.covariances(offset.forecast(ecmwf))
ncep.cov <- se.covariances(offset.forecast(ncep))
ukmo.cov <- se.covariances(offset.forecast(ukmo))

c.sigma <- se.sigma()

c.lambda <- se.lambda()

ens.mean.error <- apply(forecast.errors(superensemble())[,,,,-(1:3)], 1:4, mean)
c.eta <- apply(ens.mean.error, c(1,2,4), mean)

y.bar <- abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
               "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
               "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
               along = 0)

d <- abind("ecmwf" = ecmwf.cov / (dim(ecmwf)[[5]] - 1) + c.sigma,
           "ncep" = ncep.cov / (dim(ncep)[[5]] - 1) + c.sigma,
           "ukmo" = ukmo.cov / (dim(ukmo)[[5]] - 1) + c.sigma,
           along = 0)[,1:2,1:2,,,2:15]

c.precision <- array.solve(sweep(array.solve(apply(array.solve(d,                   # truncate to match
                                                               c(1, 4:6)),          # solve per model/day/lt
                                                   c(1:2, 4:6), sum),               # sum over models
                                             c(3:5)),                               # solve for sum
                                 c(1:3, 5), c.lambda[1:2,1:2,,2:15], "+"),          # add lambda to solution
                           3:5)                                                     # solve total

s <- array.solve(c.precision, 3:5)

c.tau <- se.tau(D = d, L = c.lambda[1:2,1:2,,2:15],
                Y = y.bar[,1:2,,,2:15], S = s, Eta = c.eta[1:2,,2:15])

spr <- sqrt(apply(apply(s, 3:5, diag), c(1,4), mean))
rmse <- se.rmse(c.tau, actual = obs[1:2,,])

c.crps <- array(crps(c(rep(obs[1:2,,], 14)),
                     cbind(c(c.tau),
                           sqrt(apply(s, 3:5, diag))))$crps, dim = dim(c.tau))
crps.mean <- apply(c.crps, c(1,4), mean)

####################################################################################################

# WORKING AREA                                                                                  ####

pdf("./Plots/Superensemble spread & RMSE - temp only.pdf"); {
    par(mfrow = c(1,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(rmse)[[1]], function(varb) {
        matplot(cbind(rmse[varb,], spr[varb,]), ylim = range(0, rmse[varb,], spr[varb,]),
                type = "l", xlab = "", ylab = "", col = c("darkred", "darkblue"), main = varb)
    }))
    
    plot.new()
    legend("left", lty = c(1,2), lwd = 1, col = c("darkred", "darkblue"), bty = "n",
           legend = c("Superensemble RMSE", "Superensemble spread"))
    
    mtext("Superensemble - RMSE and spread at each leadtime", outer = T)
}; dev.off()