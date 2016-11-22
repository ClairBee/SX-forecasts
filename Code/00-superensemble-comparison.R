
library("CB.Misc"); library("SX.weather")

# comparison of models with 0, 2, 3 principal components

    # tidy up code,, upload to GitHub.
    # further metrics: histograms etc
    # look at varying elements: what actually changes? (only S/precision, presumably)
    # write something up, including comparison of S & Tau for each model (0,2,3 PCs)
    # plot the calculated densities?

# why is perturbation for PC2 done at different leadtime to PC1? (eg. see code for 05-crps.R)

# try treating superensemble as superensemble, and fitting a single model. Difference?

####################################################################################################

# CREATE MODELS                                                                                 ####

# common matrices
{
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
               along = 0)
}

# full model: all 3 principal components
{
    precision.3 <- se.precision(d, c.lambda)
    s.3 <- array.solve(precision.3, 3:5)
    tau.3 <- se.tau(d, c.lambda, y.bar, s.3, c.eta)
    
    spr.3 <- sqrt(apply(apply(s.3, 3:5, diag), c(1,4), mean))
    rmse.3 <- se.rmse(tau.3)
    crps.3 <- se.crps(tau.3, s.3)
}

# Sichun model 1: first two principal components
{
    precision.2 <- se.precision(d, c.lambda, vars = c(1:4))
    s.2 <- array.solve(precision.2, 3:5)
    tau.2 <- se.tau(d, c.lambda, y.bar, s.2, c.eta)
    
    spr.2 <- sqrt(apply(apply(s.2, 3:5, diag), c(1,4), mean))
    rmse.2 <- se.rmse(tau.2)
    crps.2 <- se.crps(tau.2, s.2)
}

# Sichun's final model: no principal components
{
    precision.0 <- se.precision(d, c.lambda, vars = c(1:2))
    s.0 <- array.solve(precision.0, 3:5)
    tau.0 <- se.tau(d, c.lambda, y.bar, s.0, c.eta)
    
    spr.0 <- sqrt(apply(apply(s.0, 3:5, diag), c(1,4), mean))
    rmse.0 <- se.rmse(tau.0)
    crps.0 <- se.crps(tau.0, s.0)
}

####################################################################################################

# RANGE OF VALUES IDENTIFIED FOR TAU                                                            ####

# range decreases as leadtime increases. Seems odd.
tau.rng <- abind("max" = apply(tau.3, c(1,4), max),
                 "min" = apply(tau.3, c(1,4), min),
                 "q.5" = apply(tau.3, c(1,4), quantile, 0.05),
                 "q.95" = apply(tau.3, c(1,4), quantile, 0.95),
                 "mean" = apply(tau.3, c(1,4), mean), along = 0)

pdf("./Plots/Superensemble Tau.pdf"); {
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(tau.rng)[[2]][c(3:5, 1:2)], function(varb) {
        plot(0,type = "n", xlim = c(0,14), ylim = range(c(tau.rng, 0)), xlab = "", ylab = "", main = varb)
        polygon(c(0:14, 14:0), c(tau.rng["max",varb,], rev(tau.rng["min",varb,])),
                col = adjustcolor("red3", alpha = 0.1), border = NA)
        polygon(c(0:14, 14:0), c(tau.rng["q.95",varb,], rev(tau.rng["q.5",varb,])),
                col = adjustcolor("red3", alpha = 0.3), border = NA)
        lines(0:14, tau.rng["mean",varb,], col = "red3", lwd = 2)
        abline(0,0, lty = 3)
    }))
    plot.new()
    legend("center", lwd = c(NA, NA, 2), pch = c(22,22,NA), col = c(NA, NA, "red3"), bty = "n",
           pt.bg = c(adjustcolor("red3", alpha = 0.1), adjustcolor("red3", alpha = 0.3), NA), pt.cex = 2,
           legend = c("Full range", "Q.05 - Q.95", "Mean"))
    mtext(expression(paste("Range of values for ", tau, " at each leadtime")), outer = T)
}; dev.off()


####################################################################################################

# RMSE PLOTS FOR EACH VARIABLE SET                                                              ####

# bind into single array for easier plotting

res <- array(NA, dim = c(5,15,2,3), 
             dimnames = append(dimnames(spr.3), list(c("RMSE", "spread"), c("0", "2", "3"))))
res[,,,"3"] <- abind(rmse.3, spr.3, along = 3)
res[1:4,,,"2"] <- abind(rmse.2, spr.2, along = 3)
res[1:2,,,"0"] <- abind(rmse.0, spr.0, along = 3)

pdf("./Plots/Superensemble 'smart' spread & RMSE.pdf"); {
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(res)[[1]][c(3:5, 1:2)], function(varb) {
        matplot(array(res[varb,,,], dim = c(15, 6)), type = "l", cex = 1.5, lty = rep(c(1,3),3),
                col = rep(c("black", "red", "orange"), each = 2), lwd = rep(c(3,2,1), each = 2),
                pch = rep(c(15,20,3), each = 2), xlab = "", ylab = "", main = varb,
                ylim = range(c(0,res[varb,,,]), na.rm = T))
    }))
    
    plot.new()
    legend("center", lty = rep(c(1,2),3), lwd = rep(c(3,2,1), each = 2), 
           col = rep(c("black", "red", "orange"), each = 2), bty = "n",
           legend = c("RMSE - temp only", "Spread - temp only",
                      "RMSE - first 2 PCs", "Spread - first 2 PCs",
                      "RMSE - first 3 PCs", "Spread - first 3 PCs"))
    
    mtext("'Smart' superensemble - RMSE and spread at each leadtime", outer = T)
}; dev.off()

crps.res <- array(NA, dim = c(5,15,3), 
                  dimnames = append(dimnames(crps.3[,1,1,]), list(c("0", "2", "3"))))
crps.res[,,"3"] <- apply(crps.3, c(1,4), mean)
crps.res[1:4,,"2"] <- apply(crps.2, c(1,4), mean)
crps.res[1:2,,"0"] <- apply(crps.0, c(1,4), mean)

pdf("./Plots/Superensemble 'smart' CRPS.pdf"); {
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(res)[[1]][c(3:5, 1:2)], function(varb) {
        
        matplot(crps.res[varb,,], type = "l", lty = 1, col = c("black", "red", "orange"),
                lwd = c(3,2,1), xlab = "", ylab = "", main = varb,
                ylim = range(c(0,crps.res[varb,,]), na.rm = T))
    }))
    
    plot.new()
    legend("center", lty = 1, lwd = c(3,2,1), col = c("black", "red", "orange"), bty = "n",
           legend = c("CRPS - temp only", "CRPS - first 2 PCs", "CRPS - first 3 PCs"))
    
    mtext("'Smart' superensemble - CRPS at each leadtime", outer = T)
}; dev.off()

####################################################################################################

# PERFORMANCE AGAINST INDIVIDUAL MODELS                                                         ####

# run full model (3 principal components)
model.3pc <- run.model()
model.2pc <- run.model(c("temp.n", "temp.s", "pc1", "pc2"))
model.0pc <- run.model(c("temp.n", "temp.s"))

model.3pc.res <- model.performance(model.3pc)
model.2pc.res <- model.performance(model.2pc)
model.0pc.res <- model.performance(model.0pc)

crps.mat <- abind("ecmwf" = ensemble.crps(ecmwf[,,,,-1]),
                  "ncep" = ensemble.crps(ncep[,,,,-1]),
                  "ukmo" = ensemble.crps(ukmo[,,,,-1]), 
                  "super" = ensemble.crps(superensemble()[,,,,-(1:3)]),
                  "no.pc" = rbind(apply(model.0pc.res$crps, c(1,4), mean), NA, NA, NA),     # padded
                  "two.pc" = rbind(apply(model.2pc.res$crps, c(1,4), mean), NA),            # padded
                  "all.pc" = apply(model.3pc.res$crps, c(1,4), mean),
                  along = 0)

# check against figures plotted by Sichun - this matches
all(round(crps.mat["two.pc",2:5,2:15],9) == round(CRPS_Richard_average, 9))

# reproduce fig 5.7: CRPS for all variables, for all ensemble combinations
pdf("./Plots/CRPS comparison for all models.pdf"); {
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(crps.mat)[[2]][c(3:5,1:2)], function(varb) {
        matplot(t(crps.mat[,varb,]), type = "l", main = varb, xlab = "", ylab = "",
                lty = c(rep(1, 3), 2, rep(1,3)), lwd = c(rep(1,3), rep(2, 4)),
                col = c("cornflowerblue", "coral2", "aquamarine3", "darkslategrey",
                        "green", "magenta3", "black"))
    }))
    plot.new()
    legend("center", lty = c(rep(1, 3), 2, rep(1,4)), lwd = c(rep(1,3), rep(2, 5)), 
           col = c("cornflowerblue", "coral2", "aquamarine3", "darkslategrey", NA,
                   "green", "magenta3", "black"), bty = "n",
           legend = c("ECMWF", "NCEP", "UKMO", "Superensemble", NA,
                      "Fitted - no principal components", 
                      "Fitted - 2 principal components", "Fitted - 3 principal components"))
    mtext("CRPS for each model, at each leadtime", outer = T)
}; dev.off()

####################################################################################################

# WORKING AREA                                                                                  ####





