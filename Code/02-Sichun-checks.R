
library("CB.Misc"); library("SX.weather")


# STILL TO DO...


#   - don't forget CompStat homework still to do (independence sampler)

#   - then look at superensemble & comparitive plots - once everything else is fully aligned


# calculation used is for RMSE of ensemble mean - not mean RMSE of ensemble
# check RMSE of all observations (no averaging until the very end) 

# UKMO control forecast for mslp is missing - should recheck superensemble (but prob. ok)

####################################################################################################

# PLOTS REPLICATED FROM DISSERTATION                                                            ####

plot.forecast.errors(ecmwf)
plot.forecast.errors(ncep)
plot.forecast.errors(ukmo)      # no control data for MSLP

plot.forecast.rmse(ecmwf)
plot.forecast.rmse(ncep)
plot.forecast.rmse(ukmo)

# RMSE vs ensemble spread
plot.rmse.spread(ecmwf)
plot.rmse.spread(ncep)
plot.rmse.spread(ukmo)

####################################################################################################

# CHECKING AGAINST SICHUN DATA                                                                  ####

# confirm that my calculations match those taken from Sichun's code
source("./Code/sx-code-all-data.R")

# observation series
check.observations <- function() {
    c("temp.n" = all(c(obs["temp.n",,]) == TN_obs),
      "temp.s" = all(c(obs["temp.s",,]) == TS_obs),
      "pc1" = all(c(obs["pc1",,]) == PC1_obs),
      "pc2" = all(c(obs["pc2",,]) == PC2_obs))
}


# forecast series
check.forecasts <- function() {
    
    res <- list()
    
    if(dim(ecmwf)[[4]] == 15) {r <- 2:15} else {r <- 1:14}
    
    fc.pert <- apply(ecmwf[,,,r,-1], c(1,4), c)
    
    res$ctrl <- c(all(array(ecmwf["temp.n",,,r, "c"], dim = c(630, length(r))) == TN.ctrl.leadtime_ECMWF),
                  all(array(ecmwf["temp.s",,,r, "c"], dim = c(630, length(r))) == TS.ctrl.leadtime_ECMWF),
                  all(array(ecmwf["pc1",,,r, "c"], dim = c(630, length(r))) == PC1.ctrl.leadtime_ECMWF),
                  all(array(ecmwf["pc2",,,r, "c"], dim = c(630, length(r))) == PC2.ctrl.leadtime_ECMWF))
    
    res$pert <- c(all(fc.pert[,"temp.n",] == TN.pert.leadtime_ECMWF),
                  all(fc.pert[,"temp.s",] == TS.pert.leadtime_ECMWF),
                  all(fc.pert[,"pc1",] == PC1.pert.leadtime_ECMWF),
                  all(fc.pert[,"pc2",] == PC2.pert.leadtime_ECMWF))
        
    return(res)
}


# data for mean error plots
check.mean.errors <- function() {
    # checks that mean errors match (to within 9dp, which is enough to convince me)
    
    model.errors <- forecast.errors(ecmwf)
    res <- list()
    
    # make sure same leadtimes are being compared
    if(dim(model.errors)[[4]] == 15) {r <- 2:15} else {r <- 1:14}
    
    res$ctrl <- c("temp.n" = all(round(apply(model.errors["temp.n",,,r,"c"], 3, mean), 9) == round(mean.error.ctrl_TN_ECMWF, 9)),
                  "temp.s" = all(round(apply(model.errors["temp.s",,,r,"c"], 3, mean), 9) == round(mean.error.ctrl_TS_ECMWF, 9)),
                  "pc1" = all(round(apply(model.errors["pc1",,,r,"c"], 3, mean), 9) == round(mean.error.ctrl_PC1_ECMWF, 9)),
                  "pc2" = all(round(apply(model.errors["pc2",,,r,"c"], 3, mean), 9) == round(mean.error.ctrl_PC2_ECMWF, 9)))
    
    res$pert <- c("temp.n" = all(round(apply(model.errors["temp.n",,,r,-1], 3, mean), 9) == round(mean.error.pert_TN_ECMWF, 9)),
                  "temp.s" = all(round(apply(model.errors["temp.s",,,r,-1], 3, mean), 9) == round(mean.error.pert_TS_ECMWF, 9)),
                  "pc1" = all(round(apply(model.errors["pc1",,,r,-1], 3, mean), 9) == round(mean.error.pert_PC1_ECMWF, 9)),
                  "pc2" = all(round(apply(model.errors["pc2",,,r,-1], 3, mean), 9) == round(mean.error.pert_PC2_ECMWF, 9)))
    
    return(res)
}


# data for rmse plots
# this gives RMSE of ensemble mean - not mean RMSE of ensemble
check.rmse <- function() {
    # checks that rmses match (to within 9dp, which is enough to convince me)
    
    # add ensemble mean to forecast set
    ecmwf <- abind("em" = apply(ecmwf[,,,,-1], 1:4, mean), ecmwf, along = 5)
    rmse <- forecast.rmse(ecmwf)
    res <- list()
    
    # make sure same leadtimes are being compared
    if(dim(rmse)[[2]] == 15) {r <- 2:15} else {r <- 1:14}
    
    res$ctrl <- c("temp.n" = all(round(rmse["temp.n",r,"c"], 9) == round(RMSE.ctrl_TN_ECMWF, 9)),
                  "temp.s" = all(round(rmse["temp.s",r,"c"], 9) == round(RMSE.ctrl_TS_ECMWF, 9)),
                  "pc1" = all(round(rmse["pc1",r,"c"], 9) == round(RMSE.ctrl_PC1_ECMWF, 9)),
                  "pc2" = all(round(rmse["pc2",r,"c"], 9) == round(RMSE.ctrl_PC2_ECMWF, 9)))
    
    res$pert <- c("temp.n" = all(round(rmse["temp.n",r,"em"], 9) == round(RMSE.pert_TN_ECMWF, 9)),
                  "temp.s" = all(round(rmse["temp.s",r,"em"], 9) == round(RMSE.pert_TS_ECMWF, 9)),
                  "pc1" = all(round(rmse["pc1",r,"em"], 9) == round(RMSE.pert_PC1_ECMWF, 9)),
                  "pc2" = all(round(rmse["pc2",r,"em"], 9) == round(RMSE.pert_PC2_ECMWF, 9)))
    
    return(res)
}


# data for rmse-vs-spread plots
check.rmse.spread <- function(dp = 9) {
    # checks that rmse & spread match (to within 9dp, which is enough to convince me)
    
    spr <- ensemble.spread(ecmwf[,,,,-1])
    
    # add ensemble mean to forecast set
    ecmwf <- abind("em" = apply(ecmwf[,,,,-1], 1:4, mean), ecmwf, along = 5)
    rmse <- forecast.rmse(ecmwf)
    res <- list()
    
    # make sure same leadtimes are being compared
    if(dim(rmse)[[2]] == 15) {r <- 2:15} else {r <- 1:14}
    
    res$rmse <- c("temp.n" = all(round(rmse["temp.n",r,"em"], dp) == round(RMSE_TN_ECMWF, dp)),
                  "temp.s" = all(round(rmse["temp.s",r,"em"], dp) == round(RMSE_TS_ECMWF, dp)),
                  "pc1" = all(round(rmse["pc1",r,"em"], dp) == round(RMSE_PC1_ECMWF, dp)),
                  "pc2" = all(round(rmse["pc2",r,"em"], dp) == round(RMSE_PC2_ECMWF, dp)))
                  
    res$spread <- c("temp.n" = all(round(spr["temp.n",r], 9) == round(Root_TN_ECMWF, 9)),
                    "temp.s" = all(round(spr["temp.s",r], 9) == round(Root_TS_ECMWF, 9)),
                    "pc1" = all(round(spr["pc1",r], 9) == round(Root_PC1_ECMWF, 9)),
                    "pc2" = all(round(spr["pc2",r], 9) == round(Root_PC2_ECMWF, 9)))
    
    return(res)
}


check.observations()
check.forecasts()
check.mean.errors()
check.rmse()
check.rmse.spread()

# all match, except for PC2 spread - somehow, SX code produces lt 0:13 instead of 1:14.

####################################################################################################

# SICHUN'S ORIGINAL PLOT CODE                                                                   ####

# mean error
{
    par(mfrow = c(2,2))
    
    # PC1
    plot(c(1:14),mean.error.ctrl_PC1_ECMWF,xlab='lead time(day)',ylab='mean error',type='l',col='red',
         main='(c) PC1 of Pressure',ylim = c(-0.22,0))
    axis(1,at=c(1:14))
    lines(c(1:14),mean.error.pert_PC1_ECMWF,lty=2,col='blue')
    abline(h=0,col="grey",lty=2,lwd=1)
    legend(7,0,legend=c('Control Forecast','Perturbed Forecast'),col=c('red','blue'),
           lty=c(1,2),cex = 0.7,bty = 'n')
    
    # PC2
    plot(c(1:14),mean.error.ctrl_PC2_ECMWF,xlab='lead time(day)',ylab='mean error',type='l',col='red',
         main='(d) PC2 of Pressure',ylim = c(0.1,0.25))
    axis(1,at=c(1:14))
    lines(c(1:14),mean.error.pert_PC2_ECMWF,lty=2,col='blue')
    abline(h=0,col="grey",lty=2,lwd=1)
    legend(7,0.12,legend=c('Control Forecast','Perturbed Forecast'),col=c('red','blue'),
           lty=c(1,2),cex = 0.7,bty = 'n')
    
    # Temperature in North UK
    plot(c(1:14),mean.error.ctrl_TN_ECMWF,xlab='lead time(day)',ylab=expression('mean error'~(degree*C)),type='l',col='red',main='(a) Temperature in North UK',ylim = c(-2.3,-1.85))
    axis(1,at=c(1:14))
    lines(c(1:14),mean.error.pert_TN_ECMWF,lty=2,col='blue')
    abline(h=0,col="grey",lty=2,lwd=1)
    legend(7,-1.85,legend=c('Control Forecast','Perturbed Forecast'),col=c('red','blue'),
           lty=c(1,2),cex = 0.7,bty = 'n')
    
    # Temperature in South UK
    plot(c(1:14),mean.error.ctrl_TS_ECMWF,xlab='lead time(day)',ylab=expression('mean error'~(degree*C)),type='l',col='red',main='(b) Temperature in South UK')
    axis(1,at=c(1:14))
    lines(c(1:14),mean.error.pert_TS_ECMWF,lty=2,col='blue')
    abline(h=0,col="grey",lty=2,lwd=1)
    legend(7,-0.85,legend=c('Control Forecast','Perturbed Forecast'),col=c('red','blue'),
           lty=c(1,2),cex = 0.7,bty = 'n')
    
    par(mfrow = c(1,1))
}

# RMSE
{
    par(mfrow = c(2,2))
    
    # PC1
    plot(c(1:14),RMSE.ctrl_PC1_ECMWF,xlab='lead time(day)',ylab='Root Mean Square Error',type='l',col='red',main='(c) PC1 of Pressure',ylim = c(0,1.2))
    axis(1,at=c(1:14))
    lines(c(1:14),RMSE.pert_PC1_ECMWF,lty=2,col='blue')
    abline(h=0,col="grey",lty=2,lwd=1)
    legend(7,0.2,legend=c('Control Forecast','Perturbed Forecast'),col=c('red','blue'),
           lty=c(1,2),cex = 0.7,bty = 'n')
    
    # PC2
    plot(c(1:14),RMSE.ctrl_PC2_ECMWF,xlab='lead time(day)',ylab='Root Mean Square Error',type='l',col='red',main='(d) PC2 of Pressure',ylim = c(0,1.2))
    axis(1,at=c(1:14))
    lines(c(1:14),RMSE.pert_PC2_ECMWF,lty=2,col='blue')
    abline(h=0,col="grey",lty=2,lwd=1)
    legend(7,0.2,legend=c('Control Forecast','Perturbed Forecast'),col=c('red','blue'),
           lty=c(1,2),cex = 0.7,bty = 'n')
    
    # Temperature in North UK
    plot(c(1:14),RMSE.ctrl_TN_ECMWF,xlab='lead time(day)',ylab='Root Mean Square Error',type='l',col='red',main='(a) Temperature in North UK')
    axis(1,at=c(1:14))
    lines(c(1:14),RMSE.pert_TN_ECMWF,lty=2,col='blue')
    abline(h=0,col="grey",lty=2,lwd=1)
    legend(7,2.3,legend=c('Control Forecast','Perturbed Forecast'),col=c('red','blue'),
           lty=c(1,2),cex = 0.7,bty = 'n')
    
    # Temperature in South UK
    plot(c(1:14),RMSE.ctrl_TS_ECMWF,xlab='lead time(day)',ylab='Root Mean Square Error',type='l',col='red', main='(b) Temperature in South UK')
    axis(1,at=c(1:14))
    lines(c(1:14),RMSE.pert_TS_ECMWF,lty=2,col='blue')
    abline(h=0,col="grey",lty=2,lwd=1)
    legend(7,1.5,legend=c('Control Forecast','Perturbed Forecast'),col=c('red','blue'),
           lty=c(1,2),cex = 0.7,bty = 'n')
    
    par(mfrow = c(1,1))
}

####################################################################################################

# CHECK SUPERENSEMBLE RESULTS                                                                   ####

source("./Sichun code/SUPER_Mean Error and RMSE.R")
source("./Sichun code/SUPER_RMSE vs Spread.R")

check.superens.mean.errors <- function(dp = 9) {
    
    super <- abind("c" = ecmwf[,,,,"c"], "c" = ncep[,,,,"c"], "c" = ukmo[,,,,"c"], 
                   ecmwf[,,,,-1], ncep[,,,,-1], ukmo[,,,,-1], along = 5)
    
    s.err <- forecast.errors(super)
    
    s.ctrl.mean <- round(apply(s.err[,,,2:15,1:3], c(1,4), mean), dp)
    s.ens.mean <- round(apply(s.err[,,,2:15,-(1:3)], c(1,4), mean), dp)
    
    res <- list()
    
    res$crtl <- c("temp.n" = all(s.ctrl.mean["temp.n",] == round(mean.error.ctrl_TN_ALL, dp)),
                  "temp.s" = all(s.ctrl.mean["temp.s",] == round(mean.error.ctrl_TS_ALL, dp)),
                  "pc1" = all(s.ctrl.mean["pc1",] == round(mean.error.ctrl_PC1_ALL, dp)),
                  "pc2" = all(s.ctrl.mean["pc2",] == round(mean.error.ctrl_PC2_ALL, dp)))
    
    res$pert <- c("temp.n" = all(s.ens.mean["temp.n",] == round(mean.error.pert_TN_ALL, dp)),
                  "temp.s" = all(s.ens.mean["temp.s",] == round(mean.error.pert_TS_ALL, dp)),
                  "pc1" = all(s.ens.mean["pc1",] == round(mean.error.pert_PC1_ALL, dp)),
                  "pc2" = all(s.ens.mean["pc2",] == round(mean.error.pert_PC2_ALL, dp)))
    
    return(res)
}

check.superens.rmse.spread <- function(dp = 9) {
    
    super <- superensemble(list(ecmwf, ncep, ukmo))
    err <- forecast.errors(super)[,,,2:15,]
    
    ctrl.rmse <- sqrt(apply(err[,,,,1:3]^2, c(1,4), mean))
    pert.rmse <- round(sqrt(apply(apply(err[,,,,-(1:3)], 1:4, mean)^2, c(1,4), mean)), dp)
    # calculated as RMSE of mean of all ensemble members
    
    se.spread <- ensemble.spread(super[,,,,-(1:3)])[,2:15]
    
    list("ctrl.rmse" = c("temp.n" = all(ctrl.rmse["temp.n",] == RMSE.ctrl_TN_ALL),
                    "temp.s" = all(ctrl.rmse["temp.s",] == RMSE.ctrl_TS_ALL),
                    "pc1" = all(ctrl.rmse["pc1",] == RMSE.ctrl_PC1_ALL),
                    "pc2" = all(ctrl.rmse["pc2",] == RMSE.ctrl_PC2_ALL)),
         "pert.rmse" = c("temp.n" = all(pert.rmse["temp.n",] == round(RMSE.pert_TN_ALL,dp)),
                    "temp.s" = all(pert.rmse["temp.s",] == round(RMSE.pert_TS_ALL,dp)),
                    "pc1" = all(pert.rmse["pc1",] == round(RMSE.pert_PC1_ALL,dp)),
                    "pc2" = all(pert.rmse["pc2",] == round(RMSE.pert_PC2_ALL,dp))),
         "pert.spread" = c("temp.n" = all(se.spread["temp.n",] == Root_TN_ALL),
                           "temp.s" = all(se.spread["temp.s",] == Root_TS_ALL),
                           "pc1" = all(se.spread["pc1",] == Root_PC1_ALL),
                           "pc2" = all(se.spread["pc2",] == Root_PC2_ALL)))
}

check.superens.mean.errors()
check.superens.rmse.spread()

####################################################################################################

# CDF & CRPS                                                                                    ####

# example based on Sichun's code from CDF_ensemble.R
# produces same plot as in Sichun's code
{
    x <- (0:23)/23
    y <- sort(c(0,ukmo["temp.s",15,"08","1",-1]))
    p <- sort(c((0:23)/23, sum(y < obs["temp.s",1,"08"])/23))
    z <- sort(c(0, ukmo["temp.s",15,"08","1",-1], obs["temp.s",1,"08"])) < obs["temp.s",1,"08"]
    
    crps <- (p - (1-z))^2
    
    plot(y, x, type = 's', ylab = 'Cumulative Distribution Function',
         xlab = 'Temperature in south UK', col='red',lwd=3, xlim = c(7,8))
    lines(c(0,obs["temp.s",1,"08"],10),c(0,1,1),type='s',lwd=3,col='blue')
    lines(sort(c(y, obs["temp.s",1,"08"])), crps, lwd = 3, type = "s", col = "green3")
    legend(7,1,legend=c('Forecast Ensemble','Observation','CRPS'),col=c('red','blue','green3'),lty=c(1,1,1),lwd=c(2,2,2),cex = 0.8,bty = 'n')
    
}

####################################################################################################

# COVARIANCE_Method1 for Lambda and Eta                                                         ####

source("./Code/03-covariances-method-1.R")              # ~ 5mins to run

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

# check common matrices
{
    cat("Confirm that all common matrices match...")
    cat("\n", "Covariances:", 
        all(array(ecmwf.cov[1:4,1:4,,,2:15], dim = dim(C_ECMWF)) == C_ECMWF,
            array(ncep.cov[1:4,1:4,,,2:15], dim =  dim(C_NCEP)) == C_NCEP,
            array(ukmo.cov[1:4,1:4,,,2:15], dim =  dim(C_UKMO)) == C_UKMO))
    
    cat("\n", "Sigma:",
        all(array(c.sigma[1:4,1:4,,,2:15], dim = dim(sigma)) == sigma))
    
    cat("\n", "Lambda:",
        all(round(c.lambda[1:4,1:4,,2:15],9) == round(lambda,9)))
    
    cat("\n", "Eta:",
        all(round(c.eta[1:4,,2:15], 9) == round(eta,9)))
    
    cat("\n", "Y_bar:",
        all(c(y.bar["ecmwf",1:4,,,2:15]) == c(Y_bar_ECMWF),
            c(y.bar["ncep",1:4,,,2:15]) == c(Y_bar_NCEP),
            c(y.bar["ukmo",1:4,,,2:15]) == c(Y_bar_UKMO)))
    
    cat("\n", "D:",
        all(array(d["ecmwf",1:4,1:4,,,2:15], dim = dim(D_ECMWF)) == D_ECMWF,
            array(d["ncep",1:4,1:4,,,2:15], dim = dim(D_NCEP)) == D_NCEP,
            array(d["ukmo",1:4,1:4,,,2:15], dim = dim(D_UKMO)) == D_UKMO))

}

# two-PC model
{
    precision.2 <- se.precision(d, c.lambda, vars = c(1:4))
    s.2 <- array.solve(precision.2, 3:5)
    tau.2 <- se.tau(d, c.lambda, y.bar, s.2, c.eta)
    
    spr.2 <- sqrt(apply(apply(s.2, 3:5, diag), c(1,4), mean))
    rmse.2 <- model.rmse(tau.2)
    crps.2 <- model.crps(tau.2, s.2)
}

# check model with two principal components
{
    cat("Confirm that all matrices match Sichun's for model with two principal components...")
    
    cat("\n", "Precision:", all(round(array(precision.2[,,,,2:15], dim = dim(precision)), 6) == round(precision, 6)))
    
    cat("\n", "S:", all(round(array(s.2[,,,,2:15], dim = dim(S)), 9) == round(S,9)))
    
    cat("\n", "Tau:", all(round(array(tau.2[,,,2:15], dim = dim(tau)), 8) == round(tau, 8)))
    
    cat("\n", "\n")
    cat("Confirm that results match Sichun's for two-PC model...")
    
    cat("\n", "Spread:", all(round(spr.2[,2:15], 9) == round(Spread, 9)))
    cat("\n", "RMSE:", all(round(rmse.2[,2:15], 9) == round(RMSE, 9)))
    cat("\n", "CRPS:", all(round(c(crps.2[,,,2:15]), 8) == round(c(CRPS_Richard), 8)))
}

####################################################################################################

# COVARIANCE for temp only                                                                      ####

source("./Code/04-covariances-for-temp-only.R")        # ~ 3 minutes

# common matrices are run as above

# check common matrices
{
    cat("Confirm that all common matrices match...")
    cat("\n", "Covariances:", 
        all(array(ecmwf.cov[1:2,1:2,,,2:15], dim = dim(C_ECMWF)) == C_ECMWF,
            array(ncep.cov[1:2,1:2,,,2:15], dim =  dim(C_NCEP)) == C_NCEP,
            array(ukmo.cov[1:2,1:2,,,2:15], dim =  dim(C_UKMO)) == C_UKMO))
    
    cat("\n", "Sigma:",
        all(array(c.sigma[1:2,1:2,,,2:15], dim = dim(sigma)) == sigma))
    
    cat("\n", "Lambda:",
        all(round(c.lambda[1:2,1:2,,2:15],9) == round(lambda,9)))
    
    cat("\n", "Eta:",
        all(round(c.eta[1:2,,2:15], 9) == round(eta,9)))
    
    cat("\n", "Y_bar:",
        all(c(y.bar["ecmwf",1:2,,,2:15]) == c(Y_bar_ECMWF),
            c(y.bar["ncep",1:2,,,2:15]) == c(Y_bar_NCEP),
            c(y.bar["ukmo",1:2,,,2:15]) == c(Y_bar_UKMO)))
    
    cat("\n", "D:",
        all(array(d["ecmwf",1:2,1:2,,,2:15], dim = dim(D_ECMWF)) == D_ECMWF,
            array(d["ncep",1:2,1:2,,,2:15], dim = dim(D_NCEP)) == D_NCEP,
            array(d["ukmo",1:2,1:2,,,2:15], dim = dim(D_UKMO)) == D_UKMO))
    
}

# zero-PC model
{
    precision.0 <- se.precision(d, c.lambda, vars = c(1:2))
    s.0 <- array.solve(precision.0, 3:5)
    tau.0 <- se.tau(d, c.lambda, y.bar, s.0, c.eta)
    
    spr.0 <- sqrt(apply(apply(s.0, 3:5, diag), c(1,4), mean))
    rmse.0 <- model.rmse(tau.0)
    crps.0 <- model.crps(tau.0, s.0)
}

# check model with no principal components
{
    cat("Confirm that all matrices match Sichun's for model with no principal components...")
    
    cat("\n", "Precision:", all(round(array(precision.0[,,,,2:15], dim = dim(precision)), 6) == round(precision, 6)))
    
    cat("\n", "S:", all(round(array(s.0[,,,,2:15], dim = dim(S)), 9) == round(S,9)))
    
    cat("\n", "Tau:", all(round(array(tau.0[,,,2:15], dim = dim(tau)), 8) == round(tau, 8)))
    
    cat("\n", "\n")
    cat("Confirm that results match Sichun's for two-PC model...")
    
    cat("\n", "Spread:", all(round(spr.0[,2:15], 9) == round(Spread, 9)))
    cat("\n", "RMSE:", all(round(rmse.0[,2:15], 9) == round(RMSE, 9)))
    cat("\n", "CRPS:", all(round(c(crps.0[,,,2:15]), 8) == round(c(CRPS_Richard3), 8)))
}


####################################################################################################

# CRPS.R                                                                                        ####

source("./Code/05-crps.R")

ens.crps <- abind("ecmwf" = ensemble.crps(ecmwf[,,,,-1]),
                  "ncep" = ensemble.crps(ncep[,,,,-1]),
                  "ukmo" = ensemble.crps(ukmo[,,,,-1]), along = 0)

superens.crps <- ensemble.crps(superensemble()[,,,,-(1:3)])

# check CRPS of individual ensembles is calculated correctly
{
    cat("Confirm that CRPS of individual ensembles is calculated correctly...")
    cat("\n", "Temp.n:",
        all(ens.crps["ecmwf", "temp.n", 2:15] == ECMWF_TN,
            ens.crps["ncep", "temp.n", 2:15] == NCEP_TN,
            ens.crps["ukmo", "temp.n", 2:15] == UKMO_TN))
    
    cat("\n", "Temp.s:",
        all(ens.crps["ecmwf", "temp.s", 2:15] == ECMWF_TS,
            ens.crps["ncep", "temp.s", 2:15] == NCEP_TS,
            ens.crps["ukmo", "temp.s", 2:15] == UKMO_TS))
    
    cat("\n", "PC1:",
        all(ens.crps["ecmwf", "pc1", 2:15] == ECMWF_PC1,
            ens.crps["ncep", "pc1", 2:15] == NCEP_PC1,
            ens.crps["ukmo", "pc1", 2:15] == UKMO_PC1))
    
    cat("\n", "PC2:",
        all(ens.crps["ecmwf", "pc2", 2:15] == ECMWF_PC2,       # think PC2 leadtime is incorrect in SX code
            ens.crps["ncep", "pc2", 2:15] == NCEP_PC2,
            ens.crps["ukmo", "pc2", 2:15] == UKMO_PC2))
    cat("     (PC2 leadtime calculation appears to be incorrect in SX code)")
    
    cat("\n")
    cat("\n", "Superensemble CRPS for all variables:",
        all(superens.crps["temp.n",2:15] == ALL_TN,
            superens.crps["temp.s",2:15] == ALL_TS,
            superens.crps["pc1",2:15] == ALL_PC1,
            superens.crps["pc2",2:15] == ALL_PC2))
}


####################################################################################################
