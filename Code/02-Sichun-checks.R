
# suspect problem in data import is to do with forecast leadtime. Look into this....

library("CB.Misc"); library("SX.weather")


# STILL TO DO...

#   - once checks are included here, tidy & rename '01-import...' to include only data import (sourceable)

#   - don't forget CompStat homework still to do (independence sampler)

#   - then look at superensemble & comparitive plots - once everything else is fully aligned


# calculation used is for RMSE of ensemble mean - not mean RMSE of ensemble
# check RMSE of all observations (no averaging until the very end) 

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
load.all()
source("./Code/sx-code-missing-data.R")

# when remaining data files are obtained, replace with
# source("./Code/sx-code-all-data.R")

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
check.rmse.spread <- function() {
    # checks that rmse & spread match (to within 9dp, which is enough to convince me)
    
    spr <- ensemble.spread(ecmwf[,,,,-1])
    
    # add ensemble mean to forecast set
    ecmwf <- abind("em" = apply(ecmwf[,,,,-1], 1:4, mean), ecmwf, along = 5)
    rmse <- forecast.rmse(ecmwf)
    res <- list()
    
    # make sure same leadtimes are being compared
    if(dim(rmse)[[2]] == 15) {r <- 2:15} else {r <- 1:14}
    
    RMSE_TN_ECMWF <- sqrt(colMeans((TN.pert.leadtime.ensemblemean_ECMWF - TN_obs)^2))
    
    
    res$rmse <- c("temp.n" = all(round(rmse["temp.n",r,"em"], 9) == round(RMSE_TN_ECMWF, 9)),
                  "temp.s" = all(round(rmse["temp.s",r,"em"], 9) == round(RMSE_TS_ECMWF, 9)),
                  "pc1" = all(round(rmse["pc1",r,"em"], 9) == round(RMSE_PC1_ECMWF, 9)),
                  "pc2" = all(round(rmse["pc2",r,"em"], 9) == round(RMSE_PC2_ECMWF, 9)))
                  
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
