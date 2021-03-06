
library("CB.Misc"); library("SX.weather")
# update.functions("SX.weather")

# currently just replicating Sichun's objects & analysis

####################################################################################################

# NOTES                                                                                         ####

# can't do the spatial plots etc because source files unavailable.

# missing data files:
    # ECMWF_europe_starttime.rda
    # ECMWF_europe_leadtime.rda
    # ERAint_pca_espace.rda
    # ECMWF_europe_lon.rda
    # ECMWF_europe_lat.rda
    # ECMWFctrl_mslp_8yr.rda
    # ERAint_europe_starttime.rda

# Sichun code used, in order:
    # ECMWF_loading.R


####################################################################################################

# ECMWF_loading.R                                                                               ####

control_ts <- load.data("./Data/ECMWFctrl_mslp_timeseries.rda")
ensemble_ts <- load.data("./Data/ECMWFpert_timeseries.rda")

southctrltemp <- load.data("./Data/ECMWFctrl_temp2m_7yr_south24hrave.rda")
northctrltemp <- load.data("./Data/ECMWFctrl_temp2m_7yr_north24hrave.rda")

# can't load start time & lead time, will work without

{
    # check alignment of control (8yr) & ensemble (7yr) time series for PC
    
    # plot(control_ts[,1,1], type = "l")
    # lines(ensemble_ts[,1,1,1], col = "red")
    # lines(ensemble_ts[,1,1,2], col = "blue")
    
    # seem to have shared start date, rather than shared end date, so can just truncate.
}


#===========================================================================================
# 1. DATA FRAMES FOR CONTROL FORECAST

# Lead time for pressure is 61, but for temp is 58, so use shorter: 6hr intervals for 14 days.
# 1-day increments are the 1st, 5th, 9th, ... data.

ctrl_obs <- data.frame(year = rep(formatC(c(8:14), width = 2, flag = "0"), each = 105),
                       day = rep(1:105, 7),
                       day.ahead = 0,
                       temp.north = c(northctrltemp[1,,])-273.15,
                       temp.south = c(southctrltemp[1,,])-273.15,
                       PC1 = c(control_ts[1:735,1,1]),
                       PC2 = c(control_ts[1:735,2,1]),
                       PC3 = c(control_ts[1:735,3,1]))

d <- 1
ctrl_aheadday <- data.frame(year = rep(formatC(c(8:14), width = 2, flag = "0"), each = 105),
                       day = rep(1:105, 7),
                       day.ahead = d,
                       temp.north = c(northctrltemp[4 * d + 1,,])-273.15,
                       temp.south = c(southctrltemp[4 * d + 1,,])-273.15,
                       PC1 = c(control_ts[1:735,1,4 * d + 1]),
                       PC2 = c(control_ts[1:735,2,4 * d + 1]),
                       PC3 = c(control_ts[1:735,3,4 * d + 1]))

# create DF for each daily leadtime
#invisible(sapply(1:14, function(i) {
#    j <- i * 4 + 1
#    tmp <- data.frame(year = rep(formatC(c(8:14), width = 2, flag = "0"), each = 105),
#                      day = rep(1:105, 7),
#                      day.ahead = i,
#                      temp.north = c(northctrltemp[j,,])-273.15,
#                      temp.south = c(southctrltemp[j,,])-273.15,
#                      PC1 = c(control_ts[1:735,1,j]),
#                      PC2 = c(control_ts[1:735,2,j]),
#                      PC3 = c(control_ts[1:735,3,j]))
#    assign(paste0("ctrl_", i, "day"), tmp, pos = globalenv())
#}))


#===========================================================================================
# 2. ADD PERTURBED ENSEMBLES INTO DATA FRAMES

npt <- load.data("./Data/ECMWFpert_temp2m_7yr_north24hrave.rda")
spt <- load.data("./Data/ECMWFpert_temp2m_7yr_south24hrave.rda")

# temp files are in different order (leadime-day-pert-year, vs dayyear-var-leadtime-pert)
# so need to rearrange them. OF COURSE

northperttemp <- southperttemp <- array(dim = c(58, 105, 7, 50))

invisible(sapply(1:50, function(i) {
    northperttemp[,,,i] <<- npt[,,i,] 
    southperttemp[,,,i] <<- spt[,,i,] 
}))

# For lead time 0
pert_obs.all <- data.frame(year = rep(formatC(c(8:14), width = 2, flag = "0"), each = 105),
                           day = rep(1:105, 7),
                           day.ahead = 0,
                           pert = rep(1:50, each = 735),
                           temp.north = c(northperttemp[1,,,])-273.15,
                           temp.south = c(southperttemp[1,,,])-273.15,
                           PC1 = c(ensemble_ts[1:735,1,1,]),
                           PC2 = c(ensemble_ts[1:735,2,1,]),
                           PC3 = c(ensemble_ts[1:735,3,1,]))

# same principal applies to import all data. Not going to do so until purpose is clearer.

####################################################################################################

# ECMWF_Mean error and RMSE.R                                                                   ####

# variable names retained, method of obtaining them corrected

fc <- readRDS("./Data/ECMWF-forecasts.rds")[,16:105,,2:15,]

# The control forecasts for 14 lead times and 630 days (90 winter days * 7 years).
PC1.ctrl.leadtime_ECMWF <- array(fc["pc1",,,,"c"], dim = c(90 * 7, 14))
PC2.ctrl.leadtime_ECMWF <- array(fc["pc2",,,,"c"], dim = c(90 * 7, 14))
PC3.ctrl.leadtime_ECMWF <- array(fc["pc3",,,,"c"], dim = c(90 * 7, 14))
TN.ctrl.leadtime_ECMWF <- array(fc["temp.n",,,,"c"], dim = c(90 * 7, 14))
TS.ctrl.leadtime_ECMWF <- array(fc["temp.s",,,,"c"], dim = c(90 * 7, 14))

# The perturbed forecasts for 14 lead times and 630 days.
PC1.pert.leadtime_ECMWF <- array(aperm(fc, c(1,2,3,5,4))["pc1",,,2:51,], dim = c(90 * 7 * 50, 14))
PC2.pert.leadtime_ECMWF <- array(aperm(fc, c(1,2,3,5,4))["pc2",,,2:51,], dim = c(90 * 7 * 50, 14))
PC3.pert.leadtime_ECMWF <- array(aperm(fc, c(1,2,3,5,4))["pc3",,,2:51,], dim = c(90 * 7 * 50, 14))
TN.pert.leadtime_ECMWF <- array(aperm(fc, c(1,2,3,5,4))["temp.n",,,2:51,], dim = c(90 * 7 * 50, 14))
TS.pert.leadtime_ECMWF <- array(aperm(fc, c(1,2,3,5,4))["temp.s",,,2:51,], dim = c(90 * 7 * 50, 14))


#===========================================================================================
# observations from ERA-INT

timeseries <- load.data("./Data/ERAint_pca_timeseries.rda")

# 3240 - 720 + 1 = 2611, this is the start time for 2007/2008-2013/2014.
PC1_obs <- timeseries[2521:3150,1] 
PC2_obs <- timeseries[2521:3150,2] 
PC3_obs <- timeseries[2521:3150,3] 

# didn't investigate indexing for this file - this was provided by SX
TN_obs <- c(load.data("./Data/ERAint_temp2m_36yr_north24hrave.rda")[seq(3,360,4),29:35]) - 273.15
TS_obs <- c(load.data("./Data/ERAint_temp2m_36yr_south24hrave.rda")[seq(3,360,4),29:35]) - 273.15

#===========================================================================================
# mean error

# for control forecast 
{
    # mean error = sum(PC1.ctrl.leadtime-PC1.ctrl_obs)/630
    mean.error.ctrl_PC1_ECMWF <- colMeans(PC1.ctrl.leadtime_ECMWF - PC1_obs)
    mean.error.ctrl_PC2_ECMWF <- colMeans(PC2.ctrl.leadtime_ECMWF - PC2_obs)
    mean.error.ctrl_PC3_ECMWF <- colMeans(PC3.ctrl.leadtime_ECMWF - PC3_obs)
    mean.error.ctrl_TN_ECMWF <- colMeans(TN.ctrl.leadtime_ECMWF - TN_obs)
    mean.error.ctrl_TS_ECMWF <- colMeans(TS.ctrl.leadtime_ECMWF - TS_obs)
}

# for perturbed forecast
{
    # first get the ensemble mean for each variable and each lead time
    PC1.pert.leadtime.ensemblemean_ECMWF <- apply(array(PC1.pert.leadtime_ECMWF, dim = c(630, 50, 14)), c(1,3), mean)
    PC2.pert.leadtime.ensemblemean_ECMWF <- apply(array(PC2.pert.leadtime_ECMWF, dim = c(630, 50, 14)), c(1,3), mean)
    PC3.pert.leadtime.ensemblemean_ECMWF <- apply(array(PC3.pert.leadtime_ECMWF, dim = c(630, 50, 14)), c(1,3), mean)
    TN.pert.leadtime.ensemblemean_ECMWF <- apply(array(TN.pert.leadtime_ECMWF, dim = c(630, 50, 14)), c(1,3), mean)
    TS.pert.leadtime.ensemblemean_ECMWF <- apply(array(TS.pert.leadtime_ECMWF, dim = c(630, 50, 14)), c(1,3), mean)
    
    # then mean error
    mean.error.pert_PC1_ECMWF <- colMeans(PC1.pert.leadtime.ensemblemean_ECMWF - PC1_obs)
    mean.error.pert_PC2_ECMWF <- colMeans(PC2.pert.leadtime.ensemblemean_ECMWF - PC2_obs)
    mean.error.pert_PC3_ECMWF <- colMeans(PC3.pert.leadtime.ensemblemean_ECMWF - PC3_obs)
    mean.error.pert_TN_ECMWF <- colMeans(TN.pert.leadtime.ensemblemean_ECMWF - TN_obs)
    mean.error.pert_TS_ECMWF <- colMeans(TS.pert.leadtime.ensemblemean_ECMWF - TS_obs)
}

# (my addition) mean error of each perturbation
{
    PC3.error.envelope <- array(PC3.pert.leadtime_ECMWF, dim = c(630, 50, 14)) - PC3_obs
    PC2.error.envelope <- array(PC2.pert.leadtime_ECMWF, dim = c(630, 50, 14)) - PC2_obs
    PC1.error.envelope <- array(PC1.pert.leadtime_ECMWF, dim = c(630, 50, 14)) - PC1_obs
    TN.error.envelope <- array(TN.pert.leadtime_ECMWF, dim = c(630, 50, 14)) - TN_obs
    TS.error.envelope <- array(TS.pert.leadtime_ECMWF, dim = c(630, 50, 14)) - TS_obs
}

# plots of mean error for each component, averaged over 630 days at each lead time
{
    qplot <- function(element, ...) {
        
        ttl <- switch(element,
                      "TS" = "Temp (south)", "TN" = "Temp (north)",
                      "PC1" = "First PC", "PC2" = "Second PC", "PC3" = "Third PC")
        
        matplot(apply(eval(parse(text = paste0(element, ".error.envelope"))), 3:2, mean),
                type = "l", col = adjustcolor("grey", alpha = 0.4), lty = 1, main = ttl,
                xlab = "", ylab = "", ...)
        
        lines(eval(parse(text = paste0("mean.error.ctrl_", element, "_ECMWF"))), col = "blue3")
        lines(eval(parse(text = paste0("mean.error.pert_", element, "_ECMWF"))))
    }
    
    par(mfrow = c(2,3), oma = c(0.5, 0.5, 2, 0.5), mar = c(2,2,3,1))
    qplot("PC1"); qplot("PC2"); qplot("PC3")
    qplot("TN"); qplot("TS")
    
    plot.new()
    legend("left", lty = 1, col = c("blue", "black", adjustcolor("grey", alpha = 0.5)), bty = "n", cex = 1.1,
           legend = c("Control forecast", "Perturbed mean", "Perturbed members"))
    
    mtext("Mean error at each forecast lead time", outer = TRUE, cex = 1.5)
    
    par(mfrow = c(1,1), oma = c(0,0,0,0), mar = c(5.1, 4.1, 4.1, 2.1))
    
    # fix same scales for easier comparison with dissertation
    par(mfrow = c(2,2), oma = c(0.5, 0.5, 2, 0.5), mar = c(2,2,3,1))
    qplot("TN", ylim = c(-2.3,-1.85)); qplot("TS", ylim = c(-1.1, -0.85)) 
    qplot("PC1", ylim = c(-0.22,0)); qplot("PC2", ylim = c(0.1,0.25))

}


#===========================================================================================
# root mean squared error

# for control forecast
# RMSE=sum(leadtime-observation)^2/630
RMSE.ctrl_PC1_ECMWF <- sqrt(colMeans((PC1.ctrl.leadtime_ECMWF - PC1_obs)^2))
RMSE.ctrl_PC2_ECMWF <- sqrt(colMeans((PC2.ctrl.leadtime_ECMWF - PC2_obs)^2))
RMSE.ctrl_PC3_ECMWF <- sqrt(colMeans((PC3.ctrl.leadtime_ECMWF - PC3_obs)^2))
RMSE.ctrl_TN_ECMWF<- sqrt(colMeans((TN.ctrl.leadtime_ECMWF - TN_obs)^2))
RMSE.ctrl_TS_ECMWF <- sqrt(colMeans((TS.ctrl.leadtime_ECMWF - TS_obs)^2))

# for perturbed forecast
RMSE.pert_PC1_ECMWF <- sqrt(colMeans((PC1.pert.leadtime.ensemblemean_ECMWF - PC1_obs)^2))
RMSE.pert_PC2_ECMWF <- sqrt(colMeans((PC2.pert.leadtime.ensemblemean_ECMWF - PC2_obs)^2))
RMSE.pert_PC3_ECMWF <- sqrt(colMeans((PC3.pert.leadtime.ensemblemean_ECMWF - PC3_obs)^2))
RMSE.pert_TN_ECMWF <- sqrt(colMeans((TN.pert.leadtime.ensemblemean_ECMWF - TN_obs)^2))
RMSE.pert_TS_ECMWF <- sqrt(colMeans((TS.pert.leadtime.ensemblemean_ECMWF - TS_obs)^2))

# Sichun's plots
{
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
}

####################################################################################################

# ECMWF_RMSE vs spread.R                                                                        ####


# For each lead time and variable,
# Calculate the average variance of the ensemble and average variance for each lead time.
#
## PC1
{
    variance_PC1 <- matrix(nrow=630,ncol=14)
    for (i in 1:14){
        group <- rep(1:630, 50)
        variance_PC1[,i] <- tapply(PC1.pert.leadtime_ECMWF[,i],group, var)
    }
    
    averagevariance_PC1 <- colMeans(variance_PC1)
    
    Root_PC1_ECMWF <- sqrt(averagevariance_PC1)
    
    ###RMSE for mean of ensemble (should be divided into 50)
    RMSE_PC1_ECMWF <- sqrt(colMeans((PC1.pert.leadtime.ensemblemean_ECMWF - PC1_obs)^2))
}

## PC2
{
    variance_PC2 <- matrix(nrow=630,ncol=14)
    for (i in 1:14){
        group <- rep(1:630, 50)
        variance_PC2[,i] <- tapply(PC2.pert.leadtime_ECMWF[,i],group, var)
    }
    
    averagevariance_PC2 <- colMeans(variance_PC2)
    
    Root_PC2_ECMWF <- sqrt(averagevariance_PC2)
    
    ###RMSE for mean of ensemble (should be divided into 50)
    RMSE_PC2_ECMWF <- sqrt(colMeans((PC2.pert.leadtime.ensemblemean_ECMWF - PC2_obs)^2))
}

## PC3
{
    variance_PC3 <- matrix(nrow=630,ncol=14)
    for (i in 1:14){
        group <- rep(1:630, 50)
        variance_PC3[,i] <- tapply(PC3.pert.leadtime_ECMWF[,i],group, var)
    }
    
    averagevariance_PC3 <- colMeans(variance_PC3)
    
    Root_PC3_ECMWF <- sqrt(averagevariance_PC3)
    
    ###RMSE for mean of ensemble (should be divided into 50)
    RMSE_PC3_ECMWF <- sqrt(colMeans((PC3.pert.leadtime.ensemblemean_ECMWF - PC3_obs)^2))
}

# temperature requires that index file again. Will assume that array calculation works for all variables.
# temp - N
{
    northperttemp <- load.data("./Data/ECMWFpert_temp2m_7yr_north24hrave.rda")
    
    ## Temperature north
    variance_TN <- matrix(nrow=630,ncol=14)
    for (i in 1:630){
        for (j in 1:14){
            variance_TN[i,j] <- var(northperttemp[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,ceiling(i/90)]-273.15)
        }
    }
    
    averagevariance_TN <- colMeans(variance_TN)
    
    Root_TN_ECMWF <- sqrt(averagevariance_TN)
    
    ###RMSE for mean of ensemble (should be divided into 50)
    RMSE_TN_ECMWF <- sqrt(colMeans((TN.pert.leadtime.ensemblemean_ECMWF - TN_obs)^2))
    
}

# temp - S
{
    southperttemp <- load.data("./Data/ECMWFpert_temp2m_7yr_south24hrave.rda")
    
    ## Temperature south
    variance_TS <- matrix(nrow=630,ncol=14)
    for (i in 1:630){
        for (j in 1:14){
            variance_TS[i,j] <- var(southperttemp[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,ceiling(i/90)]-273.15)
        }
    }
    
    averagevariance_TS <- colMeans(variance_TS)
    
    Root_TS_ECMWF <- sqrt(averagevariance_TS)
    
    ###RMSE for mean of ensemble (should be divided into 50)
    RMSE_TS_ECMWF <- sqrt(colMeans((TS.pert.leadtime.ensemblemean_ECMWF - TS_obs)^2))
    
}

# check against array calculation
daily.var <- array(apply(array(fc[, ,,,-1], dim = c(5, 630, 14, 50)), 1:3, var),
                   dim = c(5, 90, 7, 14), 
                   dimnames = list(dimnames(fc)[[1]], 1:90, dimnames(fc)[[3]], 1:14))

average.var <- apply(daily.var, c(1, 4), mean)
root.mean.var <- sqrt(average.var)

obs <- readRDS("./Data/Observations.rds")
ens.mean <- apply(fc[,,,,-1], c(1,2,3,4), mean)
ens.mean.error <- sweep(ens.mean, 1:3, obs, "-")
ens.rmse <- sqrt(apply(ens.mean.error^2, c(1, 4), mean))

# check that everything stacks up
{
    # check that everything stacks up
    all(daily.var["pc1",,,] == array(variance_PC1, dim = dim(daily.var[1,,,])))     # T
    all(daily.var["pc2",,,] == array(variance_PC2, dim = dim(daily.var[1,,,])))     # T
    all(daily.var["pc3",,,] == array(variance_PC3, dim = dim(daily.var[1,,,])))     # T
   # all(daily.var["temp.n",,,] == array(variance_temp.n, dim = dim(daily.var[1,,,])))   # can't check
   # all(daily.var["temp.s",,,] == array(variance_temp.s, dim = dim(daily.var[1,,,])))   # can't check
    
    all(average.var["pc1",] == averagevariance_PC1)
    all(average.var["pc2",] == averagevariance_PC2)
    all(average.var["pc3",] == averagevariance_PC3)
    
    all(root.mean.var["pc1",] == Root_PC1_ECMWF)
    all(root.mean.var["pc2",] == Root_PC2_ECMWF)
    all(root.mean.var["pc3",] == Root_PC3_ECMWF)
    
    all(ens.mean["pc1",,,] == array(PC1.pert.leadtime.ensemblemean_ECMWF, dim = c(90, 7, 14)))
    all(ens.mean.error["pc1",,,] == array(PC1.pert.leadtime.ensemblemean_ECMWF - PC1_obs, dim = c(90,7, 14)))
 
    all(ens.rmse["pc1",] == RMSE_PC1_ECMWF) 
    all(ens.rmse["pc2",] == RMSE_PC2_ECMWF)   
    all(ens.rmse["pc3",] == RMSE_PC3_ECMWF)   
}

# and using functions to calculate directly...
{
    spread <- ensemble.spread(fc[,,,,-1])
        all(spread["pc1",] == Root_PC1_ECMWF)
        all(spread["pc2",] == Root_PC2_ECMWF)
        all(spread["pc3",] == Root_PC3_ECMWF)
        
    
}
plot(c(1:14),RMSE_PC1_ECMWF,xlab='lead time(day)',ylab='RMSE vs Spread',type='l',col='red',
     main='(c) PC1 of Pressure',ylim = c(0,0.9))
lines(c(1:14),Root_PC1_ECMWF,lty=2,col='blue')
axis(1,at=c(1:14))
abline(h=0,col="grey",lty=2,lwd=1)
legend(5,0.15,legend=c('RMSE of Ensemble Mean','Average Ensemble Spread'),col=c('red','blue'),
       lty=c(1,2),cex = 0.7,bty = 'n')

# plot RMSE vs spread for comparison to dissertation
{
    qplot <- function(element, ...) {
        ttl <- switch(element,
                      "temp.s" = "Temp (south)", "temp.n" = "Temp (north)",
                      "pc1" = "First PC", "pc2" = "Second PC", "pc3" = "Third PC")
        
        yl <- range(ens.rmse[element,], root.mean.var[element,]) 
        
        plot(ens.rmse[element,], type = "l", col = "red", lty = 1, main = ttl, xlab = "", ylab = "", ylim = yl)
        lines(root.mean.var[element,], col = "blue3", lty = 2)
    }
    
    pdf("./Plots/RMSE-vs-spread.pdf"); {
        par(mfrow = c(2,3), oma = c(0.5, 0.5, 2, 0.5), mar = c(2,2,3,1))
        qplot("pc1"); qplot("pc2"); qplot("pc3")
        qplot("temp.n"); qplot("temp.s")
        
        plot.new()
        legend("left", lty = c(1,2), col = c("red", "blue"), bty = "n", cex = 1.1,
               legend = c("RMSE of ensemble mean", "Ensemble spread"))
        
        mtext("RMSE vs spread at each forecast lead time", outer = TRUE, cex = 1)
    }; dev.off()

}

####################################################################################################

# NCEP_loading.R                                                                                ####

control_ts <- load.data("./Data/NCEPctrl_mslp_timeseries.rda")
ensemble_ts <- load.data("./Data/NCEPpert_mslp_timeseries.rda")

# loading the temperature of south and north UK respectively.
northctrltemp <- load.data("./Data/NCEPctrl_temp2m_7yr_north24hrave.rda")
southctrltemp <- load.data("./Data/NCEPctrl_temp2m_7yr_south24hrave.rda")

# still no files for start & lead times. Create data manually for now

####################################################################################################
