
# designed as source file to replicate Sichun's analysis.
# adapted to compensate for missing date reference files

####################################################################################################

# can't load the data as Sichun does because of the missing date files
# forecast data therefore taken from my array instead. 
# will need to check these variables against SX versions when all files available

replicate.forecast.data <- function() {
    # variable names retained, method of obtaining them corrected
    
    fc <<- readRDS("./Data/ECMWF-forecasts.rds")[,16:105,,2:15,]
    
    # The control forecasts for 14 lead times and 630 days (90 winter days * 7 years).
    PC1.ctrl.leadtime_ECMWF <<- array(fc["pc1",,,,"c"], dim = c(90 * 7, 14))
    PC2.ctrl.leadtime_ECMWF <<- array(fc["pc2",,,,"c"], dim = c(90 * 7, 14))
    PC3.ctrl.leadtime_ECMWF <<- array(fc["pc3",,,,"c"], dim = c(90 * 7, 14))
    TN.ctrl.leadtime_ECMWF <<- array(fc["temp.n",,,,"c"], dim = c(90 * 7, 14))
    TS.ctrl.leadtime_ECMWF <<- array(fc["temp.s",,,,"c"], dim = c(90 * 7, 14))
    
    # The perturbed forecasts for 14 lead times and 630 days.
    PC1.pert.leadtime_ECMWF <<- array(aperm(fc, c(1,2,3,5,4))["pc1",,,2:51,], dim = c(90 * 7 * 50, 14))
    PC2.pert.leadtime_ECMWF <<- array(aperm(fc, c(1,2,3,5,4))["pc2",,,2:51,], dim = c(90 * 7 * 50, 14))
    PC3.pert.leadtime_ECMWF <<- array(aperm(fc, c(1,2,3,5,4))["pc3",,,2:51,], dim = c(90 * 7 * 50, 14))
    TN.pert.leadtime_ECMWF <<- array(aperm(fc, c(1,2,3,5,4))["temp.n",,,2:51,], dim = c(90 * 7 * 50, 14))
    TS.pert.leadtime_ECMWF <<- array(aperm(fc, c(1,2,3,5,4))["temp.s",,,2:51,], dim = c(90 * 7 * 50, 14))
}

replicate.forecast.data()

####################################################################################
# Observation from ERAint
#load("./Data/ERAint_europe_starttime.rda")
#ERAint.starttime_obs <- model_time
#day.start = as.Date(as.vector(ERAint.starttime_obs[djf.index1(ERAint.starttime_obs),])/24,"1800-01-01")
# We can find what date we need is 2007/2008-2013/2014.

# PC1 and PC2
load("./Data/ERAint_pca_timeseries.rda")
PC1_obs <- timeseries[2521:3150,1] 
# 3240 - 720 + 1 = 2611, this is the start time for 2007/2008-2013/2014.
PC2_obs <- timeseries[2521:3150,2] 

# Temp in north UK
load("./Data/ERAint_temp2m_36yr_north24hrave.rda")
TN_obs <- as.vector(model_var[seq(3,360,4),29:35]-273.15)

# Temp in south UK
load("./Data/ERAint_temp2m_36yr_south24hrave.rda")
TS_obs <- as.vector(model_var[seq(3,360,4),29:35])-273.15

#####################################################################################
# for control forecast
# mean error = sum(PC1.ctrl.leadtime-PC1.ctrl_obs)/630
mean.error.ctrl_PC1_ECMWF <- colMeans(PC1.ctrl.leadtime_ECMWF - PC1_obs)
mean.error.ctrl_PC2_ECMWF <- colMeans(PC2.ctrl.leadtime_ECMWF - PC2_obs)
mean.error.ctrl_TN_ECMWF <- colMeans(TN.ctrl.leadtime_ECMWF - TN_obs)
mean.error.ctrl_TS_ECMWF <- colMeans(TS.ctrl.leadtime_ECMWF - TS_obs)


get.ensemble.means <- function() {
    # first get the ensemble mean for each variable and each lead time
    PC1.pert.leadtime.ensemblemean_ECMWF <<- apply(array(PC1.pert.leadtime_ECMWF, dim = c(630, 50, 14)), c(1,3), mean)
    PC2.pert.leadtime.ensemblemean_ECMWF <<- apply(array(PC2.pert.leadtime_ECMWF, dim = c(630, 50, 14)), c(1,3), mean)
    PC3.pert.leadtime.ensemblemean_ECMWF <<- apply(array(PC3.pert.leadtime_ECMWF, dim = c(630, 50, 14)), c(1,3), mean)
    TN.pert.leadtime.ensemblemean_ECMWF <<- apply(array(TN.pert.leadtime_ECMWF, dim = c(630, 50, 14)), c(1,3), mean)
    TS.pert.leadtime.ensemblemean_ECMWF <<- apply(array(TS.pert.leadtime_ECMWF, dim = c(630, 50, 14)), c(1,3), mean)
}
get.ensemble.means()


mean.error.pert_PC1_ECMWF <- colMeans(PC1.pert.leadtime.ensemblemean_ECMWF - PC1_obs)
mean.error.pert_PC2_ECMWF <- colMeans(PC2.pert.leadtime.ensemblemean_ECMWF - PC2_obs)
mean.error.pert_TN_ECMWF <- colMeans(TN.pert.leadtime.ensemblemean_ECMWF - TN_obs)
mean.error.pert_TS_ECMWF <- colMeans(TS.pert.leadtime.ensemblemean_ECMWF - TS_obs)


#=================================================================================

#                         PLOTS REMOVED

#=================================================================================



# RMSE Root Mean Square Error

# for control forecast
# RMSE=sum(leadtime-observation)^2/630
RMSE.ctrl_PC1_ECMWF <- sqrt(colMeans((PC1.ctrl.leadtime_ECMWF - PC1_obs)^2))
RMSE.ctrl_PC2_ECMWF <- sqrt(colMeans((PC2.ctrl.leadtime_ECMWF - PC2_obs)^2))
RMSE.ctrl_TN_ECMWF<- sqrt(colMeans((TN.ctrl.leadtime_ECMWF - TN_obs)^2))
RMSE.ctrl_TS_ECMWF <- sqrt(colMeans((TS.ctrl.leadtime_ECMWF - TS_obs)^2))

# for perturbed forecast
RMSE.pert_PC1_ECMWF <- sqrt(colMeans((PC1.pert.leadtime.ensemblemean_ECMWF - PC1_obs)^2))
RMSE.pert_PC2_ECMWF <- sqrt(colMeans((PC2.pert.leadtime.ensemblemean_ECMWF - PC2_obs)^2))
RMSE.pert_TN_ECMWF <- sqrt(colMeans((TN.pert.leadtime.ensemblemean_ECMWF - TN_obs)^2))
RMSE.pert_TS_ECMWF <- sqrt(colMeans((TS.pert.leadtime.ensemblemean_ECMWF - TS_obs)^2))

#=================================================================================

#                         PLOTS REMOVED

#=================================================================================

####################################################################################################

# ECMWF_RMSE vs Spread                                                                          ####

# For each lead time and variable,
# Calculate the average variance of the ensemble and average variance for each lead time.
#
## PC1
variance_PC1 <- matrix(nrow=630,ncol=14)
for (i in 1:14){
    group <- rep(1:630, 50)
    variance_PC1[,i] <- tapply(PC1.pert.leadtime_ECMWF[,i],group, var)
}

averagevariance_PC1 <- colMeans(variance_PC1)

Root_PC1_ECMWF <- sqrt(averagevariance_PC1)

###RMSE for mean of ensemble (should be divided into 50)
RMSE_PC1_ECMWF <- sqrt(colMeans((PC1.pert.leadtime.ensemblemean_ECMWF - PC1_obs)^2))

#####################################################################################
## PC2
variance_PC2 <- matrix(nrow=630,ncol=14)
for (i in 1:14){
    group <- rep(1:630, 50)
    variance_PC2[,i] <- tapply(PC2.pert.leadtime_ECMWF[,i],group, var)
}

averagevariance_PC2 <- colMeans(variance_PC2)

Root_PC2_ECMWF <- sqrt(averagevariance_PC2)

###RMSE for mean of ensemble (should be divided into 50)
RMSE_PC2_ECMWF <- sqrt(colMeans((PC2.pert.leadtime.ensemblemean_ECMWF - PC2_obs)^2))

#####################################################################################
load("./Data/ECMWFpert_temp2m_7yr_north24hrave.rda")
northperttemp <- model_var
load("./Data/ECMWFpert_temp2m_7yr_south24hrave.rda")
southperttemp <- model_var
## Temperature north
variance_TN <- matrix(nrow=630,ncol=14)

get.temp.variances <- function() {
    npt <- array(aperm(northperttemp[c(1:14) * 4 + 1, 16:105,,], c(1,3,2,4)), dim = c(14, 50, 630))
    variance_TN <<- apply(npt, c(3, 1), var)
    
    spt <- array(aperm(southperttemp[c(1:14) * 4 + 1, 16:105,,], c(1,3,2,4)), dim = c(14, 50, 630))
    variance_TS <<- apply(spt, c(3, 1), var)
}
get.temp.variances()


averagevariance_TN <- colMeans(variance_TN)

Root_TN_ECMWF <- sqrt(averagevariance_TN)

###RMSE for mean of ensemble (should be divided into 50)
RMSE_TN_ECMWF <- sqrt(colMeans((TN.pert.leadtime.ensemblemean_ECMWF - TN_obs)^2))

######################################################################################
## Temperature south


averagevariance_TS <- colMeans(variance_TS)

Root_TS_ECMWF <- sqrt(averagevariance_TS)

###RMSE for mean of ensemble (should be divided into 50)
RMSE_TS_ECMWF <- sqrt(colMeans((TS.pert.leadtime.ensemblemean_ECMWF - TS_obs)^2))