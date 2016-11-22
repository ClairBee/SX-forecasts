
####################################################################################################

# ECMWF_loading                                                                                 ####

# Loading the data

# loading two Principle Compoments(PC) of pressure field. There are 3 PCs inside, but just
# using the first two PCs.
load("./Data/ECMWFctrl_mslp_timeseries.rda")
load("./Data/ECMWFpert_timeseries.rda")

# loading the temperature of south and north UK respectively.
load("./Data/ECMWFctrl_temp2m_7yr_south24hrave.rda")
southctrltemp <- model_var
load("./Data/ECMWFctrl_temp2m_7yr_north24hrave.rda")
northctrltemp <- model_var

# loading the start time and lead time
# Valid time = start time + lead time
load("./Data/ECMWF_europe_starttime.rda")
model_time1 <- model_time
load("./Data/ECMWF_europe_leadtime.rda")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# added for compatibility
{
    model_time <- model_time[,1:7]
    model_time1 <- model_time1[,1:7]
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions

# Data just for December, January and February (DJF) in a single year.(31 + 31 + 28 = 90 days)
# Not include 29th Feb.

# produces a 1:90 index for a single year (only 90 days in first year)
# the unit of starttime is hours.
# start_date: as.Date(how many days from the origin day, origin day/begin day), only need 'year'...
#             just focus on the DJF, so valid date should start from 12-01.
# end_date:   the last number of starttime, using dimension(column and row) to find that.
# startDJF:   start_time - origin in hours
# DJFindex:   first column find round(startDJF) --> last column find round(endDJF)
# Finally, the 16th number to 105th number is the data for DJF.
djf.index1 <- function(starttime){ 
    origin <- "1800-01-01"
    start_date <- paste0(format(as.Date(starttime[1,1]/24, origin), "%Y"),"-12-01")
    end_date <- paste0(format(as.Date(starttime[dim(starttime)[1],dim(starttime)[2]]/24, 
                                      origin), "%Y"),"-02-28")
    startDJF <- as.numeric(difftime(start_date, origin, units = "hours"))+12 # Midday, 12o'clock noon
    endDJF <- as.numeric(difftime(end_date, origin, units = "hours"))+12
    DJFindex <- match(round(startDJF), starttime[,1]):match(round(endDJF),starttime[,dim(starttime)[2]])
    return(DJFindex)
}

# produces a 1:630 index for 7 years (7 years * 90 days = 630 days)
djf.index2 <- function(starttime){ 
    origin <- "1800-01-01"
    start_date <- paste0(format(as.Date(starttime[1,1]/24, origin), "%Y"),"-12-01")
    end_date <- paste0(format(as.Date(starttime[dim(starttime)[1],dim(starttime)[2]]/24, origin), "%Y"),"-02-28")
    startDJF <- as.numeric(difftime(start_date, origin, units = "hours"))+12
    endDJF <- as.numeric(difftime(end_date, origin, units = "hours"))+12
    DJFindex1 <- match(round(startDJF), starttime[,1]):match(round(endDJF),starttime[,dim(starttime)[2]])
    dim.row <- dim(starttime)[1]
    djfindex2 <- dim.row*rep(0:6,each=length(DJFindex1)) + DJFindex1
    return(djfindex2)
}

######################################################################################
# task 1: create data frames for the *control* forecast
# The lead time for pressure is 61, but the lead time for temperature is 58,
# so we need the shorter lead time which is 58, then lead time is 6hr, 12hr, 18hr and 24hr, 
# we need should be the 1st, 5th, 9th,..., 14th data. ((58-1)/4=14.25)

# data frame 1: observations
# Note: focus on the structure of each dimension.
ctrl_obs <- data.frame(
    day.start = as.Date(as.vector(model_time1[djf.index1(model_time1),])/24,"1800-01-01"),
    Date.numeric = as.vector(model_time1[djf.index1(model_time1),]), # model_time in hours
    Date.reasonable = as.Date(as.vector(model_time1[djf.index1(model_time1),])/24,"1800-01-01"), 
    # model_time in date
    PC1 = as.vector(control_ts[djf.index2(model_time1),1,1]), 
    # timeseries control_ts and ensemble_ts, the last 1 is lead time. First PC
    PC2 = as.vector(control_ts[djf.index2(model_time1),2,1]), # Second PC
    temp.north = as.vector(northctrltemp[1,djf.index1(model_time1),])-273.15,# change Kelvin to Celsius.
    temp.south = as.vector(southctrltemp[1,djf.index1(model_time1),])-273.15
)

# data frame 2: 1 day ahead. 
ctrl_1day <- data.frame(
    day.start=as.Date(as.vector(model_time1[djf.index1(model_time1)-1,])/24,"1800-01-01"),
    Date.numeric=as.vector(model_time1[djf.index1(model_time1),]),
    Date.reasonable=as.Date(as.vector(model_time1[djf.index1(model_time1),])/24,"1800-01-01"), 
    PC1=as.vector(control_ts[djf.index2(model_time1)-1,1,5]), 
    PC2=as.vector(control_ts[djf.index2(model_time1)-1,2,5]),
    temp.north=as.vector(northctrltemp[5,djf.index1(model_time1)-1,])-273.15, 
    temp.south=as.vector(southctrltemp[5,djf.index1(model_time1)-1,])-273.15
)

###### In general, data frame 5 ...15: 4...14 day ahead. 
day.ahead <-1
ctrl_aheadday <- data.frame(
    day.start=as.Date(as.vector(model_time1[djf.index1(model_time1)-day.ahead,])/24,"1800-01-01"),
    Date.numeric=as.vector(model_time1[djf.index1(model_time1),]),
    Date.reasonable=as.Date(as.vector(model_time1[djf.index1(model_time1),])/24,"1800-01-01"), 
    PC1=as.vector(control_ts[djf.index2(model_time1)-day.ahead,1,day.ahead*4+1]), 
    PC2=as.vector(control_ts[djf.index2(model_time1)-day.ahead,2,day.ahead*4+1]),
    temp.north=as.vector(northctrltemp[day.ahead*4+1,djf.index1(model_time1)-day.ahead,])-273.15,
    temp.south=as.vector(southctrltemp[day.ahead*4+1,djf.index1(model_time1)-day.ahead,])-273.15
)

##########################################################################################
# task 2: Ensembles, or *perturbed* forecasts incorporated into the data.frames
load("./Data/ECMWFpert_temp2m_7yr_north24hrave.rda")
northperttemp <- model_var
load("./Data/ECMWFpert_temp2m_7yr_south24hrave.rda")
southperttemp <- model_var

# For lead time 0
pert_obs.all <- data.frame(
    day.start = as.Date(as.vector(model_time1[djf.index1(model_time1),])/24,"1800-01-01"),
    Date.numeric = as.vector(model_time1[djf.index1(model_time1),]), # model_time in hours
    Date.reasonable = as.Date(as.vector(model_time1[djf.index1(model_time1),])/24,"1800-01-01"), 
    # model_time in date
    PC1 = as.vector(ensemble_ts[djf.index2(model_time1),1,1,]), 
    # timeseries control_ts and ensemble_ts, the last 1 is lead time. First PC
    PC2 = as.vector(ensemble_ts[djf.index2(model_time1),2,1,]), # Second PC
    temp.north = as.vector(northperttemp[1,djf.index1(model_time1),,])-273.15,
    # change Kelvin to Celsius.
    temp.south = as.vector(southperttemp[1,djf.index1(model_time1),,])-273.15
)

###### In general, day.ahead should be from 0 to 14.
day.ahead <- 1
pert_aheadday.all <- data.frame(
    day.start=as.Date(as.vector(model_time1[djf.index1(model_time1)-day.ahead,])/24,"1800-01-01"),
    Date.numeric=as.vector(model_time1[djf.index1(model_time1),]),
    Date.reasonable=as.Date(as.vector(model_time1[djf.index1(model_time1),])/24,"1800-01-01"), 
    PC1=as.vector(ensemble_ts[djf.index2(model_time1)-day.ahead,1,day.ahead*4+1,]), 
    PC2=as.vector(ensemble_ts[djf.index2(model_time1)-day.ahead,2,day.ahead*4+1,]),
    temp.north=as.vector(northperttemp[day.ahead*4+1,djf.index1(model_time1)-day.ahead,,])-273.15,
    temp.south=as.vector(southperttemp[day.ahead*4+1,djf.index1(model_time1)-day.ahead,,])-273.15
)

#=================================================================================

#                         SPATIAL PLOTTING REMOVED



####################################################################################################

# ECMWF_Mean Error and RMSE                                                                     ####

# For the control forecast 
###################################################################################
# Generally, Using function to get all variables for each lead time from 1 to 14 days.
ctrl_ECMWF <- function(day.ahead){
    PC1 <- as.vector(control_ts[djf.index2(model_time1)-day.ahead,1,day.ahead*4+1])
    PC2 <- as.vector(control_ts[djf.index2(model_time1)-day.ahead,2,day.ahead*4+1])
    tempnorth <- as.vector(northctrltemp[day.ahead*4+1,djf.index1(model_time1)-day.ahead,])-273.15
    tempsouth <- as.vector(southctrltemp[day.ahead*4+1,djf.index1(model_time1)-day.ahead,])-273.15
    results <- data.frame(PC1,PC2,tempnorth,tempsouth)
    return(results)
}

# The control forecasts for 14 lead times and 630 days (90 winter days * 7 years).
PC1.ctrl.leadtime_ECMWF <- matrix(nrow = 630,ncol=14)
for (i in 1:14) {
    PC1.ctrl.leadtime_ECMWF[,i] <- ctrl_ECMWF(i)$PC1
}

PC2.ctrl.leadtime_ECMWF <- matrix(nrow = 630,ncol=14)
for (i in 1:14) {
    PC2.ctrl.leadtime_ECMWF[,i] <- ctrl_ECMWF(i)$PC2
}

TN.ctrl.leadtime_ECMWF <- matrix(nrow = 630,ncol=14)
for (i in 1:14) {
    TN.ctrl.leadtime_ECMWF[,i] <- ctrl_ECMWF(i)$tempnorth
}

TS.ctrl.leadtime_ECMWF <- matrix(nrow = 630,ncol=14)
for (i in 1:14) {
    TS.ctrl.leadtime_ECMWF[,i] <- ctrl_ECMWF(i)$tempsouth
}

#####################################################################################
# For the perturbed forecast
pert_ECMWF <- function(day.ahead){
    PC1 <- as.vector(ensemble_ts[djf.index2(model_time)-day.ahead,1,day.ahead*4+1,]) 
    PC2 <- as.vector(ensemble_ts[djf.index2(model_time)-day.ahead,2,day.ahead*4+1,])
    tempnorth <- as.vector(northperttemp[day.ahead*4+1,djf.index1(model_time)-day.ahead,,])-273.15
    tempsouth <- as.vector(southperttemp[day.ahead*4+1,djf.index1(model_time)-day.ahead,,])-273.15
    results <- data.frame(PC1,PC2,tempnorth,tempsouth)
    return(results)
}

# The perturbed forecasts for 14 lead times and 630 days.
PC1.pert.leadtime_ECMWF <- matrix(nrow = 31500,ncol=14)
for (i in 1:14) {
    PC1.pert.leadtime_ECMWF[,i] <- pert_ECMWF(i)$PC1
}

PC2.pert.leadtime_ECMWF <- matrix(nrow = 31500,ncol=14)
for (i in 1:14) {
    PC2.pert.leadtime_ECMWF[,i] <- pert_ECMWF(i-1)$PC2
}

TN.pert.leadtime_ECMWF <- matrix(nrow = 31500,ncol=14)
for (i in 1:14) {
    TN.pert.leadtime_ECMWF[,i] <- pert_ECMWF(i)$tempnorth
}

TS.pert.leadtime_ECMWF <- matrix(nrow = 31500,ncol=14)
for (i in 1:14) {
    TS.pert.leadtime_ECMWF[,i] <- pert_ECMWF(i)$tempsouth
}

####################################################################################
# Observation from ERAint
load("./Data/ERAint_europe_starttime.rda")
ERAint.starttime_obs <- model_time
day.start = as.Date(as.vector(ERAint.starttime_obs[djf.index1(ERAint.starttime_obs),])/24,"1800-01-01")
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

# for perturbed forecast
# first get the ensemble mean for each variable and each lead time
# PC1
PC1.pert.leadtime.ensemblemean_ECMWF <- matrix(nrow=630,ncol=14)
for (i in 1:630){
    for (j in 1:14){
        PC1.pert.leadtime.ensemblemean_ECMWF[i,j] <- mean(as.vector(ensemble_ts[djf.index2(model_time1)[i]-j,1,j*4+1,]))
    }
}

#PC2
PC2.pert.leadtime.ensemblemean_ECMWF <- matrix(nrow=630,ncol=14)
for (i in 1:630){
    for (j in 1:14){
        PC2.pert.leadtime.ensemblemean_ECMWF[i,j] <- mean(as.vector(ensemble_ts[djf.index2(model_time1)[i]-j,2,j*4+1,]))
    }
}

#TN
TN.pert.leadtime.ensemblemean_ECMWF <- matrix(nrow=630,ncol=14)
for (i in 1:630){
    for (j in 1:14){
        TN.pert.leadtime.ensemblemean_ECMWF[i,j] <- mean(as.vector(northperttemp[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,ceiling(i/90)])-273.15)
    }
}

#TS
TS.pert.leadtime.ensemblemean_ECMWF <- matrix(nrow=630,ncol=14)
for (i in 1:630){
    for (j in 1:14){
        TS.pert.leadtime.ensemblemean_ECMWF[i,j] <- mean(as.vector(southperttemp[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,ceiling(i/90)])-273.15)
    }
}

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
for (i in 1:630){
    for (j in 1:14){
        variance_TN[i,j] <- var(northperttemp[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,ceiling(i/90)]-273.15)
    }
}

averagevariance_TN <- colMeans(variance_TN)

Root_TN_ECMWF <- sqrt(averagevariance_TN)

###RMSE for mean of ensemble (should be divided into 50)
RMSE_TN_ECMWF <- sqrt(colMeans((TN.pert.leadtime.ensemblemean_ECMWF - TN_obs)^2))

######################################################################################
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

#=================================================================================

#                         PLOTS REMOVED

#=================================================================================