
####################################################################################################

# LOAD PREREQUISITES                                                                            ####

# timestamp processing
{
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
    
    model_time1 <- model_time <- load.data("./Data/ECMWF_europe_starttime.rda")[,1:7]
    leadtime <- load.data("./Data/ECMWF_europe_leadtime.rda")
}

# load observations
{
    timeseries <- load.data("./Data/ERAint_pca_timeseries.rda")
    PC1_obs <- timeseries[2521:3150,1] 
    PC2_obs <- timeseries[2521:3150,2] 
    TN_obs <- c(load.data("./Data/ERAint_temp2m_36yr_north24hrave.rda")[seq(3,360,4),29:35]-273.15)
    TS_obs <- c(load.data("./Data/ERAint_temp2m_36yr_south24hrave.rda")[seq(3,360,4),29:35]-273.15)
}

# The perturbed forecasts for 14 lead times and 630 days: ECMWF
{
    # load all necessary data
    {
        control_ts <- load.data("./Data/ECMWFctrl_mslp_timeseries.rda")
        ensemble_ts <- load.data("./Data/ECMWFpert_timeseries.rda")
        
        northperttemp <- load.data("./Data/ECMWFpert_temp2m_7yr_north24hrave.rda")
        southperttemp <- load.data("./Data/ECMWFpert_temp2m_7yr_south24hrave.rda")
        

    }
    
    # For the perturbed forecast
    pert_ECMWF <- function(day.ahead){
        PC1 <- as.vector(ensemble_ts[djf.index2(model_time)-day.ahead,1,day.ahead*4+1,]) 
        PC2 <- as.vector(ensemble_ts[djf.index2(model_time)-day.ahead,2,day.ahead*4+1,])
        tempnorth <- as.vector(northperttemp[day.ahead*4+1,djf.index1(model_time)-day.ahead,,])-273.15
        tempsouth <- as.vector(southperttemp[day.ahead*4+1,djf.index1(model_time)-day.ahead,,])-273.15
        results <- data.frame(PC1,PC2,tempnorth,tempsouth)
        return(results)
    }
    
    PC1.pert.leadtime_ECMWF <- matrix(nrow = 31500,ncol=14)
    PC2.pert.leadtime_ECMWF <- matrix(nrow = 31500,ncol=14)
    TN.pert.leadtime_ECMWF <- matrix(nrow = 31500,ncol=14)
    TS.pert.leadtime_ECMWF <- matrix(nrow = 31500,ncol=14)
    
    for (i in 1:14) {
        PC1.pert.leadtime_ECMWF[,i] <- pert_ECMWF(i)$PC1
        PC2.pert.leadtime_ECMWF[,i] <- pert_ECMWF(i-1)$PC2
        TN.pert.leadtime_ECMWF[,i] <- pert_ECMWF(i)$tempnorth
        TS.pert.leadtime_ECMWF[,i] <- pert_ECMWF(i)$tempsouth
    }
}

# The perturbed forecasts for 14 lead times and 630 days: NCEP
{
    # load all necessary data
    {
        control_ts <- load.data("./Data/NCEPctrl_mslp_timeseries.rda")
        ensemble_ts <- load.data("./Data/NCEPpert_mslp_timeseries.rda")
        
        northperttemp <- load.data("./Data/NCEPpert_temp2m_7yr_north24hrave.rda")
        southperttemp <- load.data("./Data/NCEPpert_temp2m_7yr_south24hrave.rda")
   }
    
    # For the perturbed forecast
    pert_NCEP <- function(day.ahead){
        PC1 <- as.vector(ensemble_ts[djf.index2(model_time)-day.ahead,1,day.ahead*4+1,]) 
        PC2 <- as.vector(ensemble_ts[djf.index2(model_time)-day.ahead,2,day.ahead*4+1,])
        tempnorth <- as.vector(northperttemp[day.ahead*4+1,djf.index1(model_time)-day.ahead,,])-273.15
        tempsouth <- as.vector(southperttemp[day.ahead*4+1,djf.index1(model_time)-day.ahead,,])-273.15
        results <- data.frame(PC1,PC2,tempnorth,tempsouth)
        return(results)
    }
    
    PC1.pert.leadtime_NCEP <- matrix(nrow = 12600,ncol=14)
    PC2.pert.leadtime_NCEP <- matrix(nrow = 12600,ncol=14)
    TN.pert.leadtime_NCEP <- matrix(nrow = 12600,ncol=14)
    TS.pert.leadtime_NCEP <- matrix(nrow = 12600,ncol=14)
    
    for (i in 1:14) {
        PC1.pert.leadtime_NCEP[,i] <- pert_NCEP(i)$PC1
        PC2.pert.leadtime_NCEP[,i] <- pert_NCEP(i-1)$PC2
        TN.pert.leadtime_NCEP[,i] <- pert_NCEP(i)$tempnorth
        TS.pert.leadtime_NCEP[,i] <- pert_NCEP(i)$tempsouth
    }
}

# The perturbed forecasts for 14 lead times and 630 days: UKMO
{
    # load all necessary data
    {
        control_ts <- load.data("./Data/UKMOctrl_mslp_timeseries.rda")
        ensemble_ts <- load.data("./Data/UKMOpert_mslp_timeseries.rda")
        
        northperttemp <- load.data("./Data/UKMOpert_temp2m_7yr_north24hrave.rda")
        southperttemp <- load.data("./Data/UKMOpert_temp2m_7yr_south24hrave.rda")
    }
    
    # For the perturbed forecast
    pert_UKMO <- function(day.ahead){
        PC1 <- as.vector(ensemble_ts[djf.index2(model_time)-day.ahead,1,day.ahead*4+1,]) 
        PC2 <- as.vector(ensemble_ts[djf.index2(model_time)-day.ahead,2,day.ahead*4+1,])
        tempnorth <- as.vector(northperttemp[day.ahead*4+1,djf.index1(model_time)-day.ahead,,])-273.15
        tempsouth <- as.vector(southperttemp[day.ahead*4+1,djf.index1(model_time)-day.ahead,,])-273.15
        results <- data.frame(PC1,PC2,tempnorth,tempsouth)
        return(results)
    }
    
    PC1.pert.leadtime_UKMO <- matrix(nrow = 14490,ncol=14)
    PC2.pert.leadtime_UKMO <- matrix(nrow = 14490,ncol=14)
    TN.pert.leadtime_UKMO <- matrix(nrow = 14490,ncol=14)
    TS.pert.leadtime_UKMO <- matrix(nrow = 14490,ncol=14)
    
    for (i in 1:14) {
        PC1.pert.leadtime_UKMO[,i] <- pert_UKMO(i)$PC1
        PC2.pert.leadtime_UKMO[,i] <- pert_UKMO(i-1)$PC2
        TN.pert.leadtime_UKMO[,i] <- pert_UKMO(i)$tempnorth
        TS.pert.leadtime_UKMO[,i] <- pert_UKMO(i)$tempsouth
    }
}

# superensemble forecast
{
    ensemble_ECMWF <- load.data("./Data/ECMWFpert_timeseries.rda")
    ensemble_NCEP <- load.data("./Data/NCEPpert_mslp_timeseries.rda")
    ensemble_UKMO <- load.data("./Data/UKMOpert_mslp_timeseries.rda")
    
    # ECWMF
    TN_ECMWF <- load.data("./Data/ECMWFpert_temp2m_7yr_north24hrave.rda")
    TS_ECMWF <- load.data("./Data/ECMWFpert_temp2m_7yr_south24hrave.rda")

    # NCEP
    TN_NCEP <- load.data("./Data/NCEPpert_temp2m_7yr_north24hrave.rda")
    TS_NCEP <- load.data("./Data/NCEPpert_temp2m_7yr_south24hrave.rda")

    # UKMO
    TN_UKMO <- load.data("./Data/UKMOpert_temp2m_7yr_north24hrave.rda")
    TS_UKMO <- load.data("./Data/UKMOpert_temp2m_7yr_south24hrave.rda")

    # PC1 and PC2
    ensemble_ALL <- abind::abind(ensemble_ECMWF,ensemble_NCEP,ensemble_UKMO)
    TN_ALL <- abind::abind(TN_ECMWF,TN_NCEP,TN_UKMO,along=3)
    TS_ALL <- abind::abind(TS_ECMWF,TS_NCEP,TS_UKMO,along=3)
    pert_ALL<- function(day.ahead){
        PC1 <- as.vector(ensemble_ALL[djf.index2(model_time1)-day.ahead,1,day.ahead*4+1,]) 
        PC2 <- as.vector(ensemble_ALL[djf.index2(model_time1)-day.ahead,2,day.ahead*4+1,])
        TN <- as.vector(TN_ALL[day.ahead*4+1,djf.index1(model_time)-day.ahead,,])-273.15
        TS <- as.vector(TS_ALL[day.ahead*4+1,djf.index1(model_time)-day.ahead,,])-273.15
        results <- data.frame(PC1,PC2,TN,TS)
        return(results)
    }
    PC1.pert.leadtime_ALL <- matrix(nrow = 58590,ncol=14)
    for (i in 1:14) {
        PC1.pert.leadtime_ALL[,i] <- pert_ALL(i)$PC1
    }
    
    PC2.pert.leadtime_ALL <- matrix(nrow = 58590,ncol=14)
    for (i in 1:14) {
        PC2.pert.leadtime_ALL[,i] <- pert_ALL(i)$PC2
    }
    
    # Temp
    TN.pert.leadtime_ALL <- matrix(nrow = 58590,ncol=14)
    for (i in 1:14) {
        TN.pert.leadtime_ALL[,i] <- pert_ALL(i)$TN
    }
    
    TS.pert.leadtime_ALL <- matrix(nrow = 58590,ncol=14)
    for (i in 1:14) {
        TS.pert.leadtime_ALL[,i] <- pert_ALL(i)$TS
    }
}

####################################################################################################

# CRPS.R                                                                                        ####

# temp.n
{
TN_ECMWF_leadtime <- array(0,dim = c(630,50,14))
for (i in 1:630){
    for (j in 1:14)
        TN_ECMWF_leadtime[i,,j]<-TN.pert.leadtime_ECMWF[seq(4500*(ceiling(i/90)-1)+(i-1)%%90+1,4500*ceiling(i/90),90),j]
}

ECMWF_TN<- array(0,dim = c(1,14))
for (j in 1:14){
    ECMWF_TN[,j] <- verification::crpsDecomposition(TN_obs,TN_ECMWF_leadtime[,,j])$CRPS
}

TN_NCEP_leadtime <- array(0,dim = c(630,20,14))
for (i in 1:630){
    for (j in 1:14)
        TN_NCEP_leadtime[i,,j]<-TN.pert.leadtime_NCEP[seq(1800*(ceiling(i/90)-1)+(i-1)%%90+1,1800*ceiling(i/90),90),j]
}
NCEP_TN<- array(0,dim = c(1,14))
for (j in 1:14){
    NCEP_TN[,j] <- verification::crpsDecomposition(TN_obs,TN_NCEP_leadtime[,,j])$CRPS
}

TN_UKMO_leadtime <- array(0,dim = c(630,23,14))
for (i in 1:630){
    for (j in 1:14)
        TN_UKMO_leadtime[i,,j]<-TN.pert.leadtime_UKMO[seq(2070*(ceiling(i/90)-1)+(i-1)%%90+1,2070*ceiling(i/90),90),j]
}
UKMO_TN<- array(0,dim = c(1,14))
for (j in 1:14){
    UKMO_TN[,j] <- verification::crpsDecomposition(TN_obs,TN_UKMO_leadtime[,,j])$CRPS
}

TN_ALL_leadtime <- array(0,dim = c(630,93,14))
for (i in 1:630){
    for (j in 1:14)
        TN_ALL_leadtime[i,,j]<-TN.pert.leadtime_ALL[seq(8370*(ceiling(i/90)-1)+(i-1)%%90+1,8370*ceiling(i/90),90),j]
}
ALL_TN<- array(0,dim = c(1,14))
for (j in 1:14){
    ALL_TN[,j] <- verification::crpsDecomposition(TN_obs,TN_ALL_leadtime[,,j])$CRPS
}
}

# temp.s
{
    TS_ECMWF_leadtime <- array(0,dim = c(630,50,14))
    for (i in 1:630){
        for (j in 1:14)
            TS_ECMWF_leadtime[i,,j]<-TS.pert.leadtime_ECMWF[seq(4500*(ceiling(i/90)-1)+(i-1)%%90+1,4500*ceiling(i/90),90),j]
    }
    ECMWF_TS<- array(0,dim = c(1,14))
    for (j in 1:14){
        ECMWF_TS[,j] <- verification::crpsDecomposition(TS_obs,TS_ECMWF_leadtime[,,j])$CRPS
    }
    
    TS_NCEP_leadtime <- array(0,dim = c(630,20,14))
    for (i in 1:630){
        for (j in 1:14)
            TS_NCEP_leadtime[i,,j]<-TS.pert.leadtime_NCEP[seq(1800*(ceiling(i/90)-1)+(i-1)%%90+1,1800*ceiling(i/90),90),j]
    }
    NCEP_TS<- array(0,dim = c(1,14))
    for (j in 1:14){
        NCEP_TS[,j] <- verification::crpsDecomposition(TS_obs,TS_NCEP_leadtime[,,j])$CRPS
    }
    
    TS_UKMO_leadtime <- array(0,dim = c(630,23,14))
    for (i in 1:630){
        for (j in 1:14)
            TS_UKMO_leadtime[i,,j]<-TS.pert.leadtime_UKMO[seq(2070*(ceiling(i/90)-1)+(i-1)%%90+1,2070*ceiling(i/90),90),j]
    }
    UKMO_TS<- array(0,dim = c(1,14))
    for (j in 1:14){
        UKMO_TS[,j] <- verification::crpsDecomposition(TS_obs,TS_UKMO_leadtime[,,j])$CRPS
    }
    
    TS_ALL_leadtime <- array(0,dim = c(630,93,14))
    for (i in 1:630){
        for (j in 1:14)
            TS_ALL_leadtime[i,,j]<-TS.pert.leadtime_ALL[seq(8370*(ceiling(i/90)-1)+(i-1)%%90+1,8370*ceiling(i/90),90),j]
    }
    ALL_TS<- array(0,dim = c(1,14))
    for (j in 1:14){
        ALL_TS[,j] <- verification::crpsDecomposition(TS_obs,TS_ALL_leadtime[,,j])$CRPS
    }
}

# pc1
{
    PC1_ECMWF_leadtime <- array(0,dim = c(630,50,14))
    for (i in 1:630){
        for (j in 1:14)
            PC1_ECMWF_leadtime[i,,j]<-PC1.pert.leadtime_ECMWF[seq(i,31500,630),j]
    }
    ECMWF_PC1<- array(0,dim = c(1,14))
    for (j in 1:14){
        ECMWF_PC1[,j] <- verification::crpsDecomposition(PC1_obs,PC1_ECMWF_leadtime[,,j])$CRPS
    }
    
    PC1_NCEP_leadtime <- array(0,dim = c(630,20,14))
    for (i in 1:630){
        for (j in 1:14)
            PC1_NCEP_leadtime[i,,j]<-PC1.pert.leadtime_NCEP[seq(i,12600,630),j]
    }
    NCEP_PC1<- array(0,dim = c(1,14))
    for (j in 1:14){
        NCEP_PC1[,j] <- verification::crpsDecomposition(PC1_obs,PC1_NCEP_leadtime[,,j])$CRPS
    }
    
    PC1_UKMO_leadtime <- array(0,dim = c(630,23,14))
    for (i in 1:630){
        for (j in 1:14)
            PC1_UKMO_leadtime[i,,j]<-PC1.pert.leadtime_UKMO[seq(i,14490,630),j]
    }
    UKMO_PC1<- array(0,dim = c(1,14))
    for (j in 1:14){
        UKMO_PC1[,j] <- verification::crpsDecomposition(PC1_obs,PC1_UKMO_leadtime[,,j])$CRPS
    }
    
    PC1_ALL_leadtime <- array(0,dim = c(630,93,14))
    for (i in 1:630){
        for (j in 1:14)
            PC1_ALL_leadtime[i,,j]<-PC1.pert.leadtime_ALL[seq(i,58590,630),j]
    }
    ALL_PC1<- array(0,dim = c(1,14))
    for (j in 1:14){
        ALL_PC1[,j] <- verification::crpsDecomposition(PC1_obs,PC1_ALL_leadtime[,,j])$CRPS
    }
}

# pc2
{
    PC2_ECMWF_leadtime <- array(0,dim = c(630,50,14))
    for (i in 1:630){
        for (j in 1:14)
            PC2_ECMWF_leadtime[i,,j]<-PC2.pert.leadtime_ECMWF[seq(i,31500,630),j]
    }
    ECMWF_PC2<- array(0,dim = c(1,14))
    for (j in 1:14){
        ECMWF_PC2[,j] <- verification::crpsDecomposition(PC2_obs,PC2_ECMWF_leadtime[,,j])$CRPS
    }
    
    PC2_NCEP_leadtime <- array(0,dim = c(630,20,14))
    for (i in 1:630){
        for (j in 1:14)
            PC2_NCEP_leadtime[i,,j]<-PC2.pert.leadtime_NCEP[seq(i,12600,630),j]
    }
    NCEP_PC2<- array(0,dim = c(1,14))
    for (j in 1:14){
        NCEP_PC2[,j] <- verification::crpsDecomposition(PC2_obs,PC2_NCEP_leadtime[,,j])$CRPS
    }
    
    PC2_UKMO_leadtime <- array(0,dim = c(630,23,14))
    for (i in 1:630){
        for (j in 1:14)
            PC2_UKMO_leadtime[i,,j]<-PC2.pert.leadtime_UKMO[seq(i,14490,630),j]
    }
    UKMO_PC2<- array(0,dim = c(1,14))
    for (j in 1:14){
        UKMO_PC2[,j] <- verification::crpsDecomposition(PC2_obs,PC2_UKMO_leadtime[,,j])$CRPS
    }
    
    PC2_ALL_leadtime <- array(0,dim = c(630,93,14))
    for (i in 1:630){
        for (j in 1:14)
            PC2_ALL_leadtime[i,,j]<-PC2.pert.leadtime_ALL[seq(i,58590,630),j]
    }
    ALL_PC2<- array(0,dim = c(1,14))
    for (j in 1:14){
        ALL_PC2[,j] <- verification::crpsDecomposition(PC2_obs,PC2_ALL_leadtime[,,j])$CRPS
    }
}

