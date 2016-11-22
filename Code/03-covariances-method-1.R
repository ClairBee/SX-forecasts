
# replication of Sichun's 'Method 1' code, with additional files loaded where necessary

####################################################################################################

# LOAD DATA                                                                                     ####

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
    
    model_time1 <- load.data("./Data/ECMWF_europe_starttime.rda")
    model_time <- load.data("./Data/ECMWF_europe_leadtime.rda")
    
}

# load all necessary data
{
    ensemble_ECMWF <- load.data("./Data/ECMWFpert_timeseries.rda")
    ensemble_NCEP <- load.data("./Data/NCEPpert_mslp_timeseries.rda")
    ensemble_UKMO <- load.data("./Data/UKMOpert_mslp_timeseries.rda")
    
    TN_ECMWF <- load.data("./Data/ECMWFpert_temp2m_7yr_north24hrave.rda")
    TS_ECMWF <- load.data("./Data/ECMWFpert_temp2m_7yr_south24hrave.rda")
    TN_NCEP <- load.data("./Data/NCEPpert_temp2m_7yr_north24hrave.rda")
    TS_NCEP <- load.data("./Data/NCEPpert_temp2m_7yr_south24hrave.rda")
    TN_UKMO <- load.data("./Data/UKMOpert_temp2m_7yr_north24hrave.rda")
    TS_UKMO <- load.data("./Data/UKMOpert_temp2m_7yr_south24hrave.rda")
    
    timeseries <- load.data("./Data/ERAint_pca_timeseries.rda")
    PC1_obs <- timeseries[2521:3150,1] 
    PC2_obs <- timeseries[2521:3150,2] 
    TN_obs <- c(load.data("./Data/ERAint_temp2m_36yr_north24hrave.rda")[seq(3,360,4),29:35]-273.15)
    TS_obs <- c(load.data("./Data/ERAint_temp2m_36yr_south24hrave.rda")[seq(3,360,4),29:35]-273.15)
}

####################################################################################################

# COVARIANCES                                                                                   ####

# ECMWF
{
    TN_ECMWF <- load.data("./Data/ECMWFpert_temp2m_7yr_north24hrave.rda")
    TS_ECMWF <- load.data("./Data/ECMWFpert_temp2m_7yr_south24hrave.rda")
    ensemble_ECMWF <- load.data("./Data/ECMWFpert_timeseries.rda")
    
    ##1. Ci: covariance for 3 simulators on each date for each lead day.
    C_ECMWF <- array(0,dim=c(4,4,630,14))
    for (j in 1:14){
        for (i in 1:630){
            C_ECMWF[,,i,j] <- cov(cbind(TN_ECMWF[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,
                                                 ceiling(i/90)]-273.15,
                                        TS_ECMWF[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,
                                                 ceiling(i/90)]-273.15,
                                        ensemble_ECMWF[djf.index2(model_time1)[i]-j,1,
                                                       j*4+1,],
                                        ensemble_ECMWF[djf.index2(model_time1)[i]-j,2,
                                                       j*4+1,]))
        }
    }
}

# NCEP
{
    TN_NCEP <- load.data("./Data/NCEPpert_temp2m_7yr_north24hrave.rda")
    TS_NCEP <- load.data("./Data/NCEPpert_temp2m_7yr_south24hrave.rda")
    ensemble_NCEP <- load.data("./Data/NCEPpert_mslp_timeseries.rda")
    
    C_NCEP <- array(0,dim=c(4,4,630,14))
    for (j in 1:14){
        for (i in 1:630){
            C_NCEP[,,i,j] <- cov(cbind(TN_NCEP[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,
                                               ceiling(i/90)]-273.15,
                                       TS_NCEP[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,
                                               ceiling(i/90)]-273.15,
                                       ensemble_NCEP[djf.index2(model_time1)[i]-j,1,
                                                     j*4+1,],
                                       ensemble_NCEP[djf.index2(model_time1)[i]-j,2,
                                                     j*4+1,]))
        }
    }
}

# UKMO
{
    TN_UKMO <- load.data("./Data/UKMOpert_temp2m_7yr_north24hrave.rda")
    TS_UKMO <- load.data("./Data/UKMOpert_temp2m_7yr_south24hrave.rda")
    ensemble_UKMO <- load.data("./Data/UKMOpert_mslp_timeseries.rda")
    
    ##1. Ci: covariance for 3 simulators on each date for each lead day.
    C_UKMO <- array(0,dim=c(4,4,630,14))
    for (j in 1:14){
        for (i in 1:630){
            C_UKMO[,,i,j] <- cov(cbind(TN_UKMO[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,
                                                 ceiling(i/90)]-273.15,
                                        TS_UKMO[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,
                                                 ceiling(i/90)]-273.15,
                                        ensemble_UKMO[djf.index2(model_time1)[i]-j,1,
                                                       j*4+1,],
                                        ensemble_UKMO[djf.index2(model_time1)[i]-j,2,
                                                       j*4+1,]))
        }
    }
}

####################################################################################################

# SIGMA                                                                                         ####

# TN
{
    northperttemp <- load.data("./Data/ECMWFpert_temp2m_7yr_north24hrave.rda")
    #TN
    TN.pert.leadtime.ensemblemean_ECMWF <- matrix(nrow=630,ncol=14)
    for (i in 1:630){
        for (j in 1:14){
            TN.pert.leadtime.ensemblemean_ECMWF[i,j] <- mean(as.vector(northperttemp[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,ceiling(i/90)])-273.15)
        }
    }
    
    northperttemp <- load.data("./Data/NCEPpert_temp2m_7yr_north24hrave.rda")
    #TN
    TN.pert.leadtime.ensemblemean_NCEP <- matrix(nrow=630,ncol=14)
    for (i in 1:630){
        for (j in 1:14){
            TN.pert.leadtime.ensemblemean_NCEP[i,j] <- mean(as.vector(northperttemp[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,ceiling(i/90)])-273.15)
        }
    }
    
    northperttemp <- load.data("./Data/UKMOpert_temp2m_7yr_north24hrave.rda")
    #TN
    TN.pert.leadtime.ensemblemean_UKMO <- matrix(nrow=630,ncol=14)
    for (i in 1:630){
        for (j in 1:14){
            TN.pert.leadtime.ensemblemean_UKMO[i,j] <- mean(as.vector(northperttemp[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,ceiling(i/90)])-273.15)
        }
    }
    
    mean_TN <- array(0,dim=c(3,630,14))
    for (j in 1:14){
        for (i in 1:630){
            mean_TN[,i,j] <- cbind(TN.pert.leadtime.ensemblemean_ECMWF[i,j],TN.pert.leadtime.ensemblemean_NCEP[i,j],TN.pert.leadtime.ensemblemean_UKMO[i,j])
        }
    }
}

# TS
{
    southperttemp <- load.data("./Data/ECMWFpert_temp2m_7yr_south24hrave.rda")
    #TS
    TS.pert.leadtime.ensemblemean_ECMWF <- matrix(nrow=630,ncol=14)
    for (i in 1:630){
        for (j in 1:14){
            TS.pert.leadtime.ensemblemean_ECMWF[i,j] <- mean(as.vector(southperttemp[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,ceiling(i/90)])-273.15)
        }
    }
    
    southperttemp <- load.data("./Data/NCEPpert_temp2m_7yr_south24hrave.rda")
    #TS
    TS.pert.leadtime.ensemblemean_NCEP <- matrix(nrow=630,ncol=14)
    for (i in 1:630){
        for (j in 1:14){
            TS.pert.leadtime.ensemblemean_NCEP[i,j] <- mean(as.vector(southperttemp[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,ceiling(i/90)])-273.15)
        }
    }
    
    southperttemp <- load.data("./Data/UKMOpert_temp2m_7yr_south24hrave.rda")
    #TS
    TS.pert.leadtime.ensemblemean_UKMO <- matrix(nrow=630,ncol=14)
    for (i in 1:630){
        for (j in 1:14){
            TS.pert.leadtime.ensemblemean_UKMO[i,j] <- mean(as.vector(southperttemp[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,ceiling(i/90)])-273.15)
        }
    }
    
    mean_TS <- array(0,dim=c(3,630,14))
    for (j in 1:14){
        for (i in 1:630){
            mean_TS[,i,j] <- cbind(TS.pert.leadtime.ensemblemean_ECMWF[i,j],TS.pert.leadtime.ensemblemean_NCEP[i,j],TS.pert.leadtime.ensemblemean_UKMO[i,j])
        }
    }
}

# PC1 & PC2
{
    
    # ECMWF
    {
        load("./Data/ECMWFpert_timeseries.rda")
        
        PC1.pert.leadtime.ensemblemean_ECMWF <- matrix(nrow=630,ncol=14)
        PC2.pert.leadtime.ensemblemean_ECMWF <- matrix(nrow=630,ncol=14)
        for (i in 1:630){
            for (j in 1:14){
                PC1.pert.leadtime.ensemblemean_ECMWF[i,j] <- mean(as.vector(ensemble_ts[djf.index2(model_time1)[i]-j,1,j*4+1,]))
                PC2.pert.leadtime.ensemblemean_ECMWF[i,j] <- mean(as.vector(ensemble_ts[djf.index2(model_time1)[i]-j,2,j*4+1,]))
            }
        }
    }

    # NCEP
    {
        load("./Data/NCEPpert_mslp_timeseries.rda")
        
        PC1.pert.leadtime.ensemblemean_NCEP <- matrix(nrow=630,ncol=14)
        PC2.pert.leadtime.ensemblemean_NCEP <- matrix(nrow=630,ncol=14)
        for (i in 1:630){
            for (j in 1:14){
                PC1.pert.leadtime.ensemblemean_NCEP[i,j] <- mean(as.vector(ensemble_ts[djf.index2(model_time1)[i]-j,1,j*4+1,]))
                PC2.pert.leadtime.ensemblemean_NCEP[i,j] <- mean(as.vector(ensemble_ts[djf.index2(model_time1)[i]-j,2,j*4+1,]))
            }
        }
    }
    
    # UKMO
    {
        load("./Data/UKMOpert_mslp_timeseries.rda")
        
        PC1.pert.leadtime.ensemblemean_UKMO <- matrix(nrow=630,ncol=14)
        PC2.pert.leadtime.ensemblemean_UKMO <- matrix(nrow=630,ncol=14)
        for (i in 1:630){
            for (j in 1:14){
        PC1.pert.leadtime.ensemblemean_UKMO[i,j] <- mean(as.vector(ensemble_ts[djf.index2(model_time1)[i]-j,1,j*4+1,]))
                PC2.pert.leadtime.ensemblemean_UKMO[i,j] <- mean(as.vector(ensemble_ts[djf.index2(model_time1)[i]-j,2,j*4+1,]))
            }
        }
    }
    
    
    mean_PC1 <- array(0,dim=c(3,630,14))
    mean_PC2 <- array(0,dim=c(3,630,14))
    for (j in 1:14){
        for (i in 1:630){
            mean_PC1[,i,j] <- cbind(PC1.pert.leadtime.ensemblemean_ECMWF[i,j],PC1.pert.leadtime.ensemblemean_NCEP[i,j],PC1.pert.leadtime.ensemblemean_UKMO[i,j])
            mean_PC2[,i,j] <- cbind(PC2.pert.leadtime.ensemblemean_ECMWF[i,j],PC2.pert.leadtime.ensemblemean_NCEP[i,j],PC2.pert.leadtime.ensemblemean_UKMO[i,j])
        }
    }
}

sigma <- array(0,dim=c(4,4,630,14))
for (j in 1:14){
    for (i in 1:630)
        sigma[,,i,j] <- cov(cbind(mean_TN[,i,j],mean_TS[,i,j],mean_PC1[,i,j],mean_PC2[,i,j]))
}

####################################################################################################

# LAMBDA & ETA                                                                                  ####

# ensemble means for each variable & leadtime
{
    ensemble_ALL <- abind::abind(ensemble_ECMWF,ensemble_NCEP,ensemble_UKMO)
    TN_ALL <- abind::abind(TN_ECMWF,TN_NCEP,TN_UKMO,along=3)
    TS_ALL <- abind::abind(TS_ECMWF,TS_NCEP,TS_UKMO,along=3)
    
    PC1.pert.leadtime.ensemblemean_ALL<- matrix(nrow=630,ncol=14)
    PC2.pert.leadtime.ensemblemean_ALL<- matrix(nrow=630,ncol=14)
    TN.pert.leadtime.ensemblemean_ALL <- matrix(nrow=630,ncol=14)
    TS.pert.leadtime.ensemblemean_ALL <- matrix(nrow=630,ncol=14)
    for (i in 1:630){
        for (j in 1:14){
            PC1.pert.leadtime.ensemblemean_ALL[i,j] <- mean(as.vector(ensemble_ALL[djf.index2(model_time1)[i]-j,1,j*4+1,]))
            PC2.pert.leadtime.ensemblemean_ALL[i,j] <- mean(as.vector(ensemble_ALL[djf.index2(model_time1)[i]-j,2,j*4+1,]))
            TN.pert.leadtime.ensemblemean_ALL[i,j] <- mean(as.vector(TN_ALL[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,ceiling(i/90)])-273.15)
            TS.pert.leadtime.ensemblemean_ALL[i,j] <- mean(as.vector(TS_ALL[j*4+1,djf.index1(model_time1)[(i-1)%%90+1]-j,,ceiling(i/90)])-273.15)
        }
    }
}


# mean errors for each variable
{
    mean.error_PC1_ALL <- array(0,dim=c(7,90,14))
    for (i in 1:7){
        for (z in 1:90){
            mean.error_PC1_ALL[i,z,] <- PC1.pert.leadtime.ensemblemean_ALL[z+90*(i-1),] - PC1_obs[z+90*(i-1)]
        }
    }
    
    mean.error_PC2_ALL <- array(0,dim=c(7,90,14))
    for (i in 1:7){
        for (z in 1:90){
            mean.error_PC2_ALL[i,z,] <- PC2.pert.leadtime.ensemblemean_ALL[z+90*(i-1),] - PC2_obs[z+90*(i-1)]
        }
    }
    
    mean.error_TN_ALL <- array(0,dim=c(7,90,14))
    for (i in 1:7){
        for (z in 1:90){
            mean.error_TN_ALL[i,z,] <- TN.pert.leadtime.ensemblemean_ALL[z+90*(i-1),] - TN_obs[z+90*(i-1)]
        }
    }
    
    mean.error_TS_ALL <- array(0,dim=c(7,90,14))
    for (i in 1:7){
        for (z in 1:90){
            mean.error_TS_ALL[i,z,] <- TS.pert.leadtime.ensemblemean_ALL[z+90*(i-1),] - TS_obs[z+90*(i-1)]
        }
    }
}

lambda <- array(0,dim = c(4,4,90,14))
for (j in 1:14){
    for (z in 1:90)
        lambda[,,z,j] <- cov(cbind(mean.error_TN_ALL[,z,j],mean.error_TS_ALL[,z,j],mean.error_PC1_ALL[,z,j],mean.error_PC2_ALL[,z,j]))
}


eta <- array(0,dim=c(4,90,14))
for (j in 1:14){
    for (i in 1:90){
        eta[,i,j] <- c(mean(mean.error_TN_ALL[,i,j]),mean(mean.error_TS_ALL[,i,j]),mean(mean.error_PC1_ALL[,i,j]),mean(mean.error_PC2_ALL[,i,j]))
    }
}

####################################################################################################

# Y_BAR                                                                                         ####

Y_bar_ECMWF <- array(0,dim=c(4,630,14))
for (j in 1:14){
    for (i in 1:630){
        Y_bar_ECMWF[,i,j] <- cbind(TN.pert.leadtime.ensemblemean_ECMWF[i,j],TS.pert.leadtime.ensemblemean_ECMWF[i,j],PC1.pert.leadtime.ensemblemean_ECMWF[i,j],PC2.pert.leadtime.ensemblemean_ECMWF[i,j])
    }
}

Y_bar_NCEP <- array(0,dim=c(4,630,14))
for (j in 1:14){
    for (i in 1:630){
        Y_bar_NCEP[,i,j] <- cbind(TN.pert.leadtime.ensemblemean_NCEP[i,j],TS.pert.leadtime.ensemblemean_NCEP[i,j],PC1.pert.leadtime.ensemblemean_NCEP[i,j],PC2.pert.leadtime.ensemblemean_NCEP[i,j])
    }
}

Y_bar_UKMO <- array(0,dim=c(4,630,14))
for (j in 1:14){
    for (i in 1:630){
        Y_bar_UKMO[,i,j] <- cbind(TN.pert.leadtime.ensemblemean_UKMO[i,j],TS.pert.leadtime.ensemblemean_UKMO[i,j],PC1.pert.leadtime.ensemblemean_UKMO[i,j],PC2.pert.leadtime.ensemblemean_UKMO[i,j])
    }
}


####################################################################################################

# PRECISION ETC                                                                                 ####

#For 1 Dec 2007, one lead day
D_ECMWF <- array(0,dim = c(4,4,630,14))
D_NCEP <- array(0,dim = c(4,4,630,14))
D_UKMO <- array(0,dim = c(4,4,630,14))
precision <- array(0,dim = c(4,4,630,14))
S<- array(0,dim = c(4,4,630,14))
tau <- array(0,dim = c(4,630,14))
for (j in 1:14) {
    for (i in 1:630) {
        D_ECMWF[,,i,j] <- C_ECMWF[,,i,j]/length(ensemble_ECMWF[1,1,1,])+sigma[,,i,j]
        D_NCEP[,,i,j]  <- C_NCEP[,,i,j]/length(ensemble_NCEP[1,1,1,])+sigma[,,i,j]
        D_UKMO[,,i,j]  <- C_UKMO[,,i,j]/length(ensemble_UKMO[1,1,1,])+sigma[,,i,j]
        precision[,,i,j] <- solve(lambda[,,(i-1)%%90+1,j]+
                                      solve(solve(D_ECMWF[,,i,j])+solve(D_NCEP[,,i,j])+
                                                solve(D_UKMO[,,i,j])))
        S[,,i,j] <- solve(precision[,,i,j])
        tau[,i,j] <- S[,,i,j]%*%(solve(diag(4)+solve(D_ECMWF[,,i,j])%*%lambda[,,(i-1)%%90+1,j]+solve(D_NCEP[,,i,j])%*%lambda[,,(i-1)%%90+1,j]+solve(D_UKMO[,,i,j])%*%lambda[,,(i-1)%%90+1,j])%*%(solve(D_ECMWF[,,i,j])%*%Y_bar_ECMWF[,i,j]+solve(D_NCEP[,,i,j])%*%Y_bar_NCEP[,i,j]+solve(D_UKMO[,,i,j])%*%Y_bar_UKMO[,i,j]))-eta[,(i-1)%%90+1,j]
    }
}

S_diag <- array(0,dim = c(4,630,14))
Spread <- array(0,dim=c(4,14))
for (j in 1:14) {
    for (i in 1:630) {
        for (k in 1:4){
            S_diag[,i,j] <- diag(S[,,i,j])
            Spread[k,j] <- sqrt(mean(S_diag[k,,j]))
        }
    }
}

####################################################################################################

# RMSE, SPREAD, CRPS                                                                            ####

observation <- cbind(TN_obs,TS_obs,PC1_obs,PC2_obs)
RMSE <- array(0,dim=c(4,14))
for (j in 1:14){
    for(k in 1:4) {
        RMSE[k,j] <- sqrt(mean((tau[k,,j]- observation[,k])^2))
    }
}

# CRPS
# mean(verification::crps(observation[,1],c(mean(tau[1,,1]),Spread[1,1]))$crps)
# verification::crps(observation[1,],cbind(tau[,1,1],sqrt(diag(S[,,1,1]))))$crps
# colMeans(verification::crps(observation[,],cbind(tau[,,1],Spread[,1]))$crps)
CRPS_Richard <- array(0,dim = c(4,630,14))
CRPS_Richard_average <- array(0,dim=c(4,14))

for (i in 1:630){
    for (j in 1:14){
        for (k in 1:4) {
            CRPS_Richard[,i,j] <- verification::crps(observation[i,],cbind(tau[,i,j],sqrt(diag(S[,,i,j]))))$crps
            CRPS_Richard_average[k,j] <- mean(CRPS_Richard[k,,j])
        }
    }
}
