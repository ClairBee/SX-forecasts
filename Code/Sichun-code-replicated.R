
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


