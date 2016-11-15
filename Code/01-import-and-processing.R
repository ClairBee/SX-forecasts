
library("CB.Misc"); library("SX.weather")

fc <- readRDS("./Data/Forecasts.rds")

####################################################################################################

# LOAD ALL DATA INTO SINGLE ARRAY FOR EASIER MANIPULATION                                       ####

# create single array: var . day . year . leadtime . version (0 = control)
fc <- array(dim = c(5, 105, 7, 15, 51),
            dimnames = list(c("temp.n", "temp.s", "pc1", "pc2", "pc3"),
                            c(1:105), formatC(8:14, width = 2, flag = "0"), c(0:14),
                            c("c",1:50)))

# mslp control forecast
{
    fc.mslp.ctrl <- load.data("./Data/ECMWFctrl_mslp_timeseries.rda")
    invisible(sapply(0:14, function(i) {
        fc["pc1",,,toString(i),"c"] <<- fc.mslp.ctrl[1:735,1, 4 * i + 1]
        fc["pc2",,,toString(i),"c"] <<- fc.mslp.ctrl[1:735,2, 4 * i + 1]
        fc["pc3",,,toString(i),"c"] <<- fc.mslp.ctrl[1:735,3, 4 * i + 1]
    }))
}

# mslp perturbed forecast
{
    fc.mslp.pert <- load.data("./Data/ECMWFpert_timeseries.rda")
    invisible(sapply(0:14, function(i) {
        fc["pc1",,,toString(i), 2:51] <<- fc.mslp.pert[1:735,1, 4 * i + 1,]
        fc["pc2",,,toString(i), 2:51] <<- fc.mslp.pert[1:735,2, 4 * i + 1,]
        fc["pc3",,,toString(i), 2:51] <<- fc.mslp.pert[1:735,3, 4 * i + 1,]
    }))
}

# temp control forecast
{
    fc.temp.s.ctrl <- load.data("./Data/ECMWFctrl_temp2m_7yr_south24hrave.rda")
    fc.temp.n.ctrl <- load.data("./Data/ECMWFctrl_temp2m_7yr_north24hrave.rda")
    
    invisible(sapply(0:14, function(i) {
        fc["temp.n",,,toString(i), 1] <<- fc.temp.n.ctrl[4 * i + 1, ,] - 273.15
        fc["temp.s",,,toString(i), 1] <<- fc.temp.s.ctrl[4 * i + 1, ,] - 273.15
    }))
}

# temp perturbed forecast
{
    fc.temp.n.pert <- load.data("./Data/ECMWFpert_temp2m_7yr_north24hrave.rda")
    fc.temp.s.pert <- load.data("./Data/ECMWFpert_temp2m_7yr_south24hrave.rda")
    
    invisible(sapply(0:14, function(i) {
        sapply(1:50, function(j) {
            fc["temp.n",,,toString(i), j + 1] <<- fc.temp.n.pert[4 * i + 1,, j,] - 273.15
            fc["temp.s",,,toString(i), j + 1] <<- fc.temp.s.pert[4 * i + 1,, j,] - 273.15
        })
    }))
}

saveRDS(fc, "./Data/Forecasts.rds")

####################################################################################################

