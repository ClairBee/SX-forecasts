
library("CB.Misc"); library("SX.weather")

load.all()

####################################################################################################

# NOTES                                                                                         ####

# checked alignment of forecast & obs visually - data appear to be aligned (days 16:105)


####################################################################################################

# LOAD ALL FORECASTS INTO SINGLE ARRAY FOR EASIER MANIPULATION                                  ####

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

# LOAD ALL OBSERVATIONS INTO SINGLE ARRAY                                                       ####

obs <- array(dim = c(5, 90, 7), dimnames = list(dimnames(fc)[[1]], c(1:90), dimnames(fc)[[3]]))

# principal components of mean sea level pressure
{
    obs.mslp <- load.data("./Data/ERAint_pca_timeseries.rda")
    
    # 3240 - 720 + 1 = 2611, this is the start time for 2007/2008-2013/2014.
    obs["pc1",,] <- obs.mslp[2521:3150,1] 
    obs["pc2",,] <- obs.mslp[2521:3150,2] 
    obs["pc3",,] <- obs.mslp[2521:3150,3] 
}

# temperatures
{
    # didn't investigate indexing for this file - this was provided by SX
    # take every 4th, starting from 3rd (midday measurement?); 29:35 gives 7 years' data
    obs["temp.n",,] <- load.data("./Data/ERAint_temp2m_36yr_north24hrave.rda")[seq(3,360,4), 29:35] - 273.15
    obs["temp.s",,] <- load.data("./Data/ERAint_temp2m_36yr_south24hrave.rda")[seq(3,360,4), 29:35] - 273.15
    
}

# check that these are properly aligned with forecasts
{
    par(mfrow = c(5, 7), mar = c(0,0,0,0), lwd = 3, oma = c(0,0,3,0))
    
    invisible(sapply(1:5, function(varbl) {
        invisible(sapply(formatC(8:14, width = 2, flag = "0"), function(yr) {
            plot(obs[varbl,,yr], type = "l", xlab = "", ylab = "", yaxt = "none", xaxt = "none", lwd = 1)
            lines(fc[varbl,16:105,yr,"0","c"], col = "red", lwd = 1)
        }))
    }))
    mtext("Check alignment of control forecast (red) & observations (black)", outer = TRUE)
}
# control forecast from days 16:105 seems to be matched to observations. Assume dates properly aligned.

saveRDS(obs, "./Data/Observations.rds")

####################################################################################################

# CALCULATE MEAN ERROR, RMSE                                                                    ####

####################################################################################################

ens.mean <- apply(fc[,,,,2:51], c(1,2,3,4), mean)

ens.mean.error <- apply(sweep(ens.mean, 1:3, obs, "-"), c(1, 4), mean)
ctrl.mean.error <- apply(sweep(fc[,,,,1], 1:3, obs, "-"), c(1, 4), mean)

suppressWarnings(ens.member.errors <- apply(sweep(fc[,,,,2:51], c(1:3, 5), obs, "-"), c(1, 4, 5), mean))

# plot forecast errors for all variables - compare to figure 5.1 in dissertation
    # plots are similar, but not identical (particularly after 10 days?)
    # need to get date reference file in order to investigate further.
{
    # quick plotting function
    qplot <- function(element, ...) {
        
        ttl <- switch(element,
                      "temp.s" = "Temp (south)", "temp.n" = "Temp (north)",
                      "pc1" = "First PC", "pc2" = "Second PC", "pc3" = "Third PC")
        
        matplot(ens.member.errors[element,,], type = "l", col = adjustcolor("grey", alpha = 0.4),
                lty = 1, main = ttl, xlab = "", ylab = "", ...)
        
        lines(ens.mean.error[element,])
        lines(ctrl.mean.error[element,], col = "blue3")
    }
    
    par(mfrow = c(2,3), oma = c(0.5, 0.5, 2, 0.5), mar = c(2,2,3,1))
    qplot("pc1"); qplot("pc2"); qplot("pc3")
    qplot("temp.n"); qplot("temp.s")
    
    plot.new()
    legend("left", lty = 1, col = c("blue", "black", adjustcolor("grey", alpha = 0.5)), bty = "n", cex = 1.1,
           legend = c("Control forecast", "Perturbed mean", "Perturbed members"))
    
    mtext("Mean error at each forecast lead time", outer = TRUE, cex = 1)
}

ctrl.rmse <- apply(sweep(fc[,,,,1], 1:3, obs, "-"), c(1, 4), function(err) sqrt(mean(err^2)))
ens.rmse <- apply(sweep(ens.mean, 1:3, obs, "-"), c(1, 4), function(err) sqrt(mean(err^2)))
suppressWarnings(ens.member.rmse <- apply(sweep(fc[,,,,2:51], c(1:3, 5), obs, "-"), c(1, 4, 5), function(err) sqrt(mean(err^2))))

# plot forecast rmse for all variables - compare to figure 5.2 in dissertation
    # 
{
    # quick plotting function
    qplot <- function(element, ...) {
        
        ttl <- switch(element,
                      "temp.s" = "Temp (south)", "temp.n" = "Temp (north)",
                      "pc1" = "First PC", "pc2" = "Second PC", "pc3" = "Third PC")
        
        matplot(ens.member.rmse[element,,], type = "l", col = adjustcolor("grey", alpha = 0.4),
                lty = 1, main = ttl, xlab = "", ylab = "", ...)
        
        lines(ens.rmse[element,])
        lines(ctrl.rmse[element,], col = "blue3")
    }
    
    par(mfrow = c(2,3), oma = c(0.5, 0.5, 2, 0.5), mar = c(2,2,3,1))
    qplot("pc1"); qplot("pc2"); qplot("pc3")
    qplot("temp.n"); qplot("temp.s")
    
    plot.new()
    legend("left", lty = 1, col = c("blue", "black", adjustcolor("grey", alpha = 0.5)), bty = "n", cex = 1.1,
           legend = c("Control forecast", "Perturbed mean", "Perturbed members"))
    
    mtext("RMSE at each forecast lead time", outer = TRUE, cex = 1)
}

####################################################################################################

# CASE STUDY: CONTROL FORECAST FOR TEMP (NORTH), 2008                                           ####

plot(obs["temp.n",,"08"], type = "l", lwd = 2, ylab = "", xlab = "", main = "Obs vs control forecast")
lines(fc["temp.n", , "08", "0", "c"], col = "blue")

error <- fc["temp.n", , , "0", "c"] - obs["temp.n",,]

matplot(error, type = "l", lty = 1, main = "Same-day control forecast error")
apply(error, 2, mean)
mean(error)
ctrl.mean.error["temp.n", "0"]

sqrt(mean(error^2))
ctrl.rmse["temp.n", "0"]

# same result. Don't think this is a problem with my calculations.