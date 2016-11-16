
library("CB.Misc"); library("SX.weather")

kc <- 273.15

load.all()

# missing: UKMOctrl_mslp_timeseries.rda

####################################################################################################

# NOTES                                                                                         ####

# checked alignment of forecast & obs visually - data appear to be aligned (days 16:105)
# why do perturbed & control forecasts of principal components start in different places?


####################################################################################################

# LOAD ALL ECMWF FORECASTS INTO SINGLE ARRAY FOR EASIER MANIPULATION                            ####

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
        fc["temp.n",,,toString(i), 1] <<- fc.temp.n.ctrl[4 * i + 1, ,] - kc
        fc["temp.s",,,toString(i), 1] <<- fc.temp.s.ctrl[4 * i + 1, ,] - kc
    }))
}

# temp perturbed forecast
{
    fc.temp.n.pert <- load.data("./Data/ECMWFpert_temp2m_7yr_north24hrave.rda")
    fc.temp.s.pert <- load.data("./Data/ECMWFpert_temp2m_7yr_south24hrave.rda")
    
    invisible(sapply(0:14, function(i) {
        sapply(1:50, function(j) {
            fc["temp.n",,,toString(i), j + 1] <<- fc.temp.n.pert[4 * i + 1,, j,] - kc
            fc["temp.s",,,toString(i), j + 1] <<- fc.temp.s.pert[4 * i + 1,, j,] - kc
        })
    }))
}

saveRDS(fc[,,,,], "./Data/ECMWF-forecasts.rds")

# perturbed ensemble has an offset from the control forecast. Why?
{
    pc.qplot <- function(n, lt = 0) {
        diff <- sweep(fc.mslp.pert[,n, 4 * lt + 1,], 1, fc.mslp.ctrl[1:735,n,4 * lt + 1], "-")
        matplot(diff, pch = 20, cex = 0.4, col = adjustcolor("darkblue", alpha = 0.4),
                main = paste0("PC", n), ylim = c(-0.05, 0.15))
        abline(h = 0, lwd = 2)
        invisible(apply(diff, 2, function(mm) abline(line(mm), col = adjustcolor("gold", alpha = 0.6))))
    }
    
    temp.qplot <- function(a, lt = 0) {
        diff <- switch(a,
                       "n" = array(aperm(sweep(fc.temp.n.pert[4 * lt + 1,,,], c(1,3), fc.temp.n.ctrl[4 * lt + 1,,], "-"), c(1,3,2)), dim = c(735,50)),
                       "s" = array(aperm(sweep(fc.temp.s.pert[4 * lt + 1,,,], c(1,3), fc.temp.s.ctrl[4 * lt + 1,,], "-"), c(1,3,2)), dim = c(735,50)))
        
        matplot(diff, pch = 20, cex = 0.4, col = adjustcolor("darkblue", alpha = 0.4),
                main = paste0("Temp - ", toupper(a)), ylim = c(-1.5, 1.5))
        invisible(apply(diff, 2, function(mm) abline(line(mm), col = adjustcolor("gold", alpha = 0.6))))
        abline(h = 0, lwd = 2)
    }
    
    pdf("./Plots/Perturbation-offset.pdf"); {
        par(mfrow = c(2,3), oma = c(0,0,2,0), mar = c(3, 3, 3, 1))
        pc.qplot(1); pc.qplot(2); pc.qplot(3)
        temp.qplot("n"); temp.qplot("s")
        
        plot.new()
        legend("left", lty = 1, col = "gold", bty = "n", legend = "Fitted line per ensemble member")
        
        mtext("Difference between perturbed ensemble members and control forecast", outer = T)
    }; dev.off()

}


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
    obs["temp.n",,] <- load.data("./Data/ERAint_temp2m_36yr_north24hrave.rda")[seq(3,360,4), 29:35] - kc
    obs["temp.s",,] <- load.data("./Data/ERAint_temp2m_36yr_south24hrave.rda")[seq(3,360,4), 29:35] - kc
}

# check that these are properly aligned with forecasts
{
    pdf("./Plots/ECMWF-forecast-obs-alignment.pdf"); {
        par(mfrow = c(5, 7), mar = c(0,0,0,0), lwd = 3, oma = c(0,0,3,0))
        
        invisible(sapply(1:5, function(varbl) {
            invisible(sapply(formatC(8:14, width = 2, flag = "0"), function(yr) {
                plot(obs[varbl,,yr], type = "l", xlab = "", ylab = "", yaxt = "none", xaxt = "none", lwd = 1)
                lines(fc[varbl,16:105,yr,"0","c"], col = "red", lwd = 1)
            }))
        }))
        mtext("Check alignment of ECMWF control forecast (red) & observations (black)", outer = TRUE)
    }; dev.off()

}
# control forecast from days 16:105 seems to be matched to observations. Assume dates properly aligned.

saveRDS(obs, "./Data/Observations.rds")

####################################################################################################

# LOAD NCEP FORECASTS INTO SINGLE ARRAY                                                         ####

# create single array: var . day . year . leadtime . version (0 = control)
ncep <- array(dim = c(5, 105, 7, 15, 21),
            dimnames = list(c("temp.n", "temp.s", "pc1", "pc2", "pc3"),
                            c(1:105), formatC(8:14, width = 2, flag = "0"), c(0:14),
                            c("c",1:20)))

# mslp control forecast
{
    fc.mslp.ctrl <- load.data("./Data/NCEPctrl_mslp_timeseries.rda")
    invisible(sapply(0:14, function(i) {
        ncep["pc1",,,toString(i),"c"] <<- fc.mslp.ctrl[1:735,1, 4 * i + 1]
        ncep["pc2",,,toString(i),"c"] <<- fc.mslp.ctrl[1:735,2, 4 * i + 1]
        ncep["pc3",,,toString(i),"c"] <<- fc.mslp.ctrl[1:735,3, 4 * i + 1]
    }))
}

# mslp perturbed forecast
{
    fc.mslp.pert <- load.data("./Data/NCEPpert_mslp_timeseries.rda")
    invisible(sapply(0:14, function(i) {
        ncep["pc1",,,toString(i), 2:21] <<- fc.mslp.pert[1:735,1, 4 * i + 1,]
        ncep["pc2",,,toString(i), 2:21] <<- fc.mslp.pert[1:735,2, 4 * i + 1,]
        ncep["pc3",,,toString(i), 2:21] <<- fc.mslp.pert[1:735,3, 4 * i + 1,]
    }))
}

# temp control forecast
{
    fc.temp.s.ctrl <- load.data("./Data/NCEPctrl_temp2m_7yr_south24hrave.rda")
    fc.temp.n.ctrl <- load.data("./Data/NCEPctrl_temp2m_7yr_north24hrave.rda")
    
    invisible(sapply(0:14, function(i) {
        ncep["temp.n",,,toString(i), 1] <<- fc.temp.n.ctrl[4 * i + 1, ,] - kc
        ncep["temp.s",,,toString(i), 1] <<- fc.temp.s.ctrl[4 * i + 1, ,] - kc
    }))
}

# temp perturbed forecast
{
    fc.temp.n.pert <- load.data("./Data/NCEPpert_temp2m_7yr_north24hrave.rda")
    fc.temp.s.pert <- load.data("./Data/NCEPpert_temp2m_7yr_south24hrave.rda")
    
    invisible(sapply(0:14, function(i) {
        sapply(1:dim(fc.temp.s.pert)[[3]], function(j) {
            ncep["temp.n",,,toString(i), j + 1] <<- fc.temp.n.pert[4 * i + 1,, j,] - kc
            ncep["temp.s",,,toString(i), j + 1] <<- fc.temp.s.pert[4 * i + 1,, j,] - kc
        })
    }))
}

# check that forecasts & observations are properly aligned
{
    obs <- readRDS("./Data/Observations.rds")
    
    pdf("./Plots/NCEP-forecast-obs-alignment.pdf"); {
        par(mfrow = c(5, 7), mar = c(0,0,0,0), lwd = 3, oma = c(0,0,3,0))
        
        invisible(sapply(1:5, function(varbl) {
            invisible(sapply(formatC(8:14, width = 2, flag = "0"), function(yr) {
                plot(obs[varbl,,yr], type = "l", xlab = "", ylab = "", yaxt = "none", xaxt = "none", lwd = 1)
                lines(ncep[varbl,16:105,yr,"0","c"], col = "red", lwd = 1)
            }))
        }))
        mtext("Check alignment of NCEP control forecast (red) & observations (black)", outer = TRUE)
    }; dev.off()
}

####################################################################################################

# LOAD UKMO FORECASTS INTO SINGLE ARRAY                                                         ####

# create single array: var . day . year . leadtime . version (0 = control)
ukmo <- array(dim = c(5, 105, 7, 15, 24),
              dimnames = list(c("temp.n", "temp.s", "pc1", "pc2", "pc3"),
                              c(1:105), formatC(8:14, width = 2, flag = "0"), c(0:14),
                              c("c",1:23)))

# mslp control forecast - DOES NOT EXIST
{
    fc.mslp.ctrl <- load.data("./Data/UKMOctrl_mslp_timeseries.rda")
    invisible(sapply(0:14, function(i) {
        ukmo["pc1",,,toString(i),"c"] <<- fc.mslp.ctrl[1:735,1, 4 * i + 1]
        ukmo["pc2",,,toString(i),"c"] <<- fc.mslp.ctrl[1:735,2, 4 * i + 1]
        ukmo["pc3",,,toString(i),"c"] <<- fc.mslp.ctrl[1:735,3, 4 * i + 1]
    }))
}

# mslp perturbed forecast
{
    fc.mslp.pert <- load.data("./Data/UKMOpert_mslp_timeseries.rda")
    d <- dim(fc.mslp.pert)[[4]] + 1
    invisible(sapply(0:14, function(i) {
        ukmo["pc1",,,toString(i), 2:d] <<- fc.mslp.pert[1:735,1, 4 * i + 1,]
        ukmo["pc2",,,toString(i), 2:d] <<- fc.mslp.pert[1:735,2, 4 * i + 1,]
        ukmo["pc3",,,toString(i), 2:d] <<- fc.mslp.pert[1:735,3, 4 * i + 1,]
    }))
}

# temp control forecast
{
    fc.temp.s.ctrl <- load.data("./Data/UKMOctrl_temp2m_7yr_south24hrave.rda")
    fc.temp.n.ctrl <- load.data("./Data/UKMOctrl_temp2m_7yr_north24hrave.rda")
    
    invisible(sapply(0:14, function(i) {
        ukmo["temp.n",,,toString(i), 1] <<- fc.temp.n.ctrl[4 * i + 1, ,] - kc
        ukmo["temp.s",,,toString(i), 1] <<- fc.temp.s.ctrl[4 * i + 1, ,] - kc
    }))
}

# temp perturbed forecast
{
    fc.temp.n.pert <- load.data("./Data/UKMOpert_temp2m_7yr_north24hrave.rda")
    fc.temp.s.pert <- load.data("./Data/UKMOpert_temp2m_7yr_south24hrave.rda")
    
    invisible(sapply(0:14, function(i) {
        sapply(1:dim(fc.temp.s.pert)[[3]], function(j) {
            ukmo["temp.n",,,toString(i), j + 1] <<- fc.temp.n.pert[4 * i + 1,, j,] - kc
            ukmo["temp.s",,,toString(i), j + 1] <<- fc.temp.s.pert[4 * i + 1,, j,] - kc
        })
    }))
}

# check that forecasts & observations are properly aligned
{
    obs <- readRDS("./Data/Observations.rds")
    
    pdf("./Plots/UKMO-forecast-obs-alignment.pdf"); {
        par(mfrow = c(5, 7), mar = c(0,0,0,0), lwd = 3, oma = c(0,0,3,0))
        
        invisible(sapply(1:5, function(varbl) {
            invisible(sapply(formatC(8:14, width = 2, flag = "0"), function(yr) {
                plot(obs[varbl,,yr], type = "l", xlab = "", ylab = "", yaxt = "none", xaxt = "none", lwd = 1)
                lines(ukmo[varbl,16:105,yr,"0","c"], col = "red", lwd = 1)
            }))
        }))
        mtext("Check alignment of UKMO control forecast (red) & observations (black)", outer = TRUE)
    }; dev.off()
}

saveRDS(ukmo, "./Data/UKMO-forecasts.rds")
####################################################################################################

# CALCULATE MEAN ERROR, RMSE                                                                    ####

load.all()
fc <- fc[,,,2:15,]

ens.mean <- apply(fc[,,,,-1], c(1,2,3,4), mean)
ens.mean.error <- sweep(ens.mean, 1:3, obs, "-")
suppressWarnings(ens.member.errors <- apply(sweep(fc[,,,,-1], c(1:3, 5), obs, "-"), c(1, 4, 5), mean))

ens.rmse <- sqrt(apply(ens.mean.error^2, c(1, 4), mean))
ctrl.mean.error <- apply(sweep(fc[,,,,1], 1:3, obs, "-"), c(1, 4), mean)

daily.var <- array(apply(array(fc[,,,,-1], dim = c(5, 630, 14, 50)), 1:3, var),
                   dim = c(5, 90, 7, 14), 
                   dimnames = list(dimnames(fc)[[1]], 1:90, dimnames(fc)[[3]], 1:14))

average.var <- apply(daily.var, c(1, 4), mean)
root.mean.var <- sqrt(average.var)


# plot forecast errors for all variables - compare to figure 5.1 in dissertation
    # plots are similar, but not identical (particularly after 10 days?)
    # need to get date reference file in order to investigate further.
{
    # quick plotting function
    qplot <- function(element, ...) {
        
        ttl <- switch(element,
                      "temp.s" = "Temp (south)", "temp.n" = "Temp (north)",
                      "pc1" = "First PC", "pc2" = "Second PC", "pc3" = "Third PC")
        
        yl <- range(ens.member.errors[element,,], ens.mean.error[element,], ctrl.mean.error[element,])
            
        matplot(ens.member.errors[element,,], type = "l", col = adjustcolor("grey", alpha = 0.4),
                lty = 1, main = ttl, xlab = "", ylab = "", ylim = yl)
        
        lines(ens.mean.error[element,])
        lines(ctrl.mean.error[element,], col = "blue3")
    }
    
    pdf("./Plots/Mean error.pdf"); {
        par(mfrow = c(2,3), oma = c(0.5, 0.5, 2, 0.5), mar = c(2,2,3,1))
        qplot("pc1"); qplot("pc2"); qplot("pc3")
        qplot("temp.n"); qplot("temp.s")
        
        plot.new()
        legend("left", lty = 1, col = c("blue", "black", adjustcolor("grey", alpha = 0.5)), bty = "n", cex = 1.1,
               legend = c("Control forecast", "Perturbed mean", "Perturbed members"))
        
        mtext("Mean error at each forecast lead time", outer = TRUE, cex = 1)
    }; dev.off()

}

ctrl.rmse <- apply(sweep(fc[,,,,1], 1:3, obs, "-"), c(1, 4), function(err) sqrt(mean(err^2)))
suppressWarnings(ens.member.rmse <- apply(sweep(fc[,,,,2:51], c(1:3, 5), obs, "-"), c(1, 4, 5), function(err) sqrt(mean(err^2))))

# plot forecast rmse for all variables - compare to figure 5.2 in dissertation
    # 
{
    # quick plotting function
    qplot <- function(element, ...) {
        
        ttl <- switch(element,
                      "temp.s" = "Temp (south)", "temp.n" = "Temp (north)",
                      "pc1" = "First PC", "pc2" = "Second PC", "pc3" = "Third PC")
        
        yl <- range(ens.member.rmse[element,,], ens.rmse[element,], ctrl.rmse[element,])
        matplot(ens.member.rmse[element,,], type = "l", col = adjustcolor("grey", alpha = 0.4),
                lty = 1, main = ttl, xlab = "", ylab = "", ylim = yl)
        
        lines(ens.rmse[element,])
        lines(ctrl.rmse[element,], col = "blue3")
    }
    
    pdf("./Plots/RMSE.pdf"); {
        par(mfrow = c(2,3), oma = c(0.5, 0.5, 2, 0.5), mar = c(2,2,3,1))
        qplot("pc1"); qplot("pc2"); qplot("pc3")
        qplot("temp.n"); qplot("temp.s")
        
        plot.new()
        legend("left", lty = 1, col = c("blue", "black", adjustcolor("grey", alpha = 0.5)), bty = "n", cex = 1.1,
               legend = c("Control forecast", "Perturbed mean", "Perturbed members"))
        
        mtext("RMSE at each forecast lead time", outer = TRUE, cex = 1)
    }; dev.off()

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