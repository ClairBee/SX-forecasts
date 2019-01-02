
# animated plotting function to represent sequences of ensemble forecasts

####################################################################################################

# BASE PLOTTING FUNCTION                                                                        ####

# 'death and the spider' plot

ens.cols <- c("blue", "red", "green3")

tn.rng <- range(obs["temp.n",,], ecmwf["temp.n",,,,], ncep["temp.n",,,,], ukmo["temp.n",,,,])
ts.rng <- range(obs["temp.s",,], ecmwf["temp.s",,,,], ncep["temp.s",,,,], ukmo["temp.s",,,,])

ds.plot <- function(dd, yy, lt, alph = 0.2, main = paste0("Death & the spider - year ", yy, ", day ", dd, ", LT ", toString(lt))) {
    
    lt <- toString(lt)
    plot(0, type = "n", xlim = tn.rng, ylim = ts.rng, xlab = "temp.n", ylab = "temp.s", asp = T, main = main)
    
    do <- dd + 15
    
    ens.mean <- t(apply(abind(apply(ecmwf[1:2, do, yy, lt, -1], 1, mean),
                      apply(ncep[1:2, do, yy, lt, -1], 1, mean),
                      apply(ukmo[1:2, do, yy, lt, -1], 1, mean), along = 0),
                      2, mean))
    
    points(t(ecmwf[1:2, do, yy, lt, -1]), pch = 20, col = adjustcolor(ens.cols[1], alpha = alph))
    points(t(ecmwf[1:2, do, yy, lt, 1]), col = ens.cols[1])
    lines(rbind(ens.mean, t(apply(ecmwf[1:2, do, yy, lt, -1], 1, mean))),
          col = ens.cols[1], lwd = 2)
    points(t(apply(ecmwf[1:2, do, yy, lt, -1], 1, mean)), pch = 21, bg = ens.cols[1])
    
    points(t(ncep[1:2, do, yy, lt, -1]), pch = 20, col = adjustcolor(ens.cols[2], alpha = alph))
    points(t(ncep[1:2, do, yy, lt, 1]), col = ens.cols[2])
    lines(rbind(ens.mean, t(apply(ncep[1:2, do, yy, lt, -1], 1, mean))),
          col = ens.cols[2], lwd = 2)
    points(t(apply(ncep[1:2, do, yy, lt, -1], 1, mean)), pch = 21, bg = ens.cols[2])
    
    points(t(ukmo[1:2, do, yy, lt, -1]), pch = 20, col = adjustcolor(ens.cols[3], alpha = alph))
    points(t(ukmo[1:2, do, yy, lt, 1]), col = ens.cols[3])
    lines(rbind(ens.mean, t(apply(ukmo[1:2, do, yy, lt, -1], 1, mean))),
          col = ens.cols[3], lwd = 2)
    points(t(apply(ukmo[1:2, do, yy, lt, -1], 1, mean)), pch = 21, bg = ens.cols[3])
    
    lines(rbind(ens.mean, t(obs[1:2, dd, yy])), lwd = 2)
    points(ens.mean, pch = 19)
    points(t(obs[1:2, dd, yy]), pch = 23, bg = "gold", cex = 1)
}

ds.plot(19, 4, 14)

which(obs[1,,] == min(obs[1,,]), arr.ind = T)

####################################################################################################

# ANIMATE FOR SINGLE DATE, DECREASING LEADTIME                                                  ####

library(animation)

DD <- 19; YY <- 4

saveGIF({
    for(i in 14:0){
        ds.plot(DD, YY, i)
    }
}, interval = 1, ani.width = 1500, ani.height = 1500)

####################################################################################################

# ANIMATE ACROSS SEQUENTIAL DATES                                                               ####
YY <- 1

saveGIF({
    par(mfrow = c(2,2), mar = c(3,3,2,1), oma = c(0,0,2,0), cex = 2, lwd = 2)
    for(d in 1:90){
        ds.plot(d, YY, lt = 14, main = "LT = 14")
        ds.plot(d, YY, lt = 10, main = "LT = 10")
        ds.plot(d, YY, lt = 5, main = "LT = 5")
        ds.plot(d, YY, lt = 0, main = "LT = 0")
    
        mtext(paste0("Death & the spider - year ", YY, ", day ", d), cex = 2, outer = T) 
    }
}, interval = 1, ani.width = 1500, ani.height = 1500)

####################################################################################################

# ANIMATE PATHS THROUGH STATE SPACE                                                             ####

lts <- c(0,3,7,10,14)
plot(0, type = "n", xlim = tn.rng, ylim = ts.rng, xlab = "temp.n", ylab = "temp.s", asp = T,
     main = "Ensemble paths through state space")

points(t(apply(ecmwf[1:2, dd, yy, lts + 1, -1], 1:2, mean)), type = "o", pch = 20, col = sapply(1/c(5:1), function(a) adjustcolor(ens.cols[1], alpha = a)))
points(t(apply(ncep[1:2, dd, yy, lts + 1, -1], 1:2, mean)), type = "o", pch = 20, col = sapply(1/c(5:1), function(a) adjustcolor(ens.cols[2], alpha = a)))
points(t(apply(ukmo[1:2, dd, yy, lts + 1, -1], 1:2, mean)), type = "o", pch = 20, col = sapply(1/c(5:1), function(a) adjustcolor(ens.cols[3], alpha = a)))



ens.mean <- t(apply(abind(apply(ecmwf[1:2, do, yy, lt, -1], 1, mean),
                          apply(ncep[1:2, do, yy, lt, -1], 1, mean),
                          apply(ukmo[1:2, do, yy, lt, -1], 1, mean), along = 0),
                    2, mean))

points(t(ecmwf[1:2, dd, yy, lt, -1]), pch = 20, col = adjustcolor(ens.cols[1], alpha = alph))
points(t(ecmwf[1:2, do, yy, lt, 1]), col = ens.cols[1])
lines(rbind(ens.mean, t(apply(ecmwf[1:2, do, yy, lt, -1], 1, mean))),
      col = ens.cols[1], lwd = 2)
