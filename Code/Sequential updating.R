
library("SX.weather")
setwd("~/Documents/PhD/Miniprojects/03-Sichun-paper/Code")


# stop aaply complaining about duplicate levels in factors
options(warn = -1)

####################################################################################################

# EXPLORATORY: PLOTS BY LEADTIME                                                                ####

mimic <- readRDS("~/Documents/PhD/Miniprojects/03-Sichun-paper/Models/an25-mimic.rds")

fc <- abind("ec" = apply(offset.forecast(ecmwf)[1:2,,,,-1], 1:4, mean),
            "nc" = apply(offset.forecast(ncep)[1:2,,,,-1], 1:4, mean),
            "mo" = apply(offset.forecast(ukmo)[1:2,,,,-1], 1:4, mean), along = 0)
e.cols <- c("steelblue", "red3", "green3")

d <- 5; y <- 5; varb <- "temp.n"

lt.plot <- function(d, y, varb) {
    y.rng <- range(fc[,varb,d,y,], obs[varb,d,y], mimic[d,y,,varb, "Tau"])

    matplot(rev(0:14), t(fc[,varb,d,y,]), type = "l", col = e.cols, lty = 1, xlab = "", ylab = "", xaxt = "n",
            main = paste0(varb, ": ", formatC(7:14, width = 2, flag = "0")[y], " d", d), ylim = y.rng)
    axis(1, at = (0:7)*2, labels = (7:0)*2)
    abline(h = obs[varb,d,y], lty = 3)
    lines(rev(0:14), mimic[d,y,,varb, "Tau"], lty = 2)
    points(14, obs[varb,d,y], pch = 4)
}

par(mfrow = c(2,1), mar = c(2,2,3,1))

# mid-December
invisible(sapply(1:7, function(yy) {
    lt.plot(15,yy,"temp.n"); lt.plot(5,yy,"temp.s")
}))


# MIMIC USING PREVIOUS LEADTIME AS PRIOR                                                        ####
all.dat <- abind("ec" = offset.forecast(ecmwf)[,,,,-1],
                 "nc" = offset.forecast(ncep)[,,,,-1],
                 "mo" = offset.forecast(ukmo)[,,,,-1],
                 "tr" = readRDS("../Models/tr-set-default.rds"),
                 along = 5)

# initial forecast at LT 14 - no prior
mim.14 <- aaply(all.dat[,,,"14",], 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118])
})

tr.dat <- readRDS("../Models/tr-set-default.rds")
d <- 5; y <- 5
mimi.14 <- fit.mimic(fc = list("ecmwf" = ecmwf[,5, 5, "14", -1],
                               "ncep" = ncep[,5, 5, "14", -1],
                               "ukmo" = ukmo[,5, 5, "14", -1]),
                     tr = tr.dat[, d, y, 15, ])
mimi.14


mim.13 <- aaply(abind(all.dat[,,,"13",], aperm(mim.14, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

mim.12 <- aaply(abind(all.dat[,,,"12",], aperm(mim.13, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

mim.11 <- aaply(abind(all.dat[,,,"11",], aperm(mim.12, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

mim.10 <- aaply(abind(all.dat[,,,"10",], aperm(mim.11, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

mim.9 <- aaply(abind(all.dat[,,,"9",], aperm(mim.10, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

mim.8 <- aaply(abind(all.dat[,,,"8",], aperm(mim.9, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

mim.7 <- aaply(abind(all.dat[,,,"7",], aperm(mim.8, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

mim.6 <- aaply(abind(all.dat[,,,"6",], aperm(mim.7, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

mim.5 <- aaply(abind(all.dat[,,,"5",], aperm(mim.6, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

mim.4 <- aaply(abind(all.dat[,,,"4",], aperm(mim.5, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

mim.3 <- aaply(abind(all.dat[,,,"3",], aperm(mim.4, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

mim.2 <- aaply(abind(all.dat[,,,"2",], aperm(mim.3, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

mim.1 <- aaply(abind(all.dat[,,,"1",], aperm(mim.2, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})

mim.0 <- aaply(abind(all.dat[,,,"0",], aperm(mim.1, c(3,1:2,4)), along = 4), 2:3, function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})


seq.mimic <- abind(mim.0, mim.1, mim.2, mim.3, mim.4, mim.5, mim.6, mim.7, mim.8, mim.9, mim.10, mim.11, mim.12, mim.13, mim.14, along = 2.5)

seq.rmse <- sqrt(apply(sweep(seq.mimic[,,,,"Tau"], c(4,1:2), obs, "-")^2, c(3:4), mean, na.rm = T))
mimic.rmse <- sqrt(apply(sweep(mimic[,,,,1], c(4,1:2), obs, "-")^2, c(3:4), mean, na.rm = T))

plot(obs[1,,1], type = "l", lwd = 2)
lines(mim.14[,1,1,"Tau"], col = "steelblue")
lines(mim.10[,1,1,"Tau"], col = "blue")

matplot(0:14, mimic.rmse[,1:2], type = "l")
matplot(0:14, seq.rmse[,1:2], col = c("blue", "green3"), type = "l", add = T)

# well, that looks very very wrong. Try again when more time...

# select single day & plot temp.n vs temp.s for obs, 3 ensembles, prior, posterior.
# represent each distribution as contours

d <- 20; y <- 5
x.rng <- c(0,7); y.rng <- c(3,10)


# SINGLE-DATE MIMIC USING PREVIOUS LEADTIME'S POSTERIOR AS PRIOR                                ####

d <- 20; y <- 5
x.rng <- c(0,7); y.rng <- c(3,10)

all.tr <- readRDS("../Models/tr-set-default.rds")

#---------------------------------------------------------------------------------
mim.14 <- fit.mimic(fc = list("ecmwf" = offset.forecast(ecmwf)[,d,y,15,-1],
                              "ncep" = offset.forecast(ncep)[,d,y,15,-1],
                              "ukmo" = offset.forecast(ukmo)[,d,y,15,-1]),
                    tr = all.tr[,d,y,15,])

# plot model progression
{
    dat.14 <- abind("obs" = obs[1:2,d,y],
                 "ec" = apply(offset.forecast(ecmwf)[1:2,d,y,15,-1], 1, mean),
                 "nc" = apply(offset.forecast(ncep)[1:2,d,y,15,-1], 1, mean),
                 "mo" = apply(offset.forecast(ukmo)[1:2,d,y,15,-1], 1, mean),
                 "post" = mim.14[1:2, "Tau"],
                 along = 0)
    plot(dat.14, pch = c(4, rep(20, 3), 1), lwd = 2, col = c("black", "steelblue", "red3", "green3", "black"),
         xlim = x.rng, ylim = y.rng, main = "14 days ahead")
}

#---------------------------------------------------------------------------------
lt <- 14
mim.13 <- fit.mimic(fc = list("ecmwf" = offset.forecast(ecmwf)[,d,y,lt,-1],
                              "ncep" = offset.forecast(ncep)[,d,y,lt,-1],
                              "ukmo" = offset.forecast(ukmo)[,d,y,lt,-1]),
                    tr = all.tr[,d,y,lt,],
                    prior = list(alpha = mim.14[,"Tau"],
                                 Gamma = mim.14[,-1]))

# plot model progression
{
    dat.13 <- abind("obs" = obs[1:2,d,y],
                 "ec" = apply(offset.forecast(ecmwf)[1:2,d,y,lt,-1], 1, mean),
                 "nc" = apply(offset.forecast(ncep)[1:2,d,y,lt,-1], 1, mean),
                 "mo" = apply(offset.forecast(ukmo)[1:2,d,y,lt,-1], 1, mean),
                 "post" = mim.13[1:2, "Tau"],
                 "prior" = mim.14[1:2, "Tau"],
                 along = 0)
    plot(dat.13, pch = c(4, rep(20, 3), 1, 20), lwd = 2, col = c("black", "steelblue", "red3", "green3", "black", "grey"),
         xlim = x.rng, ylim = y.rng, main = paste0(lt-1, " days ahead"))
}

#---------------------------------------------------------------------------------
lt <- 13
mim.12 <- fit.mimic(fc = list("ecmwf" = offset.forecast(ecmwf)[,d,y,lt,-1],
                              "ncep" = offset.forecast(ncep)[,d,y,lt,-1],
                              "ukmo" = offset.forecast(ukmo)[,d,y,lt,-1]),
                    tr = all.tr[,d,y,lt,],
                    prior = list(alpha = mim.13[,"Tau"],
                                 Gamma = mim.13[,-1]))

# plot model progression
{
    dat.12 <- abind("obs" = obs[1:2,d,y],
                    "ec" = apply(offset.forecast(ecmwf)[1:2,d,y,lt,-1], 1, mean),
                    "nc" = apply(offset.forecast(ncep)[1:2,d,y,lt,-1], 1, mean),
                    "mo" = apply(offset.forecast(ukmo)[1:2,d,y,lt,-1], 1, mean),
                    "post" = mim.12[1:2, "Tau"],
                    "prior" = mim.13[1:2, "Tau"],
                    along = 0)
    plot(dat.12, pch = c(4, rep(20, 3), 1, 20), lwd = 2, col = c("black", "steelblue", "red3", "green3", "black", "grey"),
         xlim = x.rng, ylim = y.rng, main = paste0(lt-1, " days ahead"))
}

#---------------------------------------------------------------------------------
lt <- 13
mim.12 <- fit.mimic(fc = list("ecmwf" = offset.forecast(ecmwf)[,d,y,lt,-1],
                              "ncep" = offset.forecast(ncep)[,d,y,lt,-1],
                              "ukmo" = offset.forecast(ukmo)[,d,y,lt,-1]),
                    tr = all.tr[,d,y,lt,],
                    prior = list(alpha = mim.13[,"Tau"],
                                 Gamma = mim.13[,-1]))

# plot model progression
{
    dat.12 <- abind("obs" = obs[1:2,d,y],
                    "ec" = apply(offset.forecast(ecmwf)[1:2,d,y,lt,-1], 1, mean),
                    "nc" = apply(offset.forecast(ncep)[1:2,d,y,lt,-1], 1, mean),
                    "mo" = apply(offset.forecast(ukmo)[1:2,d,y,lt,-1], 1, mean),
                    "post" = mim.12[1:2, "Tau"],
                    "prior" = mim.13[1:2, "Tau"],
                    along = 0)
    plot(dat.12, pch = c(4, rep(20, 3), 1, 20), lwd = 2, col = c("black", "steelblue", "red3", "green3", "black", "grey"),
         xlim = x.rng, ylim = y.rng, main = paste0(lt-1, " days ahead"))
}

pdf("../Plots/Sequential.pdf", width = 7, height = 4); {
    
    plot(dat.14, pch = c(4, rep(20, 3), 1), lwd = 2, col = c("black", "steelblue", "red3", "green3", "black"),
         xlim = x.rng, ylim = y.rng, main = "14 days ahead")
    legend("topleft", cex = 0.7, bty = "n", pch = c(4, rep(20, 3), 1), col = c("black", "steelblue", "red3", "green3", "black", "grey"),
           legend = c("Observation", "ECMWF mean", "NCEP mean", "UKMO mean", "posterior", "prior"))
    plot(dat.13, pch = c(4, rep(20, 3), 1, 20), lwd = 2, col = c("black", "steelblue", "red3", "green3", "black", "grey"),
         xlim = x.rng, ylim = y.rng, main = paste0("13 days ahead"))
    legend("topleft", cex = 0.7, bty = "n", pch = c(4, rep(20, 3), 1), col = c("black", "steelblue", "red3", "green3", "black", "grey"),
           legend = c("Observation", "ECMWF mean", "NCEP mean", "UKMO mean", "posterior", "prior"))
    plot(dat.12, pch = c(4, rep(20, 3), 1, 20), lwd = 2, col = c("black", "steelblue", "red3", "green3", "black", "grey"),
         xlim = x.rng, ylim = y.rng, main = paste0("12 days ahead"))
    legend("topleft", cex = 0.7, bty = "n", pch = c(4, rep(20, 3), 1), col = c("black", "steelblue", "red3", "green3", "black", "grey"),
           legend = c("Observation", "ECMWF mean", "NCEP mean", "UKMO mean", "posterior", "prior"))
    
}; dev.off()

####################################################################################################

# PLOT OF SYSTEM FOR SINGLE DATE                                                                ####

d <- 20; y <- 5; lt <- 14
mimic <- readRDS("~/Documents/PhD/Miniprojects/03-Sichun-paper/Models/an25-mimic.rds")
tr <- readRDS("../Models/tr-set-default.rds")


# maybe add convex hull of ensemble members?
draw.system <- function(d, y, lt, mimic, tr = NULL, prior = NULL, c.lvl = c(0.75, 0.95), xl, yl, ens.detail = F) {
    
    ec <- offset.forecast(ecmwf)[1:2,d,y,lt,-1]
    nc <- offset.forecast(ncep)[1:2,d,y,lt,-1]
    mo <- offset.forecast(ukmo)[1:2,d,y,lt,-1]
    
    if(missing(xl)) xl <- range(ec[1,], nc[1,], mo[1,], obs[1,d,y]) * c(0.9, 1.1)
    if(missing(yl)) yl <- range(ec[2,], nc[2,], mo[2,], obs[2,d,y]) * c(0.9, 1.1)
    
    
    # observation
    plot(t(obs[1:2,d,y]), pch = 4, lwd = 2, xlim = xl, ylim = yl,
         main = paste0("20", formatC(7:14, width = 2, flag = "0")[y], " d", d, "; LT ", lt - 1))
    
    
    # ensembles
    kd.ec <- kde2d(x = ec[1,], y = ec[2,], n = 1000, lims = c(xl, yl))
    kd.nc <- kde2d(x = nc[1,], y = nc[2,], n = 1000, lims = c(xl, yl))
    kd.mo <- kde2d(x = mo[1,], y = mo[2,], n = 1000, lims = c(xl, yl))
    
    image(kd.ec, breaks = quantile(kd.ec$z, c(c.lvl, 1)), col = mapply(adjustcolor, "skyblue", c.lvl - 0.4), add = T)
    image(kd.nc, breaks = quantile(kd.nc$z, c(c.lvl, 1)), col = mapply(adjustcolor, "red3", c.lvl - 0.4), add = T)
    image(kd.mo, breaks = quantile(kd.mo$z, c(c.lvl, 1)), col = mapply(adjustcolor, "green3", c.lvl - 0.4), add = T)
    
    if (ens.detail) {
        points(t(ec), col = "blue", pch = 20)
        points(t(nc), col = "red3", pch = 20)
        points(t(mo), col = "green3", pch = 20)
    }
    
    points(t(apply(ec, 1, mean)), bg = adjustcolor("skyblue", alpha = 0.6), pch = 21)
    points(t(apply(nc, 1, mean)), bg = adjustcolor("red3", alpha = 0.6), pch = 21)
    points(t(apply(mo, 1, mean)), bg = adjustcolor("green3", alpha = 0.6), pch = 21)
    
    
    # posterior
    sapply(c.lvl, function(ll) {
        lines(ellipse(mimic[d,y,lt,1:2,2:3], centre = mimic[d,y,lt,1:2, 1], level = ll))})
    points(t(mimic[d,y,lt,1:2, 1]), pch = 16, col = "black")
    
    # prior
    if (!is.null(prior)) {
        sapply(c.lvl, function(ll) {
            lines(ellipse(prior$Gamma, centre = prior$alpha, level = ll), lty = 2, col = "grey")})
        points(t(prior$alpha), col = "grey")
    }
    
    # training data - should be verifying obs, NOT error
    if (!is.null(tr)) {
        points(t(all.tr[1:2,d,y,lt,]), pch = 20, col = adjustcolor("grey", alpha = 0.5))
    }
    
    # restate the observation in case it's disappeared
    points(t(obs[1:2,d,y]), pch = 4, lwd = 2)
}

draw.system(2,1,15, mimic, xl = c(-2,12), yl = c(-2,15))      # all fairly well in agreement
draw.system(2,1,14, mimic, xl = c(-2,12), yl = c(-2,15))
draw.system(2,1,11, mimic, xl = c(-2,12), yl = c(-2,15))
draw.system(2,1,8, mimic, xl = c(-2,12), yl = c(-2,15))
draw.system(2,1,5, mimic, xl = c(-2,12), yl = c(-2,15))
draw.system(2,1,3, mimic, xl = c(-2,12), yl = c(-2,15))
draw.system(2,1,2, mimic, xl = c(-2,12), yl = c(-2,15))
draw.system(2,1,1, mimic, xl = c(-2,12), yl = c(-2,15), c.lvl = 0.95, ens.detail = T)

draw.system(30,5,15, mimic)      # all fairly well in agreement
draw.system(30,5,14, mimic)      # all fairly well in agreement
draw.system(30,5,1, mimic)      # all fairly well in agreement


####################################################################################################

# COMPARISON TO OTHER PRIORS USED                                                               ####

d <- 20; y <- 5; os.d <- (y*90)+d

oo <- apply(obs, 1, cbind)
persistence.prior <- list(alpha = obs[,d-15,y], Gamma = cov(apply(obs, 1, cbind)[os.d - 14 - (25:1),]))

pers.13 <- fit.mimic(fc = list("ecmwf" = offset.forecast(ecmwf)[,d,y,lt,-1],
                              "ncep" = offset.forecast(ncep)[,d,y,lt,-1],
                              "ukmo" = offset.forecast(ukmo)[,d,y,lt,-1]),
                    tr = all.tr[,d,y,lt,],
                    prior = list(alpha = persistence.prior$alpha,
                                 Gamma = persistence.prior$Gamma))

plot(dat.13, pch = c(4, rep(20, 3), 1, 1), lwd = 2, col = c("black", "steelblue", "red3", "green3", "black", "grey"),
     xlim = x.rng, ylim = y.rng, main = paste0(lt-1, " days ahead"))

points(t(persistence.prior$alpha), col = "grey", pch = 20)
points(t(pers.13[1:2,"Tau"]), pch = 20)

eta <- obs[,d,y] - apply(all.tr[,d,y,lt,], 1, mean)
points(t(eta[1:2]), pch = 4, lwd = 2, col = "magenta3")

# FILLED PLOTTING (WILL HAVE TO WAIT FOR NOW)                                                   ####
mimic <- readRDS("~/Documents/PhD/Miniprojects/03-Sichun-paper/Models/an25-mimic.rds")

# bivariate normal plots of density
d <- 5; y <- 5
mu <- apply(offset.forecast(ecmwf)[1:2,d,y,15,-1], 1, mean)
sig <- cov(t(offset.forecast(ecmwf)[1:2,d,y,15,-1]))

bvn.ellipse <- function(mu, sig, n = 100, x.rng = c(-15,15), y.rng = c(-15,15), c.lvl = c(0.5, 0.75, 0.95), ...) {
    
    x <- seq(x.rng[1], x.rng[2], length.out = n)
    y <- seq(y.rng[1], y.rng[2], length.out = n)
    
    bvn.density <- function(s,t) {
        u <- c(s,t) - mu
        u %*% solve(sig) %*% u / 2
    }
    
    z <- mapply(bvn.density, as.vector(rep(x,n)), as.vector(outer(rep(0,n), y, "+")))
    contour(x, y, matrix(z, n, n), drawlabel = F, levels = qchisq(c.lvl, 2), ...)
}

bvn.ellipse(mu, sig, xlim = c(-5, 10), ylim = c(-5, 10), n = 500, c.lvl = 0.75, col = "steelblue")
bvn.ellipse(mu = apply(offset.forecast(ncep)[1:2,d,y,15,-1], 1, mean),
            sig = cov(t(offset.forecast(ncep)[1:2,d,y,15,-1])), add = T, n = 500, c.lvl = 0.75, col = "red3")
bvn.ellipse(mu = apply(offset.forecast(ukmo)[1:2,d,y,15,-1], 1, mean),
            sig = cov(t(offset.forecast(ukmo)[1:2,d,y,15,-1])), add = T, n = 500, c.lvl = 0.75, col = "green3")
points(t(obs[1:2,d,y]), pch = 4, lwd = 2)

bvn.ellipse(mu = mimic[d,y,15,1:2,"Tau"], sig = mimic[d,y,15,1:2,2:3], col = "grey", add = T, c.lvl = 0.75)

ens.plot <- function(d, y, lt, mimic, xlim = c(-5, 10), ylim = c(-5,10), c.lvl = 0.75, n = 100, x.rng = c(-15,15), y.rng = c(-15,15), ...) {
    
    plot(t(obs[1:2,d,y]), pch = 4, lwd = 2, xlim = xlim, ylim = ylim)
    
    bvn.ellipse(mu = mimic[d,y,lt,1:2,"Tau"], sig = mimic[d,y,lt,1:2,2:3], lwd = 2, lty = 3, add = T, c.lvl = c.lvl)
    points(t(mimic[d,y,lt,1:2,"Tau"]), pch = 20)
    
    ec <- offset.forecast(ecmwf)[1:2,d,y,lt,-1]
    nc <- offset.forecast(ncep)[1:2,d,y,lt,-1]
    mo <- offset.forecast(ukmo)[1:2,d,y,lt,-1]
    
    #    points(t(ec), pch = 20, col = adjustcolor("steelblue", alpha = 0.4))
    #    points(t(nc), pch = 20, col = adjustcolor("red3", alpha = 0.4))
    #    points(t(mo), pch = 20, col = adjustcolor("green3", alpha = 0.4))
    
    contour(kde2d(x = ec[1,], y = ec[2,], n = n, lims = c(xlim, ylim)), 
            nlevel = 6, col = mapply(adjustcolor, "steelblue", c(0:5)/10), drawlabels = F, add = T)
    contour(kde2d(x = nc[1,], y = nc[2,], n = n, lims = c(xlim, ylim)), 
            nlevel = 6, col = mapply(adjustcolor, "red3", c(0:5)/10), drawlabels = F, add = T)
    contour(kde2d(x = mo[1,], y = mo[2,], n = n, lims = c(xlim, ylim)), 
            nlevel = 6, col = mapply(adjustcolor, "green3", c(0:5)/10), drawlabels = F, add = T)
}


kd.ec <- kde2d(x = ec[1,], y = ec[2,], n = n, lims = c(xlim, ylim))
filled.contour(kd.ec$x, kd.ec$y, kd.ec$z, nlevels = 6, col = mapply(adjustcolor, "steelblue", c(0, 0.1, 0.2, 0.3, 0.4, 0.5)))

kd.nc <- kde2d(x = nc[1,], y = nc[2,], n = n, lims = c(xlim, ylim))
filled.contour(kd.nc$x, kd.nc$y, kd.nc$z, nlevels = 6, col = mapply(adjustcolor, "red3", c(0, 0.1, 0.2, 0.3, 0.4, 0.5)))


ens.plot(20, 1, 15, mimic)
ens.plot(40, 1, 13, mimic)
ens.plot(40, 1, 2, mimic)


d <- 20; y <- 1; lt <- 1
smoothScatter(t(offset.forecast(ecmwf)[1:2,d,y,lt,-1]), xlim = c(-5, 10), ylim = c(-5,10), nrpoints = 0)



# snaffled from StackExchange. This definitely works.
{
    p <- rmvnorm(1000, c(250000, 20000), matrix(c(100000^2, 22000^2, 22000^2, 6000^2),2,2))
    center <- apply(p, 2, mean)
    sigma <- cov(p)
    
    sigma.inv = solve(sigma, matrix(c(1,0,0,1),2,2))
    ellipse <- function(s,t) {u<-c(s,t)-center; u %*% sigma.inv %*% u / 2}
    
    n <- 50
    x <- (0:(n-1)) * (500000/(n-1))
    y <- (0:(n-1)) * (50000/(n-1))
    
    z <- mapply(ellipse, as.vector(rep(x,n)), as.vector(outer(rep(0,n), y, "+")))
    plot(p, pch=20, xlim=c(0,500000), ylim=c(0,50000), xlab="Packets", ylab="Flows")
    contour(x,y,matrix(z,n,n), levels=(0:10), col = terrain.colors(11), add=TRUE)
    contour(x,y,matrix(z,n,n), levels = qchisq(c(0.5, 0.75, 0.95), 2), col = c("green3", "cyan3", "purple"), 
            lwd = 2, drawlabels = F, add = T)
}

# plot as filled contours using 'raster'
library('raster')

bvn.density <- function(mu, sig, n = 100, x.rng = c(-15,15), y.rng = c(-15,15)) {
    
    x <- seq(x.rng[1], x.rng[2], length.out = n)
    y <- seq(y.rng[1], y.rng[2], length.out = n)
    
    dbvn <- function(s,t) {
        u <- c(s,t) - mu
        u %*% solve(sig) %*% u / 2
    }
    
    z <- mapply(dbvn, as.vector(rep(x,n)), as.vector(outer(rep(0,n), y, "+")))
    z <- array(z, dim = rep(sqrt(length(z)), 2))
    list("x" = x, "y" = y, "z" = z)
}

z <- bvn.density(mu = apply(offset.forecast(ncep)[1:2,d,y,15,-1], 1, mean),
                 sig = cov(t(offset.forecast(ncep)[1:2,d,y,15,-1])), n = 100)
contour(z, levels = qchisq(c(0.75, 0.95), 2), drawlabel = F)
, levels = qchisq(c(0.75, 0.95), 2), col = 1:2, drawlabel = F)

ens <- ecmwf; d <- 10; y <- 5; lt <- 15

plot.kde <- function(ens, d, y, lt, c.lvl = c(0.5, 0.75, 0.95), col = "steelblue",
                     xlim = c(-5, 10), ylim = c(-5,10), n = 100) {
    
    ec <- offset.forecast(ens)[1:2,d,y,lt,-1]
    kd.ec <- kde2d(x = ec[1,], y = ec[2,], n = n, lims = c(xlim, ylim))
    
    contour(kd.ec)
}

contour(z, levels = qchisq(.75, 2), drawlabels = F)

z.75 <- which(z$z < qchisq(.75, 2), arr.ind = T)
z.75[,1] <- z$x[z.75[,1]];  z.75[,2] <- z$y[z.75[,2]]
points(z.75, pch = 20, col = adjustcolor("steelblue", alpha = 0.3))
polygon(z.75)
polygon(z$x[chull(z.75)], z$y[chull(z.75)], col = "red")

# support function to calculate density of bivariate normal distribution
bvn.density <- function(mu, sig, n = 100, x.rng = c(-15,15), y.rng = c(-15,15)) {
    
    x <- seq(x.rng[1], x.rng[2], length.out = n)
    y <- seq(y.rng[1], y.rng[2], length.out = n)
    
    dbvn <- function(s,t) {
        u <- c(s,t) - mu
        u %*% solve(sig) %*% u / 2
    }
    
    z <- mapply(dbvn, as.vector(rep(x,n)), as.vector(outer(rep(0,n), y, "+")))
    z <- array(z, dim = rep(sqrt(length(z)), 2))
    list("x" = x, "y" = y, "z" = z)
}