
# should prior be multiplicative rather than additive?
library("SX.weather")
require(ellipse)
setwd("~/Documents/PhD/Miniprojects/03-Sichun-paper/Code")

# sequential updating (with plots) for a single day

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
        lines(ellipse(mimic[1:2,2:3], centre = mimic[1:2, 1], level = ll))})
    points(t(mimic[1:2, 1]), pch = 16, col = "black")
    
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

# check forecast alignment
all(ecmwf[,2,,"14",] == offset.forecast(ecmwf)[,1,,"14",])
all(ecmwf[,16,,"0",] == offset.forecast(ecmwf)[,1,,"0",])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# target is December 20th: so 14-day forecast made on December 6th
d <- 20; y <- 5; lt <- "14"

# check forecast alignment
all(ecmwf[,d+(15-as.numeric(lt)),y, lt,] == offset.forecast(ecmwf)[,d,y,lt,])

tr <- readRDS("../Models/tr-set-default.rds")
dimnames(tr)[[4]] <- formatC(0:14)

# LEADTIME 14                                                                                       ####

m.14 <- fit.mimic(fc = list("ecmwf" = offset.forecast(ecmwf)[,d,y,lt,-1],
                              "ncep" = offset.forecast(ncep)[,d,y,lt,-1],
                              "ukmo" = offset.forecast(ukmo)[,d,y,lt,-1]),
                    tr = tr[,d,y,lt,])
draw.system(d, y, as.numeric(lt)+1, m.14)

# LEADTIME 13                                                                                       ####

lt <- "13"
prior.13 <- list("alpha" = m.14[,1], "Gamma" = m.14[,-1])

m.13 <- fit.mimic(fc = list("ecmwf" = offset.forecast(ecmwf)[,d,y,lt,-1],
                            "ncep" = offset.forecast(ncep)[,d,y,lt,-1],
                            "ukmo" = offset.forecast(ukmo)[,d,y,lt,-1]),
                  tr = tr[,d,y,lt,],
                  prior = prior.13)

draw.system(d, y, as.numeric(lt)+1, m.13, prior = prior.13)

# LEADTIME 12                                                                                       ####

lt <- "12"
prior.12 <- list("alpha" = m.13[,1], "Gamma" = m.13[,-1])

m.12 <- fit.mimic(fc = list("ecmwf" = offset.forecast(ecmwf)[,d,y,lt,-1],
                            "ncep" = offset.forecast(ncep)[,d,y,lt,-1],
                            "ukmo" = offset.forecast(ukmo)[,d,y,lt,-1]),
                  tr = tr[,d,y,lt,],
                  prior = prior.12)

draw.system(d, y, as.numeric(lt)+1, m.12, prior = prior.12)

########################################################################################################

# look at elemtnts of tau for LT12 model

fc <- list("ecmwf" = offset.forecast(ecmwf)[1:2,d,y,lt,-1],
                "ncep" = offset.forecast(ncep)[1:2,d,y,lt,-1],
                "ukmo" = offset.forecast(ukmo)[1:2,d,y,lt,-1])

Y.bar <- abind(lapply(fc, apply, 1, mean), along = 0)
C <- abind(lapply(lapply(fc, t), cov), along = 0)
D.i <- aaply(sweep(sweep(C, 1, sapply(fc, dim)[2, ], "/"), 
                   2:3, cov(Y.bar), "+"), 1, solve)

Eta <- apply(t(tr[1:2,d,y,lt,]), 2, mean)
Lambda <- cov(t(tr[1:2,d,y,lt,]))

prior <- prior.12

# posterior without prior
S.u <- Lambda + solve(apply(D.i, 2:3, sum))
tau.u <- (S.u %*% ((solve(diag(ncol(Y.bar)) + (apply(D.i, 2:3, sum) %*% 
                                              Lambda))) %*% apply(aaply(abind(Y.bar, D.i, along = 2), 
                                                                        1, function(arr) arr[-1, ] %*% arr[1, ]), 2, sum))) - 
    Eta

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot.density <- function(ens, col = "skyblue", n = 1000, c.lvl = c(0.75, 0.95), add = T, xl = c(-10,15), yl = c(-5,15)) {

    ens.kd <- kde2d(x = ens[1,], y = ens[2,], n = n, lims = c(xl, yl))
    
    image(ens.kd, breaks = quantile(ens.kd$z, c(c.lvl, 1)), col = mapply(adjustcolor, col, c.lvl - 0.4), add = add)
    
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot(Y.bar, pch = 21, bg = c("steelblue", "green3", "red3"), xlim = c(-10,15), ylim = c(-5,15))

lines(ellipse(cov(t(fc$ecmwf)), centre = Y.bar["ecmwf",]), col = "steelblue")
lines(ellipse(cov(t(fc$ncep)), centre = Y.bar["ncep",]), col = "red3")
lines(ellipse(cov(t(fc$ukmo)), centre = Y.bar["ukmo",]), col = "green3")

# show actual ensemle density: 0.75 & 0.95-quantiles
{
    plot.density(fc$ecmwf)
    plot.density(fc$ncep, col = "coral")
    plot.density(fc$ukmo, col = "green3")
}

points(t(prior$alpha[1:2]), pch = 3, lwd = 2, col = "grey")
lines(ellipse(prior$Gamma[1:2,1:2], centre = prior$alpha[1:2]), col = "grey")

points(t(obs[1:2,d,y]), pch = 4, lwd = 3)

lines(ellipse(S.u, centre = t(tau.u)))
points(t(tau.u), pch = 20)

# elements of final output
prior.contribution <- solve(prior$Gamma) %*% prior$alpha
points(t((solve(prior$Gamma) %*% prior$alpha)[1:2]), col = "magenta")
