
# effect of choosing different priors

library("SX.weather")
setwd("~/Documents/PhD/Miniprojects/03-Sichun-paper/Code")

tr.error <- readRDS("../Models/tr-set-default.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plotting function                                                                             ####

draw.system <- function(d, y, lt, mimic, tr = NULL, prior = NULL, c.lvl = c(0.75, 0.95), xl, yl, ens.detail = F,
                        main = paste0("20", formatC(7:14, width = 2, flag = "0")[y], " d", d, "; LT ", lt - 1)) {
    
    mimic <- mimic[d,y,lt,,]
    ec <- offset.forecast(ecmwf)[1:2,d,y,lt,-1]
    nc <- offset.forecast(ncep)[1:2,d,y,lt,-1]
    mo <- offset.forecast(ukmo)[1:2,d,y,lt,-1]
    
    if(missing(xl)) xl <- range(ec[1,], nc[1,], mo[1,], obs[1,d,y]) * c(0.9, 1.1)
    if(missing(yl)) yl <- range(ec[2,], nc[2,], mo[2,], obs[2,d,y]) * c(0.9, 1.1)
    
    
    # observation
    plot(t(obs[1:2,d,y]), pch = 4, lwd = 2, xlim = xl, ylim = yl,
         main = main)
    
    
    # ensembles
    kd.ec <- kde2d(x = ec[1,], y = ec[2,], n = 1000, lims = c(xl * c(0.9, 1.1), yl * c(0.9, 1.1)))
    kd.nc <- kde2d(x = nc[1,], y = nc[2,], n = 1000, lims = c(xl * c(0.9, 1.1), yl * c(0.9, 1.1)))
    kd.mo <- kde2d(x = mo[1,], y = mo[2,], n = 1000, lims = c(xl * c(0.9, 1.1), yl * c(0.9, 1.1)))
    
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
    #if (!is.null(tr)) {
    #    points(t(tr[1:2,d,y,lt,]), pch = 4, lwd = 2, col = adjustcolor("grey", alpha = 0.5))
    #}
    
    # restate the observation in case it's disappeared
    points(t(obs[1:2,d,y]), pch = 4, lwd = 2)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Nninformative prior                                                                           ####

uninf.mimic <- readRDS("../Models/an25-mimic.rds")
uninf.rmse <- sqrt(apply(sweep(uninf.mimic[,,,,1], c(4,1:2), obs, "-")^2, c(3:4), mean, na.rm = T))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Analogue prior                                                                                ####

an.indx <- readRDS("../Models/Analogue-indices-25.rds")
tr.obs <- aaply(an.indx, 1:3, function(a) apply(obs, 1, "[", a))

# alpha & Gamma taken from analogue observations (same indices as training data)
a.data <- abind(invisible(sapply(1:15, function(lt) {
    abind("ec" = offset.forecast(ecmwf)[,,,lt,-1],
          "nc" = offset.forecast(ncep)[,,,lt,-1],
          "mo" = offset.forecast(ukmo)[,,,lt,-1], 
          "tr" = tr.error[,,,lt,], 
          "alpha" = apply(tr.obs[,,lt,,], c(4, 1:2), mean),
          "Gamma" = aperm(aaply(tr.obs[,,lt,,], c(1:2), cov), c(3, 1:2, 4)),
          along = 4)}, simplify = F)), along = 0)

analogue.mimic <- aaply(a.data, c(1,3:4), function(arr) {
    fit.mimic(fc = list("ecmwf" = arr[,1:50],
                        "ncep" = arr[,51:71],
                        "ukmo" = arr[,72:93]),
              tr = arr[,94:118],
              prior = list(alpha = arr[,119],
                           Gamma = arr[,120:124]))
})
analogue.mimic <- aperm(analogue.mimic, c(4,3,5,2,1))

analogue.rmse <- sqrt(apply(sweep(analogue.mimic[,,,,"Tau"], c(4,1:2), obs, "-")^2, 3:4, mean, na.rm = T))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Sequential prior                                                                              ####

all.dat <- abind("ec" = offset.forecast(ecmwf)[,,,,-1],
                 "nc" = offset.forecast(ncep)[,,,,-1],
                 "mo" = offset.forecast(ukmo)[,,,,-1],
                 "tr" = readRDS("../Models/tr-set-default.rds"),
                 along = 5)

seq.mimic <- list()
invisible(lapply(15:1, function(lt) {
    
    if(lt == 15) {
        seq.mimic[[16-lt]] <<- suppressWarnings(aaply(all.dat[,,,lt,], 2:3, function(arr) {
            fit.mimic(fc = list("ecmwf" = arr[,1:50],
                                                      "ncep" = arr[,51:71],
                                                      "ukmo" = arr[,72:93]),
                                            tr = arr[,94:118])
        }))
    } else {
        seq.mimic[[16-lt]] <<- suppressWarnings(aaply(abind(all.dat[,,,lt,], 
                                                           aperm(seq.mimic[[16-lt-1]], c(3,1:2,4)), along = 4), 2:3, 
                                                     function(arr) {
            fit.mimic(fc = list("ecmwf" = arr[,1:50],
                                "ncep" = arr[,51:71],
                                "ukmo" = arr[,72:93]),
                      tr = arr[,94:118],
                      prior = list(alpha = arr[,119],
                                   Gamma = arr[,120:124]))
        }))
    }
}))

seq.mimic <- abind(rev(seq.mimic), along = 2.5)
seq.rmse <- sqrt(apply(sweep(seq.mimic[,,,,"Tau"], c(4,1:2), obs, "-")^2, 3:4, mean, na.rm = T))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# comparison plots                                                                              ####

rmse <- abind("No prior" = uninf.rmse, 
              "Analogue prior" = analogue.rmse, 
              "Sequential prior" = seq.rmse, along = 0)[,,1:2]

# RMSE plots
pdf("../Plots/Prior-comparison/RMSE.pdf", height = 4); {
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(rmse)[[3]], function(varb) {
        matplot(0:14, t(rmse[,,varb]), type = "l", lty = 1, ylim = range(0, rmse), main = varb)      
        legend("bottomright", bty = "n", col = c("black", "red", "green3"), lty = 1, dimnames(rmse)[[1]])
    }))
    mtext("RMSE with different prior distributions", outer = T)
}; dev.off()


d <- 5; y <- 5; lt <- 14

# system plots
pdf("../Plots/Prior-comparison/system-lt14.pdf", height = 4); {
    par(mfrow = c(1,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    draw.system(d,y,lt, uninf.mimic, tr = tr.error, main = "Uninformative prior",
                xl = c(-5,10), yl = c(0,10))
    draw.system(d,y,lt, analogue.mimic, tr = tr.error,
                prior = list(alpha = apply(tr.obs[d,y,lt,,], 2, mean), Gamma = cov(tr.obs[d,y,lt,,])),
                main = "Analogue prior", xl = c(-5,10), yl = c(0,10))
    draw.system(d,y,lt, seq.mimic, tr = tr.error,
                prior = list(alpha = seq.mimic[d,y,lt+1,,1], Gamma = seq.mimic[d,y,lt+1,,-1]),
                main = "Sequential prior", xl = c(-5,10), yl = c(0,10))
    mtext(paste0("Example: 20", formatC(7:14, width = 2, flag = "0")[y], " d", d, ", LT ", lt - 1), outer = T)
}; dev.off()

# compare verification rank histograms - marginal PIT

uninf.pit <- aaply(uninf.mimic[,,,1:2, 1:3], 3, function(mim) {
    abind(sapply(dimnames(mim)[[3]], function(varb) {
        pnorm(obs[varb,,], mim[,,varb,"Tau"], sqrt(mim[,,varb,varb]))
    }), along = 0)
})
analogue.pit <- aaply(analogue.mimic[,,,1:2, 1:3], 3, function(mim) {
    abind(sapply(dimnames(mim)[[3]], function(varb) {
        pnorm(obs[varb,,], mim[,,varb,"Tau"], sqrt(mim[,,varb,varb]))
    }), along = 0)
})
seq.pit <- aaply(seq.mimic[,,,1:2, 1:3], 3, function(mim) {
    abind(sapply(dimnames(mim)[[3]], function(varb) {
        pnorm(obs[varb,,], mim[,,varb,"Tau"], sqrt(mim[,,varb,varb]))
    }), along = 0)
})

lt <- 14
pdf(paste0("../Plots/Prior-comparison/PIT-hists-",lt-1,".pdf"), height = 4); {
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(dimnames(uninf.pit)[[3]], function(varb) {
        hist(uninf.pit[lt,,varb], prob = T, main = varb, col = "skyblue")
        abline(h = 1, col = "red")
        mtext(paste0("Uninformative prior - leadtime ", lt-1), outer = T)
    }))
    invisible(sapply(dimnames(analogue.pit)[[3]], function(varb) {
        hist(analogue.pit[lt,,varb], prob = T, main = varb, col = "skyblue")
        abline(h = 1, col = "red")
        mtext(paste0("Analogue prior - leadtime ", lt-1), outer = T)
    }))
    invisible(sapply(dimnames(seq.pit)[[3]], function(varb) {
        hist(seq.pit[lt,,varb], prob = T, main = varb, col = "skyblue")
        abline(h = 1, col = "red")
        mtext(paste0("Sequential prior - leadtime ", lt-1), outer = T)
    }))
}; dev.off()

# divergent ensemble...
par(mfrow = c(1,1))
plot(t(ecmwf[1:2,3,"11",10,-1]), col = "red3", pch = 20)
