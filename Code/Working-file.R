
library("SX.weather")
library("CB.Misc")

# clean copy of code recreating Sichun's plots & results

    # further metrics: histograms etc
    # look at varying elements: what actually changes? (only S/precision, presumably)
    # write something up, including comparison of S & Tau for each model (0,2,3 PCs)

    # try treating superensemble as single ensemble, and fitting a single model. Difference?lunch

####################################################################################################

# SPATIAL PLOTS                                                                                 ####

# figure 1.2 - spatial plot of principal components
map.principal.components <- function() {
    espace <- load.data("./Data//ERAint_pca_espace.rda")
    lon.Europe <- load.data("./Data/ECMWF_europe_lon.rda")
    lat.Europe <- load.data("./Data/ECMWF_europe_lat.rda")
    
    org.par <- par()
    
    par(mfrow = c(1,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(1:3, function(p) {
        image.plot(lon.Europe, lat.Europe, matrix(data = espace$vectors[,p], nrow=51),
                   xlab = "", ylab = "", main = paste0("PC", p), zlim = range(espace$vectors[,1:3]),
                   xaxt = "n", yaxt = "n", asp = T)
        map("world", add = T)
    }))
    
    mtext("Map of first 3 principal components", outer = T)
    
    par(mfrow = org.par$mfrow, mar = org.par$mar, oma = org.par$oma)
}
map.principal.components()

####################################################################################################

# INDIVIDUAL/GROUPED FORECAST PERFORMANCE                                                       ####

# figure 4.1 - CDF of one-day-ahead forecast ensemble for temp.S on Dec 1st 2007.
plot.forecast.cdf <- function(model = ukmo, varb = "temp.s", yr = "08", d = 1, lt = 1) {
    
    org.par <- par()
    
    yr <- formatC(yr, width = 2, flag = "0")
    lt <- toString(lt)
    
    y <- sort(c(0, offset.forecast(model)[varb, d, yr, lt, -1]))
    l <- length(y)-1
    o <- obs[varb,d,yr]
    
    x <- (0:l)/l  
    p <- sort(c((0:l)/l, sum(y < o)/l))
    z <- sort(c(y, o)) < o
    
    crps <- (p - (1-z))^2
    
    rng <- range(offset.forecast(model)[varb, d, yr, lt, -1])
    
    par(mar = c(3,2,3,1))
    plot(y, x, type = 's', ylab = "", xlab = varb, col = "red", lwd = 2, xlim = rng,
         main = paste0("CDF & CRPS of ", varb, " for ", lt, "-day-ahead forecast (day ", d, ", 20", yr, ")"))
    lines(c(min(rng)-1, o, max(rng) + 1), c(0,1,1), type = "s", lwd = 2, col = "blue")
    lines(sort(c(y, o)), crps, lwd = 2, type = "s", col = "green3")
    legend("topleft", col = c("red", "blue", "green3"), lty = 1, lwd = 2, bty = "n", cex = 0.8,
           legend = c("Forecast Ensemble", "Observation", "CRPS"))
    par(mar = org.par$mar)
}
plot.forecast.cdf()

# figure 5.1 - mean temperature error
plot.forecast.errors <- function(model) {
    
    org.par <- par()
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    if (model == "super") {
        super <- superensemble(list(ecmwf, ncep, ukmo))
        
        s.err <- forecast.errors(super)
        s.mean.error <- apply(s.err, c(1, 4, 5), mean)
        
        ens.cols <- c("blue", "red", "green3")
        
        invisible(sapply(dimnames(s.mean.error)[[1]][c(3:5, 1:2)], function(varb) {
            
            # plot all ensemble members
            matplot(s.mean.error[varb, ,], type = "l", lty = 1,
                    col = c(rep(NA, 3), adjustcolor(ens.cols[attributes(super)$m], alpha = 0.1)),
                    xlab = "", ylab = "", main = varb)
            
            # ensemble means
            invisible(sapply(1:3, function(n) lines(apply(s.mean.error[varb,,attributes(super)$m == n], 1, mean),
                                                    col = ens.cols[n], lty = 2)))
            
            # controls
            matplot(s.mean.error[varb,,1:3], lwd = 1, col = ens.cols, add = T, type = "l", lty = 1)
            
            # overall ens. mean
            lines(apply(s.err[varb,,,,-(1:3)], 3, mean), lty = 2, lwd = 2)
            
            # overall control mean
            lines(apply(s.err[varb,,,,(1:3)], 3, mean), lty = 1, lwd = 2)
            
            abline(h = 0, col = "darkgrey", lty = 3)
        }))
        
        plot.new()
        legend("left", lty = c(1,1,1,1,NA,1,2), lwd = c(1,1,1,1,NA,2,2), col = c(ens.cols, rep("black", 4)), bty = "n",
               legend = c("ECMWF", "NCEP", "UKMO", "Superensemble", NA, "Control forecast", "Perturbations"))
        
        mtext("Superensemble - mean forecast error at each leadtime", outer = T)
    } else {
        
        # PLOT FORECAST ERRORS FOR SINGLE MODEL
        # add ensemble mean to model
        m.plus <- abind("em" = apply(model[,,,,-1], 1:4, mean), model, along = 5)
        
        errors <- forecast.errors(m.plus)
        mean.errors <- apply(errors, c(1, 4, 5), mean)
        
        invisible(sapply(dimnames(mean.errors)[[1]][c(3:5, 1:2)], function(varbl) {
            
            matplot(mean.errors[varbl,,dim(mean.errors)[[3]]:1], type = "l", lty = 1,
                    col = c(rep(adjustcolor("grey", alpha = 0.5), dim(mean.errors)[[3]]-2), "blue", "black"),
                    xlab = "", ylab = "", main = varbl)
        }))
        
        # add legend in final box
        plot.new()
        legend("left", lty = 1, col = c("blue", "black", adjustcolor("grey", alpha = 0.5)), bty = "n", cex = 1.1,
               legend = c("Control forecast", "Perturbed mean", "Perturbed members"))
        
        # add overall title
        mtext(paste0(toupper(toString(as.list(match.call())$model)), " mean error at each forecast lead time"),
              outer = TRUE, cex = 1)
        
        # reset device parameters to original values
        par(mfrow = org.par$mfrow, mar = org.par$mar, oma = org.par$oma)
    }
    
    par(mfrow = org.par$mfrow, mar = org.par$mar, oma = org.par$oma)
}

plot.forecast.errors(ecmwf)
plot.forecast.errors(ncep)
plot.forecast.errors(ukmo)
plot.forecast.errors("super")       # overlay all forecasts on single plot

# figure 5.2 - RMSE
plot.forecast.rmse <- function(model) {
    
    org.par <- par()
    
    # add ensemble mean to model
    m.plus <- abind("em" = apply(model[,,,,-1], 1:4, mean), model, along = 5)
    
    rmse <- forecast.rmse(m.plus)
    rmse <- abind("em2" = apply(rmse[,,-c(1:2)], 1:2, mean), rmse)
    
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(rmse)[[1]][c(3:5, 1:2)], function(varbl) {
        
        matplot(rmse[varbl,,dim(rmse)[[3]]:1], type = "l", lty = 1,
                col = c(rep(adjustcolor("grey", alpha = 0.5), dim(rmse)[[3]]-3), "blue", "black", "red"),
                xlab = "", ylab = "", main = varbl)
    }))
    
    # add legend in final box
    plot.new()
    legend("left", lty = 1, col = c("red", "blue", "black", adjustcolor("grey", alpha = 0.5)), bty = "n", cex = 1.1,
           legend = c("Mean of perturbed RMSE", "Control forecast", "RMSE of perturbed mean", "Perturbed members"))
    
    # add overall title
    mtext(paste0(toupper(toString(as.list(match.call())$model)), " RMSE at each forecast lead time"),
          outer = TRUE, cex = 1)
    
    # reset device parameters to original values
    par(mfrow = org.par$mfrow, mar = org.par$mar, oma = org.par$oma)
}
# not run - basically replicated with more information in 5.3

# figure 5.3 - RMSE vs spread
plot.rmse.spread <- function(model) {
    
    org.par <- par()
    
    spread <- ensemble.spread(model[,,,,-1])
    
    # add ensemble mean to model
    m.plus <- abind("em" = apply(model[,,,,-1], 1:4, mean), model, along = 5)
    rmse <- forecast.rmse(m.plus)
    rmse <- abind("em2" = apply(rmse[,,-c(1:2)], 1:2, mean), rmse)
    
    dat <- abind("spread" = spread, rmse, along = 3)
    
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(dat)[[1]][c(3:5, 1:2)], function(varbl) {
        
        matplot(dat[varbl,,dim(dat)[[3]]:1], type = "l", lty = c(rep(1, dim(dat)[[3]]-1), 2),
                col = c(rep(adjustcolor("grey", alpha = 0.5), dim(dat)[[3]]-4), "blue", "black", "red", "green3"),
                xlab = "", ylab = "", main = varbl)
    }))
    
    # add legend in final box
    plot.new()
    legend("left", lty = c(2, rep(1, dim(dat)[[3]]-1)), col = c("green3", "red", "blue", "black", adjustcolor("grey", alpha = 0.5)),
           bty = "n", cex = 1.1,
           legend = c("Ensemble spread", "RMSE of all perturbations", "Control forecast RMSE",
                      "RMSE of mean perturbation", "Perturbed member RMSE"))
    
    # add overall title
    mtext(paste0(toupper(toString(as.list(match.call())$model)), " RMSE vs spread at each forecast lead time"),
          outer = TRUE, cex = 1)
    
    # reset device parameters to original values
    par(mfrow = org.par$mfrow, mar = org.par$mar, oma = org.par$oma)
    
}
plot.rmse.spread(ecmwf)
plot.rmse.spread(ncep)
plot.rmse.spread(ukmo)

# figure 5.4 - CRPS 
# not run here - included in full evaluation of all models in 5.7

####################################################################################################

# POST-PROCESSING                                                                               ####

# models including 0, 2 and 3 principal components are created for comparison
# run full model (3 principal components)
    model.3pc <- run.model()
    model.2pc <- run.model(c("temp.n", "temp.s", "pc1", "pc2"))
    model.0pc <- run.model(c("temp.n", "temp.s"))

# save models to avoid having to re-fit later
    # saveRDS(model.3pc, "./Data/fitted-3pc.rds")
    # saveRDS(model.2pc, "./Data/fitted-2pc.rds")
    # saveRDS(model.0pc, "./Data/fitted-0pc.rds")
 
model.3pc.res <- model.performance(model.3pc)
model.2pc.res <- model.performance(model.2pc)
model.0pc.res <- model.performance(model.0pc)

# or just import fitted models & evaluate directly
model.3pc.res <- model.performance(readRDS("./Data/fitted-3pc.rds"))
model.2pc.res <- model.performance(readRDS("./Data/fitted-2pc.rds"))
model.0pc.res <- model.performance(readRDS("./Data/fitted-0pc.rds"))

crps.mat <- abind("ecmwf" = ensemble.crps(ecmwf[,,,,-1]),
                  "ncep" = ensemble.crps(ncep[,,,,-1]),
                  "ukmo" = ensemble.crps(ukmo[,,,,-1]), 
                  "super" = ensemble.crps(superensemble()[,,,,-(1:3)]),
                  "no.pc" = rbind(apply(model.0pc.res$crps, c(1,4), mean), NA, NA, NA),     # padded
                  "two.pc" = rbind(apply(model.2pc.res$crps, c(1,4), mean), NA),            # padded
                  "all.pc" = apply(model.3pc.res$crps, c(1,4), mean),
                  along = 0)

# figure 5.6 - RMSE vs spread for post-processed models
plot.postprocessed.rmse.spread <- function(ccols = c("black", "red", "orange"), weights = c(3,2,1)) {
    
    res.mat <- abind("rmse" = abind("no.pc" = rbind(model.0pc.res$rmse, NA, NA, NA),     # padded
                                    "two.pc" = rbind(model.2pc.res$rmse, NA),            # padded
                                    "all.pc" = model.3pc.res$rmse, rev.along = 0),
                     "spread" = abind("no.pc" = rbind(model.0pc.res$spread, NA, NA, NA),     # padded
                                      "two.pc" = rbind(model.2pc.res$spread, NA),            # padded
                                      "all.pc" = model.3pc.res$spread, rev.along = 0), rev.along = 0)
    
    org.par <- par()
    
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(res.mat)[[1]][c(3:5,1:2)], function(varb) {
        matplot(t(apply(res.mat[varb,,,], 1, rbind)), type = "l", xlab = "", ylab = "", main = varb,
                lty = rep(c(1,3), each = 3), col = rep(ccols, 2), 
                lwd = rep(weights, 2), ylim = range(c(0,res.mat[varb,,,]), na.rm = T))
    }))
    
    plot.new()
    legend("center", lty = rep(c(1,3), each = 3), col = rep(ccols, 2), 
           lwd = rep(weights, 2), bty = "n",
           legend = c("RMSE - temp only", "RMSE - first 2 PCs", 
                      "RMSE - first 3 PCs", "Spread - temp only",
                      "Spread - first 2 PCs","Spread - first 3 PCs"))
    
    mtext("RMSE vs spread for post-processed superensemble", outer = T)
    
    par(mfrow = org.par$mfrow, mar = org.par$mar, oma = org.par$oma)
}
plot.postprocessed.rmse.spread()

# figure 5.7 - CRPS including post-processing
plot.crps.all <- function() {
    
    org.par <- par()
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    invisible(sapply(dimnames(crps.mat)[[2]][c(3:5,1:2)], function(varb) {
        matplot(t(crps.mat[,varb,]), type = "l", main = varb, xlab = "", ylab = "",
                lty = c(rep(1, 3), 2, rep(1,3)), lwd = c(rep(1,3), rep(2, 4)),
                col = c("cornflowerblue", "coral2", "aquamarine3", "darkslategrey",
                        "green", "magenta3", "black"))
    }))
    plot.new()
    legend("center", lty = c(rep(1, 3), 2, rep(1,4)), lwd = c(rep(1,3), rep(2, 5)), 
           col = c("cornflowerblue", "coral2", "aquamarine3", "darkslategrey", NA,
                   "green", "magenta3", "black"), bty = "n",
           legend = c("ECMWF", "NCEP", "UKMO", "Superensemble", NA,
                      "Fitted - no principal components", 
                      "Fitted - 2 principal components", "Fitted - 3 principal components"))
    mtext("CRPS for each model, at each leadtime", outer = T)
    
    par(mfrow = org.par$mfrow, mar = org.par$mar, oma = org.par$oma)
}
plot.crps.all()

# figure 5.8 - RMSE & CRPS for model with no principal components

# additional: plot RMSE vs spread for all models
# including superensemble average and fitted model
plot.rmse.spread.all <- function(ccols = c("cornflowerblue", "coral2", "aquamarine3", "red", "black")) {
    
    # only run if not already publicly accessible
    if(!exists("res")) {
        res <- abind("rmse" = abind("ecmwf" = forecast.rmse(apply(ecmwf[,,,,-1], 1:4, mean)),
                                    "ncep" = forecast.rmse(apply(ncep[,,,,-1], 1:4, mean)),
                                    "ukmo" = forecast.rmse(apply(ukmo[,,,,-1], 1:4, mean)),
                                    "super" = sqrt(apply(apply(forecast.errors(superensemble())[,,,,-(1:3)], 1:4, mean)^2,
                                                         c(1,4), mean)),
                                    "pc3" = model.3pc.res$rmse, 
                                    rev.along = 0),
                     "spread" = abind("ecmwf" = ensemble.spread(ecmwf[,,,,-1]),
                                      "ncep" =  ensemble.spread(ncep[,,,,-1]),
                                      "ukmo" =  ensemble.spread(ukmo[,,,,-1]),
                                      "super" = ensemble.spread(superensemble()[,,,,-(1:3)]),
                                      "pc3" = model.3pc.res$spread,
                                      rev.along = 0),
                     rev.along = 0)
    }
    
    org.par <- par()
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(dimnames(res)[[1]][c(3:5,1:2)], function(varb) {
        
        matplot(t(apply(res[varb,,,], 1, cbind)), type = "l", ylim = range(c(0, res[varb,,,])),
                lty = rep(c(1,2), each = 5), lwd = rep(c(1,1,1,2,2),2), xlab = "", ylab = "", main = varb,
                col = rep(ccols, 2))
    }))
    plot.new()
    legend("center",  bty = "n", lty = c(rep(1, 5)), lwd = c(1,1,1,2,2), 
           col = ccols,
           legend = c("ECMWF", "NCEP", "UKMO", "Superensemble", "Fitted - 3 principal components"))
    mtext("RMSE (solid) vs spread (dashed) at each leadtime", outer = T)
    
    par(mfrow = org.par$mfrow, mar = org.par$mar, oma = org.par$oma)
}
plot.rmse.spread.all()

####################################################################################################

# CHECK ELEMENTS BY PRINCIPAL COMMPONENT                                                        ####

cov.3 <- se.covariances(offset.forecast(ecmwf))
cov.0 <- se.covariances(offset.forecast(ecmwf[1:2,,,,]))
all(cov.0 == cov.3[1:2,1:2,,,])

sig.3 <- se.sigma()
sig.0 <- se.sigma(list(ecmwf[1:2,,,,], ncep[1:2,,,,], ukmo[1:2,,,,]))
all(sig.3[1:2,1:2,,,] == sig.0)

lambda.3 <- se.lambda()
eme.0 <- apply(forecast.errors(superensemble())[1:2,,,,-(1:3)], 1:4, mean)
lambda.0 <- {array(apply(aperm(eme.0, c(3,1,2,4)), c(3:4), cov),
                  dim = c(dim(eme.0)[[1]], dim(eme.0[,,1,])),
                  dimnames = append(dimnames(eme.0[,,1,]), 
                                    list(dimnames(eme.0)[[1]]), after = 1))}
all(lambda.0 == lambda.3[1:2,1:2,,])

eta.3 <- apply(apply(forecast.errors(superensemble())[,,,,-(1:3)], 1:4, mean), c(1,2,4), mean)
eta.0 <- apply(apply(forecast.errors(superensemble())[1:2,,,,-(1:3)], 1:4, mean), c(1,2,4), mean)
all(eta.0 == eta.3[1:2,,])

d.3 <- {abind("ecmwf" = se.covariances(offset.forecast(ecmwf[,,,,])) / (dim(ecmwf)[[5]] - 1) + sig.3,
           "ncep" = se.covariances(offset.forecast(ncep[,,,,])) / (dim(ncep)[[5]] - 1) + sig.3,
           "ukmo" = se.covariances(offset.forecast(ukmo[,,,,])) / (dim(ukmo)[[5]] - 1) + sig.3,
           along = 0)}
d.0 <- {abind("ecmwf" = se.covariances(offset.forecast(ecmwf[1:2,,,,])) / (dim(ecmwf)[[5]] - 1) + sig.0,
              "ncep" = se.covariances(offset.forecast(ncep[1:2,,,,])) / (dim(ncep)[[5]] - 1) + sig.0,
              "ukmo" = se.covariances(offset.forecast(ukmo[1:2,,,,])) / (dim(ukmo)[[5]] - 1) + sig.0,
              along = 0)}
all(d.0 == d.3[,1:2,1:2,,,])

prec.3 <- se.precision(d.3, lambda.3)
prec.0 <- se.precision(d.3, lambda.3, vars = c("temp.n", "temp.s"))
all(prec.0 == prec.3[1:2,1:2,,,])
prec.0[,,1,1,1]; prec.3[1:2,1:2,1,1,1]

# inversion of D matrix is where difference begins
s.3 <- array.solve(prec.3, 3:5)
s.0 <- array.solve(prec.0, 3:5)

require(reshape)

pdf("./Plots/Components of S for full model.pdf"); {
    par(mfrow = c(5,1), mar = c(2,2,1,1), oma = c(0,0,2,0))
    invisible(sapply(1:5, function(n) {
        boxplot(value ~ X1, melt(s.3[n,,,,]), ylim = c(-5,10))
        abline(0,0,col = "red")
        abline(h = c(1,-1), col = "cyan3")
        legend("topleft", legend = dimnames(s.3)[[1]][n], bty = "n")
    }))
    mtext("Boxplots of elements of S for full model", outer = T)
}; dev.off()


####################################################################################################

# WORKING AREA                                                                                  ####

require(ensembleBMA)

varb <- "temp.s"
lt <- "0"

    verifRankHist(apply(offset.forecast(ecmwf)[varb,,,lt,-1], 3, rbind),
                  c(obs[varb,,]))
    legend("topleft", legend = paste0("ECMWF - ", varb, ", LT ", lt), bty = "n")
    # tends to underforecast temperature - cold bias
    
    verifRankHist(apply(offset.forecast(ncep)[varb,,,lt,-1], 3, rbind),
                  c(obs[varb,,]))
    legend("topleft", legend = paste0("NCEP - ", varb), bty = "n")
    # underdispersive - spread too small
    
    verifRankHist(apply(offset.forecast(ukmo)[varb,,,lt,-1], 3, rbind),
                  c(obs[varb,,]))
    legend("topleft", legend = paste0("UKMO - ", varb), bty = "n")
    # tends to underforecast temperature - cold bias
    
    verifRankHist(apply(offset.forecast(superensemble())[varb,,,lt,-1], 3, rbind),
                  c(obs[varb,,]))
    legend("topleft", legend = paste0("Pooled - ", varb), bty = "n")
    # presume os falls between two subgroups more often than not (two sets of cold biases)
}))

d.inv.sum.inv <- apply(array.solve(d.3,
                              c(1, 4:6)),          # solve per model/day/lt
                  c(1:2, 4:6), sum)# sum over models
            
solv <- array.solve(d.inv.sum.inv, solve.over = 3:5)
dim(d.inv.sum.inv)
