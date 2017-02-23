
library("SX.weather")

# 'Partial model' study: fit mimic to two of three ensembles
# to check sensitivity to model composition

# 'Perfect model' study, as suggested by Chris Brierley:
# use one model as 'truth' and others to predict results.

# uninformative prior

# load original mimic for comparison
mimic <- readRDS("../Models/an25-mimic.rds")

####################################################################################################

# COMMON SETUP: GET TRAINING SET                                                                ####

# use same training set for all models (just for speed)
an.indx <- readRDS("../Models/Analogue-indices-25.rds")

# extract ensemble means & verifying observations for each analogue
em <- abind("ecmwf" = apply(offset.forecast(ecmwf)[,,,,-1], 1:4, mean),
            "ncep" = apply(offset.forecast(ncep)[,,,,-1], 1:4, mean),
            "ukmo" = apply(offset.forecast(ukmo)[,,,,-1], 1:4, mean),
            along = 0)

# training data: ensemble mean forecast & observation for each analogue for each day
tr.dat <- abind(sapply(1:15, function(lt) {
    aaply(an.indx, 1:2, function(a) {
        abind("o" = apply(obs, 1, "[", a[lt,,]), 
              apply(em[,,,,lt], 1:2, "[", a[lt,,]),
              along = 2)
    })}, simplify = F), along = 2.5)


####################################################################################################

# FIT PARTIAL MODELS                                                                            ####

# training error is calculated against verifying observation
tr.error <- sweep(tr.dat[,,,,-1,], c(1:4, 6), tr.dat[,,,,1,], "-")

# fit partial models
{
    # omit ECMWF ensemble
    part.mimic.ec <- aaply(abind("nc" = offset.forecast(ncep)[,,,,-1],
                                 "mo" = offset.forecast(ukmo)[,,,,-1],
                                 "tr" = apply(tr.error[,,,,-1,], c(6, 1:4), mean),
                                 along = 5),
                           2:4, function(arr) {
                               fit.mimic(fc = list("ncep" = arr[,1:20],
                                                   "ukmo" = arr[,21:43]),
                                         tr = arr[,44:68])})
    
    # omit NCEP ensemble
    part.mimic.nc <- aaply(abind("ec" = offset.forecast(ecmwf)[,,,,-1],
                                 "mo" = offset.forecast(ukmo)[,,,,-1],
                                 "tr" = apply(tr.error[,,,,-2,], c(6, 1:4), mean),
                                 along = 5),
                           2:4, function(arr) {
                               fit.mimic(fc = list("ecmwf" = arr[,1:50],
                                                   "ukmo" = arr[,51:73]),
                                         tr = arr[,74:98])})
    
    # omit UKMO ensemble
    part.mimic.mo <- aaply(abind("ec" = offset.forecast(ecmwf)[,,,,-1],
                                 "nc" = offset.forecast(ncep)[,,,,-1],
                                 "tr" = apply(tr.error[,,,,-3,], c(6, 1:4), mean),
                                 along = 5),
                           2:4, function(arr) {
                               fit.mimic(fc = list("ecmwf" = arr[,1:50],
                                                   "ncep" = arr[,51:70]),
                                         tr = arr[,71:95])})
}

# RMSE
part.rmse <- sqrt(apply(sweep(abind("All three included" = mimic,
                         "ECMWF omitted" = part.mimic.ec,
                         "NCEP omitted" = part.mimic.nc,
                         "UKMO omitted" = part.mimic.mo,
                         along = 0)[,,,,1:2,"Tau"], c(5,2:3), obs[1:2,,], "-")^2, c(5,4,1), mean, na.rm = T))

invisible(sapply(dimnames(part.rmse)[[1]], function(varb) {
    matplot(part.rmse[varb,,], type = "l", lty = c(1,2,2,2), main = varb, ylim = range(0, part.rmse),
            col = c("black", "steelblue", "red3", "green3"), xlab = "", ylab = "")
    legend("bottomright", col = c("black", "steelblue", "red3", "green3"), lty = c(1,2,2,2),
           legend = dimnames(part.rmse)[[3]], bty = "n")
}))

####################################################################################################

# FIT PERFECT MODELS                                                                            ####

# create training data: calculate against target model
p.error.ec <- sweep(tr.dat[,,,,c(3:4),], c(1:4, 6), tr.dat[,,,,2,], "-")
p.error.nc <- sweep(tr.dat[,,,,c(2,4),], c(1:4, 6), tr.dat[,,,,3,], "-")
p.error.mo <- sweep(tr.dat[,,,,c(2:3),], c(1:4, 6), tr.dat[,,,,4,], "-")

# fit models
{
    # target ECMWF
    perf.mimic.ec <- aaply(abind("nc" = offset.forecast(ncep)[,,,,-1],
                                 "mo" = offset.forecast(ukmo)[,,,,-1],
                                 "tr" = apply(p.error.ec, c(6, 1:4), mean),
                                 along = 5),
                           2:4, function(arr) {
                               fit.mimic(fc = list("ncep" = arr[,1:20],
                                                   "ukmo" = arr[,21:43]),
                                         tr = arr[,44:68])})
    
    # target NCEP
    perf.mimic.nc <- aaply(abind("ec" = offset.forecast(ecmwf)[,,,,-1],
                                 "mo" = offset.forecast(ukmo)[,,,,-1],
                                 "tr" = apply(p.error.nc, c(6, 1:4), mean),
                                 along = 5),
                           2:4, function(arr) {
                               fit.mimic(fc = list("ecmwf" = arr[,1:50],
                                                   "ukmo" = arr[,51:73]),
                                         tr = arr[,74:98])})
    
    # target MO
    perf.mimic.mo <- aaply(abind("ec" = offset.forecast(ecmwf)[,,,,-1],
                                 "nc" = offset.forecast(ncep)[,,,,-1],
                                 "tr" = apply(p.error.mo, c(6, 1:4), mean),
                                 along = 5),
                           2:4, function(arr) {
                               fit.mimic(fc = list("ecmwf" = arr[,1:50],
                                                   "ncep" = arr[,51:70]),
                                         tr = arr[,71:95])})
}

# RMSE
pm.rmse <- sqrt(apply(abind("Org. model" = sweep(mimic[,,,,"Tau"], c(4,1:2), obs, "-"),
                 "Target ECMWF" = sweep(perf.mimic.ec[,,,,"Tau"], c(4,1:3), em["ecmwf",,,,], "-"),
                 "Target NCEP" = sweep(perf.mimic.nc[,,,,"Tau"], c(4,1:3), em["ncep",,,,], "-"),
                 "Target UKMO" = sweep(perf.mimic.mo[,,,,"Tau"], c(4,1:3), em["ukmo",,,,], "-"),
                 along = 0)^2, c(5,4,1), mean, na.rm = T))[1:2,,]

pdf("../Plots/Perfect-model-study.pdf", width = 7, height = 4); {
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(dimnames(pm.rmse)[[1]], function(varb) {
        matplot(0:14, pm.rmse[varb,,], type = "l", lty = c(1,2,2,2), main = varb, ylim = range(0, pm.rmse),
                col = c("black", "steelblue", "red3", "green3"), xlab = "", ylab = "")
        legend("bottomright", col = c("black", "steelblue", "red3", "green3"), lty = c(1,2,2,2),
               legend = dimnames(pm.rmse)[[3]], bty = "n", cex = 0.7)
        abline(v = 7, col = "cyan3", lty = 3)
    }))
    mtext("'Perfect model' study of each ensemble forecast", outer = T)
}; dev.off()

                                   

####################################################################################################

# PREDICTING OTHER ENSEMBLE FORECASTS AS WELL AS TRUTH                                          ####

# quick check: what is RMSE of each pair of ensembles, vs observation?

rms.ind <- sqrt(apply(abind("ecmwf : ncep" = em["ecmwf",,,,] - em["ncep",,,,],
                            "ecmwf : ukmo" = em["ecmwf",,,,] - em["ukmo",,,,],
                            "ncep : ukmo" = em["ncep",,,,] - em["ukmo",,,,],
                            "ecmwf : obs" = sweep(em["ecmwf",,,,], 1:3, obs, "-"),
                            "ncep : obs" = sweep(em["ncep",,,,], 1:3, obs, "-"),
                            "ukmo : obs" = sweep(em["ukmo",,,,], 1:3, obs, "-"),
                            along = 0)[,1:2,,,]^2, c(5,1:2), mean, na.rm = T))


invisible(sapply(dimnames(rms.ind)[[3]], function(varb) {
    matplot(0:14, rms.ind[,,varb], type = "l", lty = rep(c(2,1), each = 3),
            col = c("orange","magenta3","green3","red3", "gold", "blue"),
            xlab = "", ylab = "", main = varb, ylim = range(0, rms.ind))
    legend("topleft", cex = 0.7, legend = dimnames(rms.ind)[[2]], bty = "n",
           col = c("orange","magenta3","green3","red3", "gold", "blue"), lty = rep(c(2,1), each = 3))
}))
mtext("RMS of individual ensemble means as predictors of each other & verifying observation", outer = T)
