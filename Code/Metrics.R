
library("SX.weather"); library("CB.Misc")
require(ensembleMOS); require(ensembleBMA); require(mvtnorm)

bma.fc <- readRDS("./Data/ensBMA-fc.rds")
mos.fc <- readRDS("./Data/ensMOS-fc.rds")
fitted.bayes <- readRDS("./Data/fitted-0pc.rds")
obs606 <- rbind("temp.n" = c(obs["temp.n",,]), "temp.s" = c(obs["temp.s",,]))[,25:630]

fitted.bma <- readRDS("./Data/ensBMA-temps.rds")
fitted.mos <- readRDS("./Data/ensMOS-temps.rds")
ens.data <- readRDS("./Data/Ens-data.rds")

# function to quickly produce plots of metrics across leadtimes
comp.plot <- function(res, main = "") {
    
    matplot(res, type = "l", lty = 1, lwd = c(1,1,2), col = c("blue", "red", "black"),
            xlab = "", ylab = "", main = main)
    legend("bottomright", lty = 1, lwd = c(1,1,2), col = c("blue", "red", "black"),
           legend = c("BMA", "MOS", "Bayesian"), bty = "n")
}

####################################################################################################

# FEW EXPLORATORY PLOTS                                                                         ####
d <- 30

par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
invisible(sapply(formatC(14:0), function(lt) {
    plot(obs["temp.n",25 + (0:d),1], type = "l", lwd = 2, main = "temp.n", xlab = "", ylab = "")
    lines(bma.fc["temp.n",2,1 + (0:d),lt], col = "red")
    lines(mos.fc["temp.n",2,1 + (0:d),lt], col = "blue")
    lines(fitted.bayes$tau["temp.n",25 + (0:d),1,lt], col = "cyan2", lwd = 2)
    
    plot(obs["temp.s",25 + (0:d),1], type = "l", lwd = 2, main = "temp.s", xlab = "", ylab = "")
    lines(bma.fc["temp.s",2,1 + (0:d),lt], col = "red")
    lines(mos.fc["temp.s",2,1 + (0:d),lt], col = "blue")
    lines(fitted.bayes$tau["temp.s",25 + (0:d),1,lt], col = "cyan2", lwd = 2)
    
    mtext(paste0("First ", d, " days' forecasts at ",lt," days' lead time"), outer = T)
}))

####################################################################################################

# MEAN ABSOLUTE ERROR                                                                           ####

mae.comp <- cbind("bma" = apply(bma.fc[,"0.5",,], 3, fc.mae, obs606),
                  "mos" = apply(mos.fc[,"0.5",,], 3, fc.mae, obs606),
                  "post" = apply(fitted.bayes$tau, 4, fc.mae, obs[1:2,,]))
comp.plot(mae.comp, main = "MAE (overall)")


# check results against built-in functions
{
    bma.mae.auto <- invisible(sapply(0:14, function(lt) {
        lt <- toString(lt)
        MAE(fitted.bma[[lt]], ens.data[[lt]])["BMA"]
    }))
    all(mae.comp[,"bma"] == bma.mae.auto)
}

# More useful to compare results per variable:
comp.plot(cbind("bma" = apply(bma.fc["temp.n","0.5",,], 2, fc.mae, obs606["temp.n",]),
                "mos" = apply(mos.fc["temp.n","0.5",,], 2, fc.mae, obs606["temp.n",]),
                "post" = apply(fitted.bayes$tau["temp.n",,,], 3, fc.mae, obs["temp.n",,])),
          main = "MAE - temp.n")

comp.plot(cbind("bma" = apply(bma.fc["temp.s","0.5",,], 2, fc.mae, obs606["temp.s",]),
                "mos" = apply(mos.fc["temp.s","0.5",,], 2, fc.mae, obs606["temp.s",]),
                "post" = apply(fitted.bayes$tau["temp.s",,,], 3, fc.mae, obs["temp.s",,])),
          main = "MAE - temp.s")

####################################################################################################

# RMSE                                                                                          ####

comp.plot(cbind("bma" = apply(bma.fc[,"0.5",,], 3, fc.rmse, obs606),
                "mos" = apply(mos.fc[,"0.5",,], 3, fc.rmse, obs606),
                "post" = apply(fitted.bayes$tau, 4, fc.rmse, obs[1:2,,])),
          main = "RMSE (overall)")

# check results against built-in functions
{
    rmse.bayes <- model.performance(fitted.bayes)$rmse
    
    all(rbind(apply(fitted.bayes$tau["temp.n",,,], 3, fc.rmse, obs["temp.n",,]),
              apply(fitted.bayes$tau["temp.s",,,], 3, fc.rmse, obs["temp.s",,]))
        == model.performance(fitted.bayes)$rmse)
        }

# plot results split by region
{
    comp.plot(cbind("bma" = apply(bma.fc["temp.n","0.5",,], 2, fc.rmse, obs606["temp.n",]),
                    "mos" = apply(mos.fc["temp.n","0.5",,], 2, fc.rmse, obs606["temp.n",]),
                    "post" = apply(fitted.bayes$tau["temp.n",,,], 3, fc.rmse, obs["temp.n",,])),
              main = "RMSE - temp.n")
    
    comp.plot(cbind("bma" = apply(bma.fc["temp.s","0.5",,], 2, fc.rmse, obs606["temp.s",]),
                    "mos" = apply(mos.fc["temp.s","0.5",,], 2, fc.rmse, obs606["temp.s",]),
                    "post" = apply(fitted.bayes$tau["temp.s",,,], 3, fc.rmse, obs["temp.s",,])),
              main = "RMSE - temp.s")
}

####################################################################################################

# SPREAD                                                                                        ####

# not available automatically for BMA. Need to figure out if this can even be computed.

# built-in functions
{
    spread.bayes <- model.performance(fitted.bayes)$spread
}

####################################################################################################

# CRPS                                                                                          ####

# automated versions available for all of these. No urgent need to rework.
# Will look at reworking formula later so that each var can be split out
# (may need to refit model on only one station otherwise)

# built-in functions
    crps.bayes <- model.performance(fitted.bayes)$crps
    crps.mos <- invisible(sapply(formatC(0:14), function(lt) {
        crps(fitted.mos[[lt]], ens.data[[lt]])
    }))
    crps.bma <- invisible(sapply(formatC(0:14), function(lt) {
        crps(fitted.bma[[lt]], ens.data[[lt]])[,"BMA"]
    }))
   
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    plot(apply(crps.mos[1:606,], 2, mean), type = "l", col = "red", xlab = "", ylab = "", 
         main = "temp.n", ylim = c(0,1.6))
    lines(apply(crps.bma[1:606,], 2, mean), col = "blue")
    lines(apply(crps.bayes[1,,,], 3, mean), lwd = 2)
    legend("bottomright", lty = 1, lwd = c(1,1,2), col = c("blue", "red", "black"),
           legend = c("BMA", "MOS", "Bayesian"), bty = "n")
    
    plot(apply(crps.mos[607:1212,], 2, mean), type = "l", col = "red", xlab = "", ylab = "",
         main = "temp.s", ylim = c(0,1.6))
    lines(apply(crps.bma[607:1212,], 2, mean), col = "blue")
    lines(apply(crps.bayes[2,,,], 3, mean), lwd = 2)
    legend("bottomright", lty = 1, lwd = c(1,1,2), col = c("blue", "red", "black"),
           legend = c("BMA", "MOS", "Bayesian"), bty = "n")

    mtext("CRPS for each model at each leadtime", outer = T)

####################################################################################################

# HISTOGRAM (RANK / PIT)                                                                        ####

pit.hist <- function(model, lt = "0") {
    
    lt <- toString(lt)
    
    pit.data <- switch(tolower(model),
                       "bma" = matrix(diag(cdf(fitted.bma[[lt]], ens.data[[lt]], t(obs606))), ncol = 2),
                       "mos" = matrix(diag(cdf(fitted.mos[[lt]], ens.data[[lt]], t(obs606))), ncol = 2),
                       "mimic" = cbind(c(pnorm(obs[1,,], fitted.bayes$tau[1,,,lt], fitted.bayes$s[1,1,,,lt])),
                                       c(pnorm(obs[2,,], fitted.bayes$tau[2,,,lt], fitted.bayes$s[2,2,,,lt]))))
    
    yrng <- range(0, hist(pit.data[1:606], plot = F)$density, hist(pit.data[607:1212], plot = F)$density)
    
    org.par <- par()
    
    par(mfrow = c(1,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    hist(pit.data, breaks = c(0:10)/10, col = "skyblue", ylim = yrng,
         main = "overall", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.data[,1], breaks = c(0:10)/10, col = "skyblue", ylim = yrng,
         main = "temp.n", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.data[,2], breaks = c(0:10)/10, col = "skyblue", ylim = yrng,
         main = "temp.s", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    mtext(paste0("PIT histograms for ", model, " at leadtime ", lt), outer = T)
    par(mfrow = org.par$mfrow, mar = org.par$mar, oma = org.par$oma)
}

pit.hist("BMA", lt = "0")
pit.hist("MOS", lt = "0")
pit.hist("mimic", lt = "0")


# check against built-in functions
{
    # how to extract full PIT without using PIT?
    cdf.pit <- cdf(fitted.bma[["14"]], ens.data[["14"]], t(obs606))
    all(diag(cdf.pit) == pit(fitted.bma[["14"]], ens.data[["14"]]))
    
    # SO can use this to extract PIT for eMOS, for which no built-in PIT function exists
    bma.pit <- diag(cdf(fitted.bma[["14"]], ens.data[["14"]], t(obs606)))
    mos.pit <- diag(cdf(fitted.mos[["14"]], ens.data[["14"]], t(obs606)))
    
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    hist(bma.pit, col = "skyblue", prob = T, main = "BMA")
    hist(pit(fitted.bma[["14"]], ens.data[["14"]]), prob = T, border = "darkred", lwd = 2, add = T)
    # yup, same. Good-oh.
}

    
####################################################################################################

# BRIER SCORE                                                                                   ####

    # NB: BMA etc are forecasting each site separately (not sure if actually independent) & scored across both.
    # MVRnorm gives p(0,0), need to find marginal prob of 0 on either axis?
    
fc.brier <- function(prob, occurence) {
    mean((prob - occurence)^2)
}

# marginal probabilities of freezing (rather than joint prob of both freezing)
bayes.f.probs <- abind("temp.n" = pnorm(0, fitted.bayes$tau["temp.n",,,], fitted.bayes$s["temp.n", "temp.n",,,]),
                       "temp.s" = pnorm(0, fitted.bayes$tau["temp.s",,,], fitted.bayes$s["temp.s", "temp.s",,,]),
                         along = 0)

ens.f.probs <- invisible(lapply(formatC(0:14), function(lt) {
    cbind("bma" = cdf(fitted.bma[[lt]], ens.data[[lt]], 0),
          "mos" = cdf(fitted.mos[[lt]], ens.data[[lt]], 0))
}))

hh <- abind(ens.f.probs, along = )

bs.bma.13 <- brierScore(fitted.bma[["13"]], ens.data[["13"]], 0)$bma
all(bs.bma.13 == brier.comp["bma", "13"])

brier.comp <- invisible(sapply(formatC(0:14), function(lt) {
    c("bma" = fc.brier(cdf(fitted.bma[[lt]], ens.data[[lt]], 0), c(obs606[1,], obs606[2,]) <= 0),
      "mos" = fc.brier(cdf(fitted.mos[[lt]], ens.data[[lt]], 0), c(obs606[1,], obs606[2,]) <= 0),
      "bayes" = fc.brier(bayes.f.probs[,,,lt], obs[1:2,,] <= 0))
}))

comp.plot(t(brier.comp), main = "Brier score")

# check against built-in function for ensemble BMA
{
    brier.verif <- sapply(formatC(0:14), function(lt) brier(c(obs606[1,], obs606[2,]) <= 0, cdf(fitted.bma[[lt]], ens.data[[lt]], 0))$bs)
    brier.bma <- sapply(formatC(0:14), function(lt) brierScore(fitted.bma[[lt]], ens.data[[lt]], 0)[,"bma"])
    brier.cb <- sapply(formatC(0:14), function(lt) fc.brier(cdf(fitted.bma[[lt]], ens.data[[lt]], 0), c(obs606[1,], obs606[2,]) <= 0))
    
    all(brier.bma == brier.cb)
    all(brier.bma == brier.verif)
    # matches to result from ensembleBMA package, not for verification package. Need to investigate why.
}

# plot results split by region
{
    comp.plot(cbind("bma" = apply(bma.fc["temp.n","0.5",,], 2, fc.brier, obs606["temp.n",]),
                    "mos" = apply(mos.fc["temp.n","0.5",,], 2, fc.brier, obs606["temp.n",]),
                    "post" = apply(fitted.bayes$tau["temp.n",,,], 3, fc.brier, obs["temp.n",,])),
              main = "temp.n")
    
    comp.plot(cbind("bma" = apply(bma.fc["temp.s","0.5",,], 2, fc.brier, obs606["temp.s",]),
                    "mos" = apply(mos.fc["temp.s","0.5",,], 2, fc.brier, obs606["temp.s",]),
                    "post" = apply(fitted.bayes$tau["temp.s",,,], 3, fc.brier, obs["temp.s",,])),
              main = "temp.s")
}
    
####################################################################################################

# Q-Q PLOTS VS OBSERVATIONS                                                                     ####

# model on X, theoretical on Y

QQ <- function(lt = "0") {
    
    org.par <- par()
    
    par(mfrow = c(1,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    rng <- range(obs[1:2,,], bma.fc[,"0.5",,lt], mos.fc[,"0.5",,lt], fitted.bayes$tau[,,,lt])
    
    plot(sort(obs606), sort(bma.fc[,"0.5",,lt]), col = adjustcolor("blue", alpha = 0.5), pch = 20,
         main = "All forecasts", xlab = "Observed", ylab = "Fitted", xlim = rng, ylim = rng)
    points(sort(obs606), sort(mos.fc[,"0.5",,lt]), col = adjustcolor("red", alpha = 0.5), pch = 20)
    points(sort(obs[1:2,,]), sort(fitted.bayes$tau[,,,lt]), col = adjustcolor("black", alpha = 0.5), pch = 20)
    abline(0,1,col = "cyan3", lty = 1, lwd = 1)
    legend("bottomright", pch = 20, col = c("blue", "red", "black"), legend = c("BMA", "MOS", "Mimic"), bty = "n")
    
    plot(sort(obs606["temp.n",]), sort(bma.fc["temp.n","0.5",,lt]), col = adjustcolor("blue", alpha = 0.5), pch = 20,
         main = "temp.n", xlab = "Observed", ylab = "Fitted", xlim = rng, ylim = rng)
    points(sort(obs606["temp.n",]), sort(mos.fc["temp.n","0.5",,lt]), col = adjustcolor("red", alpha = 0.5), pch = 20)
    points(sort(obs["temp.n",,]), sort(fitted.bayes$tau["temp.n",,,lt]), col = adjustcolor("black", alpha = 0.5), pch = 20)
    abline(0,1,col = "cyan3", lty = 1, lwd = 1)
    legend("bottomright", pch = 20, col = c("blue", "red", "black"), legend = c("BMA", "MOS", "Mimic"), bty = "n")
    
    plot(sort(obs606["temp.s",]), sort(bma.fc["temp.s","0.5",,lt]), col = adjustcolor("blue", alpha = 0.5), pch = 20,
         main = "temp.s", xlab = "Observed", ylab = "Fitted", xlim = rng, ylim = rng)
    points(sort(obs606["temp.s",]), sort(mos.fc["temp.s","0.5",,lt]), col = adjustcolor("red", alpha = 0.5), pch = 20)
    points(sort(obs["temp.s",,]), sort(fitted.bayes$tau["temp.s",,,lt]), col = adjustcolor("black", alpha = 0.5), pch = 20)
    abline(0,1,col = "cyan3", lty = 1, lwd = 1)
    legend("bottomright", pch = 20, col = c("blue", "red", "black"), legend = c("BMA", "MOS", "Mimic"), bty = "n")
    
    mtext(paste0("Q-Q plots for leadtime ", lt), outer = T)
}

QQ("14")

####################################################################################################

# GENERAL PONDERINGS                                                                            ####

bma.cdf <- cdf(fitted.bma[[14]], edat, c(-10:20))
plot(bma.quants[1,], type = "l")
    
# trying to work out how to recreate PIT using only forecasts
bma.quants <- quantileForecast(fitted.bma[[14]], ens.data[[14]], c(1:19)/20)
zz <- cbind("o" = c(obs606), bma.quants)
bins <- apply(zz, 1, function(rr) findInterval(rr[1], rr[-1]))

hist(bins, prob = T)

pit.hist <- function(model, lt, nbins = 20, main = "", col = "skyblue", ...) {
    
    lt <- toString(lt)
    quants <- quantileForecast(model[[lt]], ens.data[[lt]], c(1:99) / 100)
    
    ofc <- cbind("o" = c(obs606), quants)
    bins <- apply(ofc, 1, function(rr) findInterval(rr[1], rr[-1]))
    
    hist(bins, prob = T, xlab = "", ylab = "", main = main, ...)
    abline(h = 1/nbins, lty = 2)
    
    b2 <- colnames(ofc[,-1])[bins]
    
}

pit(fitted.bma[["5"]], temps)

all(round(auto.pit,2) == b2)

pit.hist(fitted.bma, "5", nbins = 10)

zz.tf <- apply(zz[1:5,], 1, function(rr) rr[-1] < rr[1])

matplot(zz[1:5,], type = "l") 
    
    

####################################################################################################

# ANIMATED PLOT TO SHOW DEVELOPMENT OF FORECAST OVER TIME                                       ####

library(animation)

d <- 50
saveGIF({
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    for(i in 15:1){
        plot(obs["temp.n", 25 + (0:d),1], type = "l", lwd = 2, main = "temp.n", xlab = "", ylab = "",
             ylim = range(obs[1:2, 25 + (0:d),1]))
        lines(bma.fc["temp.n",2,1 + (0:d),i], col = "red")
        lines(mos.fc["temp.n",2,1 + (0:d),i], col = "blue")
        lines(fitted.bayes$tau["temp.n",25 + (0:d),1,i], col = "cyan3", lwd = 2)
        
        plot(obs["temp.s",25 + (0:d),1], type = "l", lwd = 2, main = "temp.s", xlab = "", ylab = "",
             ylim = range(obs[1:2, 25 + (0:d),1]))
        lines(bma.fc["temp.s",2,1 + (0:d),i], col = "red")
        lines(mos.fc["temp.s",2,1 + (0:d),i], col = "blue")
        lines(fitted.bayes$tau["temp.s",25 + (0:d),1,i], col = "cyan3", lwd = 2)
        
        mtext(paste0("First ", d, " days' forecasts at ", i-1 ," days' lead time"), outer = T)
    }
}, interval = 1, ani.width = 1500, ani.height = 700)

####################################################################################################

# COMBINED PLOTS                                                                                ####

pdf("./Plots/Mimic-perf/MAE.pdf", height = 4, width = 7); {
    
    par(mfrow = c(1,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    comp.plot(mae.comp, main = "All forecasts")
    
    comp.plot(cbind("bma" = apply(bma.fc["temp.n","0.5",,], 2, fc.mae, obs606["temp.n",]),
                    "mos" = apply(mos.fc["temp.n","0.5",,], 2, fc.mae, obs606["temp.n",]),
                    "post" = apply(fitted.bayes$tau["temp.n",,,], 3, fc.mae, obs["temp.n",,])),
              main = "MAE - temp.n")
    
    comp.plot(cbind("bma" = apply(bma.fc["temp.s","0.5",,], 2, fc.mae, obs606["temp.s",]),
                    "mos" = apply(mos.fc["temp.s","0.5",,], 2, fc.mae, obs606["temp.s",]),
                    "post" = apply(fitted.bayes$tau["temp.s",,,], 3, fc.mae, obs["temp.s",,])),
              main = "MAE - temp.s")
    
    mtext("Mean Absolute Error", outer = T)
}; dev.off()
    
pdf("./Plots/Mimic-perf/RMSE.pdf", height = 4, width = 7); {
    
    par(mfrow = c(1,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    comp.plot(cbind("bma" = apply(bma.fc[,"0.5",,], 3, fc.rmse, obs606),
                    "mos" = apply(mos.fc[,"0.5",,], 3, fc.rmse, obs606),
                    "post" = apply(fitted.bayes$tau, 4, fc.rmse, obs[1:2,,])),
              main = "All forecasts")
    
    comp.plot(cbind("bma" = apply(bma.fc["temp.n","0.5",,], 2, fc.rmse, obs606["temp.n",]),
                    "mos" = apply(mos.fc["temp.n","0.5",,], 2, fc.rmse, obs606["temp.n",]),
                    "post" = apply(fitted.bayes$tau["temp.n",,,], 3, fc.rmse, obs["temp.n",,])),
              main = "temp.n")
    
    comp.plot(cbind("bma" = apply(bma.fc["temp.s","0.5",,], 2, fc.rmse, obs606["temp.s",]),
                    "mos" = apply(mos.fc["temp.s","0.5",,], 2, fc.rmse, obs606["temp.s",]),
                    "post" = apply(fitted.bayes$tau["temp.s",,,], 3, fc.rmse, obs["temp.s",,])),
              main = "temp.s")
    
    mtext("Root Mean Squared Error", outer = T)
}; dev.off()
 
pdf("./Plots/Mimic-perf/CRPS.pdf", height = 4, width = 7); {
    
    par(mfrow = c(1,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    plot(apply(crps.mos, 2, mean), type = "l", col = "red", xlab = "", ylab = "", 
         main = "All forecasts", ylim = c(0,1.6))
    lines(apply(crps.bma, 2, mean), col = "blue")
    lines(apply(crps.bayes, 4, mean), lwd = 2)
    legend("bottomright", lty = 1, lwd = c(1,1,2), col = c("blue", "red", "black"),
           legend = c("BMA", "MOS", "Bayesian"), bty = "n")
    
    plot(apply(crps.mos[1:606,], 2, mean), type = "l", col = "red", xlab = "", ylab = "", 
         main = "temp.n", ylim = c(0,1.6))
    lines(apply(crps.bma[1:606,], 2, mean), col = "blue")
    lines(apply(crps.bayes[1,,,], 3, mean), lwd = 2)
    legend("bottomright", lty = 1, lwd = c(1,1,2), col = c("blue", "red", "black"),
           legend = c("BMA", "MOS", "Bayesian"), bty = "n")
    
    plot(apply(crps.mos[607:1212,], 2, mean), type = "l", col = "red", xlab = "", ylab = "", 
         main = "temp.s", ylim = c(0,1.6))
    lines(apply(crps.bma[607:1212,], 2, mean), col = "blue")
    lines(apply(crps.bayes[2,,,], 3, mean), lwd = 2)
    legend("bottomright", lty = 1, lwd = c(1,1,2), col = c("blue", "red", "black"),
           legend = c("BMA", "MOS", "Bayesian"), bty = "n")
    
    mtext("Continuous Rank Probability Score", outer = T)
}; dev.off()

pdf("./Plots/Mimic-perf/Brier.pdf", height = 4, width = 7); {
    
    par(mfrow = c(1,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    comp.plot(t(brier.comp), main = "Overall")
    
    comp.plot(t(invisible(sapply(formatC(0:14), function(lt) {
        c("bma" = fc.brier(cdf(fitted.bma[[lt]], ens.data[[lt]], 0)[1:606], obs606[1,] <= 0),
          "mos" = fc.brier(cdf(fitted.mos[[lt]], ens.data[[lt]], 0)[1:606], obs606[1,] <= 0),
          "bayes" = fc.brier(bayes.f.probs[1,,,lt], obs[1,,] <= 0))
    }))), main = "temp.n")
    
    comp.plot(t(invisible(sapply(formatC(0:14), function(lt) {
        c("bma" = fc.brier(cdf(fitted.bma[[lt]], ens.data[[lt]], 0)[607:1212], obs606[2,] <= 0),
          "mos" = fc.brier(cdf(fitted.mos[[lt]], ens.data[[lt]], 0)[607:1212], obs606[2,] <= 0),
          "bayes" = fc.brier(bayes.f.probs[2,,,lt], obs[2,,] <= 0))
    }))), main = "temp.s")
    
    mtext("Brier score for probability of freezing", outer = T)
}; dev.off()

lt <- "0"
pdf(paste0("./Plots/Mimic-perf/PIT-hists-lt-", lt, ".pdf"), height = 3, width = 7); {
    
    layout(matrix(c(1,2,3), ncol = 1), widths = 1, heights = 1/3)
    pit.hist("BMA", lt = lt)
    pit.hist("MOS", lt = lt)
    pit.hist("mimic", lt = lt)
    
}; dev.off()

pdf("./Plots/Mimic-perf/Q-Q.pdf", height = 3, width = 7); {
    layout(matrix(c(1:4), ncol = 1), widths = 1, heights = 1/4)
    QQ("0"); QQ("5"); QQ("10"); QQ("14")
}; dev.off()

####################################################################################################

# MULTIVARIATE PERFORMANCE                                                                      ####

m.hist <- readRDS("./Models/lambda-hist.rds")

# Box density ordinate transform
{
    # check calculation on single instance first
    m <- m.hist$tau[,1,1]
    s <- m.hist$s[,,1,1]
    
    o <- obs[,26,1]
    
    
    1- dchisq(t(o-m) %*% solve(s) %*% (o-m), dim(m)[[1]])
    
    # now over array
    mat <- abind("x" = aperm(apply(obs, 1, c), c(2,1))[1:2,26:630],
                 "mu" = m.hist$tau[1:2,1,],
                 "sig" = m.hist$s[1:2,1:2,1,],
                 along = 1)
    
    bot <- apply(mat, 3, function(arr) {
        1- dchisq(t(arr["x",] - arr["mu",]) %*% solve(arr[-c(1:2),]) %*% (arr["x",] - arr["mu",]), dim(arr)[[2]])
    })
    
    hist(bot, breaks = "fd")
}


