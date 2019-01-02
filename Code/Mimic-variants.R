
library("SX.weather"); library("CB.Misc")
require(ensembleMOS); require(ensembleBMA); require(mvtnorm)

# check dimensions of mimic covariances

####################################################################################################

# SUPPORT FUNCTIONS                                                                             ####

# support function - re-square an array after applying eg. cov
square.mat <- function(arr) {
    
    d <- dim(arr)

    if(sqrt(d[1]) %% 1 != 0) {
        cat("Cannot produce square array.", "\n")
        return(arr)
    } else {
        d1 <- sqrt(d[1])
        return(array(arr, dim = c(d1, d1, d[-1])))
    }
}

####################################################################################################

# FIT SINGLE VARIABLE AT A TIME                                                                 ####

#mimic.tn <- run.model(varbls = "temp.n")    # can't run - arrays collapse

# NB. everything is same until D is inverted, so can just truncate

# covariance of each forecast ensemble
C <- abind("ecmwf" = apply(offset.forecast(ecmwf)[1,,,,-1], 1:3, var),
           "ncep" = apply(offset.forecast(ncep)[1,,,,-1], 1:3, var),
           "ukmo" = apply(offset.forecast(ukmo)[1,,,,-1], 1:3, var), along = 0)
{
    all(C["ecmwf",,,] == se.covariances(offset.forecast(ecmwf))[1, 1,,,],
        C["ncep",,,] == se.covariances(offset.forecast(ncep))[1, 1,,,],
        C["ukmo",,,] == se.covariances(offset.forecast(ukmo))[1, 1,,,])
}


# covariance of mean forecast error across all 93 superensemble members
se.mean.error <- apply(forecast.errors(superensemble())[1,,,,-(1:3)], 1:3, mean)
Lam <- apply(ens.mean.error, c(1,3), var)
{
    all(Lam == se.lambda()[1,1,,])
}


# covariances of ensemble means
mean.fc <- abind(lapply(list(ecmwf, ncep, ukmo),
                        function (model) {
                            apply(offset.forecast(model)[1,,,,-1], 1:3, mean)
                        }), along = 0)
Sig <- apply(mean.fc, 2:4, var)
{
    all(Sig == se.sigma()[1,1,,,])
}


# daily average of mean FC error across all 93 superensemble members
Eta <- apply(se.mean.error, c(1,3), mean)
{
    all(Eta == apply(apply(forecast.errors(superensemble())[,,,,-(1:3)], 1:4, mean), c(1,2,4), mean)[1,,])
}


# ensemble mean forecasts
Y.bar <- abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1])[1,,,,], 1:3, mean),
               "ncep" = apply(offset.forecast(ncep[,,,,-1])[1,,,,], 1:3, mean),
               "ukmo" = apply(offset.forecast(ukmo[,,,,-1])[1,,,,], 1:3, mean),
               along = 0)
{
    all(Y.bar == abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
                       "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
                       "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
                       along = 0)[,1,,,])
}


# combined covariances of ensemble means & ensemble members
# normalised covariances per ensemble, plus covariances of ensemble means
D <- sweep(sweep(C, 1, table(attr(superensemble(), "m"))[2:4], "/"), 2:4, Sig, "+")
{
    all(D == abind("ecmwf" = se.covariances(offset.forecast(ecmwf)) / (dim(ecmwf)[[5]] - 1) + se.sigma(),
                   "ncep" = se.covariances(offset.forecast(ncep)) / (dim(ncep)[[5]] - 1) + se.sigma(),
                   "ukmo" = se.covariances(offset.forecast(ukmo)) / (dim(ukmo)[[5]] - 1) + se.sigma(),
                   along = 0)[,1,1,,,])
}


# variance at each timestep
S <- sweep(apply(D^-1, 2:4, sum)^-1, c(1,3), Lam, "+")

# mean at each timestep
Tau <- S * (1 + sweep(apply(D^-1, 2:4, sum), c(1,3), Lam, "*"))^-1 * (apply(D^-1 * Y.bar, 2:4, sum))

# repeat for temp.s only
fit.single.var <- function(varb) {
    
    # mean of each forecast ensemble
    Y.bar <- abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1])[varb,,,,], 1:3, mean),
                   "ncep" = apply(offset.forecast(ncep[,,,,-1])[varb,,,,], 1:3, mean),
                   "ukmo" = apply(offset.forecast(ukmo[,,,,-1])[varb,,,,], 1:3, mean),
                   along = 0)
    
    # covariance of each forecast ensemble
    C <- abind("ecmwf" = apply(offset.forecast(ecmwf)[varb,,,,-1], 1:3, var),
               "ncep" = apply(offset.forecast(ncep)[varb,,,,-1], 1:3, var),
               "ukmo" = apply(offset.forecast(ukmo)[varb,,,,-1], 1:3, var), 
               along = 0)

    
    # daily average of mean FC error across all 93 superensemble members
    se.mean.error <- apply(forecast.errors(superensemble())[varb,,,,-(1:3)], 1:3, mean)
    Eta <- apply(se.mean.error, c(1,3), mean)
    
    # covariance of mean forecast error across all 93 superensemble members
    Lam <- apply(se.mean.error, c(1,3), var)
    
    
    # covariances of ensemble means
    mean.fc <- abind(lapply(list(ecmwf, ncep, ukmo),
                            function (model) {
                                apply(offset.forecast(model)[varb,,,,-1], 1:3, mean)
                            }), along = 0)
    Sig <- apply(mean.fc, 2:4, var)

    
    # combined covariances of ensemble means & ensemble members
    # normalised covariances per ensemble, plus covariances of ensemble means
    D <- sweep(sweep(C, 1, table(attr(superensemble(), "m"))[2:4], "/"), 2:4, Sig, "+")

    # variance at each timestep
    S <- sweep(apply(D^-1, 2:4, sum)^-1, c(1,3), Lam, "+")
    
    # mean at each timestep
    Tau <- S * (1 + sweep(apply(D^-1, 2:4, sum), c(1,3), Lam, "*"))^-1 * (apply(D^-1 * Y.bar, 2:4, sum))
    
    list(Tau = Tau, S = S)
}

# compare to observations & original calculation: has it gone utterly mad?
model.0 <- run.model()
model.tn <- fit.single.var("temp.n")
model.ts <- fit.single.var("temp.s")

plot(obs["temp.n",,1], type = "l")
lines(model.tn$Tau[,1,1], col = "red")
lines(model.0$tau["temp.n",,1,1], col = "blue")

# metrics
res <- abind("model.0" = abind("rmse" = apply(model.0$tau["temp.n",,,], 3, fc.rmse, obs["temp.n",,]),
                               "crps" = apply(model.performance(model.0)$crps["temp.n",,,], 3, mean),
                               along = 0),
             "temp.n" = abind("rmse" = apply(model.tn$Tau, 3, fc.rmse, obs[1,,]),
                              "crps" = apply(array(verification::crps(rep(obs["temp.n",,], 15), 
                                                                      cbind(c(model.tn$Tau), c(model.tn$S)))$crps,
                                                   dim(model.tn$Tau)), 3, mean),
                              along = 0),
             "temp.s" = abind("rmse" = apply(model.ts$Tau, 3, fc.rmse, obs[2,,]),
                              "crps" = apply(array(verification::crps(rep(obs["temp.s",,], 15), 
                                                                      cbind(c(model.ts$Tau), c(model.ts$S)))$crps,
                                                   dim(model.ts$Tau)), 3, mean),
                              along = 0),
             along = 0)


matplot(t(apply(res, 3, cbind)), type = "l", col = rep(c("black", "blue", "red"), 2), lty = rep(c(1, 2), each = 3), ylim = range(0, res))
legend("topleft", lty = rep(c(1, 2), each = 3), col = rep(c("black", "blue", "red"), 2), bty = "n",
       legend = c("RMSE - model 0", "RMSE - temp.n only", "RMSE - temp.s only", 
                  "CRPS - model 0", "CRPS - temp.n only", "CRPS - temp.s only"))

# fitting a single-variable model definitely not an improvement.
# BUT useful to fully understand how things fit together (easier to visualise in 1D)

####################################################################################################

# LAMBDA & ETA PER DAY, NOT AGGREGATED OVER ALL YEARS                                           ####

# lambda & eta still calculated over all 93 superensemble members
# daily average of mean FC error across all 93 superensemble members
# not convinced on weighting here... currently using superensemble centre,
# suspect this should be ensemble average instead?
# Either way, Eta & Lambda should be calculated over same thing.
# (Here, both are 90·7·15 for all variables)
Eta <- apply(forecast.errors(superensemble())[,,,,-(1:3)], c(1:4), mean)

# covariance of mean forecast error across all 93 superensemble members
Lam <- square.mat(apply(aperm(forecast.errors(superensemble())[,,,,-(1:3)], c(5,1:4)), 3:5, cov))

# no change to remaining calculations
{
    # mean of each forecast ensemble (per day / year)
    Y.bar <- abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
                   "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
                   "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
                   along = 0)
    
    # covariance of each forecast ensemble
    C <- abind("ecmwf" = square.mat(apply(aperm(offset.forecast(ecmwf)[,,,,-1], c(5,1:4)), 3:5, cov)),
               "ncep" = square.mat(apply(aperm(offset.forecast(ncep)[,,,,-1], c(5,1:4)), 3:5, cov)),
               "ukmo" = square.mat(apply(aperm(offset.forecast(ukmo)[,,,,-1], c(5,1:4)), 3:5, cov)), 
               along = 0)
    
    # covariances of ensemble means
    Sig <- square.mat(apply(abind(lapply(list(ecmwf, ncep, ukmo),
                                         function (model) {
                                             apply(offset.forecast(model)[,,,,-1], 1:4, mean)
                                         }), along = 0), 3:5, cov))
    
    # combined covariances of ensemble means & ensemble members
    # (normalised covariances per ensemble, plus covariances of ensemble means)
    D <- sweep(sweep(C, 1, table(attr(superensemble(), "m"))[2:4], "/"), 2:6, Sig, "+")
    
    # posterior covariance matrix at each timestep
    S <- array.solve(apply(array.solve(D, c(1, 4:6)), c(1:2, 4:6), sum), c(3:5)) + Lam
    
    # elements of posterior mean
    t1 <- array.solve(sweep(square.mat(apply(abind(apply(array.solve(D, c(1, 4:6)), c(1:2, 4:6), sum),
                                                   Lam, along = 2), 
                                             3:5, function(arr) arr[,1:5] %*% arr[,-(1:5)])),
                            1:2, diag(5), "+"), solve.over = 3:5)

    t2 <- apply(apply(abind("Y" = aperm(Y.bar, c(2,1,3:5)), array.solve(D, c(1, 4:6)), along = 1), 3:6,
                      function(arr) arr[-1,] %*% arr[1,]), c(1,3:5), sum)
    
    t1.t2 <- apply(abind(t2, t1, along = 1), 3:5, function(arr) arr[-1,] %*% arr[1,])
    
    # posterior mean at each timestep
    Tau <- apply(abind(t1.t2, S, along = 1), 3:5, function(arr) arr[1,] %*% arr[-1,]) - Eta
    dimnames(Tau) <- dimnames(ecmwf[,16:105, , ,1])
}

plot(obs[1,,1], type = "l")
lines(Tau[1,,1,1], col = "red")
lines(model.0$tau[1,,1,1], col = "blue")

# start by including only temp variables (better to compare with plots already made) - 
# should re-run comparison to principal components as well later though

# sort out dimnames if necessary
{
    dn <- append(dimnames(Eta), list(dimnames(Eta)[[1]]), after = 1)
    
    dimnames(D) <- append(dn, list(dimnames(D)[[1]]), after = 0)
    dimnames(S) <- dimnames(Lam) <- dimnames(Sig) <- dn
}

####################################################################################################

# COMPARE PERFORMANCE (DAILY LAMBDA/ETA)                                                        ####

# ideally would check predictive performance - do we have another year's obs & fc available?
# in which case, use (say) December to train, Jan-Feb to predict

model.0 <- run.model()                   # Sichun's original model with all PCs retained - baseline
model.93 <- list(tau = Tau, s = S)       # revised model, fitted to all available data

saveRDS(model.93, "./Models/lambda93.rds")

model.0.res <- model.performance(model.0)
model.93.res <- model.performance(model.93)

model.crps <- abind("m.0" = apply(model.0.res$crps, c(1,4), mean),
                    "m.93" = apply(model.93.res$crps, c(1,4), mean), rev.along = 0)

model.rmse.spread <- abind("rmse.m0" = model.0.res$rmse,
                           "rmse.m93" = model.93.res$rmse,
                           "spread.m0" = model.0.res$spread,
                           "spread.m93" = model.93.res$spread, rev.along = 0)

par(mfrow = c(2,2), mar = c(2,2,3,1))
matplot(model.93rps["temp.n",,], type = "l", lty = 1, ylim = range(0, model.93rps), xlab = "", ylab = "",
        main = "CRPS - temp.n")
legend("topleft", bty = "n", lty = 1, col = c("black", "red"), legend = c("Original model", "Daily-Lambda model"))
matplot(model.93rps["temp.s",,], type = "l", lty = 1, ylim = range(0, model.93rps), xlab = "", ylab = "",
        main = "CRPS - temp.s")
legend("topleft", bty = "n", lty = 1, col = c("black", "red"), legend = c("Original model", "Daily-Lambda model"))


matplot(model.rmse.spread["temp.n",,], type = "l", lty = c(1,1,2,2), col = rep(c("black", "red"), 2), 
        ylim = range(0, model.rmse.spread), xlab = "", ylab = "", main = "RMSE / spread - temp.n")
legend("topleft", bty = "n", lty = c(1,1,2,2), col = rep(c("black", "red"), 2), 
       legend = c("RMSE - Original", "RMSE - Daily-Lambda",
                  "Spread - Original", "Spread - Daily-lambda"))

matplot(model.rmse.spread["temp.s",,], type = "l", lty = c(1,1,2,2), col = rep(c("black", "red"), 2),
        ylim = range(0, model.rmse.spread), xlab = "", ylab = "", main = "RMSE / spread - temp.s")
legend("topleft", bty = "n", lty = c(1,1,2,2), col = rep(c("black", "red"), 2), 
       legend = c("RMSE - Original", "RMSE - Daily-Lambda",
                  "Spread - Original", "Spread - Daily-lambda"))

lt <- "5"
# plot PIT histograms
{
    pit.0 <- cbind(c(pnorm(obs[1,,], model.0$tau[1,,,lt], model.0$s[1,1,,,lt])),
                   c(pnorm(obs[2,,], model.0$tau[2,,,lt], model.0$s[2,2,,,lt])))
    
    yrng.0 <- range(0, hist(pit.0[,1], plot = F)$density, hist(pit.0[,2], plot = F)$density)
    
    pit.93 <- cbind(c(pnorm(obs[1,,], model.93$tau[1,,,lt], model.93$s[1,1,,,lt])),
                   c(pnorm(obs[2,,], model.93$tau[2,,,lt], model.93$s[2,2,,,lt])))
    
    yrng.93 <- range(0, hist(pit.93[,1], plot = F)$density, hist(pit.93[,2], plot = F)$density)
    
    #-----------------------------------------------------------------------
    # plot histograms
    par(mfrow = c(2,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    hist(pit.0, breaks = "fd", col = "skyblue", ylim = yrng.0,
         main = "overall", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.0[,1], breaks = "fd", col = "skyblue", ylim = yrng.0,
         main = "temp.n", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.0[,2], breaks = "fd", col = "skyblue", ylim = yrng.0,
         main = "temp.s", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.93, breaks = "fd", col = "coral", ylim = yrng.93,
         main = "overall", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.93[,1], breaks = "fd", col = "coral", ylim = yrng.93,
         main = "temp.n", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.93[,2], breaks = "fd", col = "coral", ylim = yrng.93,
         main = "temp.s", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    mtext(paste0("PIT histograms for SX model (top row) and revised model (bottom row) - LT ", lt), outer = T)
}

pdf("./Plots/Mimic-vars/SX-vs-m93.pdf")
dev.off()

####################################################################################################

# LAMBDA & ETA DAILY ACROSS ENSEMBLE MEANS                                                      ####

# lambda & eta now calculated over 3 ensemble means
ens.mean.error <- abind("ecmwf" = apply(forecast.errors(ecmwf[,,,,-1]), 1:4, mean),
                        "ncep" = apply(forecast.errors(ncep[,,,,-1]), 1:4, mean),
                        "ukmo" = apply(forecast.errors(ukmo[,,,,-1]), 1:4, mean),
                        rev.along = 0)

Eta <- apply(ens.mean.error, c(1:4), mean)

# covariance of mean forecast error across all 93 superensemble members
Lam <- square.mat(apply(aperm(ens.mean.error, c(5,1:4)), 3:5, cov))

# no change to remaining calculations
{
    # mean of each forecast ensemble (per day / year)
    Y.bar <- abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
                   "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
                   "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
                   along = 0)
    
    # covariance of each forecast ensemble
    C <- abind("ecmwf" = square.mat(apply(aperm(offset.forecast(ecmwf)[,,,,-1], c(5,1:4)), 3:5, cov)),
               "ncep" = square.mat(apply(aperm(offset.forecast(ncep)[,,,,-1], c(5,1:4)), 3:5, cov)),
               "ukmo" = square.mat(apply(aperm(offset.forecast(ukmo)[,,,,-1], c(5,1:4)), 3:5, cov)), 
               along = 0)
    
    # covariances of ensemble means
    Sig <- square.mat(apply(abind(lapply(list(ecmwf, ncep, ukmo),
                                         function (model) {
                                             apply(offset.forecast(model)[,,,,-1], 1:4, mean)
                                         }), along = 0), 3:5, cov))
    
    # combined covariances of ensemble means & ensemble members
    # (normalised covariances per ensemble, plus covariances of ensemble means)
    D <- sweep(sweep(C, 1, table(attr(superensemble(), "m"))[2:4], "/"), 2:6, Sig, "+")
    
    # posterior covariance matrix at each timestep
    S <- array.solve(apply(array.solve(D, c(1, 4:6)), c(1:2, 4:6), sum), c(3:5)) + Lam
    
    # elements of posterior mean
    t1 <- array.solve(sweep(square.mat(apply(abind(apply(array.solve(D, c(1, 4:6)), c(1:2, 4:6), sum),
                                                   Lam, along = 2), 
                                             3:5, function(arr) arr[,1:5] %*% arr[,-(1:5)])),
                            1:2, diag(5), "+"), solve.over = 3:5)
    
    t2 <- apply(apply(abind("Y" = aperm(Y.bar, c(2,1,3:5)), array.solve(D, c(1, 4:6)), along = 1), 3:6,
                      function(arr) arr[-1,] %*% arr[1,]), c(1,3:5), sum)
    
    t1.t2 <- apply(abind(t2, t1, along = 1), 3:5, function(arr) arr[-1,] %*% arr[1,])
    
    # posterior mean at each timestep
    Tau <- apply(abind(t1.t2, S, along = 1), 3:5, function(arr) arr[1,] %*% arr[-1,]) - Eta
    dimnames(Tau) <- dimnames(ecmwf[,16:105, , ,1])
}

# assess against original model & daily lambda over all 93 superensemble members
model.3 <- list(tau = Tau, s = S)       # revised model, fitted to all available data

dimnames(model.3$s) <- dimnames(model.0$s); dimnames(model.3$tau) <- dimnames(model.0$tau)
saveRDS(model.3, "./Models/lambda3.rds")

model.3.res <- model.performance(model.3)

model.crps <- abind("m.0" = apply(model.0.res$crps, c(1,4), mean),
                    "m.93" = apply(model.93.res$crps, c(1,4), mean),
                    "m.3" = apply(model.3.res$crps, c(1,4), mean), rev.along = 0)

model.rmse.spread <- abind("rmse.m0" = model.0.res$rmse,
                           "rmse.m93" = model.93.res$rmse,
                           "rmse.m3" = model.3.res$rmse,
                           "spread.m0" = model.0.res$spread,
                           "spread.m93" = model.93.res$spread,
                           "spread.m3" = model.3.res$spread, rev.along = 0)

par(mfrow = c(2,2), mar = c(2,2,3,1))
matplot(model.crps["temp.n",,], type = "l", lty = 1, ylim = range(0, model.93rps), xlab = "", ylab = "",
        main = "CRPS - temp.n")
legend("topleft", bty = "n", lty = 1, col = c("black", "red", "green3"), 
       legend = c("Original model", "Daily-Lambda - SE", "Daily-Lambda - ens. mean"))
matplot(model.crps["temp.s",,], type = "l", lty = 1, ylim = range(0, model.93rps), xlab = "", ylab = "",
        main = "CRPS - temp.s")
legend("topleft", bty = "n", lty = 1, col = c("black", "red", "green3"), 
       legend = c("Original model", "Daily-Lambda - SE", "Daily-Lambda - ens. mean"))


matplot(model.rmse.spread["temp.n",,], type = "l", lty = rep(c(1,2), each = 3), col = rep(c("black", "red", "green3"), 2),
        ylim = range(0, model.rmse.spread), xlab = "", ylab = "", main = "RMSE / spread - temp.n")
legend("topleft", bty = "n", lty = rep(c(1,2), each = 3), col = rep(c("black", "red", "green3"), 2), 
       legend = c("RMSE - Original", "RMSE - Daily-Lambda", "RMSE - Daily-Lambda - ens. mean",
                  "Spread - Original", "Spread - Daily-lambda - SE", "Spread - Daily-Lambda - ens. mean"))

matplot(model.rmse.spread["temp.s",,], type = "l", lty = rep(c(1,2), each = 3), col = rep(c("black", "red", "green3"), 2),
        ylim = range(0, model.rmse.spread), xlab = "", ylab = "", main = "RMSE / spread - temp.s")
legend("topleft", bty = "n", lty = rep(c(1,2), each = 3), col = rep(c("black", "red", "green3"), 2), 
       legend = c("RMSE - Original", "RMSE - Daily-Lambda", "RMSE - Daily-Lambda - ens. mean",
                  "Spread - Original", "Spread - Daily-lambda - SE", "Spread - Daily-Lambda - ens. mean"))

lt <- "5"
# plot PIT histograms
{
    pit.0 <- cbind(c(pnorm(obs[1,,], model.0$tau[1,,,lt], model.0$s[1,1,,,lt])),
                   c(pnorm(obs[2,,], model.0$tau[2,,,lt], model.0$s[2,2,,,lt])))
    
    yrng.0 <- range(0, hist(pit.0[,1], plot = F)$density, hist(pit.0[,2], plot = F)$density)
    
    pit.93 <- cbind(c(pnorm(obs[1,,], model.93$tau[1,,,lt], model.93$s[1,1,,,lt])),
                    c(pnorm(obs[2,,], model.93$tau[2,,,lt], model.93$s[2,2,,,lt])))
    
    yrng.93 <- range(0, hist(pit.93[,1], plot = F)$density, hist(pit.93[,2], plot = F)$density)
    
    pit.3 <- cbind(c(pnorm(obs[1,,], model.3$tau[1,,,lt], model.3$s[1,1,,,lt])),
                    c(pnorm(obs[2,,], model.3$tau[2,,,lt], model.3$s[2,2,,,lt])))
    
    yrng.3 <- range(0, hist(pit.3[,1], plot = F)$density, hist(pit.3[,2], plot = F)$density)
    
    #-----------------------------------------------------------------------
    # plot histograms
    par(mfrow = c(3,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    hist(pit.0, breaks = "fd", col = "skyblue", ylim = yrng.0,
         main = "overall", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.0[,1], breaks = "fd", col = "skyblue", ylim = yrng.0,
         main = "temp.n", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.0[,2], breaks = "fd", col = "skyblue", ylim = yrng.0,
         main = "temp.s", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.93, breaks = "fd", col = "coral", ylim = yrng.93,
         main = "overall", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.93[,1], breaks = "fd", col = "coral", ylim = yrng.93,
         main = "temp.n", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.93[,2], breaks = "fd", col = "coral", ylim = yrng.93,
         main = "temp.s", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.3, breaks = "fd", col = "coral", ylim = yrng.3,
         main = "overall", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.3[,1], breaks = "fd", col = "coral", ylim = yrng.3,
         main = "temp.n", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.3[,2], breaks = "fd", col = "coral", ylim = yrng.3,
         main = "temp.s", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    mtext(paste0("PIT histograms for SX model (top row) and revised models (middle: SE lambda; bottom: ens. mean Lambda) - LT ", lt), outer = T)
}

pdf("./Plots/Mimic-vars/Lambda-using-all-obs.pdf")
dev.off()

####################################################################################################

# ETA, LAMBDA HISTORICAL                                                                        ####

# calculating Lambda, Eta per day is clearly overfitting.
# Use historical data to construct Lambda, Eta instead 
# (treat all obs as continuous string for this purpose)
# will use 25-day training period, for direct comparison with BMA / MOS

# may need to change dimensions of remaining calculations
ens.mean.fc <- apply(apply(abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
                                 "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
                                 "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
                                 along = 0), 2:5, mean), c(1, 4), rbind)

ens.mean.error <- sweep(ens.mean.fc, 1:2, apply(obs, 1, rbind), "-")

# get offset error: 25 prior forecast errors for 5 vars, 5 leadtimes, 605 forecasts
hist.mean.error <- abind(invisible(lapply(1:(630-25), 
                                              function(i) ens.mean.error[i+(0:24),,])),
                             rev.along = 0)
    
# get mean & covariance of prior forecasts
Eta <- apply(hist.mean.error, 2:4, mean)
Lam <- square.mat(apply(hist.mean.error, 3:4, cov))

# RESHAPE ENSEMBLE SUMMARIES TO CORRECT DIMENSIONS
# mean of each forecast ensemble (per day / year)
Y.bar <- abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
               "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
               "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
               along = 0)
Y.bar <- aperm(apply(Y.bar, c(1:2, 5), rbind), c(2:4, 1))[,,,26:630]

# covariance of each forecast ensemble
C <- abind("ecmwf" = square.mat(apply(aperm(offset.forecast(ecmwf)[,,,,-1], c(5,1:4)), 3:5, cov)),
           "ncep" = square.mat(apply(aperm(offset.forecast(ncep)[,,,,-1], c(5,1:4)), 3:5, cov)),
           "ukmo" = square.mat(apply(aperm(offset.forecast(ukmo)[,,,,-1], c(5,1:4)), 3:5, cov)), 
           along = 0)
C <- aperm(apply(C, c(1:3, 6), rbind), c(2:5, 1))[,,,,26:630]

# RESHAPE SIGMA TO CORRECT DIMENSIONS
# covariances of ensemble means
Sig <- square.mat(apply(abind(lapply(list(ecmwf, ncep, ukmo),
                                     function (model) {
                                         apply(offset.forecast(model)[,,,,-1], 1:4, mean)
                                     }), along = 0), 3:5, cov))
Sig <- aperm(apply(Sig, c(1:2, 5), rbind), c(2:4, 1))[,,,26:630]


# combined covariances of ensemble means & ensemble members
# (normalised covariances per ensemble, plus covariances of ensemble means)
D <- sweep(sweep(C, 1, table(attr(superensemble(), "m"))[2:4], "/"), 2:5, Sig, "+")

# posterior covariance matrix at each timestep
S <- array.solve(apply(array.solve(D, c(1, 4:5)), c(1:2, 4:5), sum), c(3:4)) + Lam

# elements of posterior mean
t1 <- array.solve(sweep(square.mat(apply(abind(apply(array.solve(D, c(1, 4:5)), c(1:2, 4:5), sum),
                                   Lam, along = 2),
                             3:4, function(arr) arr[,1:5] %*% arr[,-(1:5)])),
            1:2, diag(5), "+"), solve.over = c(3:4))
            
t2 <- apply(apply(abind("Y" = aperm(Y.bar, c(2,1,3:4)), "D" = array.solve(D, c(1, 4:5)), along = 1),
           3:5, function(arr) arr[-1,] %*% arr[1,]), c(1, 3:4), sum)

t1.t2 <- apply(abind(t2, t1, along = 1), 3:4, function(arr) arr[-1,] %*% arr[1,])

# posterior mean at each timestep
Tau <- apply(abind(t1.t2, S, along = 1), 3:4, function(arr) arr[1,] %*% arr[-1,]) - Eta

dimnames(Tau)[[1]] <- dimnames(ecmwf)[[1]]

model.hist <- list(tau = Tau, s = S)
saveRDS(model.hist, "./Models/lambda-hist.rds")

####################################################################################################

# LAMBDA, ETA OVER 25-DAY TRAINING PERIOD: PERFORMANCE                                          ####

m.hist <- readRDS("./Models/lambda-hist.rds")
m.0 <- readRDS("./Models/fitted-3pc.rds")

m.bma <- readRDS("./Models/ensBMA-temps.rds")
m.mos <- readRDS("./Models/ensMOS-temps.rds")
ens.data <- readRDS("./Models/Ens-data.rds")

bma.fc <- readRDS("./Models/ensBMA-fc.rds")
mos.fc <- readRDS("./Models/ensMOS-fc.rds")
obs606 <- rbind("temp.n" = c(obs["temp.n",,]), "temp.s" = c(obs["temp.s",,]))[,25:630]

res.0 <- model.performance(m.0)

#-------------------------------------------------------------------------------------------

# RMSE
rmse <- abind("bma" = sqrt(apply(sweep(bma.fc[,"0.5",,], 1:2, obs606, "-")^2, c(1,3), mean)),
              "mos" = sqrt(apply(sweep(mos.fc[,"0.5",,], 1:2, obs606, "-")^2, c(1,3), mean)),
              "sx" = res.0$rmse[1:2,],
              "hist" = sqrt(apply(sweep(m.hist$tau[1:2,,], c(1,3), obs606[,-1], "-")^2, 1:2, mean)),
              along = 0)

pdf("./Plots/Mimic-vars/25d-Lambda.pdf", height = 4); {
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    matplot(dimnames(rmse)[[3]], t(rmse[,"temp.n",]), type = "l", xlab = "", ylab = "",
            main = "temp.n", ylim = range(rmse, 0),
            lty = rep(c(2,1), each = 2), col = c("blue", "green3", "red", "black"))
    legend("bottomright", lty = rep(c(2,1), each = 2), col = c("blue", "green3", "red", "black"), bty = "n",
           legend = c("BMA", "MOS", expression(paste("Original ", Lambda)), expression(paste("Hist. ", Lambda))))
    
    matplot(dimnames(rmse)[[3]], t(rmse[,"temp.s",]), type = "l", xlab = "", ylab = "",
            main = "temp.s", ylim = range(rmse, 0),
            lty = rep(c(2,1), each = 2), col = c("blue", "green3", "red", "black"))
    legend("bottomright", lty = rep(c(2,1), each = 2), col = c("blue", "green3", "red", "black"), bty = "n",
           legend = c("BMA", "MOS", expression(paste("Original ", Lambda)), expression(paste("Hist. ", Lambda))))
    mtext("RMSE", outer = T)
}; dev.off()

#-------------------------------------------------------------------------------------------

# CRPS

crps <- array(dim = dim(rmse), dimnames = dimnames(rmse))

invisible(sapply(1:15, function(lt) {
    
    bma.crps <- crps(m.bma[[lt]], ens.data[[lt]])[,"BMA"]
    mos.crps <- crps(m.mos[[lt]], ens.data[[lt]])
    
    crps["bma", 1, lt] <<- mean(bma.crps[1:606])
    crps["bma", 2, lt] <<- mean(bma.crps[607:1212])
    
    crps["mos", 1, lt] <<- mean(mos.crps[1:606])
    crps["mos", 2, lt] <<- mean(mos.crps[607:1212])
    
    invisible(sapply(1:2, function(varb) {
        crps["hist", varb, lt] <<- verification::crps(obs = obs606[varb,-1],
                                                     pred = cbind(tau = m.hist$tau[varb,lt,],
                                                                  s = sqrt(m.hist$s[varb,varb,lt,])))$CRPS
        
        crps["sx", varb, lt] <<- verification::crps(obs = c(obs[varb,,]),
                                   pred = cbind(tau = c(m.0$tau[varb,,,lt]),
                                                s = c(sqrt(m.0$s[varb,varb,,,lt]))))$CRPS
    }))
}))

# CRPS plots
{
    par(mfrow = c(1,2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    matplot(dimnames(crps)[[3]], t(crps[,"temp.n",]), type = "l", xlab = "", ylab = "",
            main = "temp.n", ylim = range(crps, 0),
            lty = rep(c(2,1), each = 2), col = c("blue", "green3", "red", "black"))
    legend("bottomright", lty = rep(c(2,1), each = 2), col = c("blue", "green3", "red", "black"), bty = "n",
           legend = c("BMA", "MOS", expression(paste("Original ", Lambda)), expression(paste("Hist. ", Lambda))))
    
    matplot(dimnames(crps)[[3]], t(crps[,"temp.s",]), type = "l", xlab = "", ylab = "",
            main = "temp.s", ylim = range(crps, 0),
            lty = rep(c(2,1), each = 2), col = c("blue", "green3", "red", "black"))
    legend("bottomright", lty = rep(c(2,1), each = 2), col = c("blue", "green3", "red", "black"), bty = "n",
           legend = c("BMA", "MOS", expression(paste("Original ", Lambda)), expression(paste("Hist. ", Lambda))))
    mtext("CRPS", outer = T)
}

#-------------------------------------------------------------------------------------------

# PIT

lt <- 5



plot.pit.hist <- function(model, breaks = "fd", col = "skyblue", 
                          title = paste0("PIT histogram - ", model, ", LT ", lt-1), ...) {
    
    pit.model <- pit.data[[model]]
    yrng <- range(0, hist(pit.model[,1], plot = F)$density, hist(pit.model[,2], plot = F)$density)
    
    par(mfrow = c(1,3), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    hist(pit.model, breaks = breaks, col = col, ylim = yrng,
         main = "overall", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.model[,1], breaks = breaks, col = col, ylim = yrng,
         main = "temp.n", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    hist(pit.model[,2], breaks = breaks, col = col, ylim = yrng,
         main = "temp.s", xlab = "", ylab = "", prob = T)
    abline(1,0, lty = 2)
    
    mtext(title, outer = T)
}

lt <- 15

# re-run PIT data for new leadtime
{
    bma.pit <- diag(cdf(m.bma[[lt]], ens.data[[lt]], t(obs606)))
    mos.pit <- diag(cdf(m.mos[[lt]], ens.data[[lt]], t(obs606)))
    
    pit.data <- list("BMA" = cbind(bma.pit[1:606], bma.pit[607:1212]),
                     "MOS" = cbind(mos.pit[1:606], mos.pit[607:1212]),
                     "Original" = cbind(c(pnorm(obs[1,,], m.0$tau[1,,,lt], sqrt(m.0$s[1,1,,,lt]))),
                                        c(pnorm(obs[2,,], m.0$tau[2,,,lt], sqrt(m.0$s[2,2,,,lt])))),
                     "Historic" = cbind(pnorm(obs606[1,-1], m.hist$tau[1,lt,], sqrt(m.hist$s[1,1,lt,])),
                                        pnorm(obs606[2,-1], m.hist$tau[2,lt,], sqrt(m.hist$s[2,2,lt,]))))
}


pdf(paste0("./Plots/Mimic-vars/25d-Lambda-PIT-lt", lt-1, ".pdf"), height = 3); {
    plot.pit.hist("BMA")
    plot.pit.hist("MOS")
    plot.pit.hist("Original")
    plot.pit.hist("Historic")
}; dev.off()

####################################################################################################

# IDENTIFYING ANALOGUES TO THE CURRENT FORECAST                                                 ####

fc.all <- abind("ecmwf" = apply(offset.forecast(ecmwf)[1:2,,,,-1], 1:4, mean),
            "ncep" = apply(offset.forecast(ncep)[1:2,,,,-1], 1:4, mean),
            "ukmo" = apply(offset.forecast(ukmo)[1:2,,,,-1], 1:4, mean),
            along = 0)


# ideally, search for analogues that are similar in terms of current forecast & preceding weather
# (may have to reserve preceding weather for prior on alpha/Gamma?)

# consider analogues to day 5, year 1 only - longest leadtime
fc <- fc.all[,,5,1,"14"]

# look for analogues to each model's forecast, then average, or find analogue to average?
# try individual analogues first (mainly to understand similarities & code)


# Mahalanobis distance to single nearest forecast
