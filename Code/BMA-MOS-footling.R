
# replicate output of BMA & MOS models to confirm that my representation of the models is ok

library("SX.weather"); library("ensembleBMA"); library("ensembleMOS")
setwd("~/Documents/PhD/Miniprojects/03-Sichun-paper/Code")

lt <- 5; d <- 32; y <- 5 
timestamps <- gsub("-", "", as.Date(load.data("../Data/ECMWF_europe_starttime.rda")[16:105,1:7]/ 24, "1800-01-01"))
dt <- timestamps[90*(y-1) + d]

exch <- c(rep(1, 50), rep(2, 20), rep(3, 23))
names(exch) <- c(paste0("ec", 1:50), paste0("nc", 1:20), paste0("mo", 1:23))

par(pch = 20, mar = c(2,2,3,1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Create necessary data sets                                                                    ####

    fc <- cbind(apply(apply(offset.forecast(ecmwf)[1:2,,,lt,-1], c(1,4), c), 3, c),
                apply(apply(offset.forecast(ncep)[1:2,,,lt,-1], c(1,4), c), 3, c),
                apply(apply(offset.forecast(ukmo)[1:2,,,lt,-1], c(1,4), c), 3, c))
    colnames(fc) <- names(exch)
    
    o <- c(apply(obs[1:2,,], 1, c))
    
    temp.dat <- ensembleData(forecasts = fc,
                             dates = rep(timestamps, 2),
                             observations = o,
                             forecastHour = 0,
                             station = rep(c("N", "S"), each = 630),
                             latitude = rep(c(56, 52), each = 630),
                             longitude = rep(c(-3, 0), each = 630),
                             initializationTime = "00",
                             exchangeable = exch)
    
    temp.n <- temp.dat[temp.dat$station == "N",]
    temp.s <- temp.dat[temp.dat$station == "S",]
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # extract training data for the target date
    tr <- trainingData(temp.dat, trainingDays = 25, date = dt)
    tr.n <- trainingData(temp.n, trainingDays = 25, date = dt)
    tr.s <- trainingData(temp.s, trainingDays = 25, date = dt)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # datasets using ensemble mean values only
    tm.dat <- ensembleData(forecasts = cbind("ec" = apply(fc[,1:50], 1, mean),
                                          "nc" = apply(fc[,51:70], 1, mean),
                                          "mo" = apply(fc[,71:93], 1, mean)),
                           dates = rep(timestamps, 2),
                           observations = o,
                           forecastHour = 0,
                           station = rep(c("N", "S"), each = 630),
                           latitude = rep(c(56, 52), each = 630),
                           longitude = rep(c(-3, 0), each = 630),
                           initializationTime = "00")
    
    tm.n <- tm.dat[tm.dat$station == "N",]
    tm.s <- tm.dat[tm.dat$station == "S",]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# fit models using packages                                                                     ####
    
# models fitted using individual ensemble members
    bma.fit <- ensembleBMA(temp.dat, model = "normal", trainingDays = 25, dates = dt)
    bma.n <- ensembleBMA(temp.n, model = "normal", trainingDays = 25, dates = dt)
    bma.s <- ensembleBMA(temp.s, model = "normal", trainingDays = 25, dates = dt)
    
    mos.fit <- ensembleMOS(temp.dat, model = "normal", trainingDays = 25, dates = dt)
    mos.n <- ensembleMOS(temp.n, model = "normal", trainingDays = 25, dates = dt)
    mos.s <- ensembleMOS(temp.s, model = "normal", trainingDays = 25, dates = dt)
    
    
    bma.mean.fit <- ensembleBMA(tm.dat, model = "normal", trainingDays = 25, dates = dt)
    mos.mean.fit <- ensembleMOS(tm.dat, model = "normal", trainingDays = 25, dates = dt)
    
    bma.mean.n <- ensembleBMA(tm.n, model = "normal", trainingDays = 25, dates = dt)
    bma.mean.s <- ensembleBMA(tm.s, model = "normal", trainingDays = 25, dates = dt)
    
    mos.mean.n <- ensembleMOS(tm.n, model = "normal", trainingDays = 25, dates = dt)
    mos.mean.s <- ensembleMOS(tm.s, model = "normal", trainingDays = 25, dates = dt)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # compare to results using ensemble means only
    mos.mean.fit$a; mos.fit$a
    mos.mean.fit$B; ddply(data.frame(model = sapply(rownames(mos.fit$B), substr, 1, 2), b = mos.fit$B)
                          , .(model), summarise, sum(X20120101))
    mos.mean.fit$c; mos.fit$c
    mos.mean.fit$d; mos.fit$d
    
    bma.mean.fit$biasCoefs; bma.fit$biasCoefs[, c(1,51,71),]
    bma.mean.fit$sd; bma.fit$sd
    bma.mean.fit$weights; bma.fit$weights[c(1,51,71)] * c(50, 20, 23)
    
    # MOS fits same model regardless of inclusion of individual ensemble members. 
    # BMA produces quite different fit with & without individual members.
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Plot different ensembleBMA models fitted                                                      ####
    
    tr.fc <- list("ec" = as.matrix(tr[,1:50]), "nc" = as.matrix(tr[,51:70]), "mo" = as.matrix(tr[,71:93]))
  
    tr.bma.bc <- as.matrix(sweep(sweep(tr[,1:93], 2, bma.fit$biasCoefs[2,,1], "*"), 2, bma.fit$biasCoefs[1,,1], "+"))
    
    plot(0, type = "n", xlim = range(tr[,1:93]), ylim = range(tr[,1:93]), xlab = "", ylab = "")
    
    # original ensemble
    points(tr.fc$ec[1:25,], tr.fc$ec[26:50, ], pch = 20, col = adjustcolor("steelblue", alpha = 0.1))
    points(tr.fc$nc[1:25,], tr.fc$nc[26:50, ], pch = 20, col = adjustcolor("green3", alpha = 0.1))
    points(tr.fc$mo[1:25,], tr.fc$mo[26:50, ], pch = 20, col = adjustcolor("red3", alpha = 0.1))
    
    points(tr[25,"observations"], tr[50,"observations"], pch = 4, lwd = 2)
    
    # bias-corrected ensemble
    points(tr.bma.bc[1:25,1:50], tr.bma.bc[26:50,1:50], col = adjustcolor("steelblue", alpha = 0.3))
    
    points(t(quantileForecast(bma.fit, temp.dat, quantiles = 0.5)), lwd = 2)
    
    plot(tr.fc$ec[1:25,], tr.fc$ec[26:50, ], pch = 20, col = adjustcolor("steelblue", alpha = 0.1))
    points(tr.bma.bc[1:25,1:50], tr.bma.bc[26:50,1:50], col = adjustcolor("green3", alpha = 0.3))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Marginal plots of BMA modelling....                                                           ####
    
    fc.org <- abind("temp.n" = tr[1:25,1:93], "temp.s" = tr[26:50, 1:93], along = 0)
    fc.mean <- abind("ec" = apply(fc.org[,,1:50], 1:2, mean),
                     "nc" = apply(fc.org[,,51:70], 1:2, mean),
                     "mo" = apply(fc.org[,,71:93], 1:2, mean), along = 3)
        
    bma.fit.bc <- sweep(sweep(fc.org, 3, bma.fit$biasCoefs[2,,1], "*"), 3, bma.fit$biasCoefs[1,,1], "+")
    bma.n.bc <- sweep(sweep(fc.org["temp.n",,], 2, bma.n$biasCoefs[2,,1], "*"), 2, bma.n$biasCoefs[1,,1], "+")
    bma.s.bc <- sweep(sweep(fc.org["temp.s",,], 2, bma.s$biasCoefs[2,,1], "*"), 2, bma.s$biasCoefs[1,,1], "+")

    bma.mean.bc <- sweep(sweep(fc.mean, 3, bma.mean.fit$biasCoefs[2,,1], "*"), 3, bma.mean.fit$biasCoefs[1,,1], "+")
    bma.mean.n.bc <- sweep(sweep(fc.mean[1,,], 2, bma.mean.n$biasCoefs[2,,1], "*"), 2, bma.mean.n$biasCoefs[1,,1], "+")
    bma.mean.s.bc <- sweep(sweep(fc.mean[2,,], 2, bma.mean.s$biasCoefs[2,,1], "*"), 2, bma.mean.s$biasCoefs[1,,1], "+")
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
    plot(density(fc.org[1,,], bw = 0.6), type = "l", col = "black", ylim = c(0,0.2),
         ylab = "", xlab = "Temp (N)", main = "Density before & after bias correction")
    lines(density(bma.fit.bc[1,,], bw = 0.6), col = "red3")
    lines(density(bma.n.bc, bw = 0.6), col = "blue")
    
    lines(density(fc.mean[1,,], bw = 1), col = "black", lty = 2)
    lines(density(bma.mean.bc[1,,], bw = 1), col = "red3", lty = 2)
    lines(density(bma.mean.n.bc, bw = 1), col = "blue", lty = 2)
    
    abline(v = obs[1,d,y], lwd = 2, col = "orange")
    
    legend("topleft", bty = "n", col = c("black", "red3", "blue", NA, "black"), lty = c(rep(1, 4),2),
           c("Original forecasts", "BC at both stations", "BC at N station only", NA, "Ensemble mean only"))
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
    plot(density(fc.org[2,,], bw = 0.6), type = "l", col = "black", ylim = c(0,0.2),
         ylab = "", xlab = "Temp (S)", main = "Density before & after bias correction")
    lines(density(bma.fit.bc[2,,], bw = 0.6), col = "red3")
    lines(density(bma.s.bc, bw = 0.6), col = "blue")
    
    lines(density(fc.mean[2,,], bw = 1), col = "black", lty = 2)
    lines(density(bma.mean.bc[2,,], bw = 1), col = "red3", lty = 2)
    lines(density(bma.mean.s.bc, bw = 1), col = "blue", lty = 2)
    
    abline(v = obs[2,d,y], lwd = 2, col = "orange")
    
    legend("topleft", bty = "n", col = c("black", "red3", "blue", NA, "black"), lty = c(rep(1, 4),2),
           c("Original forecasts", "BC at both stations", "BC at N station only", NA, "Ensemble mean only"))
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
    plot(fc.org[1,,], fc.org[2,,], pch = 20, col = adjustcolor("grey", alpha = 0.3), 
         xlab = "temp.n", ylab = "temp.s", main = "Multivariate BMA correction")
    points(bma.fit.bc[1,,], bma.fit.bc[2,,], pch = 20, col = adjustcolor("black", 0.6))
    abline(v = obs["temp.n",d,y], h = obs["temp.s",d,y], col = "red3")
    
    plot(fc.org[1,,], fc.org[2,,], pch = 20, col = adjustcolor("grey", alpha = 0.3), 
         xlab = "temp.n", ylab = "temp.s", main = "Univariate BMA corrections")
    points(bma.n.bc, bma.s.bc, pch = 20, col = adjustcolor("black", 0.6))
    abline(v = obs["temp.n",d,y], h = obs["temp.s",d,y], col = "red3")
    }
    {
        plot(fc.mean[1,,1], fc.mean[2,,1], pch = 3, col = adjustcolor("steelblue", alpha = 0.4),
             xlim = range(fc.mean[1,,], bma.mean.bc[1,,]), ylim = range(fc.mean[2,,], bma.mean.bc[2,,]))
        points(fc.mean[1,,2], fc.mean[2,,2], pch = 3, col = adjustcolor("green3", alpha = 0.4))
        points(fc.mean[1,,3], fc.mean[2,,3], pch = 3, col = adjustcolor("red3", alpha = 0.4))
    
        points(bma.mean.bc[1,,2], bma.mean.bc[2,,2], pch = 20, col = "green3")
        
    }
    
    # final BMA model
    xl <- c(-500:1000) / 100
    
    bma.pred <- apply(sweep(apply(sweep(sweep(tr[tr$date == dt,1:93],
                                  2, bma.fit$biasCoefs[2,,1], "*"),
                            2, bma.fit$biasCoefs[1,,1], "+"),
                      1:2, function(m) dnorm(xl, m, bma.fit$sd)),
                      3, bma.fit$weights, "*"), 1:2, sum)
    

    plot(xl,bma.pred[,1], type = "l")
    
    plot(ellipse(diag(c(sig,sig)), centre = t(bma.pred[,1])), type = "l", )
    sweep(sweep(fc.org, 3, bma.fit$biasCoefs[2,,1], "*"), 3, bma.fit$biasCoefs[1,,1], "+")
        
    fc.i <- offset.forecast(ecmwf)[1:2,d,y,lt,-1]
    
#    plot(fc.org[1,,], fc.org[2,,], pch = 20, col = adjustcolor("gold", alpha = 0.3))
#    points(bma.fit.bc[1,,], bma.fit.bc[2,,], pch = 20, col = adjustcolor("red3", alpha = 0.3))
#    
#    hist(fc.org[1,,], breaks = "fd", xlab = "temp.n", ylab = "", main = "", 
#         col = adjustcolor("skyblue", alpha = 0.4), border = "skyblue")
#    hist(bma.fit.bc[1,,], breaks = "fd", add = T, col = adjustcolor("gold", alpha = 0.6), border = "gold")
#    abline(v = obs[1,d,y])
    

    
    

    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Multivariate case: estimating coefficient over both locations                                 ####

# fit a single BMA model to examine
{
    bma.fit <- ensembleBMA(temp.dat, model = "normal", trainingDays = 25,
                           dates = dt)
    
    tr <- trainingData(temp.dat, trainingDays = 25, date = dt)
}

# forecast using fitted BMA model
quantileForecast(bma.fit, temp.dat, quantiles = 0.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Plot of BMA estimated system                                                                  ####

plot(0, type = "n", xlab = "temp.n", ylab = "temp.s" , xlim = c(-0,5), ylim = c(2,10))

points(t(offset.forecast(ecmwf)[1:2,d,y,lt,-1]), pch = 20, col = adjustcolor("steelblue", alpha = 0.2))
points(t(offset.forecast(ncep)[1:2,d,y,lt,-1]), pch = 20, col = adjustcolor("green3", alpha = 0.2))
points(t(offset.forecast(ukmo)[1:2,d,y,lt,-1]), pch = 20, col = adjustcolor("red3", alpha = 0.2))

points(t(obs[1:2,d,y]), pch = 4, lwd = 2)

points(t(quantileForecast(bma.fit, temp.dat, quantiles = 0.5)), lwd = 2)

fc.0 <- abind(lapply(list("ecmwf" = ecmwf, "ncep" = ncep, "ukmo" = ukmo), function(ffc) {
    t(apply(offset.forecast(ffc)[1:2,d,y,lt,-1], 1, mean))
}), along = 1)
    
points(t(fc.0["ecmwf",]), col = "skyblue", pch = 16)
points(t(fc.0["ncep",]), col = "green2", pch = 16)
points(t(fc.0["ukmo",]), col = "coral", pch = 16)

fc.adj <- sweep(sweep(fc.0, 1, bc[2,], "*"), 1, bc[1,], "+")
points(t(fc.adj["ecmwf",]), bg = "skyblue", pch = 21)
points(t(fc.adj["ncep",]), bg = "green2", pch = 21)
points(t(fc.adj["ukmo",]), bg = "coral", pch = 21)

w.fc <- c(sum(bma.fit$weights[1:50]), sum(bma.fit$weights[51:70]), sum(bma.fit$weights[71:93]))

apply(sweep(fc.adj, 1, w.fc, "*"), 2, sum)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Manually replicate BMA post-processing                                                        ####

# ensemble BMA starts with linear regression over each ensemble's training data
bma.coefs <- lapply(list("ec" = cbind(c(obs[1:2,d-(24:0),y]), 1, c(offset.forecast(ecmwf)[1:2,d-(24:0),y,lt,-1])),
                         "nc" = cbind(c(obs[1:2,d-(24:0),y]), 1, c(offset.forecast(ncep)[1:2,d-(24:0),y,lt,-1])),
                         "mo" = cbind(c(obs[1:2,d-(24:0),y]), 1, c(offset.forecast(ukmo)[1:2,d-(24:0),y,lt,-1]))),
                    function(dat) {
                        Y.n <- as.matrix(dat[,1])
                        X.n <- as.matrix(dat[,2:3])
                        
                        solve(t(X.n) %*% X.n) %*% t(X.n) %*% Y.n
                    })

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Functions relating to EM algorithm                                                            ####

# support function
log.likelihood <- function(sig2, w, bc, fc, o) {

    l <- lapply(1:3, function(i) {
        sweep(aaply(fc[[i]], 3, function(fij) {
            log(dnorm(o, bc[1,i] + bc[2,i] * fij, sqrt(sig2)))
        }), 1, w[[i]], "*")
    })
   sum(unlist(lapply(l, sum)))
}

log.likelihood <- function(sig2, w, bc, fc, o) {
    
    # summed over all observations in training set
    l <- array(dim = dim(fc), dimnames = dimnames(fc))
    
    invisible(sapply(1:3, function(i) {
        invisible(sapply(1:25, function(t) {
            l[i,t] <<- w[i] * dnorm(o[t], bc[i,1] + bc[i,2] * fc[i,t], sqrt(sig2))
        }))}))
    sum(log(colSums(l)))
}

EM <- function(bc, fc, o, sig2 = 4, max.runs = 1000, conv = 0.0001) {
    
    n <- 0
    log.lh <- log.likelihood(sig2, w, bc, fc, o)
    new.log.lh <- abs(log.lh) + 100
    
    # initialise all weights to be equal 
    w <- lapply(fc, function(ffc) rep(1, dim(ffc)[[3]]))
    w <- lapply(w, "/", sum(unlist(w)))
    
    sqerr <- lapply(lapply(fc, sweep, 1:2, o, "-"), "^", 2)

    while ((abs(log.lh - new.log.lh) > conv) && (n < max.runs)) {
    # E step: estimate weights
    zz <- lapply(1:3, function(i) {
        sweep(aaply(fc[[i]], 3, function(fij) {
            dnorm(o, bc[1,i] + bc[2,i] * fij, sqrt(sig2))
        }), 1, w[[i]], "*")
    })
    # normalise weights
    z.ijst <- lapply(zz, function(zu) sweep(zu, 2:3, 
                                        apply(abind(lapply(zz, apply, 2:3, sum), along = 0),
                                              2:3, sum),
                                        "/"))
    
    # M step: maximise weights
    w <- mapply(rep, unlist(lapply(z.ijst, mean)), unlist(lapply(w, length)))
    
    # M step: maximise sigma
    sig2 <- sum(unlist(lapply(1:3, function(i) {
      (sum(aperm(z.ijst[[i]], c(2,3,1)) * sqerr[[i]])) / dim(sqerr[[i]])[2]
    })))
    
    # calculate new log-likelihood for convergence check
    log.lh <- new.log.lh
    new.log.lh <- log.likelihood(sig2, w, bc, fc, o)
    n <- n + 1
    }
    return(list(sig = sqrt(sig2), w = w, log.lh = new.log.lh, iter = n))
}
 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# BMA EM estimation of weights & shared variance                                                ####

# this is not calculated correctly...

o <- obs[1:2, d-(25:1), y]
fc <- list("ecmwf" = offset.forecast(ecmwf)[1:2, d-(25:1),y,lt,-1],
            "ncep" = offset.forecast(ncep)[1:2,d-(25:1),y,lt,-1],
            "ukmo" = offset.forecast(ukmo)[1:2,d-(25:1),y,lt,-1])

# linear regression over training data
bc <- cbind("ecmwf" = find.B.hat(1:50, dat = tr),
            "ncep" = find.B.hat(51:70, dat = tr),
            "ukmo" = find.B.hat(71:93, dat = tr))

em.fit <- EM(bc, fc, o, sig2 = 1, max.runs = 100)


# functions relating to EM algorithm
{
    # support function to evaluate log-likelihood
    log.likelihood <- function(sig2, w, bc, fc, o) {
        
        # summed over all observations in training set
        l <- array(dim = dim(fc), dimnames = dimnames(fc))
        
        invisible(sapply(1:3, function(i) {
            invisible(sapply(1:25, function(t) {
                l[i,t] <<- w[i] * dnorm(o[t], bc[i,1] + bc[i,2] * fc[i,t], sqrt(sig2))
            }))}))
        sum(log(colSums(l)))
    }
    
    # function to perform maximisation
    EM <- function(bc, fc, o, sig2 = 4, w = rep(1/3, 3), max.runs = 1000, conv = 0.0001) {
        
        log.lh <- log.likelihood(sig2, w, bc, fc, o)
        new.log.lh <- abs(log.lh) + 100
        n = 0
        
        sqerr <- sweep(fc, 2, o, "-")^2
        
        # create vector to store initial values
        first.10 <- c(iter = n, sig2 = sig2, w = w, log.lh = log.lh)
        
        while ((abs(log.lh - new.log.lh) > conv) && (n < max.runs)) {
            
            # E-step: estimate weights
            z <- array(dim = dim(fc), dimnames = dimnames(fc))
            
            for(i in 1:3) {
                for (t in 1:25) {
                    z[i,t] <- w[i] * dnorm(o[t], bc[i,1] + bc[i,2] * fc[i,t], sqrt(sig2))
                }
            }
            z <- sweep(z, 2, colSums(z), "/")   # normalise
            
            # M-step: maximise weights
            w <- apply(z, 1, mean)
            
            # M-step: maximise sig2
            sig2 <- mean(z * sqerr)
            
            # calculate log-likelihoods for comparison
            log.lh <- new.log.lh
            new.log.lh <- log.likelihood(sig2, w, bc, fc, o)
            n <- n + 1
            
            # save first 10 iterations
            if (n < 11) {
                next.iter <- c(iter = n, sig = sqrt(sig2), w = w, log.lh = log.lh)
                first.10 <- rbind(first.10, next.iter)
            }
        }
        
        # Output: if model hasn't converged, show error message
        #         if it has, output the parameters & first 10 iterations
        if ((abs(log.lh - new.log.lh) > conv)) {
            cat ("Data hasn't converged after", n, "iterations; \n",
                 "Difference in log-likelihoods is", 
                 round(abs(log.lh - new.log.lh),6))
        } else {
            row.names(first.10) <- first.10[,1]
            list(sig = sqrt(sig2), w = w, log.lh = new.log.lh, iter = n, 
                 first.10 = round(first.10, 2))
        }
    }
}

# simplest case: single variable, single day
o <- obs[1,1:25,1]
fc <- abind("ecmwf" = apply(offset.forecast(ecmwf)[1,1:25,1,lt,-1], 1, mean),
            "ncep" = apply(offset.forecast(ncep)[1,1:25,1,lt,-1], 1, mean),
            "ukmo" = apply(offset.forecast(ukmo)[1,1:25,1,lt,-1], 1, mean),
            along = 0)

# linear regression over training data
bc <- rbind("ecmwf" = apply(bma.fit$biasCoefs[,1:50,], 1, mean),
            "ncep" = apply(bma.fit$biasCoefs[,51:70,], 1, mean),
            "ukmo" = apply(bma.fit$biasCoefs[,71:93,], 1, mean))

em.fit <- EM(sig2 = 1)
plot.em(em.fit, xlim = c(0.6, 6))

plot(bma.temps, temp.n, dates = "20071225", ask = F, xlim = c(0,6))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Univariate case                                                                               ####

# create data in required format - northern temps
{
    fc.n <- cbind(apply(offset.forecast(ecmwf)["temp.n",,,lt,-1], 3, c),
                  apply(offset.forecast(ncep)["temp.n",,,lt,-1], 3, c),
                  apply(offset.forecast(ukmo)["temp.n",,,lt,-1], 3, c))
    colnames(fc.n) <- c(paste0("ec", 1:50), paste0("nc", 1:20), paste0("mo", 1:23))
    
    o.n <- c(obs["temp.n",,])
    
    tn.dat <- ensembleData(forecasts = fc.n,
                           dates = timestamps,
                           observations = o.n,
                           forecastHour = 0,
                           station = rep(c("N"), each = 630),
                           latitude = rep(c(56), each = 630),
                           longitude = rep(c(-3), each = 630),
                           initializationTime = "00",
                           exchangeable = exch)
}


bma.n <- ensembleBMA(tn.dat, model = "normal", trainingDays = 25, dates = dt)
tr.n <- trainingData(tn.dat, trainingDays = 25, date = dt)

bc.n <- cbind("ecmwf" = apply(bma.n$biasCoefs[,1:50,], 1, mean),
              "ncep" = apply(bma.n$biasCoefs[,51:70,], 1, mean),
              "ukmo" = apply(bma.n$biasCoefs[,71:93,], 1, mean))

# obtain single set of coefficients by treating each forecast as independent single-variable forecast
lm.n.all <- cbind("ecmwf" = lm(obs ~ fc, data.frame(obs = tr.n[,95], fc = unlist(tr.n[,1:50])))$coef,
                  "ncep" = lm(obs ~ fc, data.frame(obs = tr.n[,95], fc = unlist(tr.n[,51:70])))$coef,
                  "ukmo" = lm(obs ~ fc, data.frame(obs = tr.n[,95], fc = unlist(tr.n[,71:93])))$coef)

all(lm.n.all == bc.n)

find.B.hat <- function(c.range, c.obs = 95) {
    X.n <- as.matrix(cbind(1, unlist(tr.n[,c.range])))
    Y.n <- as.matrix(rep(tr.n[,c.obs], nrow(X.n) / nrow(tr.n)))
    
    solve(t(X.n) %*% X.n) %*% t(X.n) %*% Y.n
}

B.hat <- cbind("ecmwf" = find.B.hat(1:50), "ncep" = find.B.hat(51:70), "ukmo" = find.B.hat(71:93))

all(round(B.hat, 9) == round(bc.n, 9))

# what if we were to take the mean of all 50 ensemble members?
lm.n.mean <- cbind("ecmwf" = lm(obs ~ fc, data.frame(obs = tr.n[,95], fc = apply(tr.n[,1:50], 1, mean)))$coef,
                   "ncep" = lm(obs ~ fc, data.frame(obs = tr.n[,95], fc = apply(tr.n[,51:70], 1, mean)))$coef,
                   "ukmo" = lm(obs ~ fc, data.frame(obs = tr.n[,95], fc = apply(tr.n[,71:93], 1, mean)))$coef)

# create data in required format - southern temps
{
    fc.s <- cbind(apply(offset.forecast(ecmwf)["temp.s",,,lt,-1], 3, c),
                  apply(offset.forecast(ncep)["temp.s",,,lt,-1], 3, c),
                  apply(offset.forecast(ukmo)["temp.s",,,lt,-1], 3, c))
    colnames(fc.s) <- c(paste0("ec", 1:50), paste0("nc", 1:20), paste0("mo", 1:23))
    
    o.s <- c(obs["temp.s",,])
    
    ts.dat <- ensembleData(forecasts = fc.s,
                           dates = timestamps,
                           observations = o.s,
                           forecastHour = 0,
                           station = rep(c("S"), each = 630),
                           latitude = rep(c(56), each = 630),
                           longitude = rep(c(-3), each = 630),
                           initializationTime = "00",
                           exchangeable = exch)
}
bma.s <- ensembleBMA(ts.dat, model = "normal", trainingDays = 25, dates = dt)
tr.s <- trainingData(ts.dat, trainingDays = 25, date = dt)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# MARGINAL PLOTS - BMA vs Bayesian vs MOS                                                       ####

plot(unlist(tr.n[,1:50]))

hist(unlist(tr.n[,1:50]), col = adjustcolor("steelblue", alpha = 0.3), border = NA, main = "", xlab = "temp.n")
hist(unlist(tr.n[,51:70]), col = adjustcolor("green3", alpha = 0.3), border = NA, add = T)
hist(unlist(tr.n[,71:93]), col = adjustcolor("coral3", alpha = 0.3), border = NA, add = T)

stripchart(round(unlist(tr.n[,1:50]), 1), pch = 20, method = "stack", ylim = c(0,10), xlim = c(-7,12),
           col = adjustcolor("steelblue", alpha = 0.3))
stripchart(round(unlist(tr.n[,51:70]), 1), pch = 20, method = "stack", add = T, 
           col = adjustcolor("green3", alpha = 0.3))
stripchart(round(unlist(tr.n[,71:93]), 1), pch = 20, method = "stack", add = T, 
           col = adjustcolor("coral3", alpha = 0.3))

plot(density(unlist(tr.n[,1:50])), col = adjustcolor("steelblue", alpha = 0.3), xlim = c())
lines(density(unlist(tr.n[,51:70])), col = adjustcolor("green3", alpha = 0.3))
lines(density(unlist(tr.n[,71:93])), col = adjustcolor("coral3", alpha = 0.3))

adj.m <- sweep(sweep(tr.n[,1:93], 2, bma.n$biasCoefs[2,1,1], "*"), 2, bma.n$biasCoefs[1,1,1], "+")

xl <- c(-500:1500)/100
qq <- apply(adj.m, 1:2, function(m) dnorm(xl, m, bma.n$sd))
qq.sum <- apply(qq, 1, sum) / (25*93)

plot(0, type = "n", ylim = c(0,1), xlim = range(xl), xlab = "temp.n", ylab = "")
#invisible(apply(qq, 2:3, function(d) lines(xl, d)))
plot(xl, qq.sum, type = "l")

plot(density(unlist(tr.n[,1:50])), col = adjustcolor("steelblue", alpha = 0.3), ylim = c(0,0.2))
lines(density(unlist(tr.n[,51:70])), col = adjustcolor("green3", alpha = 0.3))
lines(density(unlist(tr.n[,71:93])), col = adjustcolor("coral3", alpha = 0.3))
lines(xl, qq.sum, lwd = 3, col = "orange")

# should show densities after bias-correction
# THEN add the weighted sum used by BMA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##


####################################################################################################

# is regression done over mean forecast?
plot(apply(tr[,1:50], 1, mean), tr[,"observations"], pch = 20, 
     col = mapply(adjustcolor, rep(c("steelblue", "green3"), each = 25), alpha = 0.4))
abline(bc[1,"ecmwf"], bc[2,"ecmwf"], col = "red3")
abline(line(apply(tr[, 1:50], 1, mean), tr[, "observations"]), col = "orange")
line(apply(tr[,1:50], 1, mean), tr[,"observations"])
lm.bc <- lm(obs ~ ., data.frame("obs" = tr$obs, tr[,1:50]))
plot(rep(tr$obs, 50), unlist(tr[,1:50]), pch = 20, col = adjustcolor("steelblue", alpha = 0.4))
line(rep(tr$obs, 50), unlist(tr[,1:50]))

abline(bc[,2], bc[,1], col = "red3")

X <- as.matrix(tr[,1:50])
Y <- as.matrix(tr[,"observations"])

B <- solve((t(X) %*% X)) %*% t(X) %*% Y

plot(train.n$ecmwf, train.n$obs, pch = 20, xlab = "forecast", ylab = "obs", main = "Bias correction - ECMWF")
abline(lm(obs ~ fc, data.frame(fc = train.n$ecmwf, obs = train.n$obs)), col = "red3")