

library("SX.weather")

tr.all <- readRDS("~/Documents/PhD/Miniprojects/03-Sichun-paper/Models/Analogue-indices-25.rds")

# check effect of Bayesian inference over training set vs simple mean/covariance
d <- 5; y <- 5; lt <- 15

tr <- tr.all[d,y,lt,,]

err <- lapply(list("ec" = aaply(tr, 1, function(an) offset.forecast(ecmwf)[1:2, an[1], an[2], lt, -1]),
                   "nc" = aaply(tr, 1, function(an) offset.forecast(ncep)[1:2, an[1], an[2], lt, -1]),
                   "mo" = aaply(tr, 1, function(an) offset.forecast(ukmo)[1:2, an[1], an[2], lt, -1])),
              sweep, 1:2, aaply(tr, 1, function(an) obs[1:2, an[1], an[2]]), "-")

par(mfrow = c(1,1), mar = c(2,2,1,1), pch = 20)
xl <- yl <- c(-7.5,7.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simple mean & covariance of ensemble means for eta, Lambda                                ####

ens.mean <- apply(abind(lapply(err, apply, 1:2, mean), along = 0), 2:3, mean)
simple.eta <- t(apply(ens.mean, 2, mean))
simple.Lambda <- cov(ens.mean)

# plot simple model
{
    plot(ens.mean, pch = 4, xlim = xl, ylim = yl, col = "blue")
    lines(ellipse(simple.Lambda, centre = simple.eta), col = "blue")
    points(simple.eta, pch = 20, col = "blue")
    
    legend("bottomright", col = "blue", pch = c(4, 20, NA), lty = c(NA, NA, 1), bty = "n",
           c("Analogue mean errors", expression(paste("Simple ", eta)), expression(paste("Simple ",Lambda))))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bayesian inference over all data (ens. members truncated)                                 ####

# ensemble means for each analogue are new ensemble members
delta.ij <- abind(lapply(err, apply, 1:2, mean), along = 0)

delta.bar <- apply(delta.ij, c(1,3), mean)
delta.C <- aaply(delta.ij, 1, cov)
delta.n <- sapply(lapply(err, dim), "[", 3)

delta.Sigma <- cov(delta.bar)

delta.D <- sweep(sweep(delta.C, 1, delta.n^{-1}, "*"), 2:3, delta.Sigma, "+")

delta.S <- solve(apply(aaply(delta.D, 1, solve), 2:3, sum))
delta.tau <- delta.S %*% apply(apply(abind("bar" = delta.bar, delta.D, along = 3), 1, function(arr) {
    solve(arr[,-1]) %*% arr[,1]  
}), 1, sum)

# plot both models for easy comparison                                                      
{
    plot(0, type = "n", xlim = xl, ylim = yl, xlab = "", ylab = "")
    
    lines(ellipse(delta.C[1,,], centre = t(delta.bar[1,])), col = adjustcolor("steelblue", alpha = 0.6))
    points(t(delta.bar[1,]), col = adjustcolor("steelblue", alpha = 0.6))
    lines(ellipse(delta.C[2,,], centre = t(delta.bar[2,])), col = adjustcolor("red3", alpha = 0.6))
    points(t(delta.bar[2,]), col = adjustcolor("red3", alpha = 0.6))
    lines(ellipse(delta.C[3,,], centre = t(delta.bar[3,])), col = adjustcolor("green3", alpha = 0.6))
    points(t(delta.bar[3,]), col = adjustcolor("green3", alpha = 0.6))
    
    lines(ellipse(delta.Sigma, centre = simple.eta), col = "grey")
    
    lines(ellipse(simple.Lambda, centre = simple.eta), col = "blue")
    points(simple.eta, pch = 20, col = "blue")
    
    lines(ellipse(delta.S, centre = t(delta.tau)), col = "red")
    points(t(delta.tau), pch = 20, col = "red")
    
    legend("bottomright", bty = "n", pch = c(rep(20, 3), NA, rep(20,2)), lty = 1, cex = 0.8,
           col = c(mapply(adjustcolor, c("steelblue", "red3", "green3"), 0.6), "grey", "red", "blue"),
           c("ECMWF", "NCEP", "UKMO", expression(paste(Sigma)), expression(paste("Bayesian ", eta)),
             expression(paste("Simple ", eta))))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# fully Bayesian (ie. no truncation)                                                        ####

# speculative model....
ik.mean <- abind(lapply(err, apply, 1:2, mean), along = 0)       # ensemble/analogue mean & cov
ik.cov <- abind(lapply(err, aaply, 1, function(ens) cov(t(ens))), along = 0)

i.mean <- apply(ik.mean, c(1,3), mean)                          # ensemble mean & cov
i.cov <- aaply(ik.mean, 1, cov)

i.sig <- cov(i.mean)

ik.n <- sapply(lapply(err, dim), "[", 3)
i.n <- sapply(lapply(err, dim), "[", 1)

i.D <- sweep(sweep(sweep(ik.cov, 1, ik.n^{-1}, "*"), 
                   c(1,3:4), sweep(i.cov, 1, i.n^{-1}, "*"), "+"),
             c(3:4), i.sig, "+")

i.S <- solve(apply(aaply(i.D, 1:2, solve), 3:4, sum))
i.tau <- t(i.S %*% apply(aaply(abind(ik.mean, i.D, along = 3), 1:2, function(arr) solve(arr[-1,]) %*% arr[1,]), 3, sum))

points(i.tau, col = "magenta3")
lines(ellipse(i.S, centre = i.tau), col = "magenta3")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bayesian inference over all analogues independently                                       ####

er.bar <- er.C <- er.D <- er.Sig <- er.S <- er.tau <- list()

for (i in 1:25) {
    fc.err <- lapply(err, "[", i,,)
    
    er.bar[[i]] <- abind(lapply(fc.err, apply, 1, mean), along = 0)
    er.C[[i]] <- abind(lapply(lapply(fc.err, t), cov), along = 0)
    er.Sig[[i]] <- cov(er.bar[[i]])
    
    er.D[[i]] <- sweep(sweep(er.C[[i]], 1, delta.n^{-1}, "*"), 2:3, er.Sig[[1]], "+")
    
    er.S[[i]] <- solve(apply(aaply(er.D[[i]], 1, solve), 2:3, sum))
    
    er.tau[[i]] <- er.S[[i]] %*% apply(aaply(abind(er.bar[[i]], er.D[[i]], along = 3), 1,
                                             function(arr) solve(arr[,-1]) %*% arr[,1]), 2, sum)
}

est.tau <- abind(er.tau, along = 0)[,,1]

# plot inferred errors
{
    plot(0, type = "n", xlim = xl, ylim = yl, xlab = "", ylab = "")
    
    points(est.tau, col = "grey")
    points(t(apply(est.tau, 2, mean)), col = "red")
    invisible(sapply(1:25, function(i) {
        lines(ellipse(er.S[[i]], centre = er.tau[[i]][,1]), col = adjustcolor("darkgrey", alpha = 0.5))
    }))
    lines(ellipse(cov(est.tau), centre = apply(est.tau, 2, mean)), col = "red")
    
    lines(ellipse(simple.Lambda, centre = simple.eta), col = "blue")
    points(simple.eta, pch = 20, col = "blue")
    
    legend("bottomright", bty = "n", col = c("grey", "red", "blue"), pch = 20, lty = 1,
           c(expression(paste("Analogue ", tau, " & ", Lambda)),
             expression(paste("Overall mean ", tau, " & ", Lambda(tau))),
             "Simple mean & covariance"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use alternative eta, Lambda in fitted model                                               ####


fit.eta <- function(fc, Eta, Lambda) {
    d <- dim(fc[[1]])[[1]]
    Y.bar <- abind(lapply(fc, apply, 1, mean), along = 0)
    C <- abind(lapply(lapply(fc, t), cov), along = 0)
    D.i <- aaply(sweep(sweep(C, 1, sapply(fc, dim)[2, ], "/"), 
                       2:3, cov(Y.bar), "+"), 1, solve)
    
    S <- Lambda + solve(apply(D.i, 2:3, sum))
    Tau <- (S %*% ((solve(diag(d) + (apply(D.i, 2:3, sum) %*% 
                                         Lambda))) %*% apply(aaply(abind(Y.bar, D.i, along = 2), 
                                                                   1, function(arr) arr[-1, ] %*% arr[1, ]), 2, sum))) - 
        Eta
    list(Tau, S)
}


simple.post <- fit.eta(fc = list("ec" = offset.forecast(ecmwf)[1:2,d,y,lt,-1],
                                   "nc" = offset.forecast(ncep)[1:2,d,y,lt,-1],
                                   "mo" = offset.forecast(ukmo)[1:2,d,y,lt,-1]),
                        Eta = simple.eta[1,], Lambda = simple.Lambda)

delta.post <- fit.eta(fc = list("ec" = offset.forecast(ecmwf)[1:2,d,y,lt,-1],
                                 "nc" = offset.forecast(ncep)[1:2,d,y,lt,-1],
                                 "mo" = offset.forecast(ukmo)[1:2,d,y,lt,-1]),
                       Eta = delta.tau[,1], Lambda = delta.S)

plot(ellipse(simple.post[[2]], centre = simple.post[[1]]), type = "l", xlim = c(-5,15), ylim = c(-5,15))
lines(ellipse(delta.post[[2]], centre = delta.post[[1]]), col = "red")
points(t(simple.post[[1]]), pch = 20)
points(t(delta.post[[1]]), pch = 20, col = "red")
points(t(obs[1:2,d,y]), pch = 4, lwd = 2)
