
library("SX.weather"); library("CB.Misc")

fitted <- readRDS("./Models/lambda-hist.rds")

####################################################################################################

# CODE MULTIVARIATE RANK VERIFICATION HISTOGRAM?

# ... AND APPLY ALL OF THE ABOVE TO THE BASE MODELS ALREADY FITTED
# (maybe print some BW plots to avoid having to carry laptop?)

####################################################################################################

# ENERGY SCORE (ensemble prediction)                                                            ####

# three variants, depending on whether ensemble, probabilistic or deterministic forecast.

# ensemble forecast:
o <- apply(obs[1:2,,], 1, c)
p <- apply(offset.forecast(ecmwf)[1:2,,,2,-1], c(1,4), cbind)

v.crps <- crpsDecomposition(obs = o, eps = p)

#-----------------------------------------------------

# for a single, univariate observation:
{
    o <- obs[1,1,1]
    p <- offset.forecast(ecmwf)[1,1,1,2,-1]
    m <- length(p)
    
    vc <- crpsDecomposition(o, array(p, dim = c(1,50)))
    
    n1 <- mean(sapply(p-o, norm, "f"))
    n2 <- sum(colSums(apply(sapply(p, function(i) p-i), 1:2, norm, "f"))) / (2 * m^2)
    
    n1 - n2 == vc$CRPS
    
    # manual calculation to confirm
    {
        nmat <- array(dim = c(50,50))
        invisible(sapply(1:50, function(i) {
            invisible(sapply(1:50, function(j) {
                nmat[i,j] <<- norm(p[i] - p[j], "f")
            }))
        }))
        all(nmat == n2)
    }
    
    # so, works in 1 dimension.
}

#-----------------------------------------------------

# for a single, multivariate observation:
{
    o <- obs[1:2,1,1]
    p <- offset.forecast(ecmwf)[1:2,1,1,2,-1]
    
    vc <- crpsDecomposition(o, p) 
    m <- ncol(p)
    
    norm.1 <- mean(apply(p-o, 2, norm, "f"))
    norm.2 <- sum(colSums(apply(p, 2, function(i) apply(p-i, 2, norm, "f")))) / (2 * m^2)
    
    norm.1 - norm.2
    
    # confirm by calculating manually
    {
        nmat <- array(dim = c(2,50,50))
        invisible(sapply(1:50, function(i) {
            invisible(sapply(1:50, function(j) {
                nmat[,i,j] <<- norm(p[,i] - p[,j], "f")
            }))
        }))
        
        all(nmat[1,,] == norm.2, nmat[2,,] == norm.2)
    }
}

#-----------------------------------------------------

# now full multivariate case over all observations (single lead-time)
{
    oo <- apply(obs[1:2,,], 1, c)
    pp <- apply(offset.forecast(ecmwf)[1:2,,,2,-1], c(1,4), cbind)
    
    m <- 50
    
    n1 <- apply(apply(sweep(pp, 1:2, oo, "-"), c(1,3), norm, "f"), 1, mean)
    norm.1 <- mean(apply(p-o, 2, norm, "f"))
    n1[1] == norm.1  # confirmed correct

    n2 <- square.mat(apply(pp, 1, function(j) apply(j, 2, function(i) apply(j-i, 2, norm, "f"))))
    norm.2 <- apply(p, 2, function(i) apply(p-i, 2, norm, "f"))
    
    norm.2 <- sum(colSums(apply(p, 2, function(i) apply(p-i, 2, norm, "f")))) / (2 * m^2)
}
# running quite slowly. Try alternative approach...

#-----------------------------------------------------

mv.crps.ens <- function(o, efc) {
    m <- ncol(efc)
    
    norm.1 <- mean(apply(efc-o, 2, function(v) sqrt(sum(v^2))))
    norm.2 <- sum(colSums(apply(efc, 2, function(i) apply(efc-i, 2, function(v) sqrt(sum(v^2)))))) / (2 * m^2)
    
    norm.1 - norm.2
} 

fudge <- abind("o" = oo, pp, along = 3)

# ~ 4s for original version over 630 observations
system.time(tt <- apply(fudge, 1, function(arr) mv.crps.ens(o = arr[,1], efc = arr[,-1])))

####################################################################################################

# ENERGY SCORE (probabilistic prediction)                                                       ####

# using MC approximation

fitted <- readRDS("./Models/lambda-hist.rds")

mv.crps.dens <- function(o, mu, sig, k = 10000) {
    
    if (length(mu) == 1) {
        
        x <- rnorm(k, mean = mu, sd = sqrt(sig))
        
        norm.1 <- mean(apply(t(x)-o, 2, function(v) sqrt(sum(v^2))))
        norm.2 <- sum(sapply(x[1:k-1] - x[2:k], function(v) sqrt(sum(v^2)))) / (2 * (k-1))
        
    } else {
        require(mvtnorm)            # needed to simulate from mimic
        
        x <- rmvnorm(k, mean = mu, sigma = sig)
        
        norm.1 <- mean(apply(t(x)-o, 2, function(v) sqrt(sum(v^2))))
        norm.2 <- sum(apply(x[1:k-1,] - x[2:k,], 1, function(v) sqrt(sum(v^2)))) / (2 * (k-1))
    }
    
    norm.1 - norm.2
}

mu <- fitted$tau[1,2,1]
sig <- fitted$s[1,1,2,1]
o <- obs[1,26,1]
k <- 1000

x <- rmvnorm(k, mean = mu, sigma = sig)
x <- rnorm(k, mean = mu, sd = sqrt(sig))

mv.crps.dens(obs[1:2,26,1], fitted$tau[1:2,2,1], fitted$s[1:2,1:2,2,1])
mv.crps.dens(obs[1,26,1], fitted$tau[1,2,1], fitted$s[1,1,2,1])

fudge <- abind("o" = array(t(apply(obs[1:2,,], 1, rbind))[,26:630], dim = c(2,1, 605)),
               "mu" = fitted$tau[1:2,2,,drop = F],
               "sig" = fitted$s[1:2,1:2,2,,drop = F],
               along = 1)

# slower to calculate all 630 observations: 
    # ~30s per leadtime with k = 10000,
    # ~ 3s with k = 1000
    # differences are quite small, so probably enough to use k=1000 for rough comparisons
    # 95% of CRPS scores are within 0.04 of each other.
system.time(qq <- apply(fudge, 4, function(arr) {
    mv.crps.dens(o = arr["o",,], mu = arr["mu",,], sig = arr[-(1:2),,], k = 1000) 
}))

boxplot(qq-zz, type = "l")

quantile(qq-zz, c(0.025, 0.5, 0.975))

####################################################################################################

# MULTIVARIATE VS UNIVARIATE CRPS                                                               ####

lt <- 2
ens.mv.crps <- apply(abind("o" = apply(obs[1:2,,], 1, c),
                           apply(offset.forecast(ecmwf)[1:2,,,lt,-1], c(1,4), cbind),
                           along = 3),
                     1, function(arr) mv.crps.ens(o = arr[,1], efc = arr[,-1]))

dens.mv.crps <- apply(abind("o" = array(t(apply(obs[1:2,,], 1, rbind))[,26:630], dim = c(2,1, 605)),
                            "mu" = fitted$tau[1:2,lt,,drop = F],
                            "sig" = fitted$s[1:2,1:2,lt,,drop = F],
                            along = 1),
                      4, function(arr) {
                          mv.crps.dens(o = arr["o",,], mu = arr["mu",,], sig = arr[-(1:2),,], k = 1000)
                      })

# now, treating as if univariate
ens.u.crps <- apply(abind("o" = apply(obs[1:2,,], 1, c),
                          apply(offset.forecast(ecmwf)[1:2,,,lt,-1], c(1,4), c),
                          along = 3),
                    1, function(arr) crpsDecomposition(arr[,1], arr[,-1])$CRPS)

dens.u.crps <- verification::crps(obs = c(apply(obs[1:2,,], 1, c)[26:630,]),
                                  pred = cbind(c(t(fitted$tau[1:2,lt,])),
                                               c(t(sqrt(apply(fitted$s[1:2,1:2, lt,], 3, diag))))))$crps


boxplot(ens.mv.crps, ens.u.crps, names = c("multivariate", "univariate"),
        main = "CRPS for ECMWF ensemble")

boxplot(dens.mv.crps, dens.u.crps, names = c("multivariate", "univariate"),
        main = "CRPS for fitted model")

####################################################################################################

# DETERMINANT SHARPNESS                                                                         ####

# suspect normalisation would be useful here?
det.sharpness <- function(sig) {
    det(sig)^(1 / (2 * dim(sig)[1]))
}

det.sharpness(fitted$s[,,1,1])
det.sharpness(fitted$s[1:2,1:2,1,1])

####################################################################################################

# 

####################################################################################################

# MULTIVARIATE RANK HISTOGRAM                                                                   ####

# consider standardizing data if including PCs?

# compare to ensemble BMA univariate version
require(ensembleBMA)

hh <- verifRankHist(apply(offset.forecast(ncep)[1,,,1,-1], 3, c),
                    c(obs[1,,]))
                    

plot(hh, type = "l")

verif.ranks <- function(o, ens) {
    
    # univariate case: don't distinguish between variables
    # collapse all dimensions into single vector (except for ensemble members)
    
    l <- length(dim(ens))

    apply(abind(o, ens, along = l), -l,
          function(v) which(sort(v) == v[1]))
}

vr.hist <- function(o, ens, breaks = c(0:10)/10, main = "", col = "skyblue", ...) {
    
    m <- dim(ens)[length(dim(ens))] + 1
    vr <- verif.ranks(o, ens)
    hist(vr/m, breaks = breaks, col = col, main = main, xlab = "", ylab = "", prob = T, ...)
    abline(h = 1, col = "red3", lty = 2)
}

vr.hist(obs[1,,], offset.forecast(ecmwf)[1,,,1,-1], main = "Verification rank hist - temp.n, ECMWF")
vr.hist(obs[1,,], offset.forecast(ncep)[1,,,1,-1], main = "Verification rank hist - temp.n, NCEP")
vr.hist(obs[1,,], offset.forecast(ukmo)[1,,,1,-1], main = "Verification rank hist - temp.n, UKMO")

vr.hist(obs[2,,], offset.forecast(ecmwf)[2,,,1,-1], main = "Verification rank hist - temp.s, ECMWF")
vr.hist(obs[2,,], offset.forecast(ncep)[2,,,1,-1], main = "Verification rank hist - temp.s, NCEP")
vr.hist(obs[2,,], offset.forecast(ukmo)[2,,,1,-1], main = "Verification rank hist - temp.s, UKMO")


# deviation from uniformity:
u.dev <- function(vr, l = max(vr)) {
    sum(sapply(1:l, function(j) abs(mean(vr == j) - 1/l)))
}

u.dev(verif.ranks(obs[1,,], offset.forecast(ecmwf)[1,,,1,-1]))
u.dev(verif.ranks(obs[1,,], offset.forecast(ncep)[1,,,1,-1]))
u.dev(verif.ranks(obs[1,,], offset.forecast(ukmo)[1,,,1,-1]))

u.dev(verif.ranks(obs[2,,], offset.forecast(ukmo)[2,,,1,-1]))

####################################################################################################

# TRENDS IN VERIFICATION RANK?                                                                  ####

# think this occurs more because of change in NCEP processing in 2010

vr <- abind("n" = abind("ecmwf" = invisible(sapply(1:15, 
                                                   function(lt) verif.ranks(obs[1,,], offset.forecast(ecmwf)[1,,,lt,-1]))),
                        "ncep" = invisible(sapply(1:15, 
                                                  function(lt) verif.ranks(obs[1,,], offset.forecast(ncep)[1,,,lt,-1]))),
                        "ukmo" = invisible(sapply(1:15, 
                                                  function(lt) verif.ranks(obs[1,,], offset.forecast(ukmo)[1,,,lt,-1]))),
                        along = 0),
            "s" = abind("ecmwf" = invisible(sapply(1:15, 
                                                   function(lt) verif.ranks(obs[2,,], offset.forecast(ecmwf)[2,,,lt,-1]))),
                        "ncep" = invisible(sapply(1:15, 
                                                  function(lt) verif.ranks(obs[2,,], offset.forecast(ncep)[2,,,lt,-1]))),
                        "ukmo" = invisible(sapply(1:15, 
                                                  function(lt) verif.ranks(obs[2,,], offset.forecast(ukmo)[2,,,lt,-1]))),
                        along = 0),
            along = 0)
            
plot(vr["n", "ecmwf", , 15], type = "l")
lines(vr["s", "ecmwf", , 15], col = "blue")
abline(line(vr["n", "ecmwf", , 15]), col = "red", lty = 2)
abline(line(vr["s", "ecmwf", , 15]), col = "red", lty = 2)

plot(vr["n", "ncep", , 15], type = "l")
lines(vr["s", "ncep", , 15], col = "blue")
abline(line(vr["n", "ncep", , 15]), col = "red", lty = 2)
abline(line(vr["s", "ncep", , 15]), col = "red", lty = 2)

plot(vr["n", "ukmo", , 15], type = "l")
lines(vr["s", "ukmo", , 15], col = "blue")
abline(line(vr["n", "ukmo", , 15]), col = "red", lty = 2)
abline(line(vr["s", "ukmo", , 15]), col = "red", lty = 2)

pdf("./Plots/Verif-hists-by-year-temp-ukmo.pdf", height = 28, width = 14); {
    par(mfrow = c(15, 7), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(1:15, function(lt) {
        invisible(sapply(1:7, function(y) {
            vr.hist(obs[1,,y], offset.forecast(ukmo)[1,,y,lt,-1], border = NA,
                    main = paste0("20", formatC(7:14, width = 2, flag = "0")[y], ", LT ", lt))
            vr.hist(obs[2,,y], offset.forecast(ukmo)[2,,y,lt,-1], add = T, col = NA, border = "darkred",
                    main = paste0("20", formatC(7:14, width = 2, flag = "0")[y], ", LT ", lt))
        }))
    }))
    title("Rank of verification against UKMO ensemble", outer = T)
}; dev.off()

# yearly mean ranks
vr.mean <- apply(array(vr, dim = c(2,3,90,7,15)), c(1:2,4:5), mean)

pdf("./Plots/Verif-rank-by-year-temp-ecmwf.pdf", height = 28, width = 14); {
    par(mfrow = c(15, 2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(1:15, function(lt) {
        plot(vr["n", "ecmwf", , lt], type = "l", xlab = "", ylab = "",
             main = paste0("temp.n, leadtime ", lt))
        abline(v = c(0:7) * 90, col = "cyan3", lty = 2)
        lines(c(rbind(90*(0:6), 90*(1:7), NA)), 
              rep(vr.mean[1, 1, , lt], each = 3),
              col = "red", lty = 2, lwd = 2)
        plot(vr["s", "ecmwf", , lt], type = "l", xlab = "", ylab = "", main = paste0("temp.s, leadtime ", lt))
        abline(v = c(0:7) * 90, col = "cyan3", lty = 2)
        lines(c(rbind(90*(0:6), 90*(1:7), NA)), 
              rep(vr.mean[2, 1, , lt], each = 3),
              col = "red", lty = 2, lwd = 2)
    }))
    title("Verification ranks for ECMWF", outer = T)
}; dev.off()
pdf("./Plots/Verif-rank-by-year-temp-ncep.pdf", height = 28, width = 14); {
    par(mfrow = c(15, 2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(1:15, function(lt) {
        plot(vr["n", "ncep", , lt], type = "l", xlab = "", ylab = "", main = paste0("temp.n for NCEP, leadtime ", lt))
        abline(v = c(0:7) * 90, col = "cyan3", lty = 2)
        lines(c(rbind(90*(0:6), 90*(1:7), NA)), 
              rep(vr.mean[1, 2, , lt], each = 3),
              col = "red", lty = 2, lwd = 2)
        plot(vr["s", "ncep", , lt], type = "l", xlab = "", ylab = "", main = paste0("temp.s for NCEP, leadtime ", lt))
        abline(v = c(0:7) * 90, col = "cyan3", lty = 2)
        lines(c(rbind(90*(0:6), 90*(1:7), NA)), 
              rep(vr.mean[2, 2, , lt], each = 3),
              col = "red", lty = 2, lwd = 2)
    }))
    title("Verification ranks for NCEP", outer = T)
}; dev.off()
pdf("./Plots/Verif-rank-by-year-temp-ukmo.pdf", height = 28, width = 14); {
    par(mfrow = c(15, 2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(1:15, function(lt) {
        plot(vr["n", "ukmo", , lt], type = "l", xlab = "", ylab = "", main = paste0("temp.n for UKMO, leadtime ", lt))
        abline(v = c(0:7) * 90, col = "cyan3", lty = 2)
        lines(c(rbind(90*(0:6), 90*(1:7), NA)), 
              rep(vr.mean[1, 3, , lt], each = 3),
              col = "red", lty = 2, lwd = 2)
        plot(vr["s", "ukmo", , lt], type = "l", xlab = "", ylab = "", main = paste0("temp.s for UKMO, leadtime ", lt))
        abline(v = c(0:7) * 90, col = "cyan3", lty = 2)
        lines(c(rbind(90*(0:6), 90*(1:7), NA)), 
              rep(vr.mean[2, 3, , lt], each = 3),
              col = "red", lty = 2, lwd = 2)
    }))
    title("Verification ranks for UKMO", outer = T)
}; dev.off()

err <- abind("ecmwf" = apply(forecast.errors(ecmwf)[,,,,-1], 1:4, mean),
             "ncep" = apply(forecast.errors(ncep)[,,,,-1], 1:4, mean),
             "ukmo" = apply(forecast.errors(ukmo)[,,,,-1], 1:4, mean),
             along = 0)[,1:2,,,]
err.mean <- apply(err, c(1:2, 4:5), mean)
err <- array(err, dim = c(3,2,630,15))
    
ens.cols <- c("red", "green3", "blue")
pdf("./Plots/FC-error-by-year.pdf", height = 4*15, width = 10); {
    
    par(mfrow = c(15, 2), mar = c(2,2,3,1), oma = c(0,0,2,0))
    
    invisible(sapply(1:15, function(lt) {
        invisible(sapply(1:2, function(varb) {
            matplot(t(err[,varb,,lt]), type = "l", lty = 1, col = adjustcolor(ens.cols, alpha = 0.5),
                    xlab = "", ylab = "", ylim = range(err),
                    main = paste0(c("temp.n", "temp.s")[varb], " - leadtime ", lt-1))
            abline(h = 0, lty = 2)
            abline(v = (0:7)*90, lty = 2)
            invisible(sapply(1:3, function(ens) {
                lines(c(rbind(90*(0:6), 90*(1:7), NA)), 
                      rep(err.mean[ens, varb, , lt], each = 3),
                      col = ens.cols[ens], lty = 2, lwd = 2)
            }))
        }))
    }))
    title("Forecast error over study period", outer = T)
}; dev.off()



####################################################################################################

####################################################################################################

# BOX DENSITY ORDINATE TRANSFORM                                                                ####

# u = 1 - P(f(x) <= f(o))                               # generally
# u = 1 - chisq(d, (o-mu)' S^-1 (o-mu))                 # MVN case
lt <- "0"

box.dot <- function(o, mu, s) {
    
    1 - pchisq(t(o - mu) %*% solve(s) %*% (o-mu), length(mu))
    
}

####################################################################################################

# TEST AUTOMATED FUNCTIONS (NOW ADDED TO PACKAGE)                                               ####

# energy score over multivariate & univariate ensembles
fitted <- readRDS("./Models/lambda-hist.rds")

uv.fitted <- es.crps(o = obs[1,90,7], mu = fitted$tau[1,1,605], sig = fitted$s[1,1,1,605])
mv.fitted <- es.crps(o = obs[1:2,90,7], mu = fitted$tau[1:2,1,605], sig = fitted$s[1:2,1:2,1,605])

uv.ens <- es.crps(o = obs[1,90,7], efc = offset.forecast(ecmwf)[1,90,7,1,-1])
mv.ens <- es.crps(o = obs[1:2,90,7], efc = offset.forecast(ecmwf)[1:2,90,7,1,-1])
