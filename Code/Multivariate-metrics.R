
library("SX.weather"); library("CB.Misc")

fitted <- readRDS("./Models/lambda-hist.rds")

####################################################################################################

# CODE MULTIVARIATE RANK VERIFICATION HISTOGRAM
# CODE BOX ORDINAL TRANSFORM?
# CODE MULTIVARIATE PIT
# CODE MST

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
    
    require(mvtnorm)            # needed to simulate from mimic
    
    x <- rmvnorm(k, mean = mu, sigma = sig)
    
    norm.1 <- mean(apply(t(x)-o, 2, function(v) sqrt(sum(v^2))))
    norm.2 <- sum(apply(x[1:k-1,] - x[2:k,], 1, function(v) sqrt(sum(v^2)))) / (2 * (k-1))
    
    norm.1 - norm.2
}

mv.crps.dens(obs[1:2,26,1], fitted$tau[1:2,2,1], fitted$s[1:2,1:2,2,1])

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