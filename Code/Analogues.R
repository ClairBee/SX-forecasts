
library("SX.weather"); library("CB.Misc")
library(plyr)           # use aaply instead of apply to avoid all that tedious reshaping

# identification of various analogue sets as training data for forecast error
# (also possibly for use in setting informative prior for alpha, Gamma)

# haven't done anything about cross-validation yet - just used everything lumped together.
# Will need to be more rigorous about using only past/CV data for the paper.

# eg. on day 2: use forecast for day 1 & day 2, observation for day 1.

####################################################################################################

# COMMON ELEMENTS                                                                               ####

# these matrices will be the same regardless of training set or variables retained.
C <- abind("ecmwf" = se.covariances(offset.forecast(ecmwf)),
           "ncep" = se.covariances(offset.forecast(ncep)),
           "ukmo" = se.covariances(offset.forecast(ukmo)),
           along = 0)

Y.bar <- abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
               "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
               "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
               along = 0)

Sig <- square.mat(apply(abind(lapply(list(ecmwf, ncep, ukmo),
                                     function (model) {
                                         apply(offset.forecast(model)[,,,,-1], 1:4, mean)
                                     }), along = 0), 3:5, cov))

D <- sweep(sweep(C, 1, table(attr(superensemble(), "m"))[2:4], "/"), 2:6, Sig, "+")

####################################################################################################

# ANALOGUES TO CURRENT FORECAST & YESTERDAY'S OBSERVATION, WITH PRINCIPAL COMPONENTS            ####

# on day t: use forecast for day t-1 & day t, observation for day t-1.
# t 2:90 for each year.

lt <- "0"

# create search space: candidate analogues (based on two days' forecast, one day's observations)
# padded for easier reference; not 'wrapping around' from year to year
cand <- abind("y.o" = obs[,1:89,],
              "y.ecmwf" = apply(offset.forecast(ecmwf)[,1:89,,lt,-1], 1:3, mean),
              "y.ncep" = apply(offset.forecast(ncep)[,1:89,,lt,-1], 1:3, mean),
              "y.ukmo" = apply(offset.forecast(ukmo)[,1:89,,lt,-1], 1:3, mean),
              "c.ecmwf" = apply(offset.forecast(ecmwf)[,2:90,,lt,-1], 1:3, mean),
              "c.ncep" = apply(offset.forecast(ncep)[,2:90,,lt,-1], 1:3, mean),
              "c.ukmo" = apply(offset.forecast(ukmo)[,2:90,,lt,-1], 1:3, mean),
              along = 0)
cand <- abind(array(NA, dim = dim(cand[,,1,])), cand, along = 3)

# standard deviation of each variable, for normalisation
cand.sd <- apply(cand, 1:2, sd, na.rm = T)

# target forecast
d <- 5; y <- 5; t <- cand[,,d,y]

# calculate distance between target & all candidates
dist <- sweep(sqrt(sweep(cand, 1:2, t, "-")^2), 1:2, cand.sd, "/")

# process is same with & without principal components up to this point.
#-------------------------------------------------------------------------------
m.dist <- apply(dist, 3:4, mean)

n <- 25

an <- which(m.dist <= sort(m.dist)[n+1], arr.ind = T)
an <- an[!(an[,1] == d & an[,2] == y),]

# training set: mean forecast errors across all 3 ensemble means
tr <- apply(apply(abind("ecmwf" = apply(forecast.errors(ecmwf)[,,,lt,-1], 1:3, mean),
                  "ncep" = apply(forecast.errors(ncep)[,,,lt,-1], 1:3, mean),
                  "ukmo" = apply(forecast.errors(ukmo)[,,,lt,-1], 1:3, mean), along = 0),
            2:4, mean), 1, "[", an)

Eta <- apply(tr, 2, mean)
Lambda <- cov(tr)

####################################################################################################

# FUNCTION TO FIND ANALOGUES & FIT RESULTING MODEL                                              ####

# add error checking to avoid that tedious padding business?
find.analogues <- function(targ, cand, n = 25) {
    
    # standard deviation for normalisation
    cand.sd <- apply(cand, 1:2, sd, na.rm = T)
    
    # distance from each point on target to equivalent in candidate space
    dist <- sweep(sqrt(sweep(cand, 1:2, targ, "-")^2), 1:2, cand.sd, "/")
    
    # mean distance across all variables
    m.dist <- apply(dist, 3:4, mean)
    
    # identify n analogues
    an <- which(m.dist <= sort(m.dist)[n+1], arr.ind = T)
    
    # find target vector and remove from analogue list
    t <- which(m.dist == 0, arr.ind = T)
    an <- an[!(an[,1] == t[1] & an[,2] == t[2]),]
    
    # return indices of analogues identified
    return(an)
}

an.dt <- abind(invisible(sapply(1:15, function(lt) {
    
    # search space for analogues: yesterday's obs + fc, today's fc
    cand <- abind("y.o" = obs[,1:89,],
                  "y.ecmwf" = apply(offset.forecast(ecmwf)[,1:89,,lt,-1], 1:3, mean),
                  "y.ncep" = apply(offset.forecast(ncep)[,1:89,,lt,-1], 1:3, mean),
                  "y.ukmo" = apply(offset.forecast(ukmo)[,1:89,,lt,-1], 1:3, mean),
                  "c.ecmwf" = apply(offset.forecast(ecmwf)[,2:90,,lt,-1], 1:3, mean),
                  "c.ncep" = apply(offset.forecast(ncep)[,2:90,,lt,-1], 1:3, mean),
                  "c.ukmo" = apply(offset.forecast(ukmo)[,2:90,,lt,-1], 1:3, mean),
                  along = 0)
    cand <- abind(array(NA, dim = dim(cand[,,1,])), cand, along = 3)
    
    all.an <- suppressWarnings(aaply(cand[,,-1,], 3:4, find.analogues, cand))
    all.an <- abind(array(NA, dim = dim(all.an[1,,,])), all.an, along = 1)          # pad
    
    all.an}, simplify = F)), along = 2.5)

# Eta & Lambda are to be trained on ens. mean errors for analogue days
em.error <- apply(abind("ecmwf" = apply(offset.forecast(ecmwf[,,,,-1]), 1:4, mean),
                        "ncep" = apply(offset.forecast(ncep[,,,,-1]), 1:4, mean),
                        "ukmo" = apply(offset.forecast(ukmo[,,,,-1]), 1:4, mean),
                        along = 0), 2:5, mean)
  
# extract training set for each day's analogues  
tr.mean.error <- aaply(an.dt, 3, function(all.an) aaply(all.an, 1:2, function(an) {aaply(em.error, 1, "[", an)}))

tr.mean.error <- suppressWarnings(aaply(all.an, 1:2, function(an) {aaply(em.error, 1, "[", an)}))

# now fit model: mean & covariance of training forecast errors
Eta <- apply(tr.mean.error, 1:3, mean)
Lambda <- aaply(tr.mean.error, 1:2, function(err) cov(t(err)))
Eta <- aperm(Eta, c(3,1:2)); Lambda <- aperm(Lambda, c(3:4, 1:2))

# S & Tau use common elements already calculated above
S <- sweep(array.solve(apply(array.solve(D, c(1,4:6)), c(1:2,4:6), sum), 3:5)[,,,,1], 1:4, Lambda, "+")

e1 <- square.mat(apply(abind("D" = apply(array.solve(D[,,,,,1], c(1,4:5)), c(1:2, 4:5), sum),
                  "L" = Lambda, 
                  along = 0), 4:5, function(arr) arr["D",,] %*% arr["L",,]))



   apply(array.solve(D[,,,,,1], c(1,4:5)), c(1:2,4:5), sum) %*% Lambda
, 1:2, diag(5), "+")


    

####################################################################################################

# ORIGINAL CODE. OFFSETTING IS NOT AS CLEAR AS IT SHOULD BE                                     ####

lt <- "0"               # search single leadtime only (repeat later)

# search space: candidate analogues (includes yesterday's observation)
cand <- abind("o" = obs[,1:89,],
              "ecmwf" = apply(offset.forecast(ecmwf)[,2:90,,lt,-1], 1:3, mean),
              "ncep" = apply(offset.forecast(ncep)[,2:90,,lt,-1], 1:3, mean),
              "ukmo" = apply(offset.forecast(ukmo)[,2:90,,lt,-1], 1:3, mean),
              along = 0)

# standard deviation of each variable, for normalisation
cand.sd <- apply(cand, 1:2, sd)

# target forecast
d <- 1; y <- 1
t <- cand[,,d,y]

dist <- sweep(sqrt(sweep(cand, 1:2, t, "-")^2), 1:2, cand.sd, "/")

# process is same with & without principal components up to this point.
#-------------------------------------------------------------------------------

m.dist <- apply(dist, 3:4, mean)
an <- which(m.dist <= sort(m.dist)[26], arr.ind = T)[-1,]

tr <- apply(cand, 1:2, "[", an)

# check arrangement of analogues around target
{
    plot(tr[,1,1:2], pch = 4, xlim = range(tr[,,1]), ylim = range(tr[,,2]))
    points(tr[,2,1:2], pch = 20, col = adjustcolor("steelblue", alpha = 0.4))
    points(tr[,3,1:2], pch = 20, col = adjustcolor("coral3", alpha = 0.4))
    points(tr[,4,1:2], pch = 20, col = adjustcolor("green3", alpha = 0.4))
    
    points(t[,1:2], pch = c(4, rep(20, 6)), col = c("black", rep(c("steelblue","coral3", "green3"), 2)), lwd = 2)

    vo <- apply(obs, 1, "[", sweep(an, 2, c(1,0), "+"))
    points(vo[,1:2], col = "orange", pch = 20)
    points(t(apply(vo, 2, mean)[1:2]), pch = 20)        # mean observation in training set
    
    points(obs[1:2,2,1], col = "red", pch = 20)
}

# compute Eta & Lambda over training data
tr.error <- apply(sweep(tr[,2:4,], c(1,3), tr[,1,], "-"), c(1,3), mean)

Eta <- apply(tr.error, 2, mean)
Lam <- cov(tr.error)

# and now, applying over all days...
tr.error <- array(dim = c(25,5,89,7, 15))
    
invisible(sapply(1:15, function(lt) {
    
    # search space: candidate analogues (includes yesterday's observation)
    cand <- abind("o" = obs[,1:89,],
                  "ecmwf" = apply(offset.forecast(ecmwf)[,2:90,,lt,-1], 1:3, mean),
                  "ncep" = apply(offset.forecast(ncep)[,2:90,,lt,-1], 1:3, mean),
                  "ukmo" = apply(offset.forecast(ukmo)[,2:90,,lt,-1], 1:3, mean),
                  along = 0)
    
    # standard deviation of each variable, for normalisation
    cand.sd <- apply(cand, 1:2, sd)

invisible(sapply(1:7, function(y) {
    invisible(sapply(1:89, function(d) {
        t <- cand[,,d,y]
        dist <- sweep(sqrt(sweep(cand, 1:2, t, "-")^2), 1:2, cand.sd, "/")
        
        m.dist <- apply(dist, 3:4, mean)
        tr.error[,,d,y,lt] <<- apply(sweep(tr[,2:4,], c(1,3), tr[,1,], "-"), c(1,3), mean)
    }))}))}))

# calculate Eta & Lambda for every day in series
Eta <- apply(tr.error, 2:5, mean)
Lambda <- square.mat(apply(tr.error, 3:5, cov))

# Possibly something wrong with Lambda calculation? Come back to this when better.

# fit model for every day in series
Tau <- array(dim = c(5,89,7,15))
S <- array(dim = c(5,5,89,7,15))

invisible(sapply(1:15, function(lt) {
    invisible(sapply(1:7, function(y) {
        invisible(sapply(2:90, function(d) {
            LL <- Lambda[,,d-1,y,lt]; DD <- D[,,,d,y,lt]; YY <- Y.bar[,,d,y,lt]
            
            S[,,d-1,y,lt] <<- solve(apply(array.solve(DD, 1), 1:2, sum)) + LL
            
            e1 <- solve((solve(apply(array.solve(DD, 1), 1:2, sum))%*% LL) + diag(5))
            e2 <- apply(apply(abind("Y" = t(YY), array.solve(DD, 1), along = 1), 3, 
                        function(arr) arr[-1,] %*% arr[1,]), 1, sum)
            
            (S[,,d-1,y,lt] %*% e1 %*% e2) - Eta[,d,y,lt]
        }))
    }))
}))
    
####################################################################################################

# ANALOGUES TO CURRENT FORECAST & YESTERDAY'S OBSERVATION, NO PRINCIPAL COMPONENTS              ####

####################################################################################################

# identical to above until mean distances are found
m.dist <- apply(dist[,1:2,,], 3:4, mean)