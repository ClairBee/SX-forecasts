
library("SX.weather")

vars <- 1
lt <- 5; d <- 32; y <- 5 

# get data
fc <- cbind(apply(apply(offset.forecast(ecmwf)[vars,,,lt,-1, drop = F], c(1,5), c), 3, c),
            apply(apply(offset.forecast(ncep)[vars,,,lt,-1, drop = F], c(1,5), c), 3, c),
            apply(apply(offset.forecast(ukmo)[vars,,,lt,-1, drop = F], c(1,5), c), 3, c))
colnames(fc) <- names(exch)

o <- c(apply(obs[vars,,], 1, c))

# get BMA bias coefficients
bc <- abind(lapply(list("ec" = cbind(c(obs[vars,d-(24:0),y]), 1, c(offset.forecast(ecmwf)[vars,d-(24:0),y,lt,-1])),
                               "nc" = cbind(c(obs[vars,d-(24:0),y]), 1, c(offset.forecast(ncep)[vars,d-(24:0),y,lt,-1])),
                               "mo" = cbind(c(obs[vars,d-(24:0),y]), 1, c(offset.forecast(ukmo)[vars,d-(24:0),y,lt,-1]))),
                          function(dat) {
                              Y.n <- as.matrix(dat[,1])
                              X.n <- as.matrix(dat[,2:3])
                              
                              solve(t(X.n) %*% X.n) %*% t(X.n) %*% Y.n
                          }), along = 0)


# log-likelihood function for EM algorithm
log.likelihood <- function(sig2, w, fc, bc, o) {
    
    # summed over all observations in training set
    l <- array(dim = dim(fc), dimnames = dimnames(fc))
    
    invisible(sapply(1:3, function(i) {
        invisible(sapply(1:25, function(t) {
            l[i,t] <<- w[i] * dnorm(o[t], bc[i,1] + bc[i,2] * fc[i,t], sqrt(sig2))
        }))}))
    sum(log(colSums(l)))
}

# EM algorithm
EM <- function(sig2 = 1, w = rep(1/3, 3), max.runs = 1000, conv = 0.0001) {
    
    log.lh <- log.likelihood(sig2, w)
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
        new.log.lh <- log.likelihood(sig2, w)
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