
# replicate output of BMA & MOS models to confirm that my representation of the models is ok

library("SX.weather"); library("ensembleBMA")

lt <- 5
dt <- "20081201"
timestamps <- gsub("-", "", as.Date(load.data("../Data/ECMWF_europe_starttime.rda")[16:105,1:7]/ 24, "1800-01-01"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Univariate case                                                                               ####

# create data in required format
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Multivariate case: estimating coefficient over both locations                                 ####

# convert data to required format
{
    fc <- cbind(apply(apply(offset.forecast(ecmwf)[1:2,,,lt,-1], c(1,4), c), 3, c),
                apply(apply(offset.forecast(ncep)[1:2,,,lt,-1], c(1,4), c), 3, c),
                apply(apply(offset.forecast(ukmo)[1:2,,,lt,-1], c(1,4), c), 3, c))
    colnames(fc) <- c(paste0("ec", 1:50), paste0("nc", 1:20), paste0("mo", 1:23))
    
    exch <- c(rep(1, 50), rep(2, 20), rep(3, 23))
    names(exch) <- colnames(fc)
    
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
}

# fit a single BMA model to examine
{
    bma.fit <- ensembleBMA(temp.dat, model = "normal", trainingDays = 25,
                           dates = dt)
    
    tr <- trainingData(temp.dat, trainingDays = 25, date = dt)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# BMA MULTIVARIATE REGRESSION COEFFICIENTS                                                      ####

# ensemble BMA starts with linear regression over each ensemble's training data
# extract bias coefficients
bc <- cbind("ecmwf" = apply(bma.fit$biasCoefs[,1:50,], 1, mean),
            "ncep" = apply(bma.fit$biasCoefs[,51:70,], 1, mean),
            "ukmo" = apply(bma.fit$biasCoefs[,71:93,], 1, mean))

lm.all <- cbind("ecmwf" = lm(obs ~ fc, data.frame(obs = tr[,95], fc = unlist(tr[,1:50])))$coef,
                "ncep" = lm(obs ~ fc, data.frame(obs = tr[,95], fc = unlist(tr[,51:70])))$coef,
                "ukmo" = lm(obs ~ fc, data.frame(obs = tr[,95], fc = unlist(tr[,71:93])))$coef)

find.B.hat <- function(c.range, c.obs = 95, dat) {
    X.n <- as.matrix(cbind(1, unlist(dat[,c.range])))
    Y.n <- as.matrix(rep(dat[,c.obs], nrow(X.n) / nrow(dat)))
    
    solve(t(X.n) %*% X.n) %*% t(X.n) %*% Y.n
}

B.hat <- cbind("ecmwf" = find.B.hat(1:50, dat = tr),
               "ncep" = find.B.hat(51:70, dat = tr),
               "ukmo" = find.B.hat(71:93, dat = tr))

all(round(B.hat, 9) == round(bc, 9))

##################################

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