
require(ensembleMOS); require(ensembleBMA); require(fields); require(maps)
library("SX.weather")

# add Brier score for probability of freezing?

# function to produce verification histograms

####################################################################################################

# PACKAGE DOCUMENTATION EXAMPLES                                                                ####

# example: surface temperature data
{
    # create ensemble data object
    data(srft)

    memberLabels <- c("CMCG","ETA","GASP","GFS","JMA","NGPS","TCWB","UKMO")
    
    srftData <- ensembleData(forecasts = srft[,memberLabels],
                             dates = srft$date, observations = srft$obs,
                             latitude = srft$lat, longitude = srft$lon,
                             forecastHour = 48, initializationTime = "00")
    
    # can specify exchangeable models by including vector of groupings, eg
    # exchangeable = c(CMCG = 1, ETA = 2, GASP = 3, GFS = 2, JMA = 4, NGPS = 5, TCWB = 6, UKMO = 7)
    
    srftFit <- ensembleBMA(srftData, dates = "2004013100",
                           model = "normal", trainingDays = 25)
    # quite slow process - 
    
    # produce a plot for the pdf of every day
    plot(srftFit, srftData, dates = "2004013100")
}


# can also be run over a grid
{
    data(srftGrid)
    memberLabels <- c("CMCG","ETA","GASP","GFS","JMA","NGPS","TCWB","UKMO")
    srftGridData <- ensembleData(forecasts = srftGrid[,memberLabels],
                                 latitude = srftGrid[,"latitude"], longitude = srftGrid[,"longitude"],
                                 date = "2004013100", forecastHour = 48, initializationTime = "00")
    gridForc <- quantileForecast( srftFit, srftGridData,
                                  quantiles = c( .1, .5, .9))
    probFreeze <- cdf( srftFit, srftGridData, date = "2004013100",
                       value = 273.15)
    
    plotProbcast( gridForc[,"0.5"], lon=srftGridData$lon,
                  lat=srftGridData$lat, type="image",
                  col=rev(rainbow(100,start=0,end=0.85)), asp = T)
    title("Median Forecast for Temperature", cex = 0.5)
    
    bluescale <- function(n)
        hsv(4/6, s = seq(from = 1/8, to = 1, length = n), v = 1)
    plotProbcast( probFreeze, lon=srftGridData$lon, lat=srftGridData$lat,
                  type="image", col=bluescale(100), asp = T)
    title("Probability of Freezing", cex = 0.5)
}

# Obtain & plot quantile forecast temperatures
{
    srftForc <- quantileForecast( rftFit, srftData, quantiles = c( .1, .5, .9))
    
    use <- as.character(srftData$dates) == "2004013100"
    lat <- srftData$latitude[use]; lon <- srftData$longitude[use]
    lonRange <- range(lon); latRange <- range(lat)
    range(srftForc[,"0.5"]) # used to determine contour levels
    
    color <- "brown"; mapColor <- "black"
 
    plotProbcast( srftForc[,"0.5"], lon, lat, interpolate = TRUE, col = color,
                  type = "contour", levels = seq(from=264, to=284, by=2), asp = T)
    title("Median Forecast")
    points(lon, lat, pch = 16, cex = 0.5, col = color) # observation locations
    
    plotProbcast( srftData$obs[use], lon, lat, interpolate = TRUE, col = color,
                  type = "contour", levels = seq(from=264, to=284, by=2), asp = T)
    title("Observed Surface Temperature")
    points(lon, lat, pch = 16, cex = 0.5, col = color)
}

# verification
{
    crps(fit = srftFit, ensembleData = srftData)
    MAE( srftFit, srftData)
}

# CRPS function not working. Run test data to compare structure
{
    data(ensBMAtest)
    
    ensMemNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo")
    
    obs <- paste("T2","obs", sep = ".")
    ens <- paste("T2", ensMemNames, sep = ".")
    
    tempTestData <- ensembleData( forecasts = ensBMAtest[,ens],
                                  dates = ensBMAtest[,"vdate"],
                                  observations = ensBMAtest[,obs],
                                  station = ensBMAtest[,"station"],
                                  forecastHour = 48,
                                  initializationTime = "00")
    
    tempTestFit <- ensembleBMAnormal( tempTestData, trainingDays = 30)
    
    # for quick run only; use more training days for forecasting
    tempTestFit <- ensembleBMAnormal( tempTestData[1:20,], trainingDays = 8)
    
    crpsValues <- crps(tempTestFit, tempTestData)
    colMeans(crpsValues)
    
    CRPS( tempTestFit, tempTestData)
}

# verification histograms
{
    use <- ensembleValidDates(srftData) >= "2004013000"
    verifRankHist( ensembleForecasts(srftData[use,]),
                   dataVerifObs(srftData[use,]))
}
####################################################################################################

# CONVERT OUR DATA TO SAME FORMAT                                                               ####

# match srft format

timestamps <- load.data("./Data/ECMWF_europe_starttime.rda")

# only a single variable at a time, and single leadtime - use temp.n to start with
# later, concatenate with temp.s

temp.n <- ensembleData(forecasts = cbind("ECMWF" = c(apply(ecmwf["temp.n",16:105,,"0",-1], 1:2, mean)),
                                         "NCEP" = c(apply(ncep["temp.n",16:105,,"0",-1], 1:2, mean)),
                                         "UKMO" = c(apply(ukmo["temp.n",16:105,,"0",-1], 1:2, mean))),
                       dates = gsub("-", "", as.Date(timestamps[16:105,1:7]/ 24, "1800-01-01")),
                       observations = c(obs["temp.n",,]),
                       forecastHour = 0,
                       latitude = rep(mean(load.data("./Data/ECMWF_europe_lat.rda")), 630),
                       longitude = rep(mean(load.data("./Data/ECMWF_europe_lon.rda")), 630),
                       initializationTime = "00")

# takes ~1m to run for all available dates
temp.n.fit <- ensembleBMA(temp.n, model = "normal", trainingDays = 25,
                          dates = gsub("-", "", as.Date(timestamps[16:105,1:7]/ 24, "1800-01-01"))[-c(1:24)])
dplot(temp.n.fit, temp.n, dates = "20071225", ask = F)
# red: ECMFW, green: NCEP, blue: UKMO

matplot(t(temp.n.fit$weights), type = "l", col = c("black", "blue", "red"))
round(temp.n.fit$weights[,1],3)

probFreeze <- cdf(temp.n.fit, temp.n, value = 0)
plot(probFreeze, type = "l")

temp.n.fc <- quantileForecast(temp.n.fit, temp.n, quantiles = c( .1, .5, .9))

plot(obs["temp.n", ,1], type = "l", lwd = 2)
lines(25:90, temp.n.fc[1:66,2], col = "red", lwd = 2)
lines(25:90, temp.n.fc[1:66,1], col = "blue")
lines(25:90, temp.n.fc[1:66,3], col = "blue")


# can't currently evaluate results; package method has an error

# PIT histogram to assess calibration of fitted model
hist(pit(temp.n.fit, temp.n), breaks = "fd", col = "skyblue",
     main = "PIT histogram of BMA model for temp.n", xlab = "", ylab = "", prob = T)
abline(1,0, lty = 2)

use <- ensembleValidDates(temp.n) >= "20071225"
verifRankHist( ensembleForecasts(temp.n[use,]),
               dataVerifObs(temp.n[use,]))

temp.n.param <- modelParameters(temp.n.fit)

####################################################################################################

# CREATE ENSEMBLE DATA FROM FORECASTS                                                           ####

tr <- 25            # length of training set
lt <- "5"             # leadtime
timestamps <- gsub("-", "", as.Date(load.data("./Data/ECMWF_europe_starttime.rda")[16:105,1:7]/ 24, "1800-01-01"))

fc.temps <- function(model, lt) {
    c(apply(apply(offset.forecast(model[1:2,,,,-1]), 1:4, mean)[,,,toString(lt)], 1, cbind))
}

# create ensemble in correct format
temps <- ensembleData(forecasts = data.frame("ecmwf" = fc.temps(ecmwf, lt),
                                            "ncep" = fc.temps(ncep, lt),
                                            "ukmo" = fc.temps(ukmo, lt)),
                       dates = rep(timestamps, 2),
                       observations = c(obs["temp.n",,], obs["temp.s",,]),
                       forecastHour = 00,
                       station = rep(c("N", "S"), each = 630),
                       latitude = rep(c(56, 52), each = 630),
                       longitude = rep(c(-3, 0), each = 630),
                       initializationTime = "00")

####################################################################################################

# ENSEMBLE BMA MODEL                                                                            ####

# fit BMA model
temps.bma <- ensembleBMA(temps, model = "normal", trainingDays = tr,
                         dates = timestamps[timestamps >= timestamps[tr]])
                             
# plot of weights per model
matplot(t(m.bma[[1]]$weights), type = "l", col = c("black", "blue", "red"),
        main = "BMA weights", xlab = "", ylab = "")

# plot probability of freezing
bmaFreeze <- cdf(temps.bma, temps, value = 0)
plot(bmaFreeze, type = "l")

# quantile forecasts
temps.bma.forecast <- quantileForecast(temps.bma, temps, quantiles = c( .1, .5, .9))
plot(obs["temp.n", ,1], type = "l", lwd = 2)
lines(25:90, temps.bma.forecast[1:66,2], col = "red", lwd = 2)
lines(25:90, temps.bma.forecast[1:66,1], col = "blue")
lines(25:90, temps.bma.forecast[1:66,3], col = "blue")


fc <- matrix(temps.bma.forecast[,"0.5"], ncol = 2)

matplot(fc, type = "l", lty = 1)

bma5 <- fitted.bma["5"]

# PIT histogram to assess calibration of fitted model
hist(pit(fitted.bma[["5"]], temps), breaks = "fd", col = "skyblue",
     main = "PIT histogram of BMA model for temps", xlab = "", ylab = "", prob = T)
abline(1,0, lty = 2)

use <- ensembleValidDates(temps) >= timestamps[tr]
verifRankHist(ensembleForecasts(temps[use,]),
              dataVerifObs(temps[use,]))

bma.brier <- brierScore(temps.bma, temps, 0)

# correlation plot
smoothScatter(temps$observations[use], temps.bma.forecast[,2])
abline(0,1, col = "darkred", lwd = 2)

# Q-Q plot
plot(sort(temps$observations[use]), sort(temps.bma.forecast[,2]), pch = 20)
abline(0,1, col = "darkred", lwd = 2)

####################################################################################################

# ENSEMBLE MOS                                                                                  ####

# fit EMOS model
temps.mos <- ensembleMOS(temps, model = "normal", trainingDays = tr,
                         dates = timestamps[timestamps >= timestamps[tr]])

# plot probability of freezing
mosFreeze <- cdf(temps.mos, temps, value = 0)
plot(mosFreeze, type = "l")

# quantile forecasts
temps.mos.forecast <- quantileForecast(temps.mos, temps, quantiles = c( .1, .5, .9))
plot(obs["temp.n", ,1], type = "l", lwd = 2)
lines(25:90, temps.mos.forecast[1:66,2], col = "red", lwd = 2)
lines(25:90, temps.mos.forecast[1:66,1], col = "blue")
lines(25:90, temps.mos.forecast[1:66,3], col = "blue")

fc <- matrix(temps.mos.forecast[,"0.5"], ncol = 2)
matplot(fc, type = "l", lty = 1)

# PIT histogram to assess calibration of fitted model
hist(pit(temps.mos, temps), breaks = "fd", col = "skyblue",
     main = "PIT histogram of EMOS model for temps", xlab = "", ylab = "", prob = T)
abline(1,0, lty = 2)

# CRPS
mos.crps <- crps(temps.mos, temps)
mean(mos.crps)

# FIT MODELS TO ALL LEADTIMES (FOR MORE DIRECT COMPARISON)                                      ####

tr <- 25            # length of training set
timestamps <- gsub("-", "", as.Date(load.data("./Data/ECMWF_europe_starttime.rda")[16:105,1:7]/ 24, "1800-01-01"))

fc.temps <- function(model, lt) {
    c(apply(apply(offset.forecast(model[1:2,,,,-1]), 1:4, mean)[,,,toString(lt)], 1, cbind))
}

fitted.mos <- fitted.bma <- list()
    
# create ensemble in correct format (~5m / 3 leadtimes)
invisible(lapply(11:14, function(lt) {
    lt <- toString(lt)
    temps <- ensembleData(forecasts = data.frame("ecmwf" = fc.temps(ecmwf, lt),
                                                 "ncep" = fc.temps(ncep, lt),
                                                 "ukmo" = fc.temps(ukmo, lt)),
                          dates = rep(timestamps, 2),
                          observations = c(obs["temp.n",,], obs["temp.s",,]),
                          forecastHour = 00,
                          station = rep(c("N", "S"), each = 630),
                          latitude = rep(c(56, 52), each = 630),
                          longitude = rep(c(-3, 0), each = 630),
                          initializationTime = "00")
    
    fitted.mos[[lt]] <<- ensembleMOS(temps, model = "normal", trainingDays = tr,
                                    dates = timestamps[timestamps >= timestamps[tr]])
    fitted.bma[[lt]] <<- ensembleBMA(temps, model = "normal", trainingDays = tr,
                                     dates = timestamps[timestamps >= timestamps[tr]])
}))

saveRDS(fitted.bma, "./Data/ensBMA-temps.rds")
saveRDS(fitted.mos, "./Data/ensMOS-temps.rds")

# full comparison of models will use:
    # RMSE
    # spread
    # CRPS
    # verification histogram (PIT / rank probability)
    # Brier score for probability of freezing
    # Q-Q plots against observations

mos.brier <- brierScore(quantileForecast(temps.bma, temps, quantiles = c( .1, .5, .9))fitted.mos[[1]], temps, 0)

#

####################################################################################################

# EXTRACT FORECASTS FROM ALL MODELS                                                             ####

tr <- 25            # length of training set
timestamps <- gsub("-", "", as.Date(load.data("./Data/ECMWF_europe_starttime.rda")[16:105,1:7]/ 24, "1800-01-01"))

fc.temps <- function(model, lt) {
    c(apply(apply(offset.forecast(model[1:2,,,,-1]), 1:4, mean)[,,,toString(lt)], 1, cbind))
}

mos.fc <- bma.fc <- ens.data <- list()

# create ensemble in correct format (~5m / 3 leadtimes)
invisible(lapply(0:14, function(lt) {
    lt <- toString(lt)
    temps <- ensembleData(forecasts = data.frame("ecmwf" = fc.temps(ecmwf, lt),
                                                 "ncep" = fc.temps(ncep, lt),
                                                 "ukmo" = fc.temps(ukmo, lt)),
                          dates = rep(timestamps, 2),
                          observations = c(obs["temp.n",,], obs["temp.s",,]),
                          forecastHour = 00,
                          station = rep(c("N", "S"), each = 630),
                          latitude = rep(c(56, 52), each = 630),
                          longitude = rep(c(-3, 0), each = 630),
                          initializationTime = "00")
    
    qf.mos <- quantileForecast(fitted.mos[[lt]], temps, c(0.05, 0.5, 0.95)) 
    qf.bma <- quantileForecast(fitted.bma[[lt]], temps, c(0.05, 0.5, 0.95)) 
    
    mos.fc[[lt]] <<- array(qf.mos, dim = c(606, 2, 3), 
                           dimnames = list(NULL, c("temp.n", "temp.s"), dimnames(qf.mos)[[2]]))
    bma.fc[[lt]] <<- array(qf.bma, dim = c(606, 2, 3), 
                           dimnames = list(NULL, c("temp.n", "temp.s"), dimnames(qf.bma)[[2]]))
    
    ens.data[[lt]] <<- temps
}))

mos.fc <- aperm(abind(mos.fc, rev.along = 0), c(2,3,1,4))
bma.fc <- aperm(abind(bma.fc, rev.along = 0), c(2,3,1,4))

saveRDS(mos.fc, "./Data/ensMOS-fc.rds")
saveRDS(bma.fc, "./Data/ensBMA-fc.rds")
saveRDS(ens.data, "./Data/Ens-data.rds")


mosmos <- aperm(abind(mos.fc, rev.along = 0), c(2,3,1,4))
