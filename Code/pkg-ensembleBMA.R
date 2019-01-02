
require(ensembleMOS); require(ensembleBMA); require(fields); require(maps)
library("SX.weather")

####################################################################################################

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
timestamps <- gsub("-", "", as.Date(load.data("../Data/ECMWF_europe_starttime.rda")[16:105,1:7]/ 24, "1800-01-01"))

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

# fit BMA model (~1.5 m to fit single model for all time points, single leadtime)
system.time(temps.bma <- ensembleBMA(temps, model = "normal", trainingDays = tr,
                         dates = timestamps[timestamps >= timestamps[tr]]))
                             
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

####################################################################################################

# ENS.BMA WITH CONTROLS & EXCHANGEABILITY                                                       ####

tr <- 25            # length of training set
lt <- "5"           # leadtime
timestamps <- gsub("-", "", as.Date(load.data("./Data/ECMWF_europe_starttime.rda")[16:105,1:7]/ 24, "1800-01-01"))
fc.mean <- function(model, lt) {
    c(apply(apply(offset.forecast(model[1:2,,,,-1]), 1:4, mean)[,,,toString(lt)], 1, cbind))
}
fc.ctrl <- function(model, lt) {
    c(apply(offset.forecast(model)[1:2,,,toString(lt),1], 1, rbind))
}

# run model across ensemble means (no controls)
{
    temps <- ensembleData(forecasts = data.frame("ecmwf" = fc.mean(ecmwf, lt),
                                                 "ncep" = fc.mean(ncep, lt),
                                                 "ukmo" = fc.mean(ukmo, lt)),
                          dates = rep(timestamps, 2),
                          observations = c(obs["temp.n",,], obs["temp.s",,]),
                          forecastHour = 00,
                          station = rep(c("N", "S"), each = 630),
                          latitude = rep(c(56, 52), each = 630),
                          longitude = rep(c(-3, 0), each = 630),
                          initializationTime = "00")
    
    # fit BMA model
    temps.bma <- ensembleBMA(temps, model = "normal", trainingDays = tr,
                                         dates = timestamps[timestamps >= timestamps[tr]])
}

# now, model with controls, no exchangeabiility (should fit same model anyway)
{
    temps.ctrl <- ensembleData(forecasts = data.frame("ecmwf.ens" = fc.mean(ecmwf, lt),
                                                      "ecmwf.ctrl" = fc.ctrl(ecmwf, lt),
                                                      "ncep.ens" = fc.mean(ncep, lt),
                                                      "ncep.ctrl" = fc.ctrl(ncep, lt),
                                                      "ukmo.ens" = fc.mean(ukmo, lt),
                                                      "ukmo.ctrl" = fc.ctrl(ukmo, lt)),
                          dates = rep(timestamps, 2),
                          observations = c(obs["temp.n",,], obs["temp.s",,]),
                          forecastHour = 00,
                          station = rep(c("N", "S"), each = 630),
                          latitude = rep(c(56, 52), each = 630),
                          longitude = rep(c(-3, 0), each = 630),
                          initializationTime = "00")
    
    system.time({
        bma.ctrl <- ensembleBMA(temps.ctrl, model = "normal", trainingDays = tr,
                                         dates = timestamps[timestamps >= timestamps[tr]])})
}

# and finally, model with controls, exchangeability specified
{
    temps.ctrl.exch <- ensembleData(forecasts = data.frame("ecmwf.ens" = fc.mean(ecmwf, lt),
                                                      "ecmwf.ctrl" = fc.ctrl(ecmwf, lt),
                                                      "ncep.ens" = fc.mean(ncep, lt),
                                                      "ncep.ctrl" = fc.ctrl(ncep, lt),
                                                      "ukmo.ens" = fc.mean(ukmo, lt),
                                                      "ukmo.ctrl" = fc.ctrl(ukmo, lt)),
                               dates = rep(timestamps, 2),
                               observations = c(obs["temp.n",,], obs["temp.s",,]),
                               forecastHour = 00,
                               station = rep(c("N", "S"), each = 630),
                               latitude = rep(c(56, 52), each = 630),
                               longitude = rep(c(-3, 0), each = 630),
                               initializationTime = "00",
                               exchangeable = c(ecmwf.ens = 1, ecmwf.ctrl = 1,
                                                ncep.ens = 2, ncep.ctrl = 2,
                                                ukmo.ens = 3, ukmo.ctrl = 3))
    system.time({
        bma.ctrl.exch <- ensembleBMA(temps.ctrl.exch, model = "normal", trainingDays = tr,
                            dates = timestamps[timestamps >= timestamps[tr]])})
}

# ~1.5min to fit single model to 3 ensemble means for all time points, single leadtime
# ~3.5min to fit single model to 3 ensemble means & 3 controls, non-exchangeable
# ~2min to fit single model to 3 ensemble means & 3 controls, exchangeable

saveRDS(temps.bma, "./Models/BMA-ens-only.rds")
saveRDS(bma.ctrl, "./Models/BMA-ctrl-nex.rds")
saveRDS(bma.ctrl.exch, "./Models/BMA-ctrl-ex.rds")

saveRDS(temps, "./Models/ensdata-ens-only.rds")
saveRDS(temps.ctrl, "./Models/ensdata-ctrl-nex.rds")
saveRDS(temps.ctrl.exch, "./Models/ensdata-ctrl-ex.rds")

# fit to all ensemble members & controls, with exchangeability
# v.slow with all ensemble members included. DON'T EVER DO THIS AGAIN.
{
#    exch.temps <- function(lt) {
#        
#        ec <- setNames(data.frame(apply(apply(offset.forecast(ecmwf)[1:2,,,toString(lt),],
#                                              c(1, 4), rbind), 3, rbind)),
#                       paste("ec", dimnames(ecmwf)[[5]], sep = "."))
#        nc <- setNames(data.frame(apply(apply(offset.forecast(ncep)[1:2,,,toString(lt),],
#                                              c(1, 4), rbind), 3, rbind)),
#                       paste("nc", dimnames(ncep)[[5]], sep = "."))
#        mo <- setNames(data.frame(apply(apply(offset.forecast(ukmo)[1:2,,,toString(lt),],
#                                              c(1, 4), rbind), 3, rbind)),
#                       paste("uk", dimnames(ukmo)[[5]], sep = "."))
#        
#        temps <- data.frame(ec, mo, nc, stringsAsFactors = F)
#        
#        exch <- c(1, rep(2, ncol(ec)-1), 
#                  3, rep(4, ncol(nc)-1), 
#                  5, rep(6, ncol(mo)-1))
#        names(exch) <- colnames(temps#
#
#        return(list(fc = temps, exch = exch))
#    }
#    
#    # create ensemble in correct format
#    exch.data <- ensembleData(forecasts = exch.temps(lt)$fc,
#                              exchangeable = exch.temps(lt)$exch,
#                              dates = rep(timestamps, 2),
#                              observations = c(obs["temp.n",,], obs["temp.s",,]),
#                              forecastHour = 00,
#                              station = rep(c("N", "S"), each = 630),
#                              latitude = rep(c(56, 52), each = 630),
#                              longitude = rep(c(-3, 0), each = 630),
#                              initializationTime = "00")
#    
#    # takes far too long to run for all dates, single leadtime. 
#    # (stopped after ~10m, had only made it to 2009...)
#    exch.bma <- ensembleBMA(exch.data, model = "normal", trainingDays = tr,
#                                dates = timestamps[timestamps >= timestamps[tr]])
}

####################################################################################################

# VERIFICATION OF ENS. MODELS WITH CONTROLS                                                     ####

# load all necessary data
{
    tr <- 25            # length of training set
    lt <- "5"           # leadtime
    timestamps <- gsub("-", "", as.Date(load.data("./Data/ECMWF_europe_starttime.rda")[16:105,1:7]/ 24, "1800-01-01"))
    
    bma.pert <- readRDS("./Models/BMA-ens-only.rds")
    bma.ctrl.nex <- readRDS("./Models/BMA-ctrl-nex.rds")
    bma.ctrl.ex <- readRDS("./Models/BMA-ctrl-ex.rds")
    
    ens.pert <- readRDS("./Models/ensdata-ens-only.rds")
    ens.ctrl.nex <- readRDS("./Models/ensdata-ctrl-nex.rds")
    ens.ctrl.ex <- readRDS("./Models/ensdata-ctrl-ex.rds")
}

# plot bias coefficients:
# same for ens-only & non-ex, but exchangeability changes bias coeffs
# sigma differs for all models
{
    matplot(t(bma.pert$biasCoefs[1,,]), type = "l", col = c("black", "red", "blue"), lty = 1,
            xlab = "", ylab = "", main = "Bias coefficient 1")
    abline(0,0,col = "grey", lty = 2)
    
    matplot(t(bma.ctrl.nex$biasCoefs[1,c(1,3,5),]), type = "l", col = c("green3", "gold", "cyan3"), lty = 2,
            xlab = "", ylab = "", main = "Bias coefficient 1", add = T)
    abline(0,0,col = "grey", lty = 2)
    
    matplot(t(bma.ctrl.ex$biasCoefs[1,c(1,3,5),]), type = "l", col = c("green3", "gold", "cyan3"), lty = 2,
            xlab = "", ylab = "", main = "Bias coefficient 1", add = T)
    abline(0,0,col = "grey", lty = 2)
    
    matplot(t(bma.pert$biasCoefs[2,,]), type = "l", col = c("black", "red", "blue"), lty = 1,
            xlab = "", ylab = "", main = "Bias coefficient 2")
    abline(1,0,col = "grey", lty = 2)
    
    matplot(t(bma.ctrl.nex$biasCoefs[2,c(1,3,5),]), type = "l", col = c("green3", "gold", "cyan3"), lty = 2,
            xlab = "", ylab = "", main = "Bias coefficient 2", add = T)
    abline(0,0,col = "grey", lty = 2)
    
    plot(bma.pert$sd, type = "l", xlab = "", ylab = "", main = "Sigma")
    lines(bma.ctrl.ex$sd, col = "blue")
    lines(bma.ctrl.nex$sd, col = "red")
}

# check weights: how do they differ?
# with exchangeability (blue), weight is split 50-50 between control & ens.mean, so max is 0.5
# black & green usually similar (ens. mean without/with ctrl) - ECMWF less so
{
    plot(bma.pert$weights["ecmwf",], type = "l", xlab = "", ylab = "", main = "ECMWF weights", ylim = c(0,1))
    matplot(t(bma.ctrl.ex$weights[c("ecmwf.ens", "ecmwf.ctrl"),]), type = "l", col = "blue", add = T)
    lines(bma.ctrl.nex$weights["ecmwf.ens",], col = "green3")
    lines(bma.ctrl.nex$weights["ecmwf.ctrl",], col = "magenta3", lty = 2, lwd = 2)

    plot(bma.pert$weights["ncep",], type = "l", xlab = "", ylab = "", main = "NCEP weights", ylim = c(0,1))
    matplot(t(bma.ctrl.ex$weights[c("ncep.ens", "ncep.ctrl"),]), type = "l", col = "blue", add = T)
    lines(bma.ctrl.nex$weights["ncep.ens",], col = "green3")
    lines(bma.ctrl.nex$weights["ncep.ctrl",], col = "magenta3", lty = 2, lwd = 2)
    
    plot(bma.pert$weights["ukmo",], type = "l", xlab = "", ylab = "", main = "UKMO weights", ylim = c(0,1))
    matplot(t(bma.ctrl.ex$weights[c("ukmo.ens", "ukmo.ctrl"),]), type = "l", col = "blue", add = T)
    lines(bma.ctrl.nex$weights["ukmo.ens",], col = "green3")
    lines(bma.ctrl.nex$weights["ukmo.ctrl",], col = "magenta3", lty = 2, lwd = 2)
}

matplot(t(bma.ctrl.nex$weights), type = "l", col = rep(c("green3", "red", "blue"), each = 2), lty = rep(c(1,2), 3))

# model performance: PIT
{
    hist(pit(bma.pert, ens.pert), breaks = c(0:10)/10, prob = T, col = "skyblue", xlab = "", ylab = "",
         main = paste0("PIT histogram: ensemble means only, LT ", lt), ylim = c(0,1.2))
    abline(1,0, col = "blue", lty = 2)
    
    hist(pit(bma.ctrl.ex, ens.ctrl.ex), breaks = c(0:10)/10, prob = T, col = "coral2", xlab = "", ylab = "",
         main = paste0("PIT histogram: controls & ensemble means exchangeable, LT ", lt), ylim = c(0,1.2))
    abline(1,0, col = "red3", lty = 2)
    
    hist(pit(bma.ctrl.nex, ens.ctrl.nex), breaks = c(0:10)/10, prob = T, col = "lightseagreen", xlab = "", ylab = "",
         main = paste0("PIT histogram: controls & ensemble means non-exchangeable, LT ", lt), ylim = c(0,1.2))
    abline(1,0, col = "seagreen3", lty = 2)
}

plot(pit(bma.ctrl.ex, ens.ctrl.ex), pit(bma.ctrl.nex, ens.ctrl.nex))

bc <- bma.pert$biasCoefs[,"ecmwf",1]
plot(apply(offset.forecast(ecmwf)[1,1:25,1,"5",-1], 1, mean), obs[1,1:25,1], pch = 20,
     xlab = "Forecast", ylab = "Observation", main = "Bias adjustment", asp = T)
points(apply(offset.forecast(ecmwf)[2,1:25,1,"5",-1], 1, mean), obs[2,1:25,1], pch = 20, col = "red")

abline(0,1, col = "green3")                                 # target
abline(bc, col = "skyblue")     # fitted line

####################################################################################################

# SINGLE-VARIABLE ENSEMBLE BMA                                                                  ####

tr <- 25            # length of training set
lt <- "5"           # leadtime
timestamps <- gsub("-", "", as.Date(load.data("./Data/ECMWF_europe_starttime.rda")[16:105,1:7]/ 24, "1800-01-01"))

# simplest model possible - should be able to replicate coefficients etc
temp.n <- ensembleData(forecasts = data.frame("ecmwf" = c(apply(offset.forecast(ecmwf)[1,,,toString(lt),-1], 1:2, mean)),
                                             "ncep" = c(apply(offset.forecast(ncep)[1,,,toString(lt),-1], 1:2, mean)),
                                             "ukmo" = c(apply(offset.forecast(ukmo)[1,,,toString(lt),-1], 1:2, mean))),
                      dates = timestamps,
                      observations = c(obs["temp.n",,]),
                      forecastHour = 00,
                      initializationTime = "00")

train.n <- trainingData(temp.n, date = "20071225", trainingDays = tr)
bma.n <- fitBMA(train.n, model = "normal")
     
# bias-correction calculation    
{
    bc <- bma.n$biasCoefs
    
    plot(train.n$ecmwf, train.n$obs, pch = 20, xlab = "forecast", ylab = "obs", main = "Bias correction - ECMWF")
    abline(lm(obs ~ fc, data.frame(fc = train.n$ecmwf, obs = train.n$obs)), col = "red3")
    points(bc[1,"ecmwf"] + bc[2,"ecmwf"] * train.n$ecmwf, train.n$obs, col = "red")
    lines(matrix(c(t(cbind(train.n$ecmwf, train.n$obs,
                           bc[1,"ecmwf"] + bc[2,"ecmwf"] * train.n$ecmwf, train.n$obs, NA, NA))),
                 ncol = 2, byrow = T), col = "red3", lty = 2)
    
    plot(train.n$ncep, train.n$obs, pch = 20, xlab = "forecast", ylab = "obs", main = "Bias correction - NCEP")
    abline(lm(obs ~ fc, data.frame(fc = train.n$ncep, obs = train.n$obs)), col = "red3")
    points(bc[1,"ncep"] + bc[2,"ncep"] * train.n$ncep, train.n$obs, col = "red")
    lines(matrix(c(t(cbind(train.n$ncep, train.n$obs,
                           bc[1,"ncep"] + bc[2,"ncep"] * train.n$ncep, train.n$obs, NA, NA))),
                 ncol = 2, byrow = T), col = "red3", lty = 2)
    
    plot(train.n$ukmo, train.n$obs, pch = 20, xlab = "forecast", ylab = "obs", main = "Bias correction - UKMO")
    abline(lm(obs ~ fc, data.frame(fc = train.n$ukmo, obs = train.n$obs)), col = "red3")
    points(bc[1,"ukmo"] + bc[2,"ukmo"] * train.n$ukmo, train.n$obs, col = "red")
    lines(matrix(c(t(cbind(train.n$ukmo, train.n$obs,
                           bc[1,"ukmo"] + bc[2,"ukmo"] * train.n$ukmo, train.n$obs, NA, NA))),
                 ncol = 2, byrow = T), col = "red3", lty = 2)
    
}

pit(bma.n, train.n)

plot((0:100)/10, bma.n$weights["ecmwf"] * dnorm((0:100)/10, bc[1,"ecmwf"] + bc[2, "ecmwf"] * train.n$"ecmwf"[25], bma.n$sd),
     type = "l", xlab = "", ylab = "", ylim = c(0,0.5), main = "Densities on 2007-12-25", col = "red")
lines((0:100)/10, bma.n$weights["ncep"] * dnorm((0:100)/10, bc[1,"ncep"] + bc[2, "ncep"] * train.n$"ncep"[25], bma.n$sd),
      col = "green3")
lines((0:100)/10, bma.n$weights["ukmo"] * dnorm((0:100)/10, bc[1,"ukmo"] + bc[2, "ukmo"] * train.n$"ukmo"[25], bma.n$sd),
      col = "blue")
lines((0:100)/10, bma.n$weights["ecmwf"] * dnorm((0:100)/10, bc[1,"ecmwf"] + bc[2, "ecmwf"] * train.n$"ecmwf"[25],
                                                 bma.n$sd) + 
          bma.n$weights["ncep"] * dnorm((0:100)/10, bc[1,"ncep"] + bc[2, "ncep"] * train.n$"ncep"[25],
                                         bma.n$sd) + 
          bma.n$weights["ukmo"] * dnorm((0:100)/10, bc[1,"ukmo"] + bc[2, "ukmo"] * train.n$"ukmo"[25],
                                                                                    bma.n$sd), lwd = 2)
abline(v = train.n$obs[25], col = "orange", lwd = 2, lty = 2)
abline(v = quantileForecast(bma.n, train.n, c(0.05, 0.95))[25,], lty = 2)

####################################################################################################

# MANUAL ENSEMBLE BMA                                                                           ####

lt <- "3"           # leadtime

# run automatic fit for comparison
{
    tr <- 25            # length of training set
    timestamps <- gsub("-", "", as.Date(load.data("./Data/ECMWF_europe_starttime.rda")[16:105,1:7]/ 24, "1800-01-01"))
    
    temp.n <- ensembleData(forecasts = data.frame("ecmwf" = c(apply(offset.forecast(ecmwf)[1,,,toString(lt),-1], 1:2, mean)),
                                                  "ncep" = c(apply(offset.forecast(ncep)[1,,,toString(lt),-1], 1:2, mean)),
                                                  "ukmo" = c(apply(offset.forecast(ukmo)[1,,,toString(lt),-1], 1:2, mean))),
                           dates = timestamps,
                           observations = c(obs["temp.n",,]),
                           forecastHour = 00,
                           initializationTime = "00")
    
    train.n <- trainingData(temp.n, date = "20071225", trainingDays = tr)
    bma.n <- fitBMA(train.n, model = "normal")
    bma.temps <- ensembleBMA(temp.n, model = "normal", trainingDays = tr,
                             dates = timestamps[timestamps >= timestamps[tr]])
}

# simplest case: single variable, single day
o <- obs[1,d-(24:0),y]
fc <- abind("ecmwf" = apply(offset.forecast(ecmwf)[1,d-(24:0),y,lt,-1], 1, mean),
               "ncep" = apply(offset.forecast(ncep)[1,d-(24:0),y,lt,-1], 1, mean),
               "ukmo" = apply(offset.forecast(ukmo)[1,d-(24:0),y,lt,-1], 1, mean),
            along = 0)

# linear regression over training data
bc <- abind(invisible(sapply(dimnames(fc)[[1]], function(model) {
    lm(obs ~ fc, data.frame(obs = o, fc = fc[model,]))$coef
}, simplify = F)), along = 0)

all(t(bma.n$biasCoefs) == bc)           # all fine so far

# sigma & weights are estimated using EM algorithm...
{
    # support function to evaluate log-likelihood
    log.likelihood <- function(sig2, w) {
        
        # summed over all observations in training set
        l <- array(dim = dim(fc), dimnames = dimnames(fc))
        
        invisible(sapply(1:3, function(i) {
            invisible(sapply(1:25, function(t) {
                l[i,t] <<- w[i] * dnorm(o[t], bc[i,1] + bc[i,2] * fc[i,t], sqrt(sig2))
            }))}))
        sum(log(colSums(l)))
    }
    
    # function to perform maximisation
    EM <- function(sig2 = 4, w = rep(1/3, 3), max.runs = 1000, conv = 0.0001) {
        
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
    
    # Support function - plot mixture model from EM algorithm 
    plot.em <- function(modelToPlot, xlab = "", ylab = "", ...) {
        
        x <- range(pretty(o))
        x <- x[1]:(x[2]*10) / 10
        
        w <- modelToPlot$w
        sig <- modelToPlot$sig
        l <- array(dim = c(3, length(x)))
        
        invisible(sapply(1:3, function(i) {
                l[i,] <<- w[i] * dnorm(x, bc[i,1] + bc[i,2] * fc[i,25], sig)
            }))
        
        matplot(x, t(l), type = "l", ylim = range(0, colSums(l)), col = c("green3", "red", "blue"), 
                lty = 1, xlab = xlab, ylab = ylab,...)
        
        # add mixture model
        lines(x, colSums(l), lwd = 2)
        
        abline(v = o[25], col = "orange", lwd = 2, lty = 2)
        
        legend("topright", bty = "n", col = c("green3", "red", "blue", NA, "black"), lwd = c(1,1,1,2, 2),
               legend = c("ECMWF", "NCEP", "UKMO", NA, "BMA mixture"))
    }
}   # EM functions

em.fit <- EM(sig2 = 1)
plot.em(em.fit, xlim = c(0.6, 6))

plot(bma.temps, temp.n, dates = "20071225", ask = F, xlim = c(0,6))

####################################################################################################

# BIAS-CORRECTION PLOTS                                                                         ####

# where does the power of ensemble BMA come from? Suspect most of the benefit is from bias correction!

fitted.bma <- readRDS("./Models/ensBMA-temps.rds")
ens.data <- readRDS("./Models/Ens-data.rds")

plot.bma.adj <- function(lt, d) {
    
    lt <- toString(lt)
    
    # find index of day specified
    i <- which(array(1:630, dim = c(90,7)) == d + 24, arr.ind = T)
    
    px <- rbind("o" = obs[1:2,i[1],i[2]],
                "ecmwf" = apply(offset.forecast(ecmwf)[1:2,i[1],i[2],lt,-1], 1, mean),
                "ncep" = apply(offset.forecast(ncep)[1:2,i[1],i[2],lt,-1], 1, mean),
                "ukmo" = apply(offset.forecast(ukmo)[1:2,i[1],i[2],lt,-1], 1, mean))
    
    bc <- rbind(o = c(0,1), t(fitted.bma[[lt]]$biasCoefs[,,d]))
    adj <-  (bc[,1] + sweep(px, 1, bc[,2], "*"))[2:4,]
   
    w <- fitted.bma[[lt]]$weights[,d]
    wm <- apply(sweep(adj, 1, w, "*"), 2, sum)
    
    qf <- quantileForecast(fitted.bma[[lt]], ens.data[[lt]], 0.5)[d + c(0,606),]
    
    plot(rbind(wm, qf, px, adj), pch = c(4, 0, 16, rep(c(16, 1), each = 3)), 
         col = c(rep("black", 3), rep(c("red", "green3", "blue"), 2)),
         main = paste0("Forecast adjustments - day ", d+24, ", LT ", lt))
    
    lines(rbind(px["ecmwf",], adj["ecmwf",], wm), col = "red", lty = 2)
    lines(rbind(px["ncep",], adj["ncep",], wm), col = "green3", lty = 2)
    lines(rbind(px["ukmo",], adj["ukmo",], wm), col = "blue", lty = 2)
    lines(rbind(wm, px["o",]))
}

plot.bma.adj(14, 60)

# how often is weighted mean closer to observation than single best BC forecast?
# weighted mean further away => ensemble biased

# how often is quantile forecast closer to observation than single best UNCORRECTED forecast?
e <- aperm(array(abind(array(NA, dim = c(24, 2, 15)), 
                       array(invisible(sapply(formatC(0:14), 
                                              function(lt) quantileForecast(fitted.bma[[lt]], ens.data[[lt]], 0.5))),
                             dim = c(606, 2, 15)), along = 1), dim = c(90,7,2, 15)), c(3,1,2,4))

err <- abind("ecmwf" = forecast.errors(apply(ecmwf[,,,,-1], 1:4, mean))[1:2,,,],
             "ncep" = forecast.errors(apply(ncep[,,,,-1], 1:4, mean))[1:2,,,],
             "ukmo" = forecast.errors(apply(ukmo[,,,,-1], 1:4, mean))[1:2,,,],
             "bma" = e, along = 0)

err.dist <- apply(err, c(1, 3:5), function(v) sqrt(sum(v^2)))

min.dist <- apply(err.dist, 2:4, which.min)
(sum(min.dist == 4) / length(min.dist)) * 100

####################################################################################################

# VERIFICATION RANKS OF EACH ENSEMBLE                                                           ####

# is there a bias trend over time? Plot TS of verification ranks to investigate. 


