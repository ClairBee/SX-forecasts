
library("CB.Misc"); library("SX.weather")

ecmwf <- readRDS("./Data/ECMWF-forecasts.rds")
obs <- readRDS("./Data/Observations.rds")

# plot improving forecast over leadtime
invisible(sapply(14:0, function(lt) {
    plot(obs["temp.n", , "09"], type = "l", lwd = 2, main = paste0("Leadtime: ", lt, " days"))
    lines(ecmwf["temp.n",(16:105)-lt,"09",toString(lt),"c"], col = "blue")
}))

# so, what needs to be done:
    # . correct load function to stop cropping the first 15 days
    # o correct error function to include offset
    # o check if anything bypasses that error calculation function

####################################################################################################

# start with manual calculation of forecast errors at each leadtime, for reference

err.os <- array(dim = dim(ecmwf[,16:105,,,]), dimnames = dimnames(ecmwf[,16:105,,,]))

err.os[,,,"14",] <- sweep(ecmwf[,(16:105) - 14,,"14",], 1:3, obs, "-")

tmp <- ecmwf["temp.n",(16:105)-14,"08","14","c"] - obs["temp.n",,"08"]

forecast.errors <- function(fc, actual = obs) {
    
    err <- array(dim = dim(fc[,16:105,,,]), dimnames = dimnames(fc[,16:105,,,]))
    
    invisible(sapply(0:14, function(lt) {
        err[,,,toString(lt),] <<- sweep(fc[,(16:105) - lt,,toString(lt),], 1:3, actual, "-")
    }))
    
    err
}

plot(err["temp.n",,"08","0","c"], type = "l")
lt <- 1

fce <- forecast.errors(ecmwf)
sum(is.na(fce))
