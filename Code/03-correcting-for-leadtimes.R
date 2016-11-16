
library("CB.Misc"); library("SX.weather")

ecmwf <- readRDS("./Data/ECMWF-forecasts.rds")

# create table of observations with leadtime offset

err.0 <- ecmwf["temp.n",,"08","0",] - obs["temp.n",,"08"]

err.os <- invisible(lapply(1:15,
                 function(lt) {
                     sweep(ecmwf[,,,lt,], 5, obs, "-")
                 }))
    
    
    
PC1.pert.leadtime.ensemblemean_ECMWF <- matrix(nrow=630,ncol=14)
for (i in 1:630){
    for (j in 1:14){
        PC1.pert.leadtime.ensemblemean_ECMWF[i,j] <- mean(as.vector(ensemble_ts[djf.index2(model_time1)[i]-j,1,j*4+1,]))
    }
}