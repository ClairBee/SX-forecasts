
# quick reference - how to obtain full file set directly from pre-processed forecasts

# ensemble mean forecasts
ens.mean <- apply(offset.forecast(superensemble())[,,,,-(1:3)], 1:4, mean)
{
    # mean of all 93 ensemble members for each day (offset-corrected) - check
    all(array(ens.mean["temp.s", , ,2:15], dim = c(630, 14)) == TS.pert.leadtime.ensemblemean_ALL,
        array(ens.mean["temp.n", , ,2:15], dim = c(630, 14)) == TN.pert.leadtime.ensemblemean_ALL,
        array(ens.mean["pc1", , ,2:15], dim = c(630, 14)) == PC1.pert.leadtime.ensemblemean_ALL,
        array(ens.mean["pc2", , ,2:15], dim = c(630, 14)) == PC2.pert.leadtime.ensemblemean_ALL)
}

# ensemble mean errors
# check that ensemble mean forecast errors match (can therefore bypass the above)
ens.mean.error <- apply(forecast.errors(superensemble())[,,,,-(1:3)], 1:4, mean)
{
    all(round(aperm(ens.mean.error["pc1",,,2:15], c(2,1,3)), 10) == round(mean.error_PC1_ALL,10),
        round(aperm(ens.mean.error["pc2",,,2:15], c(2,1,3)), 10) == round(mean.error_PC2_ALL,10),
        round(aperm(ens.mean.error["temp.n",,,2:15], c(2,1,3)), 10) == round(mean.error_TN_ALL,10),
        round(aperm(ens.mean.error["temp.s",,,2:15], c(2,1,3)), 10) == round(mean.error_TS_ALL,10))
}

TN.pert.leadtime.ensemblemean_ECMWF <- apply(apply(offset.forecast(ecmwf)["temp.n",,,2:15,-1], 1:3, mean), 3, rbind)

mean_TN <- abind(lapply(list(ecmwf, ncep, ukmo),
                        function (model) {
                            apply(apply(offset.forecast(model)["temp.n",,,2:15,-1], 1:3, mean), 3, rbind)
                        }), along = 0)