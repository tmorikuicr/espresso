# espresso 0.2.7
* `selectFeatures()` outputs NOT rejected genes.
* Modify arguments of `plotUMAP` function.
* Modify outputs format of RX-MCMC. 

# espresso 0.2.6
* Fix bugs

# espresso 0.2.4
* The default value of the parameter `maxRuns` was changed from 1000 to 500 in `selectFeatures` function.
* Rename variables `rep` to `rept` or `repl`.
* Add `max_` and `min_` for `score`, `acc`, and `ari` to MCMC results.
* Rename y-axis in MCMC plots.
* Rename `size` to `gene_size` in MCMC results.


# espresso 0.2.3
* The default value of the parameter `maxRuns` was changed from 200 to 1000 in `selectFeatures` function.
* The default value of the parameter `coef` was changed from 0 to 1.0 in `initGraphSOM` function.

# espresso 0.2.2
# Bug fixed

# espresso 0.2.1
* Bug fixed

# espresso 0.0.2

* Add `rxmcmc` function, which optimize gene sets by replica exchange MCMC approach.
* Add `plotMCMC` function, which plots RX-MCMC results.
* Add `selectFeatures` funtion, which selects feature genes by random forest based method `Boruta`.
* Add `plotDistMap`function, which plots distance matrices between domains and between samples.

# espresso 0.0.1

* A beta-version of `espresso` is released.