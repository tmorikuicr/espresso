## 1.0.0
* Add an option `version` to `graphSOM` and `rxmcmc` functions in order to enable us to perform SOM computations based on the previous version (= 0.2.17).

## 0.2.17
* Add options `size_d` and `size_s` to `plotDistMap` function.
* Add options `color_by_group`, `grpcol`, and `magnify` to `plotUMAP` function.

## 0.2.16
* Replace graphSOM and rxmcmc functions of v.0.2.15 with those of v.0.2.14

## 0.2.15
* Modify graphSOM function to guarantee the reproducibility between graphSOM and rxmcmc when the same seed value is given.
* Add a new argument `color_by_gene` to the function `plotUMAP`.


## 0.2.14
* Add `rglclose` option to `plotUMAP` function.

## 0.2.13
* The `umap` package for `plotUMAP` function was replaced with `uwot`.

## 0.2.12
* Update `plotUMAP` function. It supports 'png', 'pdf', and 'eps'.

## 0.2.11
* Add 2D plotting option to `plotUMAP`.
* Add scatterplot3d option to `plotUMAP`.

## 0.2.10
* The output file names of `writeMCMC` function are changed.
* The `gene_size` column is added to `summary.txt`.
* Bugs are fixed.

## 0.2.9
* Gibbs sampling strategy is introduced to the cluster swapping in GraphSOM.

## 0.2.8
* The `swap` option is added to `initGraphSOM` function.
* Thee `seed` options are added to `selectFeatures`, `initGraphSOM`, and `rxmcmc` functions.
* Bugs are fixed.


## 0.2.7
* `selectFeatures()` outputs NOT rejected genes.
* Modify arguments of `plotUMAP` function.
* Modify outputs format of RX-MCMC. 

## 0.2.6
* Bugs are fixed

## espresso 0.2.4
* The default value of the parameter `maxRuns` was changed from 1000 to 500 in `selectFeatures` function.
* Rename variables `rep` to `rept` or `repl`.
* Add `max_` and `min_` for `score`, `acc`, and `ari` to MCMC results.
* Rename y-axis in MCMC plots.
* Rename `size` to `gene_size` in MCMC results.


## 0.2.3
* The default value of the parameter `maxRuns` was changed from 200 to 1000 in `selectFeatures` function.
* The default value of the parameter `coef` was changed from 0 to 1.0 in `initGraphSOM` function.

## 0.2.2
* Bugs are fixed

## 0.2.1
* Bugs are fixed

## 0.0.2

* Add `rxmcmc` function, which optimize gene sets by replica exchange MCMC approach.
* Add `plotMCMC` function, which plots RX-MCMC results.
* Add `selectFeatures` funtion, which selects feature genes by random forest based method `Boruta`.
* Add `plotDistMap`function, which plots distance matrices between domains and between samples.

## 0.0.1

* A beta-version of `espresso` is released.