# =============================================================================
#' @title Initialize map by uniformly distributed random variables.
#' @description This function initializes map 
#'              by generating uniformly random values.
#' @param obj The \code{espresso} object
#' @importFrom stats runif
#' @return \code{espresso} object
#' 
.initMapByRandom <- function(obj) {
  n <- nrow(obj@topology)
  genes <- obj@gset
  newmap <- vector('list', length = obj@rept)
  names(newmap) <- paste0('rept.', 1:obj@rept)
  for (s in 1:obj@rept) {
    samples <- obj@ssets[[s]]
    exprs <- obj@exprs[match(samples, rownames(obj@exprs)), , drop = FALSE]
    exprs <- exprs[, match(genes, colnames(obj@exprs)), drop = FALSE]
    mi <- apply(exprs, 2, min)
    ma <- apply(exprs, 2, max)
    rept <- matrix(NA, n, length(genes))
    rownames(rept) <- rownames(obj@topology)
    colnames(rept) <- genes
    for (i in 1:n) {
      rept[i, ] <- mi + (ma - mi) * runif(ncol(rept))
    }
    newmap[[s]] <- rept
  }
  obj@map <- newmap
  return(obj)
}

# =============================================================================
#' @title Initialize map by input samples.
#' @description This function initializes map 
#'              by values of randomly selected samples.
#' @param obj The \code{espresso} object.
#' @return \code{espresso} object
#' 
.initMapBySamples <- function(obj) {
  n <- nrow(obj@topology)
  genes <- obj@gset
  newmap <- vector('list', length = obj@rept)
  names(newmap) <- paste0('rept.', 1:obj@rept)
  for (s in 1:obj@rept) {
    samples <- obj@ssets[[s]]
    exprs <- obj@exprs[match(samples, rownames(obj@exprs)), , drop = FALSE]
    exprs <- exprs[, match(genes, colnames(obj@exprs)), drop = FALSE]
    rept <- matrix(NA, n, length(genes))
    rownames(rept) <- rownames(obj@topology)
    colnames(rept) <- genes
    idx <- sample(1:nrow(exprs), size = n, replace = FALSE)
    rept[1:nrow(rept), ] <- exprs[idx, ]
    newmap[[s]] <- rept
  }
  obj@map <- newmap
  return(obj)
}

# =============================================================================
#' @title Initialize the best matching units (BMUs).
#' @description This function initializes BMUs.
#' @param obj The \code{espresso} object.
#' @return \code{espresso} object
#' 
.initBMU <- function(obj) {
  n <- nrow(obj@topology)
  genes <- obj@gset
  newbmu <- vector('list', length = obj@rept)
  names(newbmu) <- paste0('rept.', 1:obj@rept)
  for (s in 1:obj@rept) {
    samples <- obj@ssets[[s]]
    exprs <- obj@exprs[match(samples, rownames(obj@exprs)), , drop = FALSE]
    exprs <- exprs[, match(genes, colnames(obj@exprs)), drop = FALSE]
    bmu <- matrix(NA, nrow = nrow(exprs), ncol = obj@lsteps + 1)
    rownames(bmu) <- rownames(exprs)
    colnames(bmu) <- 0:obj@lsteps
    newbmu[[s]] <- bmu
  }
  obj@bmu <- newbmu
  return(obj)
}

# =============================================================================
#' @title Initialize an \code{espresso} object.
#' @description This function initializes an espresso object 
#'              for GraphSOM clustering.
#' @param obj The \code{espresso} object.
#' @param gset Gene set.
#' @param radius Initial value of learning radius 
#'               (default: the maximum distance of input toplology).
#' @param lsteps Number of learning steps.
#' @param stochastic If \code{stochastic = TRUE}, 
#'                   GraphSOM clustering introduces stochasticity 
#'                  for map update (default: TRUE).
#' @param rmin Minimum value of scale factor for the stochastic map update.
#' @param rmax Maximum value of scale factor for the stochastic map update.
#' @param map_method Map initialization method.
#'                   The map vector is initilaized 
#'                   by input data if 'sample' is selected. 
#'                   Otherwise, it is initialized by 
#'                   uniformly distributed random values (default: 'sample'). 
#' @param bmu_method BMU selection method.
#'                   If the unit adjacent to the BMU candidate is empty,
#'                   it is set as BMU ('fill'). 
#'                   Otherwise, the most similar unit to input data 
#'                   is selected as BMU (default: 'min').
#' @param coef Coefficient of ARI for score computation.
#' @param nmin Minimum number of genes can be analyzed.
#' @param nmax Maximum number of genes can be analyzed.
#' @param swap Logical value determins whether to swap domains or not (default: TRUE)
#' @return \code{espresso} object
#' @importFrom igraph graph.adjacency
#' @importFrom igraph shortest.paths
#' 
.initGraphSOM <- function(obj, gset, radius, lsteps, stochastic, rmin, rmax,
                          map_method, bmu_method, coef, nmin, nmax, swap) {
  obj@dist <- shortest.paths(graph.adjacency(obj@topology))
  if (is.null(radius)) {
    if (max(obj@dist) == Inf) {
      tmp <- sort(unique(as.vector(obj@dist)))
      obj@radius <- tmp[tmp != Inf][length(tmp) - 1]
    } else {
      obj@radius <- max(obj@dist)
    }
  } else {
    obj@radius <- radius
  }
  obj@lsteps <- as.integer(lsteps)
  obj@stochastic <- stochastic
  obj@rmin <- rmin
  obj@rmax <- rmax
  if (map_method == "random") {
    obj@map_method <- "random"
  } else if (map_method == "sample") {
    obj@map_method <- "sample"
  } else {
    warning(paste0("map_method = '", map_method, "' is invalid. The 'sample' is forced use."))
    obj@map_method <- "sample"
  }
  if (bmu_method == "min") {
    obj@bmu_method <- "min"
  } else if (bmu_method == "fill") {
    obj@bmu_method <- "fill"
  } else {
    warning(paste0("bmu_method = '", bmu_method, "' is invalid. The 'min' is forced use."))
    obj@map_method <- "min"
  }
  obj@coef <- coef
  obj@nmin <- as.integer(nmin)
  if (is.null(nmax)) {
    obj@nmax <- as.integer(ncol(obj@exprs))
  } else {
    obj@nmax <- as.integer(nmax)
  }
  if (swap == TRUE) {
    obj@swap <- TRUE
  } else {
    obj@swap <- FALSE
  }
  return(obj)
}

# =============================================================================
#' @title Search for the best matching units (BMU).
#' @description This function searches for the BMUs for each input sample.
#' @param obj The \code{espresso} object
#' @return BMUs given as a matrix.
#' 
.searchBMUByMin <- function(obj) {
  rept <- obj@rept_som
  bmu <- apply(obj@exprs_som, 1, 
               function(i, map = obj@map[[rept]]){
                 x <- matrix(i, nrow = nrow(obj@map[[rept]]), 
                             ncol = ncol(obj@map[[rept]]), byrow = TRUE)
                 d <- apply((x - map)^2, 1, sum)
                 return(which.min(d))
               }, 
               obj@map[[rept]])
  return(bmu)
}

# =============================================================================
#' @title Search for the best matching units (BMU) by 'fill' method.
#' @description This function searches the for BMUs for each input sample.
#'              In this function, deifferent from \code{.searchBMUByMin}, 
#'              if the adjacent units of the most similar unit are empty,
#'              one of them is selected as BMU.
#' @param obj The \code{espresso} object.
#' @param itr Iteration number.
#' @return BMUs given as a matrix.
#' 
.searchBMUByFill <- function(obj, itr) {
  rept <- obj@rept_som
  latest_bmu <- obj@bmu[[rept]][, itr-1]
  cnt <- c()
  for (i in 1:nrow(obj@topology)) {
    cnt <- c(cnt, sum(latest_bmu == i))
  }
  bmu <- c()
  for (i in 1:nrow(obj@exprs_som)) {
    map <- obj@map[[rept]]
    x <- matrix(obj@exprs_som[i, ], 
                nrow = nrow(map), 
                ncol = ncol(map), byrow = TRUE)
    d <- apply((x - map)^2, 1, sum)
    b <- which.min(d)
    adj <- (1:ncol(obj@topology))[obj@topology[b, ] == 1]
    if (cnt[b] > 0 && sum(cnt[adj] == 0) > 0) {
      tmp <- adj[cnt[adj] == 0]
      if (length(tmp) == 1) {
        b <- tmp
      } else {
        b <- sample(tmp, size = 1)
      }
    }
    bmu <- c(bmu, b)
    cnt[latest_bmu[i]] <- cnt[latest_bmu[i]] - 1
    cnt[b] <- cnt[b] + 1
  }
  return(bmu)
}

# =============================================================================
#' @title Compute neighborhood.
#' @description This function computes neighnorhood range.
#' @param obj The \code{espresso} object.
#' @param itr Iteration number.
#' @return neighborhood range
.computeNeighbor <- function(obj, itr) {
  sigma <- obj@radius * (1 - (itr - 1) / (obj@lsteps))
  d <- obj@dist
  if (obj@stochastic == TRUE) {
    r <- runif(1, min = obj@rmin, max = obj@rmax)
  }else{
    r <- 1.0
  }
  h <- exp(-r * d ^ 2 / (2 * sigma ^ 2))
  return(h)
}

# =============================================================================
#' @title Update map.
#' @description This function updates map according to 
#'              the neighborhood function and the best matching units.
#' @param obj The \code{espresso} object
#' @param h Neighborhood range.
#' @param bmu Best matching units.
#' @return updated map
#' 
.updateMap <- function(obj, h, bmu) {
  numer <- matrix(0, nrow = nrow(obj@topology), ncol = length(obj@gset))
  denom <- rep(0, nrow(obj@topology))
  for (i in 1:nrow(obj@exprs_som)) {
    exprs.i <- matrix(obj@exprs_som[i,], 
                      nrow = nrow(obj@topology), 
                      ncol = length(obj@gset), byrow = TRUE)
    numer <- numer + exprs.i * h[bmu[i], ]
    denom <- denom + h[bmu[i], ]
  }
  map_new <- numer / denom
  map_new[is.na(map_new)] <- 0
  return(map_new)
}

# =============================================================================
#' @title Swap clusters.
#' @description This function swap clusters on GraphSOM.
#' @param obj The \code{espresso} object.
#' @param bmu Best matching units.
#' @param r Repetition number.
#' @param l Learning step number.
#' @return \code{espresso} object
#' 
#' 
.swapClusters <- function(obj, bmu, r, l) {
  d <- ceiling(max(obj@dist) * (0.001 / 1) ^ (l / obj@lsteps))
  cand <- which(0 < obj@dist & obj@dist <= d, arr.ind = T)
  orig_bmu <- obj@bmu
  accs <- c()
  idx <- c()
  swapped_bmu <- list()
  acc0 <- .computeAccuracy(obj, r, l)
  for (i in 1:nrow(cand)) {
    if (cand[i, 1] > cand[i, 2]) {
      next
    }
    swap1 <- bmu == cand[i, 1]
    swap2 <- bmu == cand[i, 2]
    bmu_tmp <- bmu
    bmu_tmp[swap1] <- cand[i, 2]
    bmu_tmp[swap2] <- cand[i, 1]
    swapped_bmu[[i]] <- as.vector(bmu_tmp)
    obj@bmu[[r]][, as.character(l)] <- as.vector(bmu_tmp)
    acc <- .computeAccuracy(obj, r, l)
    accs <- c(accs, acc)
    idx <- c(idx, i)
    obj@bmu <- orig_bmu
  }
  if (var(accs) == 0) {
    p <- accs / sum(accs)
  } else {
    z <- scale(accs)
    fact = sqrt(length(z) / 2)
    score <- exp(z)^fact
    p <- score / sum(score)
  }
  selected_idx <- sample(idx, 1, prob = p)
  acc1 <- accs[which(idx == selected_idx)]
  t <- 0.001 ^ (l / (obj@lsteps - 1))
  accept_p <- .calc_accept_prob(acc0, acc1, t)
  if (runif(1) < accept_p) {
    obj@bmu[[r]][, as.character(l)] <- as.vector(swapped_bmu[[selected_idx]])
    map_tmp <- obj@map[[r]][cand[selected_idx, 1], ]
    obj@map[[r]][cand[selected_idx, 1], ] <- obj@map[[r]][cand[selected_idx, 2], ]
    obj@map[[r]][cand[selected_idx, 2], ] <- map_tmp
  }
  return(obj)
}

# =============================================================================
#' @title Sub-function of learning input data.
#' @description This sub-function trains GraphSOM with input learning data.
#' @param obj The \code{espresso} object
#' @return \code{espresso} object
#' @import progress
#' 
.trainGraphSOM <- function(obj) {
  if (obj@map_method == 'sample') {
    obj <- .initMapBySamples(obj)
  } else if (obj@map_method == 'random') {
    obj <- .initMapByRandom(obj)
  } else {
    warning(paste("'map_method' = ", obj@map_method, 
                  "is invalid. The 'random' method is forced to be used."))
    obj <- .initMapByRandom(obj)
  }
  if (obj@bmu_method != 'min' && obj@bmu_method != 'fill') {
    warning(paste("'bmu_method' = ", obj@map_method,
                  "is not allowed. The 'min' method is forced to be used."))
    obj@bmu_method <- 'min'
  }
  obj <- .initBMU(obj)
  pb <- progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = obj@rept*obj@lsteps, clear = FALSE, width = 60)
  message("GraphSOM clustering")
  genes <- obj@gset
  for (r in 1:obj@rept) {
    samples <- obj@ssets[[r]]
    exprs <- obj@exprs[match(samples, rownames(obj@exprs)), , drop = FALSE]
    exprs <- exprs[, match(genes, colnames(exprs)), drop = FALSE]
    obj@exprs_som <- exprs
    obj@rept_som <- r
    if (obj@bmu_method == 'min') {
      for (l in 1:obj@lsteps) {
        bmu <- .searchBMUByMin(obj)
        obj@bmu[[r]][,l] <- as.vector(bmu)
        h <- .computeNeighbor(obj, l)
        obj@map[[r]] <- .updateMap(obj, h, bmu)
        if (obj@swap == TRUE) {
          obj <- .swapClusters(obj, bmu, r, l - 1)
        }
        pb$tick()
      }
      bmu <- .searchBMUByMin(obj)
    } else {
      bmu <- .searchBMUByMin(obj)
      obj@bmu[[r]][, 1] <- as.vector(bmu)
      h <- .computeNeighbor(obj, 1)
      obj@map[[r]] <- .updateMap(obj, h, bmu)
      for (l in 2:obj@lsteps) {
        bmu <- .searchBMUByFill(obj, l)
        obj@bmu[[r]][,l] <- as.vector(bmu)
        h <- .computeNeighbor(obj, l)
        obj@map[[r]] <- .updateMap(obj, h, bmu)
        if (obj@swap == TRUE) {
          obj <- .swapClusters(obj, bmu, r, l - 1)
        }
        pb$tick()
      }
      bmu <- .searchBMUByFill(obj)
    }
    obj@bmu[[r]][, obj@lsteps + 1] <- as.vector(bmu)
  }
  return(obj)
}

# =============================================================================
#' @title Plot convergence curve of GraphSOM clustering.
#' @description This function outputs the convergence curve of 
#'              the GraphSOM clustering. The x-axis and y-axis indicate 
#'              the learning-steps and the number of cells 
#'              that migrated between clusters.
#' @param obj The \code{espresso} object.
#' @param rept Repetition number.
#' @param col Line color.
#' @param type Line type.
#' @param lwd Line width.
#' @importFrom graphics plot
#' @export
#' 
plotConvCurve <- function(obj, rept = NULL, 
                          col = "deepskyblue4", type = "l", lwd = 1.5) {
  x <- 1:obj@lsteps
  if (is.null(rept)) {
    rept <- 1:obj@rept
  }
  for (r in rept) {
    y <- c()
    for (i in x) {
      d <- sum(obj@bmu[[r]][, i] != obj@bmu[[r]][, i + 1])
      y <- c(y, d)
    }
    ncell <- nrow(obj@bmu[[r]])
    plot(x, y, xlab = 'Learning steps', ylab = 'Difference (#cells)',
         main = paste0('rept.', r), col = col, type = type, lwd = lwd, 
         ylim = c(0, ncell))
  }
}

# =============================================================================
#' @title Compute prediction accuracy of GraphSOM clustering results.
#' @description This function computes prediction accuracy of 
#'              GraphSOM clustering results.
#' @param obj The \code{espresso} object.
#' @param r Repetition number.
#' @param l Learning step number.
#' @return accuracy
#' 
.computeAccuracy <- function(obj, r, l) {
  score <- 0
  n_combs <- 0
  samples <- obj@ssets[[r]]
  n_samples <- length(samples)
  bmu <- obj@bmu[[r]][, as.character(l)]
  pred <- obj@domain2name[bmu, ]$name
  asgmt <- obj@asgmt[match(samples, obj@asgmt$sample), ]$domain
  score <- sum(obj@topology[pred, pred] == obj@topology[asgmt, asgmt]) - n_samples
  n_combs <- n_samples ^ 2 - n_samples
  acc <- score / n_combs
  return(acc)
}

# =============================================================================
#' @title Compute adjusted rand index (ARI).
#' @description This function computes ARI.
#' @param obj The \code{espresso} object.
#' @param r Repetition number.
#' @importFrom mclust adjustedRandIndex
#' @return adjusted rand index
#' 
.computeARI <- function(obj, r) {
  samples <- obj@ssets[[r]]
  n_samples <- length(samples)
  asgmt <- obj@asgmt[match(samples, obj@asgmt$sample), ]$domain
  asgmt <- as.integer(as.factor(asgmt))
  bmu <- obj@bmu[[r]][, as.character(obj@lsteps)]
  bmu <- bmu[match(samples, names(bmu))]
  pred <- obj@domain2name[bmu, ]$name
  pred <- as.integer(as.factor(pred))
  ari <- adjustedRandIndex(pred, asgmt)
  return(ari)
}

# =============================================================================
#' @title Summarize GraphSOM clustering results.
#' @description This function computs the mean, variance, maximum, and minimum 
#'              of scores, predicion accuracies and ARIs.
#' @param obj The \code{espresso} object.
#' @importFrom stats var
#' @return \code{espresso} object
#'
.summarizeScores <- function(obj) {
  summary <- vector("list", length = 12)
  names(summary) <- c("mean_score", "var_score", "max_score", "min_score",
                      "mean_acc", "var_acc", "max_acc", "min_acc",
                      "mean_ari", "var_ari", "max_ari", "min_ari")
  s <- obj@score
  avg <- apply(s, 2, mean)
  var <- apply(s, 2, var)
  max <- apply(s, 2, max)
  min <- apply(s, 2, min)
  summary["mean_score"] <- avg[1]
  summary["var_score"] <- var[1]
  summary["max_score"] <- max[1]
  summary["min_score"] <- min[1]
  summary["mean_acc"] <- avg[2]
  summary["var_acc"] <- var[2]
  summary["max_acc"] <- max[2]
  summary["min_acc"] <- min[2]
  summary["mean_ari"] <- avg[3]
  summary["var_ari"] <- var[3]
  summary["max_ari"] <- max[3]
  summary["min_ari"] <- min[3]
  obj@summary <- summary
  return(obj)
}

# =============================================================================
#' @title Sub-Function of evaluating GraphSOM clustering results.
#' @description This sub-function evaluates GraphSOM clustering results by 
#'              prediction accuracy and adjusted rand index (ARI).
#' @param obj The \code{espresso} object.
#' @return espresso object
#' @import progress
#' 
.evalGraphSOM <- function(obj) {
  pb <- progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = obj@rept, clear = FALSE, width= 60)
  message("Evaluation")
  score <- matrix(NA, nrow = obj@rept, ncol = 3)
  rownames(score) <- paste('rept.', 1:obj@rept, sep = '')
  colnames(score) <- c('score', 'acc', 'ari')
  for (r in 1:obj@rept) {
    acc <- .computeAccuracy(obj, r, obj@lsteps)
    ari <- .computeARI(obj, r)
    s <- acc + obj@coef * ari
    score[r, 'score'] <- s
    score[r, 'acc'] <- acc
    score[r, 'ari'] <- ari
    pb$tick()
  }
  obj@score <- as.data.frame(score)
  return(obj)
}

# =============================================================================
#' @title Select samples randomly for down sampling.
#' @description This function randomly select samples 
#'              from each domain for the purpose of down sampling.
#' @param obj The \code{espresso} object.
#' @param nsamples Number of samples to be selected from each domain.
#' @param rept Number of replications.
#' @return \code{espresso} object
#' 
.genSampleSets <- function(obj, nsamples, rept) {
  domains <- sort(unique(obj@asgmt$domain))
  obj@rept <- rept
  obj@ssets <- vector('list', length = obj@rept)
  names(obj@ssets) <- paste0('set.', 1:obj@rept)
  for (i in 1:obj@rept) {
    selected_samples <- c()
    for (d in domains) {
      sub <- subset(obj@asgmt, domain == d)
      if (nrow(sub) <= nsamples || is.null(nsamples)) {
        selected_samples <- c(selected_samples, sub$sample)
      } else {
        selected_samples <- c(selected_samples, sample(sub$sample, nsamples))
      }
    }
    obj@ssets[[i]] <- selected_samples
  }
  return(obj)
}

# =============================================================================
#' @title Initialize settings for GraohSOM clustering.
#' @description This function initializes the setting for GraphSOM clustering
#' @param obj The \code{espresso} object.
#' @param nsamples Number of samples of each domain.
#' @param rept Number of replications (default: 1).
#' @param radius Initial value of learning radius.
#' @param lsteps Number of learning steps (default: 100).
#' @param stochastic GraphSOM clustering introduces stochasticity 
#'                   for map update if \code{TRUE} (default: TRUE).
#' @param rmin Minimum value of scale factor. 
#'             for the stochastic map update (default: 0.5).
#' @param rmax Maximum value of scale factor 
#'             for the stochastic map update (default: 1.0).
#' @param nmin Minimum number of genes can be analyzed.
#' @param nmax Maximum number of genes can be analyzed.
#' @param map_method Map initialization method. 
#'                   The map vector is initilaized 
#'                   by input learning data if 'sample' is selected.
#'                   Otherwise, it is initialized uniformly random 
#'                   (default: 'sample'). 
#' @param bmu_method If the unit adjacent to the BMU candidate is empty,
#'                   it is set as BMU ('fill'). 
#'                   Otherwise, the most similar unit to input data is 
#'                   selected as BMU (default: 'min').
#' @param coef Coefficient of ARI for score computation.
#' @param swap Logical value determins whether to swap domains or not (default: TRUE)
#' @param seed Random seed.
#' @return espresso object
#' @export
#' 
initGraphSOM <- function(obj, nsamples = NULL, rept = 1, 
                         radius = NULL, lsteps = 100, 
                         stochastic = TRUE, rmin = 0.5, rmax = 1.0, nmin = 1, nmax = NULL,
                         map_method = "sample", bmu_method = "min", coef = 1.0, swap = TRUE, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(as.numeric(seed), kind = "Mersenne-Twister")
  }
  obj <- .genSampleSets(obj = obj, nsamples = nsamples, rept = as.integer(rept))
  obj <- .initGraphSOM(obj = obj, lsteps = as.integer(lsteps), radius = radius, 
                       stochastic = stochastic, rmin = rmin, rmax = rmax, 
                       nmin = nmin, nmax = nmax,
                       map_method = map_method, bmu_method = bmu_method, coef = coef, swap = swap)
  return(obj)
}

# =============================================================================
#' @title Run GraphSOM clustering.
#' @description This function executes GraohSOM clustering.
#' @param obj The \code{espresso} object.
#' @param gset Gene set.
#' @param seed Random seed.
#' @param verbose Whether to show messages.
#' @param version Character that specifies on which previous version GraphSOM computation should be performed. (e.g., "0")
#' @return \code{espresso} object
#' @export 
#' 
graphSOM <- function(obj, gset = NULL, seed = NULL, verbose = TRUE, version = NULL) {
  if (is.null(version)) {
    if (!is.null(seed)) {
      set.seed(seed, kind = "L'Ecuyer-CMRG")
    }
  } else if (version == "0") {
    if (!is.null(seed)) {
      set.seed(seed, kind = "Mersenne-Twister")
    }
  } else {
    stop(paste("version", ver, "is not available."))
    return(obj)
  }
  if (!is.null(gset)) {
    obj@gset <- gset[!is.na(match(gset, colnames(obj@exprs)))]
    if (length(obj@gset) != length(gset)) {
      ndel <- length(gset) - length(obj@gset)
      message(paste(ndel, "genes are ignored since they are not found in the expression matrix."))
    }
  }
  if (length(obj@gset) < obj@nmin || 
      (!is.null(obj@nmax) && (length(obj@gset) > obj@nmax))) {
    message(paste0("The number of genes (n = ", length(obj@gset), 
                   ") exceeds the limitation. Check '@nmin' or '@nmax' in an espresso object."))
    obj@summary <- list("mean_score" = NA, "var_score" = NA, "max_score" = NA, "min_score" = NA,
                        "mean_acc" = NA, "var_acc" = NA, "max_acc" = NA, "min_acc" = NA,
                        "mean_ari" = NA, "var_ari" = NA, "max_ari" = NA, "min_ari" = NA)
  } else {
    if (verbose == FALSE) {
      suppressMessages(obj <- .trainGraphSOM(obj = obj))
      suppressMessages(obj <- .evalGraphSOM(obj = obj))
      suppressMessages(obj <- .summarizeScores(obj = obj))
    } else {
      obj <- .trainGraphSOM(obj = obj)
      obj <- .evalGraphSOM(obj = obj)
      obj <- .summarizeScores(obj = obj)
      message("Complete")
    }
  }
  obj@exprs_som <- matrix(NA)
  return(obj)
}
