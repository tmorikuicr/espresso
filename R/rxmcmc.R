# =============================================================================
#' @title Add k genes to a gene set
#' @description This function adds k genes to a gene set to be optimized by MCMC.
#' @param gset gene set to be optimized.
#' @param bgset background gene set.
#' @param k the maximum number of genes operated in each MCMC step.
#' @importFrom utils combn
#' @return list of gene sets
.addGenes <- function(gset, bgset, k){
  cand_gsets <- list()
  bgset <- setdiff(bgset, gset)
  if (length(bgset) > 0) {
    for (i in 1:k) {
      combs <- combn(bgset, i)
      for (j in 1:ncol(combs)) {
        cand_gsets <- c(cand_gsets, list(c(gset, unlist(combs[,j]))))
      }
    }
  } else {
    gset <- sample(gset, length(gset) - 1)
  }
  return(cand_gsets)  
}

# =============================================================================
#' @title Delete k genes to a gene set
#' @description This function delete k genes to a gene set to be optimized by MCMC.
#' @param gset gene set to be optimized.
#' @param bgset background gene set.
#' @param k the maximum number of genes operated in each MCMC step.
#' @importFrom utils combn
#' @return list of gene sets
#'
.delGenes <- function(gset, bgset, k){
  cand_gsets <- list()
  for (i in 1:k) {
    if (length(gset) - k <= 1) {
      break
    }
    combs <- combn(gset, i)
    for (j in 1:ncol(combs)) {
      cand_gsets <- c(cand_gsets, list(setdiff(gset, unlist(combs[,j]))))
    }
  }
  return(cand_gsets)
}

# =============================================================================
#' @title Replace k genes
#' @description This function replace k genes
#' @param gset gene set to be optimized.
#' @param bgset background gene set.
#' @param k the maximum number of genes operated in each MCMC step.
#' @importFrom utils combn
#' @return list of gene sets
#'
.replaceGenes <- function(gset, bgset, k){
  cand_gsets <- list()
  bgset <- setdiff(bgset,gset)
  for(i in 1:k){
    del_idx <- sample(1:length(gset), i)
    del_gset <- gset[-del_idx]
    combs <- combn(bgset, i)
    for(j in 1:ncol(combs)){
      cand_gsets <- c(cand_gsets, list(c(del_gset, unlist(combs[,j]))))
    }
  }
  return(cand_gsets)
}

# =============================================================================
#' @title Define a prior distribution for MCMC sampling
#' @description This function defines a prior distribution for MCMC sampling.
#' @param tbl summary table of a GraphSOM clustering result.
#' @param fact scaling factor for computing selection probability
#' @return selection probabilities given as a matrix.
#'
.definePriorDistr <- function(tbl, fact){
  s <- tbl$mean_score
  s[is.na(s)] <- 0
  if (is.na(var(s))) {
    p <- matrix(1, ncol = 1)
  } else if (var(s) == 0) {
    p <- matrix(s / sum(s), ncol = 1)
  } else {
    z <- scale(s)
    if (is.null(fact)) {
      fact = sqrt(length(z) / 2)
    }
    score <- exp(z)^fact
    p <- score / sum(score)
  }
  rownames(p) <- rownames(tbl)
  colnames(p) <- 'prob'
  return(p)
}

# =============================================================================
#' @title Sample gene set from candidate gene sets
#' @description This function samples a gene set from candidate gene sets 
#' under the defined prior probability distribution.
#' @param probs vector of probability distribution.
#' @return index number of selected sample given as an integer.
#'
.sampleGeneSet <- function(probs){
  r <- runif(1)
  idx <- which(cumsum(probs) > r)[1]
  return(idx)
}

# =============================================================================
#' @title Calcurate the acceptance probability
#' @description This function calculates the acceptance probability for MCMC sampling.
#' The condition equation is based on the traditional Simulated Annealing approach.
#' @param acc0 accuracy of the previouse sample.
#' @param acc accuracy of the candidate sample.
#' @param t temperature parameter.
#' @return probability given as numeric value
.calc_accept_prob <- function(acc0, acc, t){
  if (is.na(acc)) {
    p <- 0
    return(p)
  } else {
    e1 <- -acc0
    e2 <- -acc
    if (e1 >= e2) {
      p <- 1
    } else {
      p <- exp(-(e2 - e1) / t)
    }
    return(p)
  }
}

# =============================================================================
#' @title Make index vector for exchange of replicas
#' @description This function makes a index vector for exchange of replicas
#' @param n_repl the number of replicas
#' @param e exchange time
#' @return vector
.idx4ex <- function(n_repl, e) {
  if (e %% 2 == 0) {
    1:floor(n_repl / 2) * 2 - 1
  } else {
    1:(floor(n_repl / 2) - 1) * 2
  }
}

# =============================================================================
#' @title Generate a data frame for replica exchange MCMC results
#' @description This function generates a data frame for storing MCMC results
#' @param n_repl the number of replicas
#' @param n the number of rows
#' @return a data frame
#' 
.generateDataFrameForRXMCMC <- function(n_repl, n) {
  df <- vector("list", length = n_repl)
  for (i in 1:n_repl) {
    df[[i]] <- data.frame(mean_score = numeric(n), var_score = numeric(n), 
                          max_score = numeric(n), min_score = numeric(n),
                          mean_acc = numeric(n), var_acc = numeric(n),
                          max_acc = numeric(n), min_acc = numeric(n),
                          mean_ari = numeric(n), var_ari = numeric(n),
                          max_ari = numeric(n), min_ari = numeric(n),
                          gene_size = integer(n), stringsAsFactors = FALSE)
    rownames(df[[i]]) <- paste0("step.", 0:(n - 1))
    #df[[i]][1, ] <- c(-1.0, 0, 0, 0, -1.0, 0, 0)
    df[[i]][1, ] <- c(-1.0, 0, -Inf, Inf,
                      0, 0, -Inf, Inf,
                      -1.0, 0, -Inf, Inf,
                      0)
  }
  return(df)
}

# =============================================================================
#' @title Delete gene sets already sampled
#' @description This function deletes gene sets which are already sampled.
#' @param glist1 list of candidate gene sets
#' @param glist2 list of gene sets already sampled
#' @return a list
#' 
.delGeneSets <- function(glist1, glist2) {
  del_list <- c()
  for (i in 1:length(glist1)) {
    for (j in 1:length(glist2)) {
      if (identical(sort(glist1[[i]]), sort(glist2[[j]]))) {
        del_list <- c(del_list, i)
        break
      }
    }
  }
  glist1[del_list] <- NULL
  return(glist1)
}

# =============================================================================
#' @title Optimize genes by replica exchange MCMC
#' @description This function performs optimization of genes included in a selected gene set.
#'              The optimization is done by replica exchange Markov Chain Monte Carlo (MCMC) approach.
#' @param obj espresso object
#' @param gset geneset 
#' @param temp maximum temperature
#' @param itr the number of iterations until exchange of replicas
#' @param k the maximum number of genes operated in each MCMC step
#' @param seed the seed
#' @param n_cl the number of clusters for parallel computing
#' @param n_ex the number of exchanges
#' @param n_repl the number of replicas
#' @param n_ig the number of initial genes randomly selected for MCMC.
#' @param fact scaling factor for computing selection probability
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import progress
#' @return espresso object
#' @export 
#'
rxmcmc <- function(obj, gset = NULL, temp = 1.0, itr = 10, k = 1, seed = 0, 
                   n_cl = detectCores(), n_ex = 10, n_repl = detectCores(), n_ig = 3, fact = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(gset)) {
    warning("Input `gset` to be optimized.")
    return(obj)
  }
  t <- 2^seq(log2(temp), log2(0.001), len = n_repl)
  repl <- 1:n_repl
  history_sampling <- .generateDataFrameForRXMCMC(n_repl, itr * n_ex + 1)
  names(history_sampling) <- paste0("repl.", 1:n_repl)
  history_max <- .generateDataFrameForRXMCMC(n_repl, itr * n_ex + 1)
  names(history_max) <- paste0("repl.", 1:n_repl)
  history_ex <- matrix(0, ncol = n_repl, nrow = itr * n_ex + 1)
  colnames(history_ex) <- paste0("repl.", 1:n_repl)
  rownames(history_ex) <- paste0("step.", 0:(itr * n_ex))
  genes <- vector("list", length = n_repl)
  best_genes <- vector("list", length = n_repl)
  names(best_genes) <- paste0("repl.", 1:n_repl)
  bg_genes <- gset
  history_cand <- vector("list", length = n_repl)
  best_map <- vector("list", length = n_repl)
  best_bmu <- vector("list", length = n_repl)
  best_score <- vector("list", length = n_repl)
  best_summary <- vector("list", length = n_repl)
  ex <- 0
  pb <- progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = itr * n_ex, clear = FALSE, width = 60)
  message("Optimize genes by Replica exchange MCMC.")
  for (i in 2:(itr * n_ex + 1)) {
    cl <- makeCluster(n_cl, setup_strategy = "sequential")
    clusterSetRNGStream(cl, seed)
    registerDoParallel(cl)
    res_rx <- foreach(j = 1:n_repl, .packages = 'espresso') %dopar% {
      map_tmp <- vector("list", length = n_repl)
      bmu_tmp <- vector("list", length = n_repl)
      score_tmp <- vector("list", length = n_repl)
      summary_tmp <- vector("list", length = n_repl)
      res_each <- vector("list", length = 7)
      names(res_each) <- c('genes', 'summary', 'map', 'bmu', 'score', 'summary0', 'cand')
      if (i == 2 && n_ig > 0) {
        cand_sets <- list(sample(gset, size = n_ig))
      } else {
        add <- .addGenes(genes[[j]], bg_genes, k)
        del <- .delGenes(genes[[j]], bg_genes, k)
        cand_sets <- c(add, del)
      }
      cand_sets <- .delGeneSets(cand_sets, history_cand[[j]])
      if (length(cand_sets) == 0) {
        cand_sets <- .replaceGenes(genes[[j]], bg_genes, k)
      }
      n_sets <- length(cand_sets)
      res <- .generateDataFrameForRXMCMC(1, n_sets)[[1]]
      for (l in 1:n_sets) {
        res_tmp <- graphSOM(obj, gset = cand_sets[[l]], verbose = FALSE)
        res[l, ] <- c(as.numeric(unlist(res_tmp@summary)), length(cand_sets[[l]]))
        map_tmp[[l]] <- res_tmp@map
        bmu_tmp[[l]] <- res_tmp@bmu
        score_tmp[[l]] <- res_tmp@score
        summary_tmp[[l]] <- res_tmp@summary
      }
      probs <- .definePriorDistr(res, fact)
      idx <- .sampleGeneSet(probs)
      acc0 <- history_sampling[[j]]$mean_score[i - 1]
      acc <- res[idx, ]$mean_score
      p <- .calc_accept_prob(acc0, acc, t[match(j, repl)])
      if (runif(1) < p) {
        res_each$genes <- cand_sets[[idx]]
        res_each$summary <- res[idx, ]
        res_each$map <- map_tmp[[idx]]
        res_each$bmu <- bmu_tmp[[idx]]
        res_each$score <- score_tmp[[idx]]
        res_each$summary0 <- summary_tmp[[idx]]
        res_each$cand <- cand_sets[[idx]]
      } else {
        res_each$genes <- genes[[j]]
        res_each$summary <- history_sampling[[j]][i - 1, ]
        res_each$map <- best_map[[j]]
        res_each$bmu <- best_bmu[[j]]
        res_each$score <- best_score[[j]]
        res_each$summary0 <- best_summary[[j]]
        res_each$cand <- cand_sets[[idx]]
      }
      res_each
    }
    stopCluster(cl)
    for (j in 1:n_repl) {
      genes[[j]] <- res_rx[[j]]$genes
      history_sampling[[j]][i, ] <- res_rx[[j]]$summary
      history_cand[[j]] <- c(history_cand[[j]], list(res_rx[[j]]$cand))
    }
    for (j in 1:n_repl) {
      if (history_sampling[[j]][i, ]$mean_score > history_max[[j]][i - 1, ]$mean_score) {
        history_max[[j]][i, ] <- history_sampling[[j]][i, ]
        best_genes[[j]] <- res_rx[[j]]$genes
        best_map[[j]] <- res_rx[[j]]$map
        best_bmu[[j]] <- res_rx[[j]]$bmu
        best_score[[j]] <- res_rx[[j]]$score
        best_summary[[j]] <- res_rx[[j]]$summary0
      } else {
        history_max[[j]][i, ] <- history_max[[j]][i - 1, ]
      }
    }
    history_ex[i, ] <- t[match(1:n_repl, repl)]
    if (!is.nan(i %% itr)) {
      if (i %% itr == 0) {
        ex <- ex + 1
        for (r in .idx4ex(n_repl, ex)) {
          inv_t1 <- 1 / t[r]
          inv_t2 <- 1 / t[r + 1]
          e1 <- res_rx[[repl[r]]]$summary$mean_score
          e2 <- res_rx[[repl[r + 1]]]$summary$mean_score
          w <- exp((inv_t1 - inv_t2) * (- e1 + e2 ))
          if (runif(1, 0, 1) < w) {
            tmp <- repl[r]
            repl[r] <- repl[r + 1]
            repl[r + 1] <- tmp
          }
        }
      }
    }
    pb$tick()
  }
  for (i in 1:n_repl) {
    history_sampling[[i]] <- history_sampling[[i]][-1, ]
    history_max[[i]] <- history_max[[i]][-1, ]
  }
  obj@mcmc[["sampling"]] <- history_sampling
  obj@mcmc[["max"]] <- history_max
  obj@mcmc[["exchange"]] <- history_ex[-1, ]
  
  best_mean_score <- -1000000
  idx <- c()
  for (i in 1:n_repl) {
    s <- obj@mcmc[["max"]][[i]]$mean_score
    if (s[itr * n_ex] > best_mean_score) {
      best_mean_score <- s[itr * n_ex]
      idx <- c(i)
    } else if (s[itr * n_ex] == best_mean_score) {
      idx <- c(idx, i)
    }
  }
  steps <- c()
  for (i in idx) {
    steps <- c(steps, which.max(obj@mcmc[["sampling"]][[i]]$mean_score))
  }
  best_repl <- as.data.frame(matrix(c(idx, steps), ncol = 2))
  colnames(best_repl) <- c("repl", "step")
  obj@mcmc[["final_genes"]] <- best_genes
  obj@mcmc[["best_genes"]] <- best_genes[idx]
  obj@mcmc[["best_repl"]] <- best_repl
  obj@map <- best_map
  obj@bmu <- best_bmu
  obj@score <- best_score
  obj@summary <- best_summary
  names(obj@map) <- paste0("repl.", 1:n_repl)
  names(obj@bmu) <- paste0("repl.", 1:n_repl)
  names(obj@score) <- paste0("repl.", 1:n_repl)
  names(obj@summary) <- paste0("repl.", 1:n_repl)
  return(obj)
}

# =============================================================================
#' @title Plot MCMC history
#' @description This function plots the sampling history and 
#'              replica exchange history.
#' @param obj espresso object
#' @import RColorBrewer
#' @importFrom graphics legend lines par
#' @export
plotMCMC <- function(obj) {
  n_steps <- nrow(obj@mcmc[['max']][[1]])
  n_repl <- length(obj@mcmc[['max']])
  col_pal <- colorRampPalette(brewer.pal(11, "Spectral"))
  colors <- col_pal(n_repl)
  replicas <- paste0('repl.', 1:n_repl)
  
  # get max score
  scores <- c()
  accs <- c()
  aris <- c()
  temps <- c()
  for (i in 1:n_repl) {
    scores <- c(scores, obj@mcmc[['sampling']][[i]]$mean_score)
    accs <- c(accs, obj@mcmc[['sampling']][[i]]$mean_acc)
    aris <- c(aris, obj@mcmc[['sampling']][[i]]$mean_ari)
  }
  
  # Plotting parameters
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar)) 
  par(mar = c(5,4,4,6))
  
  # Sampling (scores) --------------------------------
  plot(0, 0, type = 'n', 
       xlim = range(1:n_steps), ylim = range(scores), 
       xlab = 'Step', ylab = 'Mean score')
  for (i in 1:n_repl) {
    lines(1:n_steps, obj@mcmc[['sampling']][[i]]$mean_score, lty = 1, col = colors[i], lwd = 1.5)
  }
  par(xpd=T)
  legend(par()$usr[2], par()$usr[4], legend = replicas, col = colors, lty = 1, lwd = 1.5, bty = 'n')
  
  # Best (scores) ------------------------------------
  plot(0, 0, type = 'n', 
       xlim = range(1:n_steps), ylim = range(scores), 
       xlab = 'Step', ylab = 'Maximum mean score')
  for (i in 1:n_repl) {
    lines(1:n_steps, obj@mcmc[['max']][[i]]$mean_score, lty = 1, col = colors[i], lwd = 1.5)
  }
  par(xpd=T)
  legend(par()$usr[2], par()$usr[4], legend = replicas, col = colors, lty = 1, lwd = 1.5, bty = 'n')
  
  # Sampling (acc) -----------------------------------
  plot(0, 0, type = 'n', 
       xlim = range(1:n_steps), ylim = range(accs), 
       xlab = 'Step', ylab = 'Mean accuracy ')
  for (i in 1:n_repl) {
    lines(1:n_steps, obj@mcmc[['sampling']][[i]]$mean_acc, lty = 1, col = colors[i], lwd = 1.5)
  }
  par(xpd=T)
  legend(par()$usr[2], par()$usr[4], legend = replicas, col = colors, lty = 1, lwd = 1.5, bty = 'n')
  
  # Sampling (ari) -----------------------------------
  plot(0, 0, type = 'n', 
       xlim = range(1:n_steps), ylim = range(aris), 
       xlab = 'Step', ylab = 'Mean ARI')
  for (i in 1:n_repl) {
    lines(1:n_steps, obj@mcmc[['sampling']][[i]]$mean_ari, lty = 1, col = colors[i], lwd = 1.5)
  }
  par(xpd=T)
  legend(par()$usr[2], par()$usr[4], legend = replicas, col = colors, lty = 1, lwd = 1.5, bty = 'n')

  # Reeplica exchange ---------------------------------
  suppressWarnings(
    plot(0, 0, type = 'n', 
         xlim = range(1:n_steps), ylim = range(obj@mcmc[['exchange']]), 
         xlab = 'Step', ylab = 'Temperature', log = 'y')
  )
  for (i in 1:n_repl) {
    lines(1:n_steps, obj@mcmc[['exchange']][, i], lty = 1, col = colors[i], lwd = 1.5)
  }
  par(xpd=T)
  legend(par()$usr[2], par()$usr[4], legend = replicas, col = colors, lty = 1, lwd = 1.5, bty = 'n')
}
