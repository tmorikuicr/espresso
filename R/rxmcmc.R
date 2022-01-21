# =============================================================================
#' @title Add genes to the input gene set.
#' @description This function adds genes to the gene set to be optimized by MCMC.
#' @param gset Gene set to be optimized.
#' @param bgset Background gene set.
#' @param k Maximum number of genes operated in each MCMC step.
#' @importFrom utils combn
#' @return list of gene sets
#' 
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
#' @title Delete genes from the input gene set.
#' @description This function delete genes from the gene set to be optimized by MCMC.
#' @param gset Gene set to be optimized.
#' @param bgset Background gene set.
#' @param k Maximum number of genes operated in each MCMC step.
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
#' @title Replace genes of the input gene set.
#' @description This function replace genes of the input gene set.
#' @param gset Gene set to be optimized.
#' @param bgset Background gene set.
#' @param k Maximum number of genes operated in each MCMC step.
#' @importFrom utils combn
#' @return list of gene sets
#'
.replaceGenes <- function(gset, bgset, k){
  cand_gsets <- list()
  bgset <- setdiff(bgset,gset)
  for (i in 1:k) {
    del_idx <- sample(1:length(gset), i)
    del_gset <- gset[-del_idx]
    combs <- combn(bgset, i)
    for (j in 1:ncol(combs)) {
      cand_gsets <- c(cand_gsets, list(c(del_gset, unlist(combs[,j]))))
    }
  }
  return(cand_gsets)
}

# =============================================================================
#' @title Define a probability distribution for MCMC sampling.
#' @description This function defines a probability distribution for MCMC sampling.
#' @param tbl Summary table of a GraphSOM clustering result.
#' @param fact Scaling factor for computing selection probability
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
#' @title Sample a gene set from candidate gene sets.
#' @description This function samples a gene set from the candidate gene sets 
#'              under the defined probability distribution.
#' @param probs Probability distribution.
#' @return Index number of selected sample given as an integer.
#'
.sampleGeneSet <- function(probs){
  r <- runif(1)
  idx <- which(cumsum(probs) > r)[1]
  return(idx)
}

# =============================================================================
#' @title Calcurate the acceptance probability.
#' @description This function calculates the acceptance probability for MCMC sampling.
#'              The condition equation is based on the traditional Simulated Annealing approach.
#' @param acc0 Accuracy of the previouse sample.
#' @param acc Accuracy of the candidate sample.
#' @param t Temperature.
#' @return Probability given as numeric value
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
#' @title Make index vector for exchange of replicas.
#' @description This function makes a index vector for exchange of replicas.
#' @param n_repl Number of replicas
#' @param e Exchange time
#'
.idx4ex <- function(n_repl, e) {
  if (e %% 2 == 0) {
    1:floor(n_repl / 2) * 2 - 1
  } else {
    1:(floor(n_repl / 2) - 1) * 2
  }
}

# =============================================================================
#' @title Generate a data frame for storing replica exchange MCMC results.
#' @description This function generates a data frame for storing MCMC results.
#' @param n_repl Number of replicas
#' @param n Number of rows
#' @return Data frame
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
    df[[i]][1, ] <- c(-1.0, 0, -Inf, Inf, 0, 0, -Inf, Inf,
                      -1.0, 0, -Inf, Inf, 0)
  }
  return(df)
}

# =============================================================================
#' @title Delete gene sets which are already sampled.
#' @description This function deletes gene sets which are already sampled.
#' @param glist1 List of candidate gene sets.
#' @param glist2 List of gene sets already sampled.
#' @return List
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
#' @title Optimize genes by replica exchange MCMC.
#' @description This function performs optimization of genes included in a selected gene set.
#'              The optimization is done by replica exchange Markov Chain Monte Carlo (MCMC) approach.
#' @param obj The \code{espresso} object.
#' @param gset Gene set. 
#' @param temp Maximum temperature (default: 1.0).
#' @param itr Number of iterations until exchange of replicas (default: 10).
#' @param k Maximum number of genes operated in each MCMC step (default: 1).
#' @param seed Random seed.
#' @param n_cl Number of clusters for parallel computing (default: \code{detectCores()}).
#' @param n_ex Number of exchanges (default: 10)
#' @param n_repl Number of replicas (default: \code{detectCores()}).
#' @param n_ig Number of initial genes randomly selected for MCMC.
#' @param fact Scaling factor for computing selection probability.
#' @param version Character that specifies on which previous version GraphSOM computation should be performed. (e.g., "0.2.17")
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import progress
#' @return \code{espresso} object
#' @export 
#'
rxmcmc <- function(obj, gset = NULL, temp = 1.0, itr = 10, k = 1, seed = NULL, 
                   n_cl = detectCores(), n_ex = 10, n_repl = detectCores(), n_ig = 3, fact = NULL, version = NULL) {
  if (!is.null(seed)) {
    set.seed(seed, kind = "Mersenne-Twister")
  } else {
    seed <- 0
  }
  if (is.null(gset)) {
    warning("Input `gset` to be optimized.")
    return(obj)
  }
  if (version != "0.2.17" && !is.null(version)) {
    stop(paste("version", version, "is not available."))
    return(obj)
  }
  gset0 <- gset
  exist <- match(gset, colnames(obj@exprs))
  gset <- colnames(obj@exprs)[exist[!is.na(exist)]]
  diff <- setdiff(gset0, gset)
  if (length(diff) == 1) {
    warning(paste("A gene (", paste(diff, collapse = ", "), ") is removed 
                  from `gset` since it does not exist in @exprs.", sep=""))
  } else if (length(diff) > 1) {
    warning(paste("Genes (", paste(diff, collapse = ", "), ") are removed 
                  from `gset` since they do not exist in @exprs.", sep=""))
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
      if (is.null(version)) {
        for (l in 1:n_sets) {
          res_tmp <- graphSOM(obj, gset = cand_sets[[l]], verbose = FALSE, seed = seed)
          res[l, ] <- c(as.numeric(unlist(res_tmp@summary)), length(cand_sets[[l]]))
          map_tmp[[l]] <- res_tmp@map
          bmu_tmp[[l]] <- res_tmp@bmu
          score_tmp[[l]] <- res_tmp@score
          summary_tmp[[l]] <- res_tmp@summary
        }
      } else if (version == "0.2.17"){
        for (l in 1:n_sets) {
          res_tmp <- graphSOM(obj, gset = cand_sets[[l]], verbose = FALSE, version = "0.2.17")
          res[l, ] <- c(as.numeric(unlist(res_tmp@summary)), length(cand_sets[[l]]))
          map_tmp[[l]] <- res_tmp@map
          bmu_tmp[[l]] <- res_tmp@bmu
          score_tmp[[l]] <- res_tmp@score
          summary_tmp[[l]] <- res_tmp@summary
        }
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
          if (r == 0 || is.na(t[r + 1])) {
            next
          }
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
#' @title Plot MCMC history.
#' @description This function plots the sampling history and replica exchange history.
#' @param obj The \code{espresso} object.
#' @import RColorBrewer
#' @importFrom graphics legend lines par
#' @export
#' 
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

# =============================================================================
#' @title Write MCMC results.
#' @description This function outputs the results of replica exchange MCMC.
#' @param obj The \code{espresso} object.
#' @param dir Path to the output directory.
#' @export
#' 
writeMCMC <- function(obj, dir = NULL) {
  if (is.null(obj@mcmc)) {
    warning("This function becomes valid after running rxmcmc.")
  } else {
    if (is.null(dir)) {
      dt <- strsplit(as.character(Sys.time()), " ")[[1]]
      d <- dt[1]
      t <- gsub(":", "", dt[2])
      dir <- paste0("writeMCMC_", d, "-", t)
    }
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    write.table(obj@mcmc[["best_repl"]],
                file = paste(dir, "best_repl.txt", sep = "/"), row.names = FALSE, quote = FALSE)
    gene_sizes <- c()
    for (i in 1:length(obj@mcmc[["final_genes"]])) {
      gene_sizes <- c(gene_sizes, length(obj@mcmc[["final_genes"]][[i]]))
      write.table(data.frame(sort(obj@mcmc[["final_genes"]][[i]])),
                  file = paste(dir, paste0("final_genes_repl.", i, ".txt"), sep = "/"), 
                  col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
    summary <- as.data.frame(matrix(unlist(obj@summary), ncol = 12, byrow = TRUE))
    colnames(summary) <- c("mean_score", "var_score","max_score", "min_score", 
                           "mean_acc", "var_acc", "max_acc", "min_acc",
                           "mean_ari", "var_ari", "max_ari", "min_ari")
    rownames(summary) <- paste0("repl.", 1:nrow(summary))
    summary <- transform(summary, gene_size = gene_sizes)
    summary <- transform(summary, repl = rownames(summary))
    summary <- summary[, c(14, 1:13)]
    write.table(summary, file = paste(dir, "summary.txt", sep = "/"), row.names= F, quote = F, sep = "\t")
  }
}
