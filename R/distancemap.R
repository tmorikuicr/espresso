# =============================================================================
#' @title Plot distance maps
#' @description This function plots distance matrices between domain-domain and 
#'              before/after GraphSOM clustering.
#' @param obj espresso object
#' @param type type of map ("d": domain, "s": sample, "b": both)
#' @importFrom igraph shortest.paths
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @import RColorBrewer
#' @export
#' 
plotDistMap <- function(obj, type = "b") {
  if (type != "d" && type != "s" && type != "b") {
    warning("'type' option is invalid.")
  } else {
    distmat <- obj@dist
    distmat[which(distmat == Inf, arr.ind =TRUE)] <- NA
    if (type == "d" || type == "b") {
      
      # Plot distance maps of domains before GraphSOM ------------------------
      pheatmap(distmat, cluster_cols = F, cluster_rows = F, 
               color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100), na_col = "black",
               main = "Distance of domains")
    } 
    if (type == "s" || type == "b") {
      for (i in 1:length(obj@ssets)) {
        # Plot distance maps of samples before GraphSOM ------------------------
        asgmt <- obj@asgmt[match(obj@ssets[[i]], obj@asgmt$sample), ]
        rownames(asgmt) <- asgmt$sample
        asgmt <- asgmt[-1]
        tmp <- NULL
        for (j in obj@domain2name$name) {
          tmp <- rbind(tmp, subset(asgmt, domain == j))
        }

        asgmt <- tmp
        sample_dist = distmat[unlist(as.matrix(asgmt)), unlist(as.matrix(asgmt))]
        rownames(sample_dist) = as.character(rownames(asgmt))
        colnames(sample_dist) = as.character(rownames(asgmt))
        if (length(obj@ssets) == 1) {
          main_title <- "Distance of sample cells"
        } else {
          main_title <- paste("Distance of sample cells\n(replication ", i, ")", sep ="")
        }
        pheatmap(sample_dist, 
                 cluster_cols = F, cluster_rows = F, 
                 fontsize_row = 5, fontsize_col = 5, 
                 cellheight = 5, cellwidth = 5, 
                 color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100), na_col = "black",
                 main = main_title
        )
        
        # Plot distance maps of samples after GraphSOM -------------------------
        if (!is.null(obj@bmu)) {
          domains <- obj@domain2name[obj@bmu[[i]][, ncol(obj@bmu[[i]])], ]$name
          asgmt2 <- data.frame(domain = domains)
          rownames(asgmt2) <- names(obj@bmu[[i]][, ncol(obj@bmu[[i]])])
          sample_dist2 = distmat[unlist(as.matrix(asgmt2)), unlist(as.matrix(asgmt2))]
          rownames(sample_dist2) = as.character(rownames(asgmt2))
          colnames(sample_dist2) = as.character(rownames(asgmt2))
          if (length(obj@ssets) == 1) {
            main_title <- "Distance of sample cells after GraphSOM"
          } else {
            main_title <- paste("Distance of sample cells after GraphSOM\n(replication ", i, ")", sep = "")
          }
          vec <- as.vector(sample_dist2)
          if (length(unique(vec[!is.na(vec)])) != 1) {
            pheatmap(sample_dist2,
                     cluster_cols = F, cluster_rows = F, 
                     fontsize_row = 5, fontsize_col = 5, 
                     cellheight = 5, cellwidth = 5, 
                     color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100), na_col = "black",
                     main = main_title
            )
            # Plot distance maps of samples before/after GraphSOM -------------------
            ds <- sample_dist
            ds[upper.tri(ds)] <- sample_dist2[upper.tri(sample_dist2)]
            if (length(obj@ssets) == 1) {
              main_title <- "Distance of sample cells in original (lower) and after (upper) GraphSOM"
            } else {
              main_title <- paste("Distance of sample cells in original (lower) and after (upper) GraphSOM\n(replication ", i, ")", sep = "")
            }
            pheatmap(ds, 
                     cluster_cols = F, cluster_rows = F, 
                     fontsize_row = 5, fontsize_col = 5, 
                     cellheight = 5, cellwidth = 5, 
                     color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100), na_col = "black",
                     main = main_title
            )
          } else {
            warning("Distance map after graph-SOM clustering cannot be drown since only one cluseter was generated.")
          }
        }
      }
    }
  }
}