# =============================================================================
#' @title Plot distance map as heatmap.
#' @description This function plots distance maps between domain-domain and 
#'              before/after GraphSOM clustering.
#' @param obj The \code{espresso} object.
#' @param type Type of the plotted maps ("d": domain, "s": sample, "b": both).
#' @param size_d Font and cell sizes for row and column for distance map of domains (default: NULL).
#' @param size_s Font and cell sizes for row and column for distance map of cell samples (default: 5).
#' @importFrom igraph shortest.paths
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @import RColorBrewer
#' @export
#' 
plotDistMap <- function(obj, type = "b", size_d = NULL, size_s = 5) {
  if (type != "d" && type != "s" && type != "b") {
    warning("'type' option is invalid.")
  } else {
    distmat <- obj@dist
    distmat[which(distmat == Inf, arr.ind =TRUE)] <- NA
    if (type == "d" || type == "b") {
      
      # Plot distance maps of domains before GraphSOM ------------------------
      if (!is.null(size_d)) {
        pheatmap(distmat, cluster_cols = F, cluster_rows = F, 
                 color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100), na_col = "black",
                 main = "Distance of domains", 
                 fontsize_row = size_d, fontsize_col = size_d,
                 cellheight = size_d, cellwidth = size_d
                 )
      } else {
        pheatmap(distmat, cluster_cols = F, cluster_rows = F, 
                 color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100), na_col = "black",
                 main = "Distance of domains")
      }
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
          main_title <- paste("Distance of sample cells\n(rept.", i, ")", sep ="")
        }
        pheatmap(sample_dist, 
                 cluster_cols = F, cluster_rows = F, 
                 fontsize_row = size_s, fontsize_col = size_s, 
                 cellheight = size_s, cellwidth = size_s, 
                 color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100), na_col = "black",
                 main = main_title
        )
        for (repl in names(obj@bmu)) {
          # Plot distance maps of samples after GraphSOM -------------------------
          if (!is.null(obj@bmu)) {
            if ( length(grep("repl", names(obj@bmu))) > 0 ) {
              domains <- obj@domain2name[obj@bmu[[repl]][[i]][, ncol(obj@bmu[[repl]][[i]])], ]$name
              asgmt2 <- data.frame(domain = domains)
              rownames(asgmt2) <- names(obj@bmu[[repl]][[i]][, ncol(obj@bmu[[repl]][[i]])])
              sample_dist2 = distmat[unlist(as.matrix(asgmt2)), unlist(as.matrix(asgmt2))]
              rownames(sample_dist2) = as.character(rownames(asgmt2))
              colnames(sample_dist2) = as.character(rownames(asgmt2))
              if (length(obj@ssets) == 1 && length(obj@bmu) == 1) {
                main_title <- "Distance of sample cells after GraphSOM"
              } else if (length(obj@ssets) == 1 && length(obj@bmu) > 1) {
                main_title <- paste("Distance of sample cells after GraphSOM\n(", repl, ")", sep ="")
              } else if (length(obj@ssets) > 1 && length(obj@bmu) == 1) {
                main_title <- paste("Distance of sample cells after GraphSOM\n(rept.", i, ")", sep ="")
              } else {
                main_title <- paste("Distance of sample cells after GraphSOM\n(", repl, ", rept.", i, ")", sep ="")
              }
            } else {
              domains <- obj@domain2name[obj@bmu[[i]][, ncol(obj@bmu[[i]])], ]$name
              asgmt2 <- data.frame(domain = domains)
              rownames(asgmt2) <- names(obj@bmu[[i]][, ncol(obj@bmu[[i]])])
              sample_dist2 = distmat[unlist(as.matrix(asgmt2)), unlist(as.matrix(asgmt2))]
              rownames(sample_dist2) = as.character(rownames(asgmt2))
              colnames(sample_dist2) = as.character(rownames(asgmt2))
              if (length(obj@ssets) == 1) {
                main_title <- "Distance of sample cells after GraphSOM"
              } else {
                main_title <- paste("Distance of sample cells after GraphSOM\n(rept.", i, ")", sep ="")
              }
              
            }
            vec <- as.vector(sample_dist2)
            if (length(unique(vec[!is.na(vec)])) != 1) {
              pheatmap(sample_dist2,
                       cluster_cols = F, cluster_rows = F, 
                       fontsize_row = size_s, fontsize_col = size_s, 
                       cellheight = size_s, cellwidth = size_s,  
                       color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100), na_col = "black",
                       main = main_title
              )
              # Plot distance maps of samples before/after GraphSOM -------------------
              ds <- sample_dist
              ds[upper.tri(ds)] <- sample_dist2[upper.tri(sample_dist2)]
              if ( length(grep("repl", names(obj@bmu))) > 0 ) {
                if (length(obj@ssets) == 1 && length(obj@bmu) == 1) {
                  main_title <- "Distance of sample cells \nin original (lower) and after GraphSOM (upper)"
                } else if (length(obj@ssets) == 1 && length(obj@bmu) > 1) {
                  main_title <- paste("Distance of sample cells \nin original (lower) and after GraphSOM (upper)\n(", repl, ")", sep = "")
                } else if (length(obj@ssets) > 1 && length(obj@bmu) == 1) {
                  main_title <- paste("Distance of sample cells \nin original (lower) and after GraphSOM (upper)\n(rept.", i, ")", sep = "")
                } else {
                  main_title <- paste("Distance of sample cells \nin original (lower) and after GraphSOM (upper)\n(", repl, ", rept.", i, ")", sep = "")
                }
              } else {
                if (length(obj@ssets) == 1) {
                  main_title <- "Distance of sample cells \nin original (lower) and after GraphSOM (upper)"
                } else {
                  main_title <- paste("Distance of sample cells \nin original (lower) and after GraphSOM (upper)\n(rept.", i, ")", sep = "")
                }
              }
              pheatmap(ds, 
                       cluster_cols = F, cluster_rows = F, 
                       fontsize_row = size_s, fontsize_col = size_s, 
                       cellheight = size_s, cellwidth = size_s, 
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
}