# =============================================================================
#' @title Scale gene expression data by log function
#' @description This function scales gene expression data 
#'              by log function with bases 2, 10, and e.
#' @param obj \code{espresso} object.
#' @param base base of a logarithm specified 
#'             as character ('2', '10', or 'e') (default: '10').
#' @return \code{espresso} object
#' @export
#' 
logScale <- function(obj, base = "10") {
  if (base == "2") {
    obj@exprs <- log2(obj@exprs + 1)
  } else if (base == "10") {
    obj@exprs <- log10(obj@exprs + 1)
  } else if (base == "e") {
    obj@exprs <- log(obj@exprs + 1)
  } else {
    warning(paste("Base =", base, "is invalid. Forcibly converted to log10."))
    obj@exprs <- log10(obj@exprs + 1)
  }
  return(obj)
}

# =============================================================================
#' @title Filter out low variable genes
#' @description This function filters out low expressed and low variable genes.
#' @param obj \code{espresso} object.
#' @param ncell the minimum number of expressed cells 
#'              for each gene (default: 2)
#' @param expressed threshold of expression level 
#'                  for gene filtering (default: 1.0)
#' @param sd standard deviation for gene filtering (default: 0.05)
#' @return \code{espresso} object
#' @export
#' 
filterGenes <- function(obj, ncell = 2, expressed = 1.0, sd = 0.05) {
  message("Filtering genes...")
  n_before <- dim(obj@exprs)[2]
  condition1 <- apply(obj@exprs, 2, function(x) {
    if (sum(x > expressed) >= ncell) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  condition2 <- apply(obj@exprs, 2, function(x) {
    if (sd(x) > sd) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  pass <- condition1 & condition2
  if (sum(pass) == 0) {
    warning("Filtering was not performed since no gene passed the filtering.")
    return(obj)
  }
  obj@exprs <- obj@exprs[, pass]
  obj@gset <- colnames(obj@exprs)
  n_after <- dim(obj@exprs)[2]
  message(paste0(n_before - n_after, 
                 " genes are filtered out (", 
                 n_before, " -> ", n_after, " genes)."))
  return(obj)
}






