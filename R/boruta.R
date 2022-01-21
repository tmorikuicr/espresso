# =============================================================================
#' @title Select important features by Boruta.
#' @description This function gets important features by \code{Boruta} package.
#' @param obj The \code{espresso} object.
#' @param maxRuns The maximum number of importance source runs.
#' @param decision Decision type of Boruta 
#'                 ("c": Confirmed, "nr": Not Rejected (i.e., Confirmed & Tentative)
#'                 (default: "nr").
#' @param seed Random seed.
#' @import Boruta
#' @import rFerns
#' @return Data Frame of the result of feature selection.
#' @export 
#' 
selectFeatures <- function(obj, maxRuns = 500, decision = "nr", seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed, kind = "Mersenne-Twister")
  }
  X <- obj@exprs
  Y <- as.factor(obj@asgmt[, "domain"])
  res_boruta <- Boruta(X, Y, getImp = getImpFerns, maxRuns = maxRuns)
  res <- attStats(res_boruta)
  rownames(res) <- colnames(X)
  res <- res[order(res$normHits, decreasing = T), ]
  if (decision == "c") {
    conf_genes <- rownames(subset(res, decision == "Confirmed"))
  } else if (decision == "nr") {
    conf_genes <- rownames(subset(res, decision != "Rejected"))
  } else {
    warning(paste("'decision = ", decision, "' is invalid. 'nr' is used insted of it."))
    conf_genes <- rownames(subset(res, decision != "Rejected"))
  }
  return(conf_genes)
}