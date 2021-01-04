# =============================================================================
#' @title Select important features by Boruta
#' @description This function get feature importance by Boruta package
#' @param obj espresso object.
#' @param maxRuns the maximum number of importance source runs.
#' @import Boruta
#' @import rFerns
#' @return a data.frame of the result of feature selection
#' @export 
#' 
selectFeatures <- function(obj, maxRuns = 500) {
  X <- obj@exprs
  Y <- as.factor(obj@asgmt[, "domain"])
  res_boruta <- Boruta(X, Y, getImp = getImpFerns, maxRuns = maxRuns)
  res <- attStats(res_boruta)
  res <- res[order(res$normHits, decreasing = T), ]
  conf_genes <- rownames(subset(res, decision != "Rejected"))
  return(conf_genes)
}