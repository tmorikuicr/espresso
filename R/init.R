# =============================================================================
#' @title Initialize an \code{espresso} object.
#' @description Initialize an \code{espresso} object 
#'              by \code{data}, \code{topology}, and \code{asgmt}.
#' @param exprs Expression matrix (row: samples, col: features).
#' @param topology Adjacent matrix representing a domain toplology.
#' @param asgmt Domain assignment.
#' @return The \code{espresso} object.
#' @importFrom methods new
#' @export
#' 
initEsp <- function(exprs, topology, asgmt = NULL) {
  obj <- new('Espresso')
  obj@exprs <- as.matrix(exprs)
  obj@gset <- colnames(obj@exprs)
  obj@topology <- as.matrix(topology)
  domain_id <- 1:nrow(obj@topology)
  domain_name <- rownames(obj@topology)
  obj@domain2name <- data.frame("domain" = domain_id, 
                                "name" = domain_name, 
                                stringsAsFactors = FALSE)
  if (!is.null(asgmt)) {
      obj@asgmt <- asgmt
  }
  return(obj)
}
