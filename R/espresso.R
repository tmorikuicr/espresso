#' Class to define a GraphSOM object
#' 
#' Define a GraphSOM object to be trained.
#' 
#' @slot exprs Expression matrix (row: samples, col: features).
#' @slot exprs_som Data matrix for graphSOM computation (row: samples, col: features).
#' @slot rept_som Repeat number for GraphSOM clustering.
#' @slot topology Adjacent matrix representing the toplology of domains.
#' @slot dist Distance matrix corresponding to \code{topology}.
#' @slot map Values of the map's vector given by a matrix.
#' @slot map_method Method of map initialization.
#' @slot radius Initial value of learning radius.
#' @slot lsteps Number of learning steps.
#' @slot rept Number of repeats for GraphSOM clustering.
#' @slot bmu Best matchiing units (BMU) through the learnig steps.
#' @slot bmu_method BMU selection method.
#' @slot stochastic If \code{stochastic = TRUE}, the stochastic strategy is introduced to map update.
#' @slot rmin Minimum value of scale factor for the stochastic map update. 
#'            If \code{stochastic = FALSE}, \code{rmin} is ignored.
#' @slot rmax Maximum value of scale factor for the stochastic map update. 
#'            If \code{stochastic = FALSE}, \code{rmax} is ignored.
#' @slot domain2name Data frame indicating the domain IDs and their names.
#' @slot gset Feature gene set.
#' @slot ssets Sample sets.
#' @slot asgmt Domain assignment.
#' @slot score Score of GraphSOM clustering results.
#' @slot summary Summary of GraphSOM clustering results.
#' @slot mcmc History of replica exchange MCMC results.
#' @slot seed Random seed.
#' @slot coef Coefficient of ARI for score computation.
#' @slot nmin Minimum number of genes can be analyzed.
#' @slot nmax Maximum number of genes can be analyzed.
#' @slot swap Logical value determins whether to swap domains or not.
#' @slot umap UMAP result.
#' 
Espresso <- setClass("Espresso",
         slots = list(
           exprs = 'matrix',
           exprs_som = 'matrix',
           rept_som = 'integer',
           topology = 'matrix', 
           dist = 'matrix', 
           map = 'list',
           map_method = 'character',
           radius = 'numeric',
           lsteps = 'integer',
           rept = 'integer',
           bmu = 'list',
           bmu_method = 'character',
           stochastic = 'logical', 
           rmin = 'numeric',
           rmax = 'numeric',
           domain2name = 'data.frame',
           gset = 'vector',
           ssets = 'list',
           asgmt = 'data.frame',
           score = 'list',
           summary = 'list',
           mcmc = 'list',
           seed = 'numeric',
           coef = 'numeric',
           nmin = 'integer',
           nmax = 'integer',
           swap = 'logical',
           umap = 'matrix'
         ), 
         prototype = list(
           exprs = NULL,
           exprs_som = NULL,
           rept_som = NULL,
           topology = NULL,
           dist = NULL, 
           map = NULL,
           map_method = NULL,
           radius = NULL,
           lsteps = NULL,
           rept = NULL,
           bmu = NULL,
           bmu_method = NULL,
           stochastic = NULL,
           rmin = NULL,
           rmax = NULL,
           domain2name = NULL, 
           gset = NULL,
           ssets = NULL,
           asgmt = NULL,
           score = NULL,
           summary = NULL,
           mcmc = NULL,
           seed = NULL,
           coef = NULL,
           nmin = NULL,
           nmax = NULL,
           swap = NULL,
           umap = NULL
         )
)