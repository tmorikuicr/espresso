#' Class to define a GraphSOM object
#' 
#' Define a GraphSOM object to be trained.
#' 
#' @slot exprs expression matrix (row: samples, col: features).
#' @slot exprs_som data matrix for graphSOM computation (row: samples, col: features).
#' @slot rep_som replication number for GraphSOM.
#' @slot topology adjacent matrix representing a domain toplology.
#' @slot dist distance matrix corresponding to \code{topology}.
#' @slot map map vectors given as a matrix.
#' @slot map_method map initialization method.
#' @slot radius initial value of learning radius (default: the maximum distance of an input toplology).
#' @slot lstep the number of learning steps (default: #sample * 2).
#' @slot rep the number of replications for GraphSOM clustering.
#' @slot bmu the best matchiing units through the learnig steps.
#' @slot bmu_method bmu selection method.
#' @slot stochastic if \code{stochastic = TRUE}, the stochastic strategy is introduced to map update.
#' @slot rmin the minimum value of scale factor for the stochastic map update (default: 0.5). 
#'            If \code{stochastic = FALSE}, \code{rmin} is ignored.
#' @slot rmax the maximum value of scale factor for the stochastic map update (default: 1.0). 
#'            If \code{stochastic = FALSE}, \code{rmax} is ignored.
#' @slot domain2name corresponding table between domain IDs and their nemaes.
#' @slot gset feature gene set.
#' @slot ssets sample sets.
#' @slot asgmt domain assignment.
#' @slot score GraphSOM clustering results.
#' @slot summary summary of GraphSOM clustering results.
#' @slot mcmc history of MCMC results.
#' @slot seed random seed.
#' @slot coef coefficient of ARI for score computation.
#' @slot nmin the minimum number of genes can be analyzed.
#' @slot nmax the maximum number of genes can be analyzed.
#' 
Espresso <- setClass("Espresso",
         slots = list(
           exprs = 'matrix',
           exprs_som = 'matrix',
           rep_som = 'integer',
           topology = 'matrix', 
           dist = 'matrix', 
           map = 'list',
           map_method = 'character',
           radius = 'numeric',
           lsteps = 'integer',
           rep = 'integer',
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
           nmax = 'integer'
         ), 
         prototype = list(
           exprs = NULL,
           exprs_som = NULL,
           rep_som = NULL,
           topology = NULL,
           dist = NULL, 
           map = NULL,
           map_method = NULL,
           radius = NULL,
           lsteps = NULL,
           rep = NULL,
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
           nmax = NULL
         )
)