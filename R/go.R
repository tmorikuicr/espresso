# =============================================================================
#' @title Get Gene Ontology terms.
#' @description This function gets GO terms via BioMart database.
#' @param ds Dataset to be accessed via biomaRt.
#' @importFrom biomaRt getBM
#' @return Data frame whose the first and second columns indicate 
#'         GO ID and terms, respectively.
#' 
.getGO2Term <- function(ds){
  message("Downloading GO IDs and terms...")
  go2term <- getBM(attributes = c("go_id", "name_1006"), mart = ds)
  go2term <- go2term[order(go2term$go_id), ]
  go2term <- subset(go2term, subset = go_id != "")
  go2term <- unique(go2term)
  rownames(go2term) <- 1:nrow(go2term)
  return(go2term)
}

# =============================================================================
#' @title Get offsprings of each Gene Ontology.
#' @description This function gets GO and 
#'              thier offspring information via \code{GO.db}.
#' @param terms Data frame of GO terms.
#' @return List of GOs and their offsprings. 
#' @importFrom GO.db GOBPOFFSPRING GOMFOFFSPRING GOCCOFFSPRING
#' @import progress
#' 
.getGO2Offsprings <- function(terms){
  message("Searching GO offsprings...")
  go2offsprings <- vector("list", length = nrow(terms))
  names(go2offsprings) <- terms$go_id
  cnt_end <- nrow(terms)
  pb <- progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = nrow(terms), clear = FALSE, width= 60)
  for (g in terms$go_id) {
    offsprings <- unique(c(GOBPOFFSPRING[[g]], 
                           GOMFOFFSPRING[[g]], 
                           GOCCOFFSPRING[[g]]))
    go2offsprings[[g]] <- offsprings
    pb$tick()
  }
  return(go2offsprings)
}

# =============================================================================
#' @title Get gene information.
#' @description This function gets gene information for each GO via BioMart database.
#' @param obj The \code{espresso} object.
#' @param ds Dataset to be accessed via biomaRt.
#' @return List of GO and genes
#' @importFrom biomaRt getBM
#' @importFrom stats na.omit
#' @import progress
#' 
.getGO2Genes <- function(obj, ds){
  go2genes <- NULL
  terms <- .getGO2Term(ds)
  offsprings <- .getGO2Offsprings(terms)
  message("Downloading GO genes...")
  go2genes_tmp <- getBM(attributes = c("go_id", obj$gid), mart = ds)
  go2genes_tmp <- go2genes_tmp[order(go2genes_tmp$go_id), ]
  go2genes_tmp <- subset(go2genes_tmp, subset = go_id != "")
  rownames(go2genes_tmp) <- 1:nrow(go2genes_tmp)
  colnames(go2genes_tmp) <- c("go_id", "gene")
  go2genes <- vector("list", length = nrow(terms))
  names(go2genes) <- terms$go_id
  cnt_end <- nrow(terms)
  message("Generating gene sets...")
  pb <- progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = length(terms$go_id), clear = FALSE, width= 60)
  for (g in terms$go_id) {
    gos <- unique(c(g, na.omit(offsprings[[g]])))
    genes <- subset(go2genes_tmp, subset = go_id %in% gos)$gene
    genes <- genes[!is.na(genes)]
    go2genes[[g]] <- list(term = terms[match(g, terms$go_id),]$name_1006, 
                          genes = sort(unique(genes)))
    pb$tick()
  }
  return(go2genes)
}

# =============================================================================
#' @title Get Gene Ontology information.
#' @description This function gets GO information via 
#'              BioMart and Ensembl databases.
#' @param species Species name used in biomaRt (e.g., mmusculus, hsapiens).
#' @param gid Gene identifier (e.g., mgi_symbol, hgnc_symbol).
#' @param version Ensembl release version. 
#'                If \code{version = NULL}, the latest version is downloaded.
#' @return List of GO information.
#' @importFrom biomaRt useEnsembl getBM
#' @export
#' 
getGO <- function(species, gid, version = NULL){
  obj <- list(species = species, gid = gid)
  if (is.null(version)) {
    message("Use Ensembl database (latest release) via biomaRt.")
  } else {
    message(paste0("Use Ensembl database (release", version, ") via biomaRt."))
  }
  if (is.null(version)) {
    ds <- useEnsembl(biomart = "ensembl", 
                     dataset = paste0(species, "_gene_ensembl"))
  } else {
    ds <- useEnsembl(biomart = "ensembl", 
                     dataset = paste0(species, "_gene_ensembl"), 
                     version = version)
  }
  go2genes <- .getGO2Genes(obj, ds)
  return(go2genes)
}