## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  # Installation (X.X.X is a version number)
#  install.packages("espresso_X.X.X.tar.gz", repos = NULL, type = "source")

## ---- eval = FALSE------------------------------------------------------------
#  # Installation example of CRAN packages
#  install.packages(c("igraph", "mclust", "progress",
#                     "umap", "rgl", "rFerns", "Boruta",
#                     "doParallel", "foreach", "pheatmap", "RColorBrewer"),
#                   dependencies = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  # Installation example of Bioconductor packages
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install(c("biomaRt","GO.db"))

## -----------------------------------------------------------------------------
version
packageVersion("igraph")
packageVersion("mclust")
packageVersion("progress")
packageVersion("umap")
packageVersion("rgl")
packageVersion("rFerns")
packageVersion("Boruta")
packageVersion("doParallel")
packageVersion("foreach")
packageVersion("pheatmap")
packageVersion("RColorBrewer")
packageVersion("biomaRt")
packageVersion("GO.db")

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(espresso)
set.seed(0)

## ---- eval = FALSE------------------------------------------------------------
#  library(espresso)
#  
#  # Initialize espresso object
#  eobj <- initEsp(exprs, topology, asgmt)
#  
#  # GraphSOM clustering
#  eobj <- initGraphSOM(eobj)
#  eobj <- graphSOM(eobj)
#  
#  # Show prediction scores
#  eobj@summary

## -----------------------------------------------------------------------------
# Get paths of input files
path_exprs <- system.file("extdata", "exprs_soysa19_e825.txt", package="espresso")
path_topology <- system.file("extdata", "topology_soysa19_e825.txt", package="espresso")
path_asgmt <- system.file("extdata", "asgmt_soysa19_e825.txt", package="espresso")

# Load data
exprs <- t(read.table(path_exprs, header = TRUE, row.names = 1, sep = "\t", 
                      as.is = TRUE, check.names = FALSE))
topology <- as.matrix(read.table(path_topology, header = T, row.names = 1))
asgmt <- read.table(path_asgmt, header = T, as.is = TRUE)

# Show the contents of the objects
exprs[1:5, 1:5]
topology
head(asgmt)

## -----------------------------------------------------------------------------
eobj <- initEsp(exprs, topology, asgmt)

## -----------------------------------------------------------------------------
str(eobj)

## ---- eval = FALSE------------------------------------------------------------
#  # Log-scaling of the expression profile
#  # In this example, log-scalling is unnecessary since the expression values is already log-scalled.
#  # eobj <- logScale(eobj)

## -----------------------------------------------------------------------------
# Filter out genes with low expresssion and low variance
eobj <- filterGenes(eobj)

## ---- eval = FALSE------------------------------------------------------------
#  # Filter out genes with low expresssion and low variance
#  eobj <- filterGenes(eobj, ncell = 2, expressed = 1.0, sd = 0.05)

## ---- eval = FALSE------------------------------------------------------------
#  eobj <- initGraphSOM(eobj)

## -----------------------------------------------------------------------------
eobj <- initGraphSOM(eobj, nsamples = 10, coef = 0.5, rept = 3)

## ---- eval = FALSE------------------------------------------------------------
#  eobj <- graphSOM(eobj)

## -----------------------------------------------------------------------------
eobj <- graphSOM(eobj, gset = c("Actn2", "Hand2", "Vegfb", "Wnt5a"))

## -----------------------------------------------------------------------------
eobj@summary

## -----------------------------------------------------------------------------
eobj@score

## -----------------------------------------------------------------------------
plotConvCurve(eobj, rept = 1:2)

## ---- eval = FALSE------------------------------------------------------------
#  par(mfrow = c(2,2))
#  plotConvCurve(eobj)
#  dev.off()

## -----------------------------------------------------------------------------
fgenes <- selectFeatures(eobj, maxRuns = 100)

fgenes

## -----------------------------------------------------------------------------
# Mouse gene sets (Ensembl release version 101)
data(mm_symbol_v101)

# Human gene sets (Ensembl release version 101)
data(hs_symbol_v101)

## ---- eval = FALSE------------------------------------------------------------
#  # Get GO infomation
#  hs_symbol_v98 <- getGO(species = 'hsapiens', gid = 'hgnc_symbol', version = '98')
#  
#  # Save 'hs_symbol_v98' as an R data file.
#  save(hs_symbol_v98, file = "hs_symbol_v98")
#  
#  # Load 'hs_symbol_v98'
#  load(hs_symbol_v98)

## -----------------------------------------------------------------------------
mm_symbol_v101[["GO:0003209"]][["term"]]

mm_symbol_v101[["GO:0003209"]][["genes"]]

## ---- eval = FALSE------------------------------------------------------------
#  eobj <- graphSOM(eobj, gset = mm_symbol_v101[["GO:0003209"]][["genes"]],
#                   nsamples = 10, rept = 3)

## ---- message = FALSE---------------------------------------------------------
eobj <- initGraphSOM(eobj, nsamples = 10, rept = 1, coef = 0.5)
eobj <- rxmcmc(eobj, gset = fgenes, itr = 2, n_ex = 5, n_repl = 4)

## ---- eval = FALSE------------------------------------------------------------
#  eobj <- rxmcmc(eobj, gset = fgenes, itr = 2, n_ex = 5, n_repl = 4, n_cl = 4, n_ig = 3)

## -----------------------------------------------------------------------------
eobj@mcmc[["sampling"]][["repl.1"]]
eobj@mcmc[["max"]][["repl.1"]]
eobj@mcmc[["exchange"]]
eobj@mcmc[["best_genes"]]

## ---- fig.width = 8, fig.height = 8, out.width = "45%"------------------------
plotMCMC(eobj)

## ---- fig.width = 8, fig.height = 8, out.width = "45%"------------------------
plotDistMap(eobj)

## ---- eval = FALSE------------------------------------------------------------
#  plotUMAP(eobj)

## -----------------------------------------------------------------------------
plotUMAP(eobj, gset = eobj@mcmc[["best_genes"]],
         n_neighbors = 5, n_components = 3, n_epochs = 1000)

## ---- echo=FALSE, fig.cap= "UMAP", out.width = '50%'--------------------------
knitr::include_graphics("umap.png")

## ---- eval = FALSE------------------------------------------------------------
#  domain_colors <- list(AHF = "orange",
#                        Atrial = "pink",
#                        LV = "green",
#                        OFT = "lightblue",
#                        pSHF = "purple",
#                        RV = "blue",
#                        SV = "plum"
#                        )
#  plotUMAP(eobj, gset = eobj@mcmc[["best_genes"]],
#           n_neighbors = 5, n_components = 3, n_epochs = 1000,
#           movie = TRUE, domcol = domain_colors)
#  

## ---- eval = FALSE------------------------------------------------------------
#  eobj.tmp <- eobj
#  samples <- eobj.tmp@ssets[["set.1"]]
#  eobj.tmp@exprs <- eobj.tmp@exprs[match(samples, rownames(eobj.tmp@exprs)), , drop = FALSE]
#  plotUMAP(eobj.tmp, gset = eobj.tmp@mcmc[["best_genes"]],
#           n_neighbors = 5, n_components = 3, n_epochs = 1000,
#           movie = TRUE, domcol = domain_colors)

## ---- echo=FALSE, fig.cap= "UMAP", out.width = '50%'--------------------------
knitr::include_graphics("umap.all.png")

## -----------------------------------------------------------------------------
library(biomaRt)
fgenes <- c("Actn2", "Hand2", "Vegfb", "Wnt5a")
ds <- useEnsembl(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
ann <- getBM(attributes = c("mgi_symbol", "description"), 
             mart = ds, filters = "mgi_symbol", values = fgenes)
ann

## ---- eval=FALSE--------------------------------------------------------------
#  ann <- getBM(attributes = c("mgi_symbol", "description"),
#               mart = ds, filters = "mgi_symbol", values = fgenes, useCache = FALSE)

