\name{mm_symbol_v101}
\alias{mm_symbol_v101}
\docType{data}
\title{
Gene Ontology (GO) data of Mus musculus
}
\description{
This data contains GO terms and genes of Mus musculus genes, 
which are created via R \code{biomaRt} package.
The data corresponds to Ensemble Release version 101 (August 2020).
}
\usage{data("mm_symbol_v101")}
\format{
  A list of 18,583 elements.
  \describe{
    \item{\code{term}}{a character vector of GO term}
    \item{\code{genes}}{a character vector of genes}
  }
}
\examples{
data(mm_symbol_v101)
}
