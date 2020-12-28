\name{hs_symbol_v101}
\alias{hs_symbol_v101}
\docType{data}
\title{
Gene Ontology (GO) data of Homo sapiens
}
\description{
This data contains GO terms and genes of Homo sapiens genes, 
which are created via R \code{biomaRt} package.
The data corresponds to Ensemble Release version 101 (August 2020).
}
\usage{data("hs_symbol_v101")}
\format{
  A list of 18,763 elements.
  \describe{
    \item{\code{term}}{a character vector of GO term}
    \item{\code{genes}}{a character vector of genes}
  }
}
\examples{
data(hs_symbol_v101)
}
