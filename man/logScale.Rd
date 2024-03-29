% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{logScale}
\alias{logScale}
\title{Scale gene expression data by log function.}
\usage{
logScale(obj, base = "10")
}
\arguments{
\item{obj}{The \code{espresso} object.}

\item{base}{Base of a logarithm specified 
as character ('2', '10', or 'e') (default: '10').}
}
\value{
\code{espresso} object
}
\description{
This function scales gene expression data 
             by log function with bases 2, 10, and e.
}
