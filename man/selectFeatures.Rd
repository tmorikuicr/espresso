% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boruta.R
\name{selectFeatures}
\alias{selectFeatures}
\title{Select important features by Boruta.}
\usage{
selectFeatures(obj, maxRuns = 500, decision = "nr", seed = NULL)
}
\arguments{
\item{obj}{The \code{espresso} object.}

\item{maxRuns}{The maximum number of importance source runs.}

\item{decision}{Decision type of Boruta 
("c": Confirmed, "nr": Not Rejected (i.e., Confirmed & Tentative)
(default: "nr").}

\item{seed}{Random seed.}
}
\value{
Data Frame of the result of feature selection.
}
\description{
This function gets important features by \code{Boruta} package.
}
