% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init.R
\name{initEsp}
\alias{initEsp}
\title{Initialize an \code{espresso} object.}
\usage{
initEsp(exprs, topology, asgmt = NULL)
}
\arguments{
\item{exprs}{Expression matrix (row: samples, col: features).}

\item{topology}{Adjacent matrix representing a domain toplology.}

\item{asgmt}{Domain assignment.}
}
\value{
The \code{espresso} object.
}
\description{
Initialize an \code{espresso} object 
             by \code{data}, \code{topology}, and \code{asgmt}.
}
