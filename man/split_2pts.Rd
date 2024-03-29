% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{split_2pts}
\alias{split_2pts}
\title{Split rf_2pts object by markers}
\usage{
split_2pts(twopts.obj, mks)
}
\arguments{
\item{twopts.obj}{object of class rf_2pts}

\item{mks}{markers names (vector of characters) or number (vector of integers) to be removed and added to a new rf_2pts object}
}
\value{
An object of class \code{rf_2pts} with only the selected markers, which is a list containing the
following components:  \item{n.mar}{total number of markers.} \item{LOD}{minimum LOD Score to declare
linkage.} \item{max.rf}{maximum recombination fraction to declare linkage.}
}
\description{
Split rf_2pts object by markers
}
\author{
Cristiane Taniguti, \email{chtaniguti@tamu.edu}
}
