% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_mapmaker.R
\name{read_mapmaker}
\alias{read_mapmaker}
\title{Read data from a Mapmaker raw file}
\usage{
read_mapmaker(file = NULL, dir = NULL, verbose = TRUE)
}
\arguments{
\item{file}{the name of the input file which contains the data to be read.}

\item{dir}{directory where the input file is located.}

\item{verbose}{A logical, if TRUE it output progress status
information.}
}
\value{
An object of class \code{onemap}, i.e., a list with the following
components: \item{geno}{a matrix with integers indicating the genotypes
read for each marker in \code{onemap} fashion. Each column contains data
for a marker and each row represents an individual.}

\code{MAPMAKER/EXP} fashion, i.e., 1, 2, 3: AA, AB, BB, respectively; 3, 4:
BB, not BB, respectively; 1, 5: AA, not AA, respectively. Each column
contains data for a marker and each row represents an individual.

\item{n.ind}{number of individuals.} \item{n.mar}{number of markers.}
\item{segr.type}{a vector with the segregation type of each marker, as
\code{strings}. Segregation types were adapted from outcross segregation
types, using the same notation. For details see \link{read_onemap}.}
\item{segr.type.num}{a vector with the segregation type of each marker,
represented in a simplified manner as integers. Segregation types were
adapted from outcross segregation types. For details see
\link{read_onemap}.} \item{input}{the name of the input file.}
\item{n.phe}{number of phenotypes.} \item{pheno}{a matrix with phenotypic
values.  Each column contains data for a trait and each row represents an
individual. Currently ignored.} \item{error}{matrix containing HMM emission probabilities}
}
\description{
Imports data from a Mapmaker raw file.
}
\details{
For details about MAPMAKER files see \cite{Lincoln et al.} (1993). The
current version supports backcross, F2s and RIL populations. The file
can contain phenotypic data, but it will not be used in the analysis.
}
\examples{
\donttest{
 map_data <-read_mapmaker(file=system.file("extdata/mapmaker_example_f2.raw", package = "onemap"))
 #Checking 'mapmaker_example_f2'
 data(mapmaker_example_f2)
 names(mapmaker_example_f2)
}
}
\references{
Broman, K. W., Wu, H., Churchill, G., Sen, S., Yandell, B.
(2008) \emph{qtl: Tools for analyzing QTL experiments} R package version
1.09-43

Lincoln, S. E., Daly, M. J. and Lander, E. S. (1993) Constructing genetic
linkage maps with MAPMAKER/EXP Version 3.0: a tutorial and reference
manual. \emph{A Whitehead Institute for Biomedical Research Technical
Report}.
}
\seealso{
\code{mapmaker_example_bc} and \code{mapmaker_example_f2} raw files in the
package source.
}
\author{
Adapted from Karl Broman (package \pkg{qtl}) by Marcelo Mollinari,
\email{mmollina@usp.br}
}
\keyword{IO}
