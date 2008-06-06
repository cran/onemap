\name{cr3pts}
\alias{cr3pts}
\title{Call to the C routine for three-point analysis}
\description{
  An interface between \R and the C routine that performs the
  three-point analysis. Given the type of segregation of the three
  markers, this function creates the corresponding vectors and matrices
  that are necessary to the C routine and makes a call to it.
}
\usage{
cr3pts(mrk1, mrk2, mrk3, segr.type1, segr.type2, segr.type3)
}
\arguments{
  \item{mrk1}{a vector containing the genotypes of all individuals for
    the first marker.}
  \item{mrk2}{a vector containing the genotypes for the second marker.}
  \item{mrk3}{a vector containing the genotypes for the third marker.}
  \item{segr.type1}{the type of segregation of the first marker
    (\code{mrk1}).}
  \item{segr.type2}{the type of segregation of the second marker
    (\code{mrk2}).}
  \item{segr.type3}{the type of segregation of the third marker
    (\code{mrk3}).}
}
\value{
  Returns a list with the following components:
  \item{theta12}{a vector with the maximum likelihood estimates (MLEs)
    of the recombination fraction between markers \code{1} and \code{2}
    (\code{mrk1} and \code{mrk2}, respectively), under the \code{16}
    possible assignments.}
  \item{theta23}{a vector with the MLEs of the recombination fraction
    between markers \code{2} and \code{3} (\code{mrk2} and \code{mrk3},
    respectively), under the \code{16} assignments.}
  \item{theta13}{a vector with the MLEs of the recombination fraction
    between markers \code{1} and \code{3} (\code{mrk1} and \code{mrk3},
    respectively), under the \code{16} assignments.}
  \item{like}{a vector with the likelihood of the data under each of the
    \code{16} assignments.}
  \item{posterior}{a vector with the posterior probability of each of
    the \code{16} assignments.}
  \item{LOD}{a vector with the LOD Score for the \code{16} assignments -
    the base 10 logarithm of the ratio of the maximum likelihood and the
    likelihood under the null hypothesis of no linkage.}
}
\references{
  Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous
  maximum likelihood estimation of linkage and linkage phases in
  outcrossing species. \emph{Theoretical Population Biology} 61:
  349-363.
}
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br}}
\seealso{
  \code{\link[onemap]{est.rf.3pts}}; for details concerning the types of
  segregation, see \code{\link[onemap]{mrktype}}.
}
\examples{
  data(example_out)
  
  # getting values
  marker1 <- example_out$geno[,1]
  marker2 <- example_out$geno[,2]
  marker3 <- example_out$geno[,14]
  segr.type1 <- example_out$segr.type[1]
  segr.type2 <- example_out$segr.type[2]
  segr.type3 <- example_out$segr.type[14]
  
  # three-point analysis
  cr3pts(marker1,marker2,marker3,segr.type1,segr.type2,segr.type3)
}
\keyword{manip}
\keyword{interface}
