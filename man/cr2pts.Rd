\name{cr2pts}
\alias{cr2pts}
\title{Call to the C routine for two-point analysis}
\description{
  An interface between \R and the C routine that performs the two-point
  analysis. Given the type of segregation of the two markers, this
  function creates the corresponding vectors and matrices that are
  necessary to the C routine and makes a call to it.
}
\usage{
cr2pts(mrk1, mrk2, segr.type1, segr.type2)
}
\arguments{
  \item{mrk1}{a vector containing the genotypes of all individuals for
    the first marker.}
  \item{mrk2}{a vector containing the genotypes for the second marker.}
  \item{segr.type1}{the type of segregation of the first marker
    (\code{mrk1}).}
  \item{segr.type2}{the type of segregation of the second marker
    (\code{mrk2}).}
}
\value{
  Returns a list with the following components:
  \item{r}{a vector with the maximum likelihood estimates of the
    recombination fraction between both markers, under the \code{4}
    possible assignments.}
  \item{like}{a vector with the likelihood of the data under each of the
    \code{4} assignments.}
  \item{posterior}{a vector with the posterior probability of each of
    the \code{4} assignments.}
  \item{LOD}{a vector with the LOD Score for the \code{4} assignments -
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
  \code{\link[onemap]{est.rf.2pts}}; for details concerning the types of
    segregation, see \code{\link[onemap]{mrktype}}.
}
\examples{
  data(example_out)

  # getting values
  marker1 <- example_out$geno[,1]
  marker2 <- example_out$geno[,2]
  segr.type1 <- example_out$segr.type[1]
  segr.type2 <- example_out$segr.type[2]
  
  # two-point analysis
  cr2pts(marker1,marker2,segr.type1,segr.type2)
}
\keyword{manip}
\keyword{interface}
