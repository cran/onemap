\name{modify.rf.2pts}
\alias{modify.rf.2pts}
\title{Update two-point analyses with new criteria}
\description{
  Takes an object of class \code{rf.2pts} and updates the analyses with
  new criteria, i.e., new thresholds for LOD Score and maximum
  recombination fraction.
}
\usage{
modify.rf.2pts(w, LOD = NULL, max.rf = NULL, verbose = TRUE)
}
\arguments{
  \item{w}{an object of class \code{rf.2pts}.}
  \item{LOD}{minimum LOD Score to declare linkage.
    \code{NULL} (default) indicates that the LOD Score threshold will
    not be changed.}
  \item{max.rf}{maximum recombination fraction to declare
    linkage. \code{NULL} (default) indicates that the maximum
    recombination fraction threshold will not be changed.}
  \item{verbose}{logical. If \code{TRUE}, current progress is shown; if
    \code{FALSE}, no output is produced.}
}
\details{
  The two-point analyses are not re-done: the results of the former
  analyses are re-interpreted in light of the new criteria and a new
  object of class \code{rf.2pts} is produced for further investigation.
}
\value{
  An object of class \code{rf.2pts}, which is a list containing the
  following components:
  \item{n.mar}{total number of markers.}
  \item{LOD}{minimum LOD Score to declare linkage.}
  \item{max.rf}{maximum recombination fraction to declare linkage.}
  \item{recomb}{a (symmetric) matrix with two-point estimates of the
    recombination fraction between markers (for all pairs of markers),
    under the most probable assignment for each pair.}
  \item{phases}{a (symmetric) matrix with the most probable linkage
    phase (assignment) between markers, for each pair of markers.}
  \item{analysis}{an array with the complete results of the two-point
    analysis for each pair of markers.}
  \item{flags}{a (symmetric) matrix indicating if there is more than one
    equally probable linkage phase for each pair of markers. Possible
    values are: \code{1} if positive, \code{0} if negative and \code{NA}
    if linkage is not significant.}
  \item{arbitr}{a (symmetric) matrix indicating the most probable
    linkage phase for each pair of markers, when arbitrarily defined by
    the user. This function will actually return a matrix filled with
    zeroes, which indicates that the linkage phases were determined by
    the software (but may be arbitrarily changed afterwards).}
  \item{segr.type}{a vector with the segregation type of each marker.}
}
\references{
  Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous
  maximum likelihood estimation of linkage and linkage phases in
  outcrossing species. \emph{Theoretical Population Biology} 61:
  349-363.
}
\note{
  The thresholds used for \code{LOD} and \code{max.rf} will be used in
  subsequent analyses.
}
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br}}
\seealso{
  \code{\link[onemap]{est.rf.2pts}},
  \code{\link[onemap]{arbitr.rf.2pts}} and \code{\link[onemap]{cr2pts}}
}
\examples{
  data(example_out)

  twopts <- est.rf.2pts(example_out,LOD=3,max.rf=0.35) # perform two-point analyses
  twopts

  twopts <- modify.rf.2pts(twopts,max.rf=0.5) # same LOD Score threshold, new maximum recombination fraction criterion
  twopts
}
\keyword{utilities}
