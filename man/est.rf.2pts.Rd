\name{est.rf.2pts}
\alias{est.rf.2pts}
\alias{print.rf.2pts}
\title{Two-point analysis between genetic markers}
\description{
  Performs the two-point (pairwise) analysis proposed by \cite{Wu et al.
    (2002)} between all pairs of markers.
}
\usage{
est.rf.2pts(w, LOD = 3, max.rf = 0.35, verbose = TRUE)

\method{print}{rf.2pts}(x, mrk1 = NULL, mrk2 = NULL, \dots)
}
\arguments{
  \item{w}{an object of class \code{outcross}.}
  \item{LOD}{minimum LOD Score to declare linkage (defaults to
    \code{3}).}
  \item{max.rf}{maximum recombination fraction to declare linkage
    (defaults to \code{0.35}).}
  \item{verbose}{logical. If \code{TRUE}, current progress is shown; if
    \code{FALSE}, no output is produced.}
  \item{x}{an object of class \code{rf.2pts}.}
  \item{mrk1, mrk2}{optionally, two markers can be specified. If so,
    detailed results of the two-point analysis will be printed for this
    pair. Both arguments can be numeric or character strings indicating
    the numbers/names corresponding to any markers in the input file.}
  \item{\dots}{further arguments, passed to other methods. Currently
    ignored.}
}
\details{
  For \code{n} markers, there are \deqn{\frac{n(n-1)}{2}}{n*(n-1)/2}
  pairs of markers to be analyzed. Therefore, completion of the
  two-point analyses can take quite a long time.
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
  \code{\link[onemap]{modify.rf.2pts}} to modify criteria used in the
  analyses; \code{\link[onemap]{arbitr.rf.2pts}} and
  \code{\link[onemap]{cr2pts}}
}
\examples{
  data(example_out)

  twopts <- est.rf.2pts(example_out,LOD=3,max.rf=0.5) # perform two-point analyses
  twopts

  print(twopts,"M1","M2") # detailed results for markers 1 and 2
}
\keyword{utilities}
