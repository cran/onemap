\name{arbitr.rf.2pts}
\alias{arbitr.rf.2pts}
\title{Arbitrarily define the linkage phase between two markers}
\description{
  Superimposes the most probable linkage phase determined by the
  software with an \emph{ad hoc} value chosen by the user, for any two
  markers in an object of class \code{rf.2pts}.
}
\usage{
arbitr.rf.2pts(w, mrk1name = NULL, mrk2name = NULL, new.phase = NULL)
}
\arguments{
  \item{w}{an object of class \code{rf.2pts}.}
  \item{mrk1name}{a character string indicating the name of the first
    marker, corresponding to any marker in the input file used for the
    two-point analyses stored in object \code{w}.}
  \item{mrk2name}{a character string indicating the name of the second
    marker.}
  \item{new.phase}{a character string indicating the new linkage phase
    to be assigned between the two markers. Possible values are
    (\cite{Wu et al., 2002}): \code{"A1"} for coupling/coupling;
    \code{"A2"} for coupling/repulsion; \code{"A3"} for
    repulsion/coupling; \code{"A4"} for repulsion/repulsion; \code{"ns"}
    for not significant.}
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
    the user.}
  \item{segr.type}{a vector with the segregation type of each marker.}
}
\references{
  Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous
  maximum likelihood estimation of linkage and linkage phases in
  outcrossing species. \emph{Theoretical Population Biology} 61:
  349-363.
}
\note{
  Intended to be used if other sources of information are
  available. \emph{Use this function very carefully}.
}
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br}}
\seealso{
  \code{\link[onemap]{est.rf.2pts}}, \code{\link[onemap]{print.rf.2pts}},
  \code{\link[onemap]{modify.rf.2pts}} and \code{\link[onemap]{cr2pts}}
}
\examples{
  data(example_out)

  twopts <- est.rf.2pts(example_out)
  print(twopts,"M1","M2")
  
  # these arbitrary definitions are likely to be unnecessary
  twopts <- arbitr.rf.2pts(twopts,"M1","M2","A3") # repulsion/coupling
  print(twopts,"M1","M2")

  twopts <- arbitr.rf.2pts(twopts,"M1","M2","ns") # non-significant
  print(twopts,"M1","M2")

  twopts <- arbitr.rf.2pts(twopts,"M1","M2","A1") # back to coupling/coupling
  print(twopts,"M1","M2")
}
\keyword{manip}
\keyword{utilities}
