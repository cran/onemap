\name{make.map}
\alias{make.map}
\title{Manually define a genetic map}
\description{
  Allows the user to manually define a genetic map, i.e., a set of
  ordered markers on a linkage group. Genetic distances are drawn from
  two-point analyses.
}
\usage{
make.map(w, ordered.group = NULL)
}
\arguments{
  \item{w}{an object of class \code{rf.2pts}, which contains the results
    of two-point analyses.}
  \item{ordered.group}{a vector of character strings indicating the
    names of the \emph{ordered} markers, according to the original input
    file. The order of input of markers corresponds to the genetic map.}
}
\value{
  An object of class \code{map}, which is a list with the following
  components:
  \item{order}{a numeric vector with the arbitrarily given order of
    genetic markers.}
  \item{recomb}{a (symmetric) matrix with two-point estimates of the
    recombination fraction between markers, for all pairs of markers in
    the linkage group, under the most probable assignment for each
    pair.}
  \item{marnames}{names of markers present in the linkage group,
    according to the original input file.}
  \item{number}{number of the linkage group. In this case, the value is
    \code{-1} to indicate that the genetic map was arbitrarily defined.}
  \item{name}{name of the object of class \code{group} from where the
    linkage group was extracted. In this case, the value is \code{""}.}
  \item{LOD}{minimum LOD Score to declare linkage. In this  case, the
    value is \code{-1} to indicate that the genetic map was arbitrarily
    defined.}
  \item{max.rf}{maximum recombination fraction to declare linkage. In
    this case, the value is \code{-1} to indicate that the genetic map
    was arbitrarily defined.}
  \item{phases}{a (symmetric) matrix with the most probable linkage
    phase (assignment) between markers, for each pair of markers in the
    linkage group.}
}
\references{
  Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous
  maximum likelihood estimation of linkage and linkage phases in
  outcrossing species. \emph{Theoretical Population Biology} 61:
  349-363.
}
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br}}
\seealso{
  \code{\link[onemap]{map}}, \code{\link[onemap]{print.map}}
}
\examples{
  data(example_out)
  twopts <- est.rf.2pts(example_out)
  
  my_map <- c("M18","M8","M13","M22","M7") # map for the third linkage group (details in tutorial)
  user_map <- make.map(twopts,my_map)
  user_map # must be complemented with three-point analyses
}
\keyword{utilities}
