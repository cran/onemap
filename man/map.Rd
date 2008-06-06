\name{map}
\alias{map}
\alias{print.map}
\title{Estimate the genetic map}
\description{
  Does soma data manipulation and calls \code{\link[onemap]{rcd}} to
  order genetic markers on a linkage group (object of class5
  \code{extracted.group}).
}
\usage{
map(w, y, LOD = NULL, max.rf = NULL)

\method{print}{map}(x, cumulative=FALSE, \dots)
}
\arguments{
  \item{w}{an object of class \code{extracted group} (the group to be
    ordered).}
  \item{y}{the object of class \code{rf.2pts} that was used to generate
    object \code{w}, i.e., that contains the results of two-point
    analyses.}
  \item{LOD}{minimum LOD Score to declare linkage. If \code{NULL}
    (default) the threshold already present in object \code{y} is used.}
  \item{max.rf}{maximum recombination fraction to declare linkage. If
    \code{NULL} (default) the threshold already present in object
    \code{y} is used.}
  \item{x}{an object of class \code{map}.}
  \item{cumulative}{logical. If \code{FALSE} (default), distances
    between markers are printed independently for each interval. If
    \code{TRUE}, cumulative distances are displayed (for compatibility
    purposes with softwares used to draw chromosomes).}
  \item{\dots}{further arguments, passed to other methods. Currently
    ignored.}
}
\value{
  An object of class \code{map}, which is a list with the following
  components:
  \item{order}{a numeric vector with the estimated order of genetic
    markers.}
  \item{recomb}{a (symmetric) matrix with two-point estimates of the
    recombination fraction between markers, for all pairs of markers in
    the linkage group, under the most probable assignment for each
    pair.}
  \item{marnames}{names of markers present in the linkage group,
    according to the original input file.}
  \item{number}{number of the linkage group, according to the numeration
    in the object of class \code{group} from where the linkage group was
    extracted.}
  \item{name}{name of the object of class \code{group} from where the
    linkage group was extracted.}
  \item{LOD}{minimum LOD Score to declare linkage.}
  \item{max.rf}{maximum recombination fraction to declare linkage.}
  \item{phases}{a (symmetric) matrix with the most probable linkage
    phase (assignment) between markers, for each pair of markers in the
    linkage group.}
}
\note{
  Currently only \emph{Rapid Chain Delineation} is available
  (\cite{Doerge, 1996}). Other ordering algorithms will be made
  available soon.
}
\references{
  Doerge, R. W. (1996) Constructing genetic maps by rapid chain
  delineation. \emph{Journal of Quantitative Trait Loci} 2: 121-132.
  
  Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous
  maximum likelihood estimation of linkage and linkage phases in
  outcrossing species. \emph{Theoretical Population Biology} 61:
  349-363.
}
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br}}
\seealso{
  \code{\link[onemap]{group}}, \code{\link[onemap]{extract.group}},
  \code{\link[onemap]{rcd}}, \code{\link[onemap]{make.map}}
}
\examples{
  data(example_out)
  
  twopts <- est.rf.2pts(example_out)  
  link_gr <- group(twopts)
  LG1 <- extract.group(link_gr,1)

  LG1_map <- map(LG1,twopts) # order makers
  LG1_map
}
\keyword{manip}
\keyword{utilities}

