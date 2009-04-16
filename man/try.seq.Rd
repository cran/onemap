\name{try.seq}
\alias{try.seq}
\alias{print.try}
\title{Try to map a marker into every possible position between markers in a given sequence}
\description{
  For a given linkage map and an additional marker, this function estimates
  parameters for all possible maps including the new marker, while keeping the
  base linkage map inaltered.
}
\usage{
try.seq(w,mrk,tol=10E-5,verbose=FALSE)

\method{print}{try}(x,j=NULL, \dots)
}
\arguments{
  \item{w}{an object of class \code{sequence} with a predefined order.}
  \item{mrk}{the index of the marker to be tried, according to the input file.}
  \item{tol}{tolerance for the C routine, i.e., the value used to evaluate
    convergence.}
  \item{verbose}{if \code{FALSE} (default), simplified output is displayed. if
    \code{TRUE}, detailed output is displayed.}
  \item{x}{an object of class \code{try}.}
  \item{j}{if \code{NULL} (default), output is a summary of the results for
    all possible positions of the additional marker. Otherwise, an integer
    makes detailed output to be printed for the corresponding position. This
    integer must be less than or equal to the length of the original sequence plus
    1.}
  \item{\dots}{further arguments, passed to other methods. Currently
    ignored.}
}
\value{
  An object of class \code{try}, which is a list containing the following
  components:
  \item{ord}{a \code{list} containing results for every linkage map estimated.
    These results include linkage phases, recombination frequencies and
    log-likelihoods.}
  \item{LOD}{a \code{vector} with LOD-Scores for each position where the
    additional marker is placed. This Score is based on the best combination of
    linkage phases for each map.}
  \item{try.ord}{a \code{matrix} with the orders of all linkage maps.}
  \item{data.name}{name of the object of class \code{outcross} with the
    raw data.}
  \item{twopt}{name of the object of class \code{rf.2pts} with the 2-point
    analyses.}
}
\references{ 
  Broman, K. W., Wu, H., Churchill, G., Sen, S., Yandell, B. (2008)
  \emph{qtl: Tools for analyzing QTL experiments} R package version 1.09-43
  
  Jiang, C. and Zeng, Z.-B. (1997). Mapping quantitative trait loci with
  dominant and missing markers in various crosses from two inbred lines.
  \emph{Genetica} 101: 47-58.

  Lander, E. S., Green, P., Abrahamson, J., Barlow, A., Daly, M. J., Lincoln,
  S. E. and Newburg, L. (1987) MAPMAKER: An interactive computer package for
  constructing primary genetic linkage maps of experimental and natural
  populations. \emph{Genomics} 1: 174-181.

  Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002a) Simultaneous
  maximum likelihood estimation of linkage and linkage phases in
  outcrossing species. \emph{Theoretical Population Biology} 61:
  349-363.
  
  Wu, R., Ma, C.-X., Wu, S. S. and Zeng, Z.-B. (2002b). Linkage
  mapping of sex-specific differences. \emph{Genetical Research}
  79: 85-96
}
\author{Marcelo Mollinari, \email{mmollina@esalq.usp.br}}
\seealso{
  \code{\link[onemap]{make.seq}} and \code{\link[onemap]{compare}}.
}
\examples{
\dontrun{
  data(example_out)
  twopt <- rf.2pts(example_out)

  markers <- make.seq(twopt,c(2,3,12,14))
  markers.comp <- compare(markers)
  base.map <- make.seq(markers.comp,1)

  extend.map <- try.seq(base.map,30)
  extend.map
  print(extend.map,5) # best position
  print(extend.map,4) # second best position
}
}
\keyword{utilities}
