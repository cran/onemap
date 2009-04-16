\name{order.seq}
\alias{order.seq}
\title{Search for the best order of markers combining compare and try.seq
       functions}
\description{
  For a given sequence of markers, this function first uses the \code{compare}
  function to create a framework for a subset of informative  markers. Then, it tries to map
  remaining ones using the \code{try.seq} function.
}
\usage{
order.seq(w,n.init=5,THRES=3,touchdown=FALSE,tol=10E-5)
}
\arguments{
  \item{w}{an object of class \code{sequence}.}
  \item{n.init}{the number of markers to be used in the \code{compare} step
    (defaults to 5).}
  \item{THRES}{threshold to be used when positioning markers in the \code{try.seq}
    step.}
  \item{touchdown}{logical. If \code{FALSE} (default), the \code{try.seq} step is
    run only once, with the value of \code{THRES}. If \code{TRUE}, \code{try.seq}
    runs with \code{THRES} and then once more, with \code{THRES-1}. The latter
    calculations take longer, but usually are able to map more markers.}
  \item{tol}{tolerance number for the C routine, i.e., the value used to evaluate
    convergence.}
}
\details{
  After running the \code{compare} and \code{try.seq} steps, which result in a
  "safe" order, markers that could not be mapped are "forced" into the map,
  resulting in a map with all markers positioned.
}
\value{
  An object of class \code{order}, which is a list containing the following
  components:
  \item{ord}{an object of class \code{sequence} containing the "safe" order.}
  \item{mrk.unpos}{a \code{vector} with (possibly) unpositioned markers.}
  \item{LOD.unpos}{a \code{matrix} with LOD-Scores for unmapped markers, if any,
    for each position in the "safe" order.}
  \item{THRES}{the same as the input value, just for printing.}
  \item{ord.all}{an object of class \code{sequence} containing the "forced"
    order, i.e., the best order with all markers.}
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
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br} and
        Marcelo Mollinari, \email{mmollina@esalq.usp.br}}
\seealso{
  \code{\link[onemap]{make.seq}}, \code{\link[onemap]{compare}} and
  \code{\link[onemap]{try.seq}}.
}
\examples{
\dontrun{
  data(example_out)
  twopt <- rf.2pts(example_out)
  all.mark <- make.seq(twopt,"all")
  groups <- group(all.mark)
  LG2 <- make.seq(groups,2)
  LG2.ord <- order.seq(LG2,touchdown=TRUE)
  LG2.ord
  make.seq(LG2.ord) # get safe sequence
  make.seq(LG2.ord,"force") # get forced sequence
}
}
\keyword{utilities}
