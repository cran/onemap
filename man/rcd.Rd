\name{rcd}
\alias{rcd}
\title{Rapid Chain Delineation}
\description{
  Implements the marker ordering algorithm \emph{Rapid Chain
    Delineation} (\cite{Doerge, 1996}).
}
\usage{
rcd(r)
}
\arguments{
  \item{r}{A (symmetric) matrix containing the results of
    two-point analyses, i.e., pairwise recombination fractions between
    genetic markers.}
}
\details{
  \emph{Rapid Chain Delineation} (\emph{RCD}) is an algorithm for marker
  ordering in linkage groups. It is not an exhaustive search method and,
  therefore, is not computationally intensive. However, it does not
  guarantee that the best order is always found. The only requirement is
  a matrix with recombination fractions between markers.
  Next is an excerpt from QTL Cartographer Version 1.17 Manual
  describring the \emph{RCD} algorithm (\cite{Basten et al., 2005}):
  
  \emph{The linkage group is initiated with the pair of markers having the
  smallest recombination fraction. The remaining markers are placed in a
  \dQuote{pool} awaiting placement on the map. The linkage group is
  extended by adding markers from the pool of unlinked markers. Each
  terminal marker of the linkage group is a candidate for extension of
  the chain: The unlinked marker that has the smallest recombination
  fraction with either is added to the chain subject to the provision
  that the recombination fraction is statistically significant at a
  prespecified level. This process is repeated as long as markers can be
  added to the chain.}
}
\value{
  A vector with the estimated order of markers on the linkage
  group. Numbers are based on the position of each marker on the
  recombination fraction matrix given as input, i.e., on the
  \code{row}/\code{column} each marker occupies on that
  matrix. No names are displayed (for portability purposes).
}
\references{
  Basten, C. J., Weir, B. S. and Zeng, Z.-B. (2005) \emph{QTL
    Cartographer Version 1.17: A Reference Manual and Tutorial for QTL
    Mapping}.
  
  Doerge, R. W. (1996) Constructing genetic maps by rapid chain
  delineation. \emph{Journal of Quantitative Trait Loci} 2: 121-132.
}
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br}}
\seealso{
  \code{\link[onemap]{map}}, \code{\link[onemap]{print.map}}
}
\examples{
  rcmb <- c(NA,0.10,0.05,0.10,NA,0.06,0.05,0.06,NA)
  rcmb_mat <- matrix(rcmb,nrow=3,ncol=3,byrow=TRUE) # hypothetical recombination fraction matrix
  colnames(rcmb_mat) <- rownames(rcmb_mat) <- paste("M",1:3,sep="")
  rcmb_mat # not monotonic

  order <- rcd(rcmb_mat)
  order
  
  rcmb_mat[order,order] # monotonic
}
\keyword{utilities}
