\name{rf.graph.table}
\alias{rf.graph.table}
\title{Plots a pairwise recombination fraction matrix using a color scale.}
\description{
  Plots a pairwise recombination fraction matrix using a color scale.
  Values of interest can be accessed using an interactive Tcl-Tk interface.
}
\usage{
rf.graph.table(w, scale=1, axis.cex=1, inter=TRUE)
}
\arguments{
  \item{w}{an object of class \code{sequence} with a predefined order.}
  \item{scale}{controls the plot size. If \code{inter == FALSE} this value is not used.}
  \item{axis.cex}{the magnification to be used for axis annotation.}
  \item{inter}{logical. If \code{TRUE}, an interactive graphic is plotted. Otherwise, a default graphic device is used.}
}
\details{
  The color scale varies from red (small distances) to
  dark blue. When clicking on a cell, a dialog box is displayed with some
  information about corresponding markers for that cell (line
  \eqn{\times} column). The informations are: \eqn{i}) the name of the markers;
  \eqn{ii}) the number of the markers; \eqn{iii}) the segregation
  types; \eqn{iv}) the recombination fraction between the  markers and
  \eqn{v}) the LOD-Score for each possible linkage phase calculated via
  two-point analysis. For neighbor
  markers, the multipoint recombination fraction is printed; otherwise,
  the two-point recombination fraction is printed. For markers of type \code{D1}
  and \code{D2}, it's impossible to calculate recombination
  fraction via two-point, therefore the corresponding cell will be
  empty. For diagonal cells, the name, the number and the type of the marker
  are printed, as well as the percentage of missing data for that marker.
}

\author{Marcelo Mollinari, \email{mmollina@esalq.usp.br}}
\examples{
  data(example_out)
  twopt <- rf.2pts(example_out)
  all.mark <- make.seq(twopt,"all")
  groups <- group(all.mark)
  LG1 <- make.seq(groups,1)
  LG1.rcd <- rcd(LG1)
  rf.graph.table(LG1.rcd, inter=FALSE)
\dontrun{
  ##Now, using interactive Tcl-Tk
  rf.graph.table(LG1.rcd, scale=1.5, inter=TRUE)
  }
}
\keyword{utilities}
