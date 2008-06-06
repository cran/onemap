\name{remove.marker}
\alias{remove.marker}
\title{Drop a marker from a linkage group}
\description{
  Removes a genetic marker from a linkage group, i.e., from an object of
  class \code{extracted.group}.
}
\usage{
remove.marker(w, mname = NULL)
}
\arguments{
  \item{w}{an object of class \code{extracted.group}.}
  \item{mname}{a character string indicating the name of the marker to be
    removed.}
}
\value{
  Returns an object of class \code{extracted.group}, which is a list
  containing the following components:
  \item{marnames}{names of markers present in linkage group \code{i},
    according to the original input file.}
  \item{number}{number of the linkage group, according to the numeration
    in the object of class \code{group} from where the linkage group was
    extracted.}
  \item{oriname}{name of the original object of class \code{rf.2pts}
    used  to generate the object of class \code{group}, i.e., containing
    the results of two-point analyses.}
  \item{name}{name of the object of class \code{group} from where the
    linkage group was extracted.}
}
\note{
  Intended for dealing with \sQuote{problematic} markers, that is,
  markers whose position on the map is hard to estimate, or when linkage
  is uncertain.
}
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br}}
\seealso{
  \code{\link[onemap]{group}}, \code{\link[onemap]{print.group}} and
  \code{\link[onemap]{extract.group}}}
\examples{
  data(example_out)
  twopts <- est.rf.2pts(example_out)
   
  link_gr <- group(twopts)
  LG1 <- extract.group(link_gr, 1) # extracts the first linkage group
  LG1

  LG1 <- remove.marker(LG1,"M1") # linkage group 1 without first marker
  LG1
}
\keyword{manip}
