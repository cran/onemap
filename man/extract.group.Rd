\name{extract.group}
\alias{extract.group}
\title{Extracts a linkage group for further analysis}
\description{
  Takes as input an object of class \code{group} and \dQuote{extracts} a
  linkage group that is going to be ordered.
}
\usage{
extract.group(w, i = NULL)
}
\arguments{
  \item{w}{an object of class \code{group}.}
  \item{i}{number of the linkage group to be extracted, according to the
    numeration found in object \code{w}.}
}
\value{
  Returns an object of class \code{extracted.group}, which is a list
  containing the following components:
  \item{marnames}{names of markers present in linkage group \code{i},
    according to the original input file.}
  \item{number}{number of the linkage group extracted, according to the
    numeration in object \code{w}.}
  \item{oriname}{name of the original object of class \code{rf.2pts}
    used  to generate the object of class \code{group}, i.e., containing
    the results of two-point analyses.}
  \item{name}{name of the object of class \code{group} used as input,
    i.e., containing the information on linkage groups.}
}
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br}}
\seealso{
  \code{\link[onemap]{group}}, \code{\link[onemap]{print.group}} and
  \code{\link[onemap]{remove.marker}}
}
\examples{
  data(example_out)
  twopts <- est.rf.2pts(example_out)

  link_gr <- group(twopts)
  LG1 <- extract.group(link_gr, 1) # extracts the first linkage group
  LG1
}
\keyword{manip}
