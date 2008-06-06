\name{mrktype}
\alias{mrktype}
\title{Type of segregation of a marker}
\description{
  Informs the type of segregation of a marker from an object of class
  \code{rf.2pts}, according to the notation by \cite{Wu et al., 2002}.
}
\usage{
mrktype(w, mrkname = NULL)
}
\arguments{
  \item{w}{an object of class \code{rf.2pts}.}
  \item{mrkname}{the name of a marker, corresponding to any marker in the
    input file used for the two-point analyses stored in object \code{w}.}
}
\details{
  The segregation types are (\cite{Wu et al., 2002}):
  \tabular{lcc}{
    Type \tab Cross \tab Segregation \cr
    A.1 \tab ab x cd \tab 1:1:1:1 \cr
    A.2 \tab ab x ac \tab 1:1:1:1 \cr
    A.3 \tab ab x co \tab 1:1:1:1 \cr
    A.4 \tab ao x bo \tab 1:1:1:1 \cr
    B1.5 \tab ab x ao \tab 1:2:1 \cr
    B2.6 \tab ao x ab \tab 1:2:1 \cr
    B3.7 \tab ab x ab \tab 1:2:1 \cr
    C8 \tab ao x ao \tab 3:1 \cr
    D1.9 \tab ab x cc \tab 1:1 \cr
    D1.10 \tab ab x aa \tab 1:1 \cr
    D1.11 \tab ab x oo \tab 1:1 \cr
    D1.12 \tab bo x aa \tab 1:1 \cr
    D1.13 \tab ao x oo \tab 1:1 \cr
    D2.14 \tab cc x ab \tab 1:1 \cr
    D2.15 \tab aa x ab \tab 1:1 \cr
    D2.16 \tab oo x ab \tab 1:1 \cr
    D2.17 \tab aa x bo \tab 1:1 \cr
    D2.18 \tab oo x ao \tab 1:1
  }
}
\value{
  A character string giving the type of segregation of marker
  \code{mrkname}, according the above table.
}
\references{
  Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous
  maximum likelihood estimation of linkage and linkage phases in
  outcrossing species. \emph{Theoretical Population Biology} 61: 349-363.
}
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br}}
\seealso{\code{\link[onemap]{est.rf.2pts}}}
\examples{
  data(example_out)
  twopts <- est.rf.2pts(example_out)

  mrktype(twopts,"M1") # segregation type for marker 1
}
\keyword{manip}
\keyword{utilities}
