\name{read.outcross}
\alias{read.outcross}
\title{Read data from a segregating full-sib population}
\description{
  Imports data from a full-sib family derived from two outbred parents
  and creates an object of class \code{outcross}. Currently only one
  format is supported.
}
\usage{
read.outcross(dir, file)
}
\arguments{
  \item{dir}{directory where the input file is located.}
  \item{file}{the name of the input file which contains the data to be
    read.}
}
\details{
  The only format currently supported is quite similar to that of
  \code{MAPMAKER/EXP} (\cite{Lincoln et al.}, 1993). The first line
  contains two integers: the number of individuals and the number of
  markers.

  Next comes the genotype data for all markers. Each new marker is
  initiated with a \dQuote{*} (without the quotes) followed by the
  marker name, without any space between them. Each marker name is
  followed by the corresponding segregation type, which may be:
  \code{"A.1"}, \code{"A.2"}, \code{"A.3"}, \code{"A.4"}, \code{"B1.5"},
  \code{"B2.6"}, \code{"B3.7"}, \code{"C.8"}, \code{"D1.9"},
  \code{"D1.10"}, \code{"D1.11"}, \code{"D1.12"}, \code{"D1.13"},
  \code{"D2.14"}, \code{"D2.15"}, \code{"D2.16"}, \code{"D2.17"} or
  \code{"D2.18"} (without quotes) [see \code{\link[onemap]{mrktype}} and
  \cite{Wu et al.} (2002) for details].

  Finally, after the segregation type comes the genotype data for the
  corresponding marker. Depending on the segregation type, genotypes may
  be denoted by \code{ac}, \code{ad}, \code{bc}, \code{bd}, \code{a},
  \code{ba}, \code{b}, \code{bc}, \code{ab} and \code{o}, in various
  combinations. To make things easier, we have followed \strong{exactly}
  the notation used by \cite{Wu et al.} (2002). Genotypes \emph{must} be
  separated by commas. Missing values are denoted by \code{"-"}

  The \code{example} directory in the package distribution contains
  an example data file to be read with this function. Further
  instructions can be found at the tutorial distributed along with this
  package.
}
\value{
  An object of class \code{outcross}, i.e., a list with the following
  components:
  \item{geno}{a matrix with character strings indicating the genotypes
    read for each marker. Each column contains data for a marker and
    each row represents an individual.}
  \item{n.ind}{number of individuals.}
  \item{n.mar}{number of markers.}
  \item{segr.type}{a vector with the segregation type of each marker.}
}
\references{
  Lincoln, S. E., Daly, M. J. and Lander, E. S. (1993) Constructing
  genetic linkage maps with MAPMAKER/EXP Version 3.0: a tutorial and
  reference manual. \emph{A Whitehead Institute for Biomedical Research
    Technical Report}.

  Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous
  maximum likelihood estimation of linkage and linkage phases in
  outcrossing species. \emph{Theoretical Population Biology} 61:
  349-363.
}
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br}}
\seealso{
  \code{example} directory in the package source.
}
\examples{
  \dontrun{
    outcr_data <-
read.outcross(dir="work_directory",file="data_file.txt")
  }
}
\keyword{IO}
