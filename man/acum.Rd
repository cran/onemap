\name{acum}
\alias{acum}
\title{Cumulative sum}
\description{
  Calculates the cumulative sum from \code{1} to \code{w}.
}
\usage{
acum(w)
}
\arguments{
  \item{w}{a non-negative integer}
}
\value{
  Returns the cumulative sum \code{1+2+\ldots+w}, that is,
  \deqn{\frac{w \left( w+1 \right)}{2},}{w*(w+1)/2,}
  for \eqn{w \geq 0}.
}
\note{
  Internally used by \pkg{onemap} - you probably do not need to use it
  directly.
}
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br}}
\examples{
acum(0) # 0

acum(4) # 1+2+3+4
}
\keyword{arith}
