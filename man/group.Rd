\name{group}
\alias{group}
\alias{print.group}
\title{Assign markers to linkage groups}
\description{
  Identifies linkage groups of markers, using results from two-point
  (pairwise) analysis and the \emph{transitive} property of linkage.
}
\usage{
group(w)

\method{print}{group}(x, detailed=TRUE, \dots)
}
\arguments{
  \item{w}{an object of class \code{rf.2pts}.}
  \item{x}{an object of class \code{group}.}
  \item{detailed}{logical. If \code{FALSE}, only a small summary of the
    linkage groups is printed. If \code{TRUE} (default), the names of
    markers in each linkage group are also displayed.}
  \item{\dots}{further arguments, passed to other methods. Currently
    ignored.}
}
\details{
   Thresholds used to group markers, i.e., minimum LOD Score and maximum
   recombination fraction, are defined in object \code{w}.
}
\value{
  Returns an object of class \code{group}, which is a list containing
  the following components:
  \item{marnames}{marker names, according to the input file.}
  \item{n.mar}{total number of markers.}
  \item{LOD}{minimum LOD Score to declare linkage.}
  \item{max.rf}{maximum recombination fraction to declare linkage.}
  \item{n.groups}{number of groups found.}
  \item{groups}{number of the linkage group to which each marker is
    assigned.}
  \item{name}{name of the object of class \code{rf.2pts} used as input,
    i.e., containing the information used to assign markers to linkage
    groups.}
}
\references{
  Lincoln, S. E., Daly, M. J. and Lander, E. S. (1993) Constructing
  genetic linkage maps with MAPMAKER/EXP Version 3.0: a tutorial and
  reference manual. \emph{A Whitehead Institute for Biomedical Research
    Technical Report}.
  % NAO SEI SE ESSA E A MELHOR REFERENCIA
}
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br}}
\seealso{
  \code{\link[onemap]{est.rf.2pts}},
  \code{\link[onemap]{modify.rf.2pts}} and
  \code{\link[onemap]{arbitr.rf.2pts}}
}
\examples{
  data(example_out)
  twopts <- est.rf.2pts(example_out)

  link_gr <- group(twopts)
  link_gr
}
\keyword{misc}
