\name{onemap-package}
\alias{onemap-package}
\alias{onemap}
\docType{package}
\title{
  Software for constructing genetic maps in outcrossing species
}
\description{
  Implements the methodology proposed by \cite{Wu et al.} (2002), which
  uses the maximum likelihood approach to estimate linkage and linkage
  phases for a mixed set of different marker types containing various
  segregation patterns. The software was developed to construct genetic
  maps using full-sib families derived from two outbred parent plants,
  such as done for example by \cite{Garcia et al.} (2006) and
  \cite{Oliveira et al.} (2007).
}
\details{
  \tabular{ll}{
    Package: \tab onemap\cr
    Type: \tab Package\cr
    Version: \tab 0.1-0\cr
    Date: \tab 2007-09-25\cr
    License: \tab GNU GENERAL PUBLIC LICENSE (Version 2 or higher)\cr
  }
  
  Usage of \pkg{onemap} is completely described in a tutorial
  distributed with the package. You can look for it in directory
  \code{/inst/doc} of the package distribution.

  The most important functions are:
  \itemize{
    \item \code{\link[onemap]{read.outcross}} for data importing.

    \item \code{\link[onemap]{est.rf.2pts}} to perform two-point
    analyses between all pairs of markers.

    \item \code{\link[onemap]{modify.rf.2pts}} to update two-point
    analyses with different criteria.

    \item \code{\link[onemap]{arbitr.rf.2pts}} to arbitrarily (\emph{ad
      hoc}) define the most probable linkage phase between two markers.
    
    \item \code{\link[onemap]{mrktype}} to check the segregation type of
    a marker.
     
    \item \code{\link[onemap]{group}} to assign markers to linkage
    groups.

    \item \code{\link[onemap]{extract.group}} to \sQuote{extract} a
    linkage group for further analysis.

    \item \code{\link[onemap]{remove.marker}} to drop a marker from a
    linkage group.

    \item \code{\link[onemap]{map}} to estimate the genetic map (marker
    order and distances) for a given linkage group.

    \item \code{\link[onemap]{make.map}} to arbitrarily define a genetic
    map.

    \item \code{\link[onemap]{est.rf.3pts}} to refine a genetic map
    using three-point analyses for a triplet of markers.
  }
}
\note{
  See the complete tutorial distributed along with this package for
  complete examples (located in directory /inst/doc of the package
  distribution).
}
\author{Gabriel R A Margarido, \email{gramarga@esalq.usp.br}}
\references{
  Basten, C. J., Weir, B. S. and Zeng, Z.-B. (2005) \emph{QTL
    Cartographer Version 1.17: A Reference Manual and Tutorial for QTL
    Mapping}.
  
  Doerge, R. W. (1996) Constructing genetic maps by rapid chain
  delineation. \emph{Journal of Quantitative Trait Loci} 2: 121-132.
  
  Garcia, A. A. F., Kido, E. A., Meza, A. N., Souza, H. M. B., Pinto,
  L. R., Pastina, M. M., Leite, C. S., Silva, J. A. G., Ulian, E. C.,
  Figueira, A. V. O. and Souza, A. P. (2006) Development of an integrated
  genetic map of a sugarcane (\emph{Saccharum} spp.) commercial cross,
  based on a maximum-likelihood approach for estimation of linkage and
  linkage phases. \emph{Theoretical and Applied Genetics} 112: 298-314.
  
  Haldane, J. B. S. (1919) The combination of linkage values and the
  calculation of distance between the loci of linked
  factors. \emph{Journal of Genetics} 8: 299-309.
  
  Kosambi, D. D. (1944) The estimation of map distance from
  recombination values. \emph{Annuaire of Eugenetics} 12: 172-175.
  
  Lincoln, S. E., Daly, M. J. and Lander, E. S. (1993) Constructing
  genetic linkage maps with MAPMAKER/EXP Version 3.0: a tutorial and
  reference manual. \emph{A Whitehead Institute for Biomedical Research
    Technical Report}.
  
  Margarido, G. R. A., Souza, A.P. and Garcia, A. A. F. (2007) OneMap:
  software for genetic mapping in outcrossing species. \emph{Hereditas}
  144: 78-79.
  
  Oliveira, K. M., Pinto, L. R., Marconi, T. G., Margarido, G. R. A.,
  Pastina, M. M., Teixeira, L. H. M., Figueira, A. V. O., Ulian, E. C.,
  Garcia, A. A. F. and Souza, A. P. (2007) Functional genetic linkage map
  based on EST-markers for a sugarcane (\emph{Saccharum} spp.)
  commercial cross. \emph{Molecular Breeding} 20: 189-208.
  
  Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous
  maximum likelihood estimation of linkage and linkage phases in
  outcrossing species. \emph{Theoretical Population Biology} 61:
  349-363.
}
\keyword{package}

