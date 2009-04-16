\name{onemap-internal}
\alias{acum}
\alias{cr2pts}
\alias{cr3pts}
\alias{codif.data}
\alias{comb.ger}
\alias{comb}
\alias{est.map.c}
\alias{draw.rf.inter} 
\alias{plotFunction.out}
\alias{perm.tot}
\alias{perm.pars}
\alias{phases}
\alias{return.geno}
\alias{print.compare}
\alias{print.order}
\alias{print.outcross}
\alias{print.rf.3pts}
\alias{print.sequence}
\title{Internal onemap functions}
\description{
  Functions from \pkg{onemap} not to be directly called by the user.
}
\usage{
acum(w)
cr2pts(mrk1, mrk2, segr.type1, segr.type2)
cr3pts(mrk1, mrk2, mrk3, segr.type1, segr.type2, segr.type3)
codif.data(geno.in,segr.type.in)
comb.ger(f)
comb(x,y)
est.map.c(geno,type,phase,rec,verbose,tol=1e-6)
perm.tot(v)
perm.pars(v)
draw.rf.inter(info.graph, scale, cex)
plotFunction.out(info.graph, cex)
phases(w, LOD=0, max.rf=0.50)
return.geno(segr.type, link.phases)
\method{print}{compare}(x, \dots)
\method{print}{order}(x, \dots)
\method{print}{outcross}(x, \dots)
\method{print}{rf.3pts}(x, \dots)
\method{print}{sequence}(x, \dots)
}
\author{Gabriel Rodrigues Alves Margarido}
\keyword{internal}
