#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: write_map.R                                                   #
# Contains: write_map                                                 #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2010, Marcelo Mollinari                               #
#                                                                     #
# First version: 10/10/2010                                           #
# Last update: 11/30/2010                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

##' Write a genetic map to a file
##'
##' Write a genetic map to a file, base on a given map, or a list of maps. The
##' output file can be used as an input to perform QTL mapping using the package
##' R/qtl. It is also possible to create an output to be used with
##' QTLCartographer program.
##'
##' This function is avaliable only for backcross, F2 and RILs.
##'
##' @param map.list a map, i.e. an object of class \code{sequence} with a
##' predefined order, linkage phases, recombination fraction and likelihood or
##' a list of maps.
##' @param file.out output map file.
##' @author Marcelo Mollinari, \email{mmollina@@usp.br}
##' @references Broman, K. W., Wu, H., Churchill, G., Sen, S., Yandell, B.
##' (2008) \emph{qtl: Tools for analyzing QTL experiments} R package version
##' 1.09-43
##'
##' Wang S., Basten, C. J. and Zeng Z.-B. (2010) Windows QTL Cartographer 2.5.
##' Department of Statistics, North Carolina State University, Raleigh, NC.
##' @keywords rqtl
##' @examples
##'
##' \dontrun{
##' data(mapmaker_example_f2)
##' twopt<-rf_2pts(mapmaker_example_f2)
##' lg<-group(make_seq(twopt, "all"))
##'
##' ##"pre-allocate" an empty list of length lg$n.groups (3, in this case)
##'   maps.list<-vector("list", lg$n.groups)
##'
##'   for(i in 1:lg$n.groups){
##'     ##create linkage group i
##'     LG.cur <- make_seq(lg,i)
##'     ##ordering
##'     map.cur<-order_seq(LG.cur, subset.search = "sample")
##'     ##assign the map of the i-th group to the maps.list
##'     maps.list[[i]]<-make_seq(map.cur, "force")
##'   }
##'
##' ##write maps.list to "mapmaker_example_f2.map" file
##' write_map(map.list, "mapmaker_example_f2.map")
##'
##' ##Using R/qtl
##' ##you must install the package  'qtl'
##' ##install.packages("qtl")
##'
##' require(qtl)
##' file<-paste(system.file("example",package="onemap"),"mapmaker_example_f2.raw", sep="/")
##' dat1 <- read.cross("mm", file=file, mapfile="mapmaker_example_f2.map")
##' newmap <- est.map(dat1, tol=1e-6, map.function="kosambi")
##'
##' (logliks <- sapply(newmap, attr, "loglik"))
##' plot.map(dat1, newmap)
##'
##' ##Using R/qtl to generate QTL Cartographer input files (.map and .cro)
##' write.cross(dat1, format="qtlcart", filestem="mapmaker_example_f2")
##'
##' }
##'@export
write_map<-function(map.list,file.out){
   # checking for correct object
  if(!(is(map.list,"list") | is(map.list,"sequence"))) stop(deparse(substitute(map.list))," is not an object of class 'list' or 'sequnece'")
  if(!is(file.out,"character")) stop(deparse(substitute(file.out))," is an invalid output file name")

  # if map.list is just a single chormosome, convert it  into a list
  if(is(map.list,"sequence")) map.list<-list(map.list)

  write(x="",file=file.out)
  for(i in 1:length(map.list)){
    if(!is(map.list[[i]],"sequence")) stop("Object ", i , " in map.list is not an object of class 'sequnece'")
    if(is.null(map.list[[i]]$seq.like))  stop("Parameters are not estimated for object ", i, " in map.list")
    map<-cumsum(c(0,get(get(".map.fun", envir=.onemapEnv))(map.list[[i]]$seq.rf)))
    marnames<-colnames(get(map.list[[i]]$data.name, pos=1)$geno)[map.list[[i]]$seq.num]
    write(t(cbind(i,marnames,map)), file=file.out, ncolumns =3, append=TRUE)
  }
}

