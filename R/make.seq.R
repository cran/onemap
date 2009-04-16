#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: make.seq.R                                                    #
# Contains: make.seq, print.sequence                                  #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2009, Gabriel R A Margarido                           #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# Function to create sequences based on other object types
make.seq <- 
function(w, arg=NULL, phase=NULL, twopt=NULL) {
  # checking for correct object
  if(all(is.na(match(class(w),c("rf.2pts","group","compare","try","order")))))
    stop(deparse(substitute(w))," is not an object of classes 'rf.2pts', 'group', 'compare', 'try' or 'order'")

  switch(EXPR=class(w),
         'rf.2pts' = {
           if (length(arg) == 1 && arg == "all") seq.num <- 1:w$n.mar # generally used for grouping markers
		   else if(is.vector(arg) && is.numeric(arg)) seq.num <- arg
		   else stop("for an object of class 'rf.2pts', \"arg\" must be a vector of integers or the string 'all'")
		   ### CHECK IF MARKERS REALLY EXIST
           if (is.null(phase)) seq.phases <- -1 # no predefined linkage phases
           else if(length(phase) == (length(seq.num)-1)) seq.phases <- phase
           else stop("the length of 'phase' must be equal to the length of the sequence minus 1")
           seq.rf <- -1
           seq.like <- NULL
           if(is.null(twopt)) twopt <- deparse(substitute(w))
         },
         'group' = {
           if(length(arg) == 1 && is.numeric(arg) && arg <= w$n.groups) seq.num <- which(w$groups == arg)
		   else stop("for this object of class 'group', \"arg\" must be an integer less than or equal to ",w$n.groups)
           seq.phases <- -1
           seq.rf <- -1
           seq.like <- NULL
           twopt <- w$twopt
         },
		 'compare' = {
           n.ord <- max(which(head(w$best.ord.LOD,-1) != -Inf))
           unique.orders <- unique(w$best.ord[1:n.ord,])
		   if(is.null(arg)) seq.num <- unique.orders[1,] # NULL = 1 is the best order
           else if(length(arg) == 1 && is.numeric(arg) && arg <= nrow(unique.orders)) seq.num <- unique.orders[arg,]
		   else stop("for this object of class 'compare', \"arg\" must be an integer less than or equal to ",nrow(unique.orders))
           if (is.null(phase)) phase <- 1 # NULL = 1 is the best combination of phases
           chosen <- which(apply(w$best.ord[1:n.ord,],1,function(x) all(x==seq.num)))[phase]
           seq.phases <- w$best.ord.phase[chosen,]
           seq.rf <- w$best.ord.rf[chosen,]
           seq.like <- w$best.ord.like[chosen]
           twopt <- w$twopt
         },
		 'try' = {
		   if(length(arg) != 1 || !is.numeric(arg) || arg > length(w$ord))
		     stop("for this object of class 'try', \"arg\" must be an integer less than or equal to ",length(w$ord))
           if (is.null(phase)) phase <- 1 # NULL = 1 is the best combination of phases
           seq.num <- w$try.ord[arg,]
           seq.phases <- w$ord[[arg]]$phase[phase,]
           seq.rf <- w$ord[[arg]]$rf[phase,]
           seq.like <- w$ord[[arg]]$like[phase]
           twopt <- w$twopt
         },
		 'order' = {
           arg <- match.arg(arg,c("safe","force"))
           if (arg == "safe") {
		     # order with safely mapped markers
             seq.num <- w$ord$seq.num
             seq.phases <- w$ord$seq.phases
             seq.rf <- w$ord$seq.rf
             seq.like <- w$ord$seq.like
           }
           else {
		     # order with all markers
             seq.num <- w$ord.all$seq.num
             seq.phases <- w$ord.all$seq.phases
             seq.rf <- w$ord.all$seq.rf
             seq.like <- w$ord.all$seq.like
           }
           twopt <- w$twopt
         }
		 )

  # check if any marker appears more than once in the sequence
  if(length(seq.num) != length(unique(seq.num))) stop("there are duplicated markers in the sequence")
  
  structure(list(seq.num=seq.num, seq.phases=seq.phases, seq.rf=seq.rf, seq.like=seq.like,
                 data.name=w$data.name, twopt=twopt), class = "sequence")
}



# print method for object class 'sequence'
print.sequence <- function(x,...) {
  marnames <- colnames(get(x$data.name)$geno)[x$seq.num]
  if(is.null(x$seq.like)) {
    # no information available for the order
    cat("\nMarkers in the sequence:\n")
    cat(marnames,fill=TRUE)
    cat("\nParameters not estimated.\n\n")
  }
  else {
    # convert numerical linkage phases to strings
    link.phases <- matrix(NA,length(x$seq.num),2)
    link.phases[1,] <- rep(1,2)
    for (i in 1:length(x$seq.phases)) {
      switch(EXPR=x$seq.phases[i],
             link.phases[i+1,] <- link.phases[i,]*c(1,1),
             link.phases[i+1,] <- link.phases[i,]*c(1,-1),
             link.phases[i+1,] <- link.phases[i,]*c(-1,1),
             link.phases[i+1,] <- link.phases[i,]*c(-1,-1),
             )
    }
	# create diplotypes from segregation types and linkage phases
    link.phases <- apply(link.phases,1,function(x) paste(as.character(x),collapse="."))
    parents <- matrix("",length(x$seq.num),4)
    for (i in 1:length(x$seq.num))
      parents[i,] <- return.geno(get(x$data.name)$segr.type[x$seq.num[i]],link.phases[i])

    # display results
    longest.name <- max(nchar(marnames))
    marnames <- formatC(marnames,flag="-")
    distances <- formatC(c(0,cumsum(get(.map.fun)(x$seq.rf))),format="f",digits=2,width=7)
    cat("\nPrinting map:\n\n")
    cat("Markers",rep("",max(longest.name-7,0)+9),"Position",rep("",10),"Parent 1","     ","Parent 2\n\n")
    for (i in 1:length(x$seq.num)) {
      cat(marnames[i],rep("",max(7-longest.name,0)+10),distances[i],rep("",10),parents[i,1],"|  |",parents[i,2],"     ",parents[i,3],"|  |",parents[i,4],"\n")
    }
    cat("\nlog-likelihood:",x$seq.like,"\n\n")
  }
}
