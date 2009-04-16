#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: order.seq.R                                                   #
# Contains: order.seq, print.order                                    #
#                                                                     #
# Written by Gabriel R A Margarido & Marcelo Mollinari                #
# copyright (c) 2009, Gabriel R A Margarido & Marcelo Mollinari       #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function automates linkage map construction in two steps:
# first, it applies the 'compare' algorithm to a subset of markers;
# second, it adds markers sequentially with the 'try' function
order.seq <- function(w,n.init=5,THRES=3,touchdown=FALSE,tol=10E-5) {
  # checking for correct objects
  if(!any(class(w)=="sequence")) stop(deparse(substitute(w))," is not an object of class 'sequence'")
  if(n.init < 2) stop("n.init must be greater than or equal to 2")
  if(THRES <= 10E-10) stop("Threshold must be greater than 0")
  
  if(length(w$seq.num) <= n.init) {
    # in this case, only the 'compare' function is used
    cat("   Length of sequence w is less than n.init \n   Returning the best order using compare function:\n")
    ifelse(length(w$seq.num) == 2, seq.ord <- map(w,tol=tol), seq.ord <- make.seq(compare(w=w,tol=tol),1))
    structure(list(ord=seq.ord, mrk.unpos=NULL, LOD.unpos=NULL, THRES=THRES,
			       ord.all=seq.ord, data.name=w$data.name, twopt=w$twopt), class = "order")
  }
  else {
    # here, the complete algorithm will be applied
	
	# select the order in which markers will be added
	segregation.types <- get(w$data.name)$segr.type.num[w$seq.num]
	if(sum(segregation.types == 7) > sum(segregation.types == 6)) segregation.types[segregation.types == 6] <- 8 # if there are more markers of type D2 than D1, try to map those first
    seq.work <- order(segregation.types)
	
	# apply the 'compare' step to the subset of initial markers
    seq.init <- w$seq.num[seq.work[1:n.init]]
    seq.ord <- compare(w=make.seq(get(w$twopt),seq.init,twopt=w$twopt),n.best=50,tol=tol)
	
	# 'try' to map remaining markers
    w2 <- make.seq(seq.ord,1)
	cat ("\n\nRunning try algorithm\n")
    for (i in (n.init+1):length(w$seq.num)){
      seq.ord <- try.seq(w2,w$seq.num[seq.work[i]],tol=tol)
      if(all(seq.ord$LOD[-which(seq.ord$LOD==max(seq.ord$LOD))[1]] < -THRES))
		w2 <- make.seq(seq.ord,which.max(seq.ord$LOD))
    }
	
	# markers that do not meet the threshold remain unpositioned
    mrk.unpos <- w$seq.num[which(is.na(match(w$seq.num, w2$seq.num)))]
	LOD.unpos <- NULL
	cat("\nLOD threshold =",THRES,"\n\nPositioned markers:", w2$seq.num, "\n\n")
    cat("Markers not placed on the map:", mrk.unpos, "\n")
	
	if(touchdown && length(mrk.unpos) > 0) {
	  # here, a second round of the 'try' algorithm is performed, if requested
      cat("\n\n\nTrying to map remaining markers with LOD threshold ",THRES-1,"\n")
      for (i in mrk.unpos) {
        seq.ord <- try.seq(w2,i,tol=tol)
        if(all(seq.ord$LOD[-which(seq.ord$LOD==max(seq.ord$LOD))[1]] < -THRES+1))
          w2 <- make.seq(seq.ord,which.max(seq.ord$LOD))
      }
	  
	  # markers that do not meet this second threshold still remain unpositioned 
      mrk.unpos <- w$seq.num[which(is.na(match(w$seq.num, w2$seq.num)))]
      cat("\nLOD threshold =",THRES-1,"\n\nPositioned markers:", w2$seq.num, "\n\n")
      cat("Markers not placed on the map:", mrk.unpos, "\n")
    }
	
    if(length(mrk.unpos) > 0) {
	  # LOD-Scores are calculated for each position, for each unmapped marker, if any
      LOD.unpos <- matrix(NA,length(mrk.unpos),(length(w2$seq.num)+1))
      j <- 1
      cat("\n\nCalculating LOD-Scores\n")
      for (i in mrk.unpos){
        LOD.unpos[j,] <- try.seq(w=w2,mrk=i,tol=tol)$LOD
        j <- j+1
      }
    }
    else mrk.unpos <- NULL

    # to end the algorithm, possibly remaining markers are 'forced' into the map
    w3 <- w2
    if(!is.null(mrk.unpos)) {
	  cat("\n\nPlacing remaining marker(s) at most likely position\n")
	  
	  # these markers are added from the least to the most doubtful
      which.order <- order(apply(LOD.unpos,1,function(x) max(x[-which(x==0)[1]])))
	  
      for (i in mrk.unpos[which.order]) {
        seq.ord <- try.seq(w3,i,tol)
        w3 <- make.seq(seq.ord,which(seq.ord$LOD==0)[sample(sum(seq.ord$LOD==0))[1]])
      }
    }
    structure(list(ord=w2, mrk.unpos=mrk.unpos, LOD.unpos=LOD.unpos, THRES=THRES,
                   ord.all=w3, data.name=w$data.name, twopt=w$twopt), class = "order")
  }
}

print.order <- function(x,...) {
  cat("\nBest sequence found.")
  # print the 'safe' order
  print(x$ord)
  if(!is.null(x$mrk.unpos)) {
    # print LOD-Score information for unpositioned markers
    cat("\n\nThe following markers could not be uniquely placed!\n")
    cat("Printing most likely positions for each unplaced marker:\n")
    
    size1 <- max(3,max(nchar(x$mrk.unpos)))
    mrk.unpos.pr <- format(x$mrk.unpos,width=size1)
    size2 <- max(nchar(x$ord$seq.num))
    seq.pr <- format(x$ord$seq.num,width=size2)
    
    limit <- (x$THRES-2)/2
    
    cat("\n")
    cat(paste(rep("-",size2+4+length(mrk.unpos.pr)*(size1+3)),collapse=""),"\n")
    cat("| ",rep("",size2),"|")
    #### MAYBE WE SHOULD PUT A LIMIT TO THE NUMBER OF UNPOSITIONED MARKERS
    for(j in 1:length(mrk.unpos.pr)) {
      cat(rep("",max(0,3-size1)+1),mrk.unpos.pr[j],"|")
    }
    cat("\n")
    cat(paste("|",paste(rep("-",size2+2),collapse=""),"|",sep=""))
    cat(paste(rep(paste(paste(rep("-",size1+2),collapse=""),"|",sep=""),length(mrk.unpos.pr)),collapse=""),"\n")
    cat("| ",rep("",size2),"|")
    for(j in 1:length(x$mrk.unpos)) {
      if(x$LOD.unpos[j,1] > -limit) cat(rep("",max(0,3-size1)+1),"*** |")
      else if(x$LOD.unpos[j,1] > -2*limit) cat(rep("",max(0,3-size1)+1)," *  |")
      else cat(rep("",max(0,3-size1)+1),"    |")
    }
    cat("\n")
    for(i in 1:length(seq.pr)) {
      cat("|",seq.pr[i],"|")
      cat(paste(rep(paste(paste(rep(" ",size1+2),collapse=""),"|",sep=""),length(mrk.unpos.pr)),collapse=""),"\n")
      cat(paste("|",paste(rep(" ",size2+2),collapse=""),"|",sep=""))
      for(j in 1:length(x$mrk.unpos)) {
        if(x$LOD.unpos[j,i+1] > -limit) cat(rep("",max(0,3-size1)+1),"*** |")
        else if(x$LOD.unpos[j,i+1] > -2*limit) cat(rep("",max(0,3-size1)+1)," *  |")
        else cat(rep("",max(0,3-size1)+1),"    |")
      }
      cat("\n")
    }
    cat(paste(rep("-",size2+4+length(mrk.unpos.pr)*(size1+3)),collapse=""),"\n")
    cat("\n")
    cat("'***' indicates the most likely position(s)\n\n")
    cat("'*' indicates less likely positions\n\n")
  }
}

# end of file
