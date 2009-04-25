#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: try.seq.R                                                     #
# Contains: try.seq, print.try                                        #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2009, Marcelo Mollinari                               #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function tries to position a marker in a given linkage map
## CHECK IF THE ord$rf MATRIX HAS SOME NaN VALUE AND ISSUE WARNING (CHANGE FOR 'BAD' VALUE)
try.seq <-
function(w,mrk,tol=10E-5,verbose=FALSE) {
  # checking for correct objects
  if(!any(class(w)=="sequence")) stop(deparse(substitute(w))," is not an object of class 'sequence'")
  if(w$seq.phases[1] == -1) stop("You must run 'comapre' or 'map' before the 'try.seq' function")
  if(mrk > get(w$data.name)$n.mar) stop(deparse(substitute(mrk	))," exceeds the number of markers in object ", w$data.name)

  # allocate variables
  rf.init <- vector("list",length(w$seq.num))
  phase.init <- vector("list",length(w$seq.num))
  ord <- list(list(matrix(NA,16,length(w$seq.num)),
                 matrix(NA,16,length(w$seq.num)),
                 rep(-Inf,16)))
  names(ord[[1]]) <- c("rf","phase","like")
  ord <- rep(ord,length(w$seq.num)+1)
  best.seq <- rep(-Inf,length(w$seq.num)+1)
      
  ### positioning before the given sequence
  # get two-point information
  list.init <- phases(make.seq(get(w$twopt),c(mrk,w$seq.num[1]),twopt=w$twopt))
  rf.init[[1]] <- list.init$rf.init[[1]]
  for(j in 1:(length(w$seq.num)-1)) rf.init[[j+1]] <- w$seq.rf[j]
  phase.init[[1]] <- list.init$phase.init[[1]]
  for(j in 1:(length(w$seq.num)-1)) phase.init[[j+1]] <- w$seq.phases[j]
  Ph.Init <- comb.ger(phase.init)
  Rf.Init <- comb.ger(rf.init)
  mark.max<-max(nchar(colnames(get(w$data.name)$geno)))
  num.max<-nchar(ncol(get(w$data.name)$geno))
  
  # create first order
  try.ord <- c(mrk,w$seq.num)
  if(verbose) cat("TRY", 1,": ", c(mrk,w$seq.num),"\n")
  else cat(format(mrk,width=num.max) , "-->", format(colnames(get(w$data.name)$geno)[mrk], width=mark.max), ": .")
  flush.console()
  
  if(nrow(Ph.Init)>1){
    ##Removing ambigous phases
    rm.ab<-rem.amb.ph(M=Ph.Init, w=w, seq.num=c(mrk,w$seq.num))
    Ph.Init <- Ph.Init[rm.ab,]
    Rf.Init <- Rf.Init[rm.ab,]
    if(class(Ph.Init)=="integer"){
      Ph.Init<-matrix(Ph.Init,nrow=1)
      Rf.Init<-matrix(Rf.Init,nrow=1)
    }
  }
  # estimate parameters for all possible linkage phases for this order
  for(j in 1:nrow(Ph.Init)) {
    final.map <- est.map.c(geno=get(w$data.name)$geno[,c(mrk,w$seq.num)],
                           type=get(w$data.name)$segr.type.num[c(mrk,w$seq.num)],
                           phase=Ph.Init[j,],
                           rec=Rf.Init[j,],
                           verbose=FALSE,
                           tol=tol)
    ord[[1]]$rf[j,] <- final.map$rf
    ord[[1]]$phase[j,] <- Ph.Init[j,]
    ord[[1]]$like[j] <- final.map$loglike
    best.seq[1] <- max(best.seq[1],final.map$loglike)
  }
  # sort linkage phases by log-likelihood
  ord.ind <- order(ord[[1]]$like, decreasing=TRUE)
  ord[[1]]$rf <- ord[[1]]$rf[ord.ind,]
  ord[[1]]$phase <- ord[[1]]$phase[ord.ind,]
  ord[[1]]$like <- ord[[1]]$like[ord.ind]
  
  ### positioning between markers of the given sequence
  for(i in 1:(length(w$seq.num)-1)) {
    # get two-point information
    list.init <- phases(make.seq(get(w$twopt),c(w$seq.num[i],mrk,w$seq.num[i+1]),twopt=w$twopt))
    if(i!=1) {
      for(k in 1:(i-1)) {
        rf.init[[k]] <- w$seq.rf[k]
        phase.init[[k]] <- w$seq.phases[k]
      }
    }
    rf.init[[i]] <- list.init$rf.init[[1]]
    phase.init[[i]] <- list.init$phase.init[[1]]
    rf.init[[i+1]] <- list.init$rf.init[[3]]
    phase.init[[i+1]] <- list.init$phase.init[[3]]
    if(i!=(length(w$seq.num)-1)) {
      for(k in (i+2):length(w$seq.num)) {
        rf.init[[k]] <- w$seq.rf[k-1]
        phase.init[[k]] <- w$seq.phases[k-1]
      }
    }   
    Ph.Init <- comb.ger(phase.init)
    Rf.Init <- comb.ger(rf.init)
 
	# create intermediate orders
    try.ord <- rbind(try.ord,c(w$seq.num[1:i], mrk, w$seq.num[(i+1):length(w$seq.num)]))
    if(verbose) cat("TRY", i+1,": ",c(w$seq.num[1:i], mrk, w$seq.num[(i+1):length(w$seq.num)]) ,"\n")
    else cat(".")
    flush.console()
    
    if(nrow(Ph.Init)>1){
      ##Removing ambigous phases
      rm.ab<-rem.amb.ph(M=Ph.Init, w=w, seq.num=c(w$seq.num[1:i], mrk, w$seq.num[(i+1):length(w$seq.num)]))
      Ph.Init <- Ph.Init[rm.ab,]
      Rf.Init <- Rf.Init[rm.ab,]
      if(class(Ph.Init)=="integer"){
        Ph.Init<-matrix(Ph.Init,nrow=1)
        Rf.Init<-matrix(Rf.Init,nrow=1)
      }
    }
    ## estimate parameters for all possible linkage phases for the current order
    for(j in 1:nrow(Ph.Init)) {
      final.map <- est.map.c(geno=get(w$data.name)$geno[,c(w$seq.num[1:i], mrk, w$seq.num[(i+1):length(w$seq.num)])],
							 type=get(w$data.name)$segr.type.num[c(w$seq.num[1:i], mrk, w$seq.num[(i+1):length(w$seq.num)])],
                             phase=Ph.Init[j,],
                             rec=Rf.Init[j,],
                             verbose=FALSE,
                             tol=tol)
      ord[[i+1]]$rf[j,] <- final.map$rf
      ord[[i+1]]$phase[j,] <- Ph.Init[j,]
      ord[[i+1]]$like[j] <- final.map$loglike
      best.seq[i+1] <- max(best.seq[i+1],final.map$loglike)
    }
	# sort linkage phases by log-likelihood
    ord.ind <- order(ord[[i+1]]$like, decreasing=TRUE)
    ord[[i+1]]$rf <- ord[[i+1]]$rf[ord.ind,]
    ord[[i+1]]$phase <- ord[[i+1]]$phase[ord.ind,]
    ord[[i+1]]$like <- ord[[i+1]]$like[ord.ind] 
  }
  
  ### positioning after the given sequence
  # get two-point information
  list.init <- phases(make.seq(get(w$twopt),c(w$seq.num[length(w$seq.num)],mrk),twopt=w$twopt))
  rf.init[[(length(w$seq.num))]] <- list.init$rf.init[[1]]
  for(j in 1:(length(w$seq.num)-1)) rf.init[[j]] <- w$seq.rf[j]
  phase.init[[(length(w$seq.num))]] <- list.init$phase.init[[1]]
  for(j in 1:(length(w$seq.num)-1)) phase.init[[j]] <- w$seq.phases[j]
  Ph.Init <- comb.ger(phase.init)
  Rf.Init <- comb.ger(rf.init)
  
  # create last order
  try.ord <- rbind(try.ord,c(w$seq.num,mrk))
  if(verbose) cat("TRY",length(w$seq.num)+1,": ", c(w$seq.num,mrk) ,"\n")
  else cat(".\n")
  flush.console()
  if(nrow(Ph.Init)>1){
    ##Removing ambigous phases
    rm.ab<-rem.amb.ph(M=Ph.Init, w=w, seq.num=c(w$seq.num,mrk))
    Ph.Init <- Ph.Init[rm.ab,]
    Rf.Init <- Rf.Init[rm.ab,]
    if(class(Ph.Init)=="integer"){
      Ph.Init<-matrix(Ph.Init,nrow=1)
      Rf.Init<-matrix(Rf.Init,nrow=1)
    }
  }
  # estimate parameters for all possible linkage phases for this order
  for(j in 1:nrow(Ph.Init)) {
    final.map <- est.map.c(geno=get(w$data.name)$geno[,c(w$seq.num,mrk)],
                           type=get(w$data.name)$segr.type.num[c(w$seq.num,mrk)],
                           phase=Ph.Init[j,],
                           rec=Rf.Init[j,],
                           verbose=FALSE,
                           tol=tol)
    ord[[length(w$seq.num)+1]]$rf[j,] <- final.map$rf
    ord[[length(w$seq.num)+1]]$phase[j,] <- Ph.Init[j,]
    ord[[length(w$seq.num)+1]]$like[j] <- final.map$loglike
    best.seq[length(w$seq.num)+1] <- max(best.seq[length(w$seq.num)+1],final.map$loglike)
  }
  # sort linkage phases by log-likelihood
  ord.ind <- order(ord[[length(w$seq.num)+1]]$like, decreasing=TRUE)
  ord[[length(w$seq.num)+1]]$rf <- ord[[length(w$seq.num)+1]]$rf[ord.ind,]
  ord[[length(w$seq.num)+1]]$phase <- ord[[length(w$seq.num)+1]]$phase[ord.ind,]
  ord[[length(w$seq.num)+1]]$like <- ord[[length(w$seq.num)+1]]$like[ord.ind]
  
  # calculate LOD-Scores (best linkage phase combination for each position)
  LOD <- (best.seq-max(best.seq))/log(10)
  structure(list(ord=ord, LOD=LOD, try.ord=try.ord, data.name=w$data.name, twopt=w$twopt), class = "try")
}



# print method for object class 'try'
print.try <- function(x,j=NULL,...) {
  phases.char <- c("CC","CR","RC","RR")
  marker <- x$try.ord[1,1]
  
  if(is.null(j)) {
    # general summary
    seq <- x$try.ord[1,-1]
    size1 <- max(nchar(seq))
    seq.pr <- format(seq,width=size1)
    size2 <- max(nchar(formatC(x$LOD,format="f",digits=2)))
    LOD.pr <- formatC(round(x$LOD,2),format="f",digits=2,width=size2)
    LOD.pr[which(LOD.pr=="  -0.0")] <- "   0.0"
    
    cat("\nLOD scores correspond to the best linkage phase combination\nfor each position\n")
    cat("\nThe symbol \"*\" outside the box indicates that more than one\nlinkage phase is possible for the corresponding position\n")
    cat(paste("\n\n\t\t  Marker tested: ",marker,"\n\n",sep=""))
    cat("\t\t  Markers",rep("",size1+size2-4),"LOD\n")
    cat(paste("\t\t",paste(rep("=",size1+size2+13),collapse=""),"\n",sep=""))
    cat("\t\t|",rep("",size1+size2+10),"|\n")
    cat("\t\t|",rep("",size1+8),LOD.pr[1]," |")
    ifelse(max(which(x$ord[[1]]$like != -Inf)) != 1,pr <- paste("  ",1,"  *\n",sep=""),pr <- paste("  ",1,"  \n",sep=""))
    cat(pr)
    for(i in 1:length(seq.pr)) {
      cat("\t\t| ",seq.pr[i],rep("",size2+8),"|","\n")
      cat("\t\t|",rep("",size1+8),LOD.pr[i+1]," |")
      ifelse(max(which(x$ord[[i+1]]$like != -Inf)) != 1,pr <- paste("  ",i+1,"  *\n",sep=""),pr <- paste("  ",i+1,"  \n",sep=""))
      cat(pr)
    }
    cat("\t\t|",rep("",size1+size2+10),"|\n")
    cat(paste("\t\t",paste(rep("=",size1+size2+13),collapse=""),"\n",sep=""))
  }
  else {
    # detailed output for a given position
    seq <- x$try.ord[j,]
    size1 <- max(nchar(seq))
    seq.pr <- format(seq,width=size1)
    n.phase <- max(which(x$ord[[j]]$like != -Inf))
    max.like <- -Inf
    for(i in 1:length(x$ord)) max.like <- max(max.like,max(x$ord[[i]]$like[1:n.phase]))
    LOD <- round((x$ord[[j]]$like[1:n.phase]-max.like)/log(10),1)
    LOD.pr <- formatC(LOD,format="f",digits=1,width=6)
    LOD.pr[which(LOD.pr=="  -0.0")] <- "   0.0"
    nest.LOD <- round((x$ord[[j]]$like[1:n.phase]-max(x$ord[[j]]$like[1:n.phase]))/log(10),1)
    nest.LOD.pr <- formatC(nest.LOD,format="f",digits=1,width=6)
    nest.LOD.pr[which(nest.LOD.pr=="  -0.0")] <- "   0.0"
    
    cat("\nLOD is the overall LOD score (among all orders)\n")
    cat("\nNEST.LOD is the LOD score within this order\n")
    cat(paste("\nMarker tested: ",marker,"\n",sep=""))
    
    cat(paste(rep("-",max(2,size1)+5+7*(n.phase)),collapse=""),"\n")
    cat("|",rep("",max(2,size1)+2),rep("|     ",n.phase),"|\n")
    cat("|",seq.pr[1],rep("",max(3-size1,0)),rep("|     ",n.phase),"|\n|",rep("",size1))
    for(i in 2:length(seq.pr)) {
      cat(rep("",max(3-size1,0)+1))
      for(k in 1:n.phase) {
        cat("  | ",phases.char[x$ord[[j]]$phase[k,i-1]])
      }
      cat("  |",rep("",max(3-size1,0)),"\n")
      cat("|",seq.pr[i],rep("",max(3-size1,0)),rep("|     ",n.phase),"|\n|",rep("",size1))
    }
    cat(" ",rep("",max(3-size1,0)),rep("|     ",n.phase),"|\n")
    cat(paste("|",paste(rep("-",max(3,size1)+2+7*(n.phase)),collapse=""),"|",sep=""),"\n")
    cat(paste("| LOD ",rep("",max(size1-3,0)),sep=""))
    for(k in 1:n.phase) {
      cat(paste("|",LOD.pr[k],sep=""))
    }
    cat("|\n")
    cat(paste("|",paste(rep("-",max(3,size1)+2+7*(n.phase)),collapse=""),"|",sep=""),"\n")
    cat(paste("|NEST.",rep("",max(size1-3,0)),sep=""))
    cat(rep("|     ",n.phase),"|\n")
    cat(paste("| LOD ",rep("",max(size1-3,0)),sep=""))
    for(k in 1:n.phase) {
      cat(paste("|",nest.LOD.pr[k],sep=""))
    }
    cat("|\n")
    cat(paste(rep("-",max(3,size1)+4+7*(n.phase)),collapse=""),"\n")
  }
}

# end of file
