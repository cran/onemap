#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: ripple.R                                                      #
# Contains: ripple                                                    #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2009, Gabriel R A Margarido                           #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function searches for alternative orders, by comparing all possible
# orders of subsets of markers
ripple <-
function(w,ws=4,LOD=3,tol=10E-5) {
  # checking for correct objects
  if(!any(class(w)=="sequence")) stop(deparse(substitute(w))," is not an object of class 'sequence'")
  if(ws < 2) stop("ws must be greater than or equal to 2")
  if(ws > 5) cat("WARNING: this operation may take a VERY long time\n\n")
  flush.console()

  len <- length(w$seq.num)
  # computations unnecessary in this case
  if (len <= ws) stop("Length of sequence ", deparse(substitute(w))," is smaller than ws. You can use the compare function instead")
  
  # convert numerical linkage phases to strings, to facilitate rearrangements
  link.phases <- matrix(NA,len,2)
  link.phases[1,] <- rep(1,2)
  for (i in 1:length(w$seq.phases)) {
    switch(EXPR=w$seq.phases[i],
           link.phases[i+1,] <- link.phases[i,]*c(1,1),
           link.phases[i+1,] <- link.phases[i,]*c(1,-1),
           link.phases[i+1,] <- link.phases[i,]*c(-1,1),
           link.phases[i+1,] <- link.phases[i,]*c(-1,-1),
           )
  }
  
  # allocate variables
  rf.init <- rep(NA,len-1)
  phase <- rep(NA,len-1)
  tot <- prod(1:ws)
  best.ord.phase <- matrix(NA,tot,len-1)
  best.ord.like <- best.ord.LOD <- rep(-Inf,tot)
  
  # gather two-point information
  list.init <- phases(w)
  
  ### first position
  cat(w$seq.num[1:ws],"...")
  
  # create all possible alternative orders for the first subset
  all.ord <- t(apply(perm.tot(head(w$seq.num,ws)),1,function(x) c(x,tail(w$seq.num,-ws))))
  for(i in 1:nrow(all.ord)){
    all.match <- match(all.ord[i,],w$seq.num)
	
	# rearrange linkage phases according to the current order
    for(j in 1:ws) {
      temp <- paste(as.character(link.phases[all.match[j],]*link.phases[all.match[j+1],]),collapse=".")
      phase[j] <- switch(EXPR=temp,
                         '1.1'   = 1,
                         '1.-1'  = 2,
                         '-1.1'  = 3,
                         '-1.-1' = 4
                         )
    }
    if(len > (ws+1)) phase[(ws+1):length(phase)] <- tail(w$seq.phase,-ws)
	
	# get initial values for recombination fractions
    for(j in 1:(len-1)){
      if(all.match[j] > all.match[j+1]) ind <- acum(all.match[j]-2)+all.match[j+1]
      else ind <- acum(all.match[j+1]-2)+all.match[j]
	  if(length(which(list.init$phase.init[[ind]] == phase[j])) == 0) rf.init[j] <- 0.49 # for safety reasons
	  else rf.init[j] <- list.init$rf.init[[ind]][which(list.init$phase.init[[ind]] == phase[j])]
    }
	# estimate parameters
    final.map <- est.map.c(geno=get(w$data.name)$geno[,all.ord[i,]],
                           type=get(w$data.name)$segr.type.num[all.ord[i,]],
                           phase=phase,
                           rec=rf.init,
                           verbose=FALSE,
                           tol=tol)
    best.ord.phase[i,] <- phase
    best.ord.like[i] <- final.map$loglike
  }
  # calculate LOD-Scores for alternative orders
  best.ord.LOD <- round((best.ord.like-max(best.ord.like))/log(10),2)
  # which orders will be printed
  which.LOD <- which(best.ord.LOD > -LOD)

  if(length(which.LOD) > 1) {
    # if any order to print, sort by LOD-Score
    order.print <- order(best.ord.LOD,decreasing=TRUE)
    all.ord <- all.ord[order.print,]
    best.ord.phase <- best.ord.phase[order.print,]
    best.ord.LOD <- best.ord.LOD[order.print]
	
	# display results
	which.LOD <- which(best.ord.LOD > -LOD)
	LOD.print <- format(best.ord.LOD,digits=2,nsmall=2)
    cat("\n  Alternative orders:\n")
    for(j in which.LOD) {
      cat("  ",all.ord[j,1:(ws+1)],ifelse(len > (ws+1),"... : ",": "),LOD.print[j],"( linkage phases:",best.ord.phase[j,1:ws],ifelse(len > (ws+1),"... )\n",")\n"))
    }
    cat("\n")
  }
  else cat(" OK\n\n")
  
  ### middle positions
  if (len > (ws+1)) {
    for (p in 2:(len-ws)) {
      cat(w$seq.num[p:(p+ws-1)],"...")
      
	  # create all possible alternative orders for the first subset
      all.ord <- t(apply(perm.tot(w$seq.num[p:(p+ws-1)]),1,function(x) c(head(w$seq.num,p-1),x,tail(w$seq.num,-p-ws+1))))
      for(i in 1:nrow(all.ord)){
        all.match <- match(all.ord[i,],w$seq.num)
		
		# rearrange linkage phases according to the current order
        if(p > 2) phase[1:(p-2)] <- head(w$seq.phase,p-2)
        for(j in (p-1):(p+ws-1)) {
          temp <- paste(as.character(link.phases[all.match[j],]*link.phases[all.match[j+1],]),collapse=".")
          phase[j] <- switch(EXPR=temp,
                             '1.1'   = 1,
                             '1.-1'  = 2,
                             '-1.1'  = 3,
                             '-1.-1' = 4
                             )
        }
        if(p < (len-ws)) phase[(p+ws):length(phase)] <- tail(w$seq.phase,len-p-ws)
		
		# get initial values for recombination fractions
        for(j in 1:(len-1)){
          if(all.match[j] > all.match[j+1]) ind <- acum(all.match[j]-2)+all.match[j+1]
          else ind <- acum(all.match[j+1]-2)+all.match[j]
          if(length(which(list.init$phase.init[[ind]] == phase[j])) == 0) rf.init[j] <- 0.49 # for safety reasons
          else rf.init[j] <- list.init$rf.init[[ind]][which(list.init$phase.init[[ind]] == phase[j])]
        }
		# estimate parameters
        final.map <- est.map.c(geno=get(w$data.name)$geno[,all.ord[i,]],
                               type=get(w$data.name)$segr.type.num[all.ord[i,]],
                               phase=phase,
                               rec=rf.init,
                               verbose=FALSE,
                               tol=tol)
        best.ord.phase[i,] <- phase
        best.ord.like[i] <- final.map$loglike
      }
	  # calculate LOD-Scores for alternative orders
      best.ord.LOD <- round((best.ord.like-max(best.ord.like))/log(10),2)
	  # which orders will be printed
      which.LOD <- which(best.ord.LOD > -LOD)
      
      if(length(which.LOD) > 1) {
        # if any order to print, sort by LOD-Score
        order.print <- order(best.ord.LOD,decreasing=TRUE)
        all.ord <- all.ord[order.print,]
        best.ord.phase <- best.ord.phase[order.print,]
        best.ord.LOD <- best.ord.LOD[order.print]
		
		# display results
        which.LOD <- which(best.ord.LOD > -LOD)
        LOD.print <- format(best.ord.LOD,digits=2,nsmall=2)
	    cat("\n  Alternative orders:\n")
        for(j in which.LOD) {
          cat(ifelse(p>2,"  ...","  "),all.ord[j,(p-1):(p+ws)],ifelse((p+ws)<len,"... : ",": "),LOD.print[j],"( linkage phases:",ifelse(p>2,"...","\b"),best.ord.phase[j,(p-1):(p+ws-1)],ifelse((p+ws)<len,"... )\n",")\n"))
        }
        cat("\n")
      }
      else cat(" OK\n\n")
    }
  }
  
  ### last position
  cat(tail(w$seq.num,ws),"...")
  
  # create all possible alternative orders for the first subset
  all.ord <- t(apply(perm.tot(tail(w$seq.num,ws)),1,function(x) c(head(w$seq.num,-ws),x)))
  for(i in 1:nrow(all.ord)){
    all.match <- match(all.ord[i,],w$seq.num)
	
	# rearrange linkage phases according to the current order
    if(len > (ws+1)) phase[1:(len-ws-1)] <- head(w$seq.phase,-ws)
    for(j in (len-ws):(len-1)) {
      temp <- paste(as.character(link.phases[all.match[j],]*link.phases[all.match[j+1],]),collapse=".")
      phase[j] <- switch(EXPR=temp,
                         '1.1'   = 1,
                         '1.-1'  = 2,
                         '-1.1'  = 3,
                         '-1.-1' = 4
                         )
    }
    
	# get initial values for recombination fractions
    for(j in 1:(len-1)){
      if(all.match[j] > all.match[j+1]) ind <- acum(all.match[j]-2)+all.match[j+1]
      else ind <- acum(all.match[j+1]-2)+all.match[j]
	  if(length(which(list.init$phase.init[[ind]] == phase[j])) == 0) rf.init[j] <- 0.49 # for safety reasons
	  else rf.init[j] <- list.init$rf.init[[ind]][which(list.init$phase.init[[ind]] == phase[j])]
    }
	# estimate parameters
    final.map <- est.map.c(geno=get(w$data.name)$geno[,all.ord[i,]],
                           type=get(w$data.name)$segr.type.num[all.ord[i,]],
                           phase=phase,
                           rec=rf.init,
                           verbose=FALSE,
                           tol=tol)
    best.ord.phase[i,] <- phase
    best.ord.like[i] <- final.map$loglike
  }
  # calculate LOD-Scores for alternative orders
  best.ord.LOD <- round((best.ord.like-max(best.ord.like))/log(10),2)
  # which orders will be printed
  which.LOD <- which(best.ord.LOD > -LOD)
  
  if(length(which.LOD) > 1) {
    # if any order to print, sort by LOD-Score
    order.print <- order(best.ord.LOD,decreasing=TRUE)
    all.ord <- all.ord[order.print,]
    best.ord.phase <- best.ord.phase[order.print,]
    best.ord.LOD <- best.ord.LOD[order.print]
	
	# display results
	which.LOD <- which(best.ord.LOD > -LOD)
	LOD.print <- format(best.ord.LOD,digits=2,nsmall=2)
    cat("\n  Alternative orders:\n")
    for(j in which.LOD) {
      cat(ifelse(len > (ws+1),"  ...","  "),all.ord[j,(len-ws):len],": ",LOD.print[j],"( linkage phases:",ifelse(len > (ws+1),"...","\b"),best.ord.phase[j,(len-ws):(len-1)],")\n")
    }
    cat("\n")
  }
  else cat(" OK\n\n") 
}

# end of file
