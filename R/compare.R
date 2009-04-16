#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: compare.R                                                     #
# Contains: compare, print.compare                                    #
#                                                                     #
# Written by Gabriel R A Margarido & Marcelo Mollinari                #
# copyright (c) 2009, Gabriel R A Margarido & Marcelo Mollinari       #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function evaluates all n!/2 possible orders for n markers
compare <- function(w,n.best=50,tol=10E-6,verbose=FALSE) {
  # checking for correct objects
  if(!any(class(w)=="sequence")) stop(sQuote(deparse(substitute(w)))," is not an object of class 'sequence'")
  if(length(w$seq.num) > 5) cat("WARNING: this operation may take a VERY long time\n")
  flush.console()
  
  if(length(w$seq.num) == 2) return(map(w, tol=tol)) # nothing to be done for 2 markers
  else {
    # allocating variables
    rf.init <- vector("list",length(w$seq.num)-1)
    phase.init <- vector("list",length(w$seq.num)-1)
    best.ord <- matrix(NA,(n.best+1),length(w$seq.num))
    best.ord.rf <- matrix(NA,(n.best+1),length(w$seq.num)-1)
    best.ord.phase <- matrix(NA,(n.best+1),length(w$seq.num)-1)
    best.ord.like <- best.ord.LOD <- rep(-Inf,(n.best+1))
	
	# 'phases' gathers information from two-point analyses
    list.init <- phases(w)
	
	# 'perm.pars' generates all n!/2 orders
    all.ord <- perm.pars(w$seq.num)
	cat("\nComparing",nrow(all.ord),"orders:     \n\n")
    if (verbose){ 
      for(i in 1:nrow(all.ord)){
                                        # print output for each order
        cat("Order", i, ":", all.ord[i,], "\n")
        flush.console()  
                                        # get initial values for the HMM
        all.match <- match(all.ord[i,],w$seq.num)
        for(j in 1:(length(w$seq.num)-1)){
          if(all.match[j] > all.match[j+1]){
            rf.init[[j]] <- list.init$rf.init[[acum(all.match[j]-2)+all.match[j+1]]]
            phase.init[[j]] <- list.init$phase.init[[acum(all.match[j]-2)+all.match[j+1]]]
          }
          else {
            rf.init[[j]] <- list.init$rf.init[[acum(all.match[j+1]-2)+all.match[j]]]
            phase.init[[j]] <- list.init$phase.init[[acum(all.match[j+1]-2)+all.match[j]]]
          }
        }    
        Ph.Init <- comb.ger(phase.init)
        Rf.Init <- comb.ger(rf.init)
        
        for(j in 1:nrow(Ph.Init)){
                                        # estimate parameters
          final.map <- est.map.c(geno=get(w$data.name)$geno[,all.ord[i,]],
                                 type=get(w$data.name)$segr.type.num[all.ord[i,]],
                                 phase=Ph.Init[j,],
                                 rec=Rf.Init[j,],
                                 verbose=FALSE,
                                 tol=tol)
          best.ord[(n.best+1),] <- all.ord[i,]
          best.ord.rf[(n.best+1),] <- final.map$rf
          best.ord.phase[(n.best+1),] <- Ph.Init[j,]
          best.ord.like[(n.best+1)] <- final.map$loglike
          
                                        # arrange orders according to the likelihood
          like.order <- order(best.ord.like, decreasing=TRUE)
          best.ord <- best.ord[like.order,]
          best.ord.rf <- best.ord.rf[like.order,]
          best.ord.phase <- best.ord.phase[like.order,]
          best.ord.like <- sort(best.ord.like, decreasing=TRUE)
        }
      }
    }
    else{
      nc<-NA
      out.pr <- seq(from=1,to=nrow(all.ord), length.out=20)
      cat("    ")
      for(i in 1:nrow(all.ord)){
                                        # print output for each order
        if (sum(i == round(out.pr))){
          cat(rep("\b",nchar(nc)+1),sep="")
          nc<-round(i*100/nrow(all.ord))
          cat(nc,"%", sep="")
          flush.console()
        }
                                        # get initial values for the HMM
        all.match <- match(all.ord[i,],w$seq.num)
        for(j in 1:(length(w$seq.num)-1)){
          if(all.match[j] > all.match[j+1]){
            rf.init[[j]] <- list.init$rf.init[[acum(all.match[j]-2)+all.match[j+1]]]
            phase.init[[j]] <- list.init$phase.init[[acum(all.match[j]-2)+all.match[j+1]]]
          }
          else {
            rf.init[[j]] <- list.init$rf.init[[acum(all.match[j+1]-2)+all.match[j]]]
            phase.init[[j]] <- list.init$phase.init[[acum(all.match[j+1]-2)+all.match[j]]]
          }
        }    
        Ph.Init <- comb.ger(phase.init)
        Rf.Init <- comb.ger(rf.init)
        
        for(j in 1:nrow(Ph.Init)){
                                        # estimate parameters
          final.map <- est.map.c(geno=get(w$data.name)$geno[,all.ord[i,]],
                                 type=get(w$data.name)$segr.type.num[all.ord[i,]],
                                 phase=Ph.Init[j,],
                                 rec=Rf.Init[j,],
                                 verbose=FALSE,
                                 tol=tol)
          best.ord[(n.best+1),] <- all.ord[i,]
          best.ord.rf[(n.best+1),] <- final.map$rf
          best.ord.phase[(n.best+1),] <- Ph.Init[j,]
          best.ord.like[(n.best+1)] <- final.map$loglike
          
                                        # arrange orders according to the likelihood
          like.order <- order(best.ord.like, decreasing=TRUE)
          best.ord <- best.ord[like.order,]
          best.ord.rf <- best.ord.rf[like.order,]
          best.ord.phase <- best.ord.phase[like.order,]
          best.ord.like <- sort(best.ord.like, decreasing=TRUE)
        }
      }
    }
    cat("\n\nFinished\n\n")
    best.ord.LOD <- round((best.ord.like-max(best.ord.like))/log(10),4)
    
    structure(list(best.ord = best.ord, best.ord.rf = best.ord.rf,
                   best.ord.phase = best.ord.phase, best.ord.like = best.ord.like,
                   best.ord.LOD = best.ord.LOD, data.name=w$data.name, twopt=w$twopt), class = "compare")
  }
}

# print method for object class 'compare'
print.compare <-
function(x,...) {
  phases.char <- c("CC","CR","RC","RR")
  n.ord <- max(which(head(x$best.ord.LOD,-1) != -Inf))

  unique.orders <- unique(x$best.ord[1:n.ord,])
  n.ord.nest <- dim(unique.orders)[1]
  phases.nested <- vector("list",n.ord.nest)
  LOD <- vector("list",n.ord.nest)

  for (i in 1:n.ord.nest) {
    same.order <- which(apply(x$best.ord[1:n.ord,],1,function(x) all(x==unique.orders[i,])))
    ifelse(length(same.order)==1,phases.nested[[i]] <- t(as.matrix(x$best.ord.phase[same.order,])),phases.nested[[i]] <- x$best.ord.phase[same.order,])
    LOD[[i]] <- x$best.ord.LOD[same.order]
  }

  skip <- c(nchar(n.ord.nest),max(nchar(unique.orders[1,])+2))
  leng.print <- nchar(paste("order ",format(n.ord.nest,width=skip[1]),":  ",paste(format(unique.orders[1,],width=skip[2]),collapse=""),"     ",format(11.11,digits=2,format="f",width=6),"     ",format(11.11,digits=2,format="f",width=6),"\n",sep=""))
  cat("\nNumber of orders:",n.ord,"\n")
  cat(paste("Best ",n.ord.nest," unique orders",paste(rep(" ",leng.print-37),collapse=""),"LOD    Nested LOD","\n",sep=""))
  cat(paste(rep("-",leng.print),collapse=""),"\n")

  for (i in 1:n.ord.nest) {
    cat(paste("order ",format(i,width=skip[1]),":  ",paste(format(unique.orders[i,],width=skip[2]),collapse=""),"\n",sep=""))
    for (j in 1:dim(phases.nested[[i]])[1]) {
      cat(paste("\t",paste(rep(" ",1+skip[1]+skip[2]),collapse=""),paste(format(phases.char[phases.nested[[i]][j,]],width=skip[2]),collapse=""),"     ",formatC(round(LOD[[i]][j],2),digits=2,format="f",width=6),"     ",formatC(round(LOD[[i]][j]-LOD[[i]][1],2),digits=2,format="f",width=6),"\n",sep=""))
    }
    cat(paste(rep("-",leng.print),collapse=""))
    cat("\n")
  }
}

# end of file
