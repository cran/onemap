#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: map.R                                                         #
# Contains: map                                                       #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2009, Gabriel R A Margarido                           #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function constructs the linkage map for a set of markers in a given order
map <-
function(w,tol=10E-6) {
  # checking for correct object
  if(!any(class(w)=="sequence")) stop(deparse(substitute(w))," is not an object of class 'sequnece'")
  if(length(w$seq.num) < 2) stop("The sequence must have at least 2 markers")

  if((w$seq.phases != -1) && (w$seq.rf != -1)  && !is.null(w$seq.like)) {
    # if the linkage phases and recombination fractions have already been
    # estimated, nothings needs to be done
    cat("\nNo action required.\n")
    cat("Returning sequence.\n\n")
    w
  }
  else if((w$seq.phases != -1) && (w$seq.rf == -1) && is.null(w$seq.like)) {
    # if the linkage phases are provided but the recombination fractions have
    # not yet been estimated, this is done here
	
	# gather two-point information
    rf.init <- numeric(length(w$seq.num)-1)
    for(i in 1:(length(w$seq.num)-1)) {
      if(w$seq.num[i] > w$seq.num[i+1])
        rf.init[i] <- get(w$twopt)$analysis[acum(w$seq.num[i]-2)+w$seq.num[i+1],w$seq.phases[i],1]
      else
        rf.init[i] <- get(w$twopt)$analysis[acum(w$seq.num[i+1]-2)+w$seq.num[i],w$seq.phases[i],1]
    }
	# estimate parameters
    final.map <- est.map.c(geno=get(w$data.name)$geno[,w$seq.num],
                           type=get(w$data.name)$segr.type.num[w$seq.num],
                           phase=w$seq.phases,
                           rec=rf.init,
                           verbose=FALSE,
                           tol=tol)
    structure(list(seq.num=w$seq.num, seq.phases=w$seq.phases, seq.rf=final.map$rf,
                   seq.like=final.map$loglike, data.name=w$data.name, twopt=w$twopt), class = "sequence")
  }
  else if((w$seq.phases == -1) && (w$seq.rf == -1) && is.null(w$seq.like)) {
    # if only the marker order is provided, without predefined linkage phases,
	# a search for the best combination of phases is performed and recombination
	# fractions are estimated
    seq.phase <- numeric(length(w$seq.num)-1)
    results <- list(rep(NA,4),rep(-Inf,4))
    
	# linkage map is started with the first two markers in the sequence
	# gather two-point information for this pair
    phase.init <- vector("list",1)
	list.init <- phases(make.seq(get(w$twopt),c(w$seq.num[1],w$seq.num[2]),twopt=w$twopt))
    phase.init[[1]] <- list.init$phase.init[[1]]
	Ph.Init <- comb.ger(phase.init)
    for(j in 1:nrow(Ph.Init)) {
	  # call to 'map' function with predefined linkage phase
	  temp <- map(make.seq(get(w$twopt),w$seq.num[1:2],phase=Ph.Init[j],twopt=w$twopt))
      results[[1]][j] <- temp$seq.phases
      results[[2]][j] <- temp$seq.like
    }
	seq.phase[1] <- results[[1]][which.max(results[[2]])] # best linkage phase is chosen
    
	if(length(w$seq.num) > 2) {
	  # for sequences with three or more markers, these are added sequentially
	  for(mrk in 2:(length(w$seq.num)-1)) {
	    results <- list(rep(NA,4),rep(-Inf,4))
		
		# gather two-point information
	    phase.init <- vector("list",mrk)
        list.init <- phases(make.seq(get(w$twopt),c(w$seq.num[mrk],w$seq.num[mrk+1]),twopt=w$twopt))
            phase.init[[mrk]] <- list.init$phase.init[[1]]
            for(j in 1:(mrk-1)) phase.init[[j]] <- seq.phase[j]
            Ph.Init <- comb.ger(phase.init)
            for(j in 1:nrow(Ph.Init)) {
		  # call to 'map' function with predefined linkage phases
	      temp <- map(make.seq(get(w$twopt),w$seq.num[1:(mrk+1)],phase=Ph.Init[j,],twopt=w$twopt))
              results[[1]][j] <- temp$seq.phases[mrk]
              results[[2]][j] <- temp$seq.like
            }
	    seq.phase[mrk] <- results[[1]][which.max(results[[2]])] # best combination of phases is chosen
	  }
	}
	# one last call to map function, with the final map
    map(make.seq(get(w$twopt),w$seq.num,phase=seq.phase,twopt=w$twopt))
  }
  #else SHOULD NOT GET HERE
}

# end of file
