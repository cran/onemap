#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: phases.R                                                      #
# Contains: phases                                                    #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2009, Marcelo Mollinari                               #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function gathers information from two-point analyses and
# makes a list with 'plausible' linkage phases, with the respective
# recombination fractions
phases <- 
function(w, LOD=0, max.rf=0.50) {
  # checking for correct object
  if(!any(class(w)=="sequence")) stop(deparse(substitute(w))," is not an object of class 'sequence'")
  
  tot <- choose(length(w$seq.num),2)
  # initial values for the recombination fraction
  rf.init <- vector("list",tot)
  # initially relevant linkage phases
  phase.init <- vector("list",tot)
  
  temp <- matrix(NA,4,2)
  
  for (i in 2:length(w$seq.num)) {
    for (j in 1:(i-1)) {
	  # recover values from two-point analyses
      big <- pmax.int(w$seq.num[i],w$seq.num[j])
      small <- pmin.int(w$seq.num[i],w$seq.num[j])
      temp <- get(w$twopt)$analysis[acum(big-2)+small,,]
      
      # check which assignments meet the criteria
      relevant <- which(temp[,2]>(max(temp[,2])-0.005)) # maximum LOD scores
      phase.init[[acum(i-2)+j]] <- relevant[which((temp[relevant,1] <= max.rf & temp[relevant,2] >= LOD))]
      if (length(phase.init[[acum(i-2)+j]]) == 0) phase.init[[acum(i-2)+j]] <- 1:4
      rf.init[[acum(i-2)+j]] <- temp[phase.init[[acum(i-2)+j]],1]
    }
  }
  
  list(phase.init=phase.init, rf.init=rf.init)
}

# end of file
