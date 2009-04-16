#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: comb.ger.R                                                    #
# Contains: comb, comb.ger                                            #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2009, Marcelo Mollinari                               #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function combines two linkage phase vectors 
comb <-
function(x,y) {
  count <- 0
  M <- matrix(NA, nrow(x)*nrow(y), 1+ncol(y))
  for(i in 1:nrow(x)) {
    for(j in 1:nrow(y)) {
      count <- count+1
      M[count,] <- c(x[i,],y[j,])
    }
  }
  return(M)
}

# This function makes all possible combinations of a list
# containing the linkage phase vectors for each interval
comb.ger <- function(f){
  M <- as.matrix(f[[length(f)]])
  if (length(f)==1) return(M)
  else{
    for(i in (length(f)-1):1){
      M <- comb(as.matrix(f[[i]]),M)
    }
    return(M)
  }
}

# end of file