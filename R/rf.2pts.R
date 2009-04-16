#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: rf.2pts.R                                                     #
# Contains: rf.2pts, print.rf.2pts                                    #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007-9, Gabriel R A Margarido                         #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# Function to perform two-point analyses for all markers in a data set
rf.2pts <- 
function(w, LOD=3, max.rf=0.50, verbose = TRUE) {
  # checking for correct object
  if(!any(class(w)=="outcross")) stop(deparse(substitute(w))," is not an object of class 'outcross'")
  if (w$n.mar<2) stop("there must be at least two markers to proceed with analysis")

  # creating variables (result storage and progress output)
  count <- 0
  prtd <- 0
  progr.past <- 0
  tot <- choose(w$n.mar,2)
  analysis <- array(dim=c(tot,4,2))
  
  if (verbose==TRUE) {
    cat("--Progress: 0%")
	
	for (i in 2:w$n.mar) {
      for (j in 1:(i-1)) {
        count <- count+1

        # indirect call to C routine
        current <- cr2pts(w$geno[,i],w$geno[,j],w$segr.type.num[i],w$segr.type.num[j])
        analysis[acum(i-2)+j,,] <- current[,c(1,4)]

        # current progress output
        if (count==tot) cat("\n--Finished\n")
        else {
          progr <- round(100*(count/tot),0)
          if ((progr - progr.past) >= 5) {
            if (prtd==10) {
              cat("\n            ")
              prtd <- 0
            }
            cat(paste("....",progr,"%",sep=""))
            flush.console()
            progr.past <- progr
            prtd <- prtd+1
          }
        }
      }
    }
  }
  else {
    # two-point analysis for each pair of markers
    for (i in 2:w$n.mar) {
      for (j in 1:(i-1)) {
        count <- count+1

        # indirect call to C routine
        current <- cr2pts(w$geno[,i],w$geno[,j],w$segr.type.num[i],w$segr.type.num[j])
        analysis[acum(i-2)+j,,] <- current[,c(1,4)]
      }
    }
  }
  
  # results
  dimnames(analysis) <- list(NULL,c("1","2","3","4"),
                             c("Theta","LODs"))
  structure(list(data.name=as.character(sys.call())[2], n.mar=w$n.mar,marnames=colnames(w$geno),
                 LOD=LOD, max.rf=max.rf, input=w$input, analysis=analysis), class = "rf.2pts")
}



# print method for object class 'rf.2pts'
print.rf.2pts <-
function(x, mrk1=NULL, mrk2=NULL,...) {
  # checking for correct object
  if(!any(class(x)=="rf.2pts")) stop(deparse(substitute(x))," is not an object of class 'rf.2pts'")

  if (is.null(mrk1) || is.null(mrk2)) {
    # printing a brief summary
    cat("  This is an object of class 'rf.2pts'\n")
    cat("\n  Criteria: LOD =", x$LOD, ", Maximum recombination fraction =",
        x$max.rf, "\n")
    cat("\n  This object is too complex to print\n")
    cat("  Type 'print(object,mrk1=marker,mrk2=marker)' to see the analysis for two markers\n")
    cat("    mrk1 and mrk2 can be the names or numbers of both markers\n")
  }
  else {
    # printing detailed results for two markers
	
    # checking if markers exist and converting character to numeric
    if (is.character(mrk1) && is.character(mrk2)) {
      mrk1name <- mrk1
      mrk2name <- mrk2
      mrk1 <- which(x$marnames==mrk1)
      mrk2 <- which(x$marnames==mrk2)
      if (length(mrk1)==0) stop("marker ", mrk1name, " not found")
      if (length(mrk2)==0) stop("marker ", mrk2name, " not found")
    }
    if (is.numeric(mrk1) && is.numeric(mrk2)) {
      cat("  Results of the 2-point analysis for markers:", x$marnames[mrk1],
          "and", x$marnames[mrk2], "\n")

        # results found
        cat("  Criteria: LOD = ", x$LOD, ", Maximum recombination fraction = ",
            x$max.rf, "\n\n")
        
        # results of the two-point analysis
        if (mrk1 > mrk2) print(x$analysis[acum(mrk1-2)+mrk2,,])
        else print(x$analysis[acum(mrk2-2)+mrk1,,])
	}
    else stop("'mrk1' and 'mrk2' must be of the same type \"numeric\" or \"character\"")
  }
}

# end of file
