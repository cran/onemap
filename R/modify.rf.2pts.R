#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: modify.rf.2pts.R                                              #
# Contains: modify.rf.2pts                                            #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 11/07/2007                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

modify.rf.2pts <-
function(w, LOD=NULL, max.rf=NULL, verbose = TRUE) {
  # checking for correct object
  if (any(is.na(match(c("n.mar", "LOD", "max.rf", "recomb", "phases",
                        "analysis", "flags", "arbitr", "segr.type"),
                      names(w))))) 
    stop("this is not an object of class 'rf.2pts'")
  if (is.null(LOD) && is.null(max.rf)) w  # no action required
  else {
    if (is.null(LOD)) LOD <- w$LOD
    if (is.null(max.rf)) max.rf <- w$max.rf
    
    # creating variables (result storage and progress output)
    count <- 0
    prtd <- 0
    progr.past <- 0
    tot <- (w$n.mar*(w$n.mar-1))/2
    phases <- matrix(NA,w$n.mar,w$n.mar)
    recomb <- matrix(NA,w$n.mar,w$n.mar)
    flags <- matrix(NA,w$n.mar,w$n.mar)
    arbitr <- matrix(0,w$n.mar,w$n.mar)
    diag(arbitr) <- NA
    goodness <- character(4)
    
    if (verbose==TRUE) cat("--Progress: 0%")

    # updating two-point analysis for each pair of markers
    for (i in 2:w$n.mar) {
      for (j in 1:(i-1)) {
        count <- count+1
        phase <- NA

        # analyses are not performed again
        # getting previous results
        current <- w$analysis[acum(i-2)+j,,]
        
        # choosing most probable assignments
        probab <- which(current[,4]>(max(current[,4]-0.005)) &
                        current[,4]<=max(current[,4]))
        for (a in probab) {
          if (current[a,1] <= max.rf) goodness[a] <- "**"
          else goodness[a] <- "-"
        }
        goodness[-probab] <- "-"
        phase <- which(goodness=="**")
        
        # no assignment meets 'max.rf' criterion
        if (length(phase)==0) {
          phases[i,j] <- phases[j,i] <- "ns"
          recomb[i,j] <- recomb[j,i] <- 0.5
        }
        else {  # one or more assignments meet 'max.rf' criterion
          
          # if more than one probable phase, choose the first one
          if (length(phase)>1) {
            phase <- phase[1]
            flags[i,j] <- flags[j,i] <- 1
          }
          else flags[i,j] <- flags[j,i] <- 0
          
          # checking if phase meets 'LOD' criterion
          if (current[phase,4] >= LOD) {
            recomb[i,j] <- recomb[j,i] <- current[phase,1]
            if (phase==1) phases[i,j] <- phases[j,i] <- "C/C"
            else if (phase==2) phases[i,j] <- phases[j,i] <- "C/R"
            else if (phase==3) phases[i,j] <- phases[j,i] <- "R/C"
            else if (phase==4) phases[i,j] <- phases[j,i] <- "R/R"
          }
          else {  # 'LOD' criterion not met
            phases[i,j] <- phases[j,i] <- "ns"
            recomb[i,j] <- recomb[j,i] <- 0.5
            flags[i,j] <- flags[j,i] <- NA
          }
        }

        # current progress output
        if (verbose) {
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
    
    # results
    dimnames(recomb) <- dimnames(phases) <- dimnames(flags) <- dimnames(w$recomb)
    rf.2pts <- list(n.mar=w$n.mar, LOD=LOD, max.rf=max.rf, recomb=recomb,
                    phases=phases, analysis=w$analysis, flags=flags,
                    arbitr=arbitr, segr.type=w$segr.type)
    class(rf.2pts) <- "rf.2pts"
    rf.2pts
  }
}

