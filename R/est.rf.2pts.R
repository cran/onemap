#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: est.rf.2pts.R                                                 #
# Contains: est.rf.2pts, print.rf.2pts                                #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 11/07/2007                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

est.rf.2pts <-
function(w, LOD=3, max.rf=0.35, verbose = TRUE) {
  # checking for correct objects
  if (any(is.na(match(c("geno", "n.ind", "n.mar", "segr.type"),
                      names(w))))) 
    stop("this is not an object of class 'outcross'")
  if (w$n.mar<2)
    stop("there must be at least two markers to proceed with analysis")

  # creating variables (result storage and progress output)
  count <- 0
  prtd <- 0
  progr.past <- 0
  tot <- (w$n.mar*(w$n.mar-1))/2
  analysis <- array(dim=c(tot,4,4))
  phases <- matrix(NA,w$n.mar,w$n.mar)
  recomb <- matrix(NA,w$n.mar,w$n.mar)
  flags <- matrix(NA,w$n.mar,w$n.mar)
  arbitr <- matrix(0,w$n.mar,w$n.mar)
  diag(arbitr) <- NA
  goodness <- character(4)

  if (verbose==TRUE) cat("--Progress: 0%")

  # two-point analysis for each pair of markers
  for (i in 2:w$n.mar) {
    for (j in 1:(i-1)) {
      count <- count+1
      phase <- NA

      # indirect call to C routine
      current <- cr2pts(w$geno[,i],w$geno[,j],w$segr.type[i],w$segr.type[j])
      analysis[acum(i-2)+j,,] <- current

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
  dimnames(recomb) <- dimnames(phases) <- dimnames(flags) <-
    list(colnames(w$geno),colnames(w$geno))
  dimnames(analysis) <- list(NULL,c("A1","A2","A3","A4"),
                             c("Theta","log-Like","Posterior","LODs"))
  rf.2pts <- list(n.mar=w$n.mar, LOD=LOD, max.rf=max.rf, recomb=recomb,
                  phases=phases, analysis=analysis, flags=flags,
                  arbitr=arbitr, segr.type=w$segr.type)
  class(rf.2pts) <- "rf.2pts"
  rf.2pts
}



# print method for object class 'rf.2pts'
print.rf.2pts <-
function(x, mrk1=NULL, mrk2=NULL,...) {
  # checking for correct objects
  if (any(is.na(match(c("n.mar", "LOD", "max.rf", "recomb", "phases",
                        "analysis", "flags", "arbitr", "segr.type"),
                      names(x))))) 
    stop("this is not an object of class 'rf.2pts'")

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
    marnames <- colnames(x$recomb)

    # checking if markers exist and converting character to numeric
    if (is.character(mrk1) && is.character(mrk2)) {
      mrk1name <- mrk1
      mrk2name <- mrk2
      mrk1 <- which(marnames==mrk1)
      mrk2 <- which(marnames==mrk2)
      if (length(mrk1)==0) stop("marker ", mrk1name, " not found")
      if (length(mrk2)==0) stop("marker ", mrk2name, " not found")
    }
    if (is.numeric(mrk1) && is.numeric(mrk2)) {
      cat("  Results of the 2-point analysis for markers:", marnames[mrk1],
          "and", marnames[mrk2], "\n")
      
      if (is.na(x$phases[mrk1,mrk2]) || is.na(x$recomb[mrk1,mrk2]))
        cat("    Analysis unavailable for these markers.\n")
      else {
        # results found
        cat("  Criteria: LOD = ", x$LOD, ", Maximum recombination fraction = ",
            x$max.rf, "\n\n")
        
        # results of the two-point analysis
        if (mrk1 > mrk2) print(x$analysis[acum(mrk1-2)+mrk2,,])
        else print(x$analysis[acum(mrk2-2)+mrk1,,])
        
        cat("\n  Most probable linkage phase (assignment): ")
        if (x$phases[mrk1,mrk2]=="C/C")
          cat("A1 - coupling / coupling\n")
        else if (x$phases[mrk1,mrk2]=="C/R")
          cat("A2 - coupling / repulsion\n")
        else if (x$phases[mrk1,mrk2]=="R/C")
          cat("A3 - repulsion / coupling \n")
        else if (x$phases[mrk1,mrk2]=="R/R")
          cat("A4 - repulsion / repulsion \n")
        else if (x$phases[mrk1,mrk2]=="ns")
          cat("Markers not linked \n")
        cat("  Map distance in cM (KOSAMBI): ",
            kosambi(x$recomb[mrk1,mrk2]), "\n")
        cat("  Map distance in cM (HALDANE): ",
            haldane(x$recomb[mrk1,mrk2]), "\n") 
        if (!is.na(x$flags[mrk1,mrk2]) && x$flags[mrk1,mrk2]==1)
          cat("\n  Warning: more than 1 probable assignment!\n")
        if (!is.na(x$arbitr[mrk1,mrk2]) && x$arbitr[mrk1,mrk2]==1)
          cat("\n  Warning: assignment was arbitrarily chosen!\n")
      }
    }
    else stop("'mrk1' and 'mrk2' must be of the same type (",
              dQuote("numeric"), " or ", dQuote("character"), ")")
  }
}


