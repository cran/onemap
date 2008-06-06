#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: arbitr.rf.2pts.R                                              #
# Contains: arbitr.rf.2pts                                            #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 11/07/2007                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

arbitr.rf.2pts <-
function(w, mrk1name=NULL, mrk2name=NULL, new.phase=NULL) {
  # checking for correct objects
  if (any(is.na(match(c("n.mar", "LOD", "max.rf", "recomb", "phases",
                        "analysis", "flags", "arbitr", "segr.type"),
                      names(w)))))
    stop("this is not an object of class 'rf.2pts'")
  if (is.null(mrk1name) || is.null(mrk2name)) stop("you must specify 2 markers")
  if (is.null(new.phase)) stop("you must specify a new linkage phase")
  
  # checking if markers really exist and analysis is available
  if (is.character(mrk1name) && is.character(mrk2name)) {
    marnames <- colnames(w$recomb)
    mrk1 <- which(marnames==mrk1name)
    mrk2 <- which(marnames==mrk2name)
    if (length(mrk1)==0) stop("marker ", mrk1name, " not found")
    if (length(mrk2)==0) stop("marker ", mrk2name, " not found")
  }
  else stop("'mrk1name' and 'mrk2name' must be of the type character")
  if (is.na(w$phases[mrk1,mrk2]) || is.na(w$recomb[mrk1,mrk2]))
    stop("analysis unavailable for these markers")

  # assuring mrk1 > mrk2
  if (mrk1 < mrk2) {
    temp <- mrk1
    mrk1 <- mrk2
    mrk2 <- temp
    rm(temp)
  }

  # copying objects
  recomb <- w$recomb
  phases <- w$phases
  arbitr <- w$arbitr

  # updating recombination fraction and linkage phase
  # TODO: remove if-elses and use switch
  if (new.phase=="A1") {
    phases[mrk1,mrk2] <- phases[mrk2,mrk1] <- "C/C"
    recomb[mrk1,mrk2] <- recomb[mrk2,mrk1] <- w$analysis[acum(mrk1-2)+mrk2,1,1]
    arbitr[mrk1,mrk2] <- arbitr[mrk2,mrk1] <- 1
  }
  else if (new.phase=="A2") {
    phases[mrk1,mrk2] <- phases[mrk2,mrk1] <- "C/R"
    recomb[mrk1,mrk2] <- recomb[mrk2,mrk1] <- w$analysis[acum(mrk1-2)+mrk2,2,1]
    arbitr[mrk1,mrk2] <- arbitr[mrk2,mrk1] <- 1
  }
  else if (new.phase=="A3") {
    phases[mrk1,mrk2] <- phases[mrk2,mrk1] <- "R/C"
    recomb[mrk1,mrk2] <- recomb[mrk2,mrk1] <- w$analysis[acum(mrk1-2)+mrk2,3,1]
    arbitr[mrk1,mrk2] <- arbitr[mrk2,mrk1] <- 1
  }
  else if (new.phase=="A4") {
    phases[mrk1,mrk2] <- phases[mrk2,mrk1] <- "R/R"
    recomb[mrk1,mrk2] <- recomb[mrk2,mrk1] <- w$analysis[acum(mrk1-2)+mrk2,4,1]
    arbitr[mrk1,mrk2] <- arbitr[mrk2,mrk1] <- 1
  }
  else if (new.phase=="ns") {
    phases[mrk1,mrk2] <- phases[mrk2,mrk1] <- "ns"
    recomb[mrk1,mrk2] <- recomb[mrk2,mrk1] <- 0.5
    arbitr[mrk1,mrk2] <- arbitr[mrk2,mrk1] <- 1
  }
  else stop("unknown linkage phase")

  # results
  rf.2pts <- list(n.mar=w$n.mar, LOD=w$LOD, max.rf=w$max.rf, recomb=recomb,
                  phases=phases, analysis=w$analysis, flags=w$flags, arbitr=arbitr,
                  segr.type=w$segr.type)
  class(rf.2pts) <- "rf.2pts"
  rf.2pts
}

