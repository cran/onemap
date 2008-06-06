#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: mrktype.R                                                     #
# Contains: mrktype                                                   #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 11/07/2007                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

mrktype <-
function(w, mrkname=NULL) {
  # checking for correct objects
  if (any(is.na(match(c("n.mar", "LOD", "max.rf", "recomb", "phases",
                        "analysis", "flags", "arbitr", "segr.type"),
                      names(w))))) 
    stop("this is not an object of class 'rf.2pts'")
  if (is.null(mrkname)) stop("you must specify a marker")

  if (is.character(mrkname)) {
    # converting "character" to "numeric"
    mrk <- which(colnames(w$recomb)==mrkname)
    
    # checking if marker really exists
    if (length(mrk)==0) stop("marker ", mrkname, " not found")

    # printing marker type
    cat("  Marker", mrkname, "is of type", w$segr.type[mrk], "\n")
  }
  else stop("'mrkname' must be of type ", dQuote("character"))
}

