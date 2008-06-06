#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: remove.marker.R                                               #
# Contains: remove.marker                                             #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 11/07/2007                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

remove.marker <-
function(w, mname=NULL) {
  # checking for correct objects
  if (any(is.na(match(c("marnames", "number", "oriname", "name"),
                      names(w))))) 
    stop("this is not an object of class 'extracted.group'")
  if (is.null(mname)) stop("you must specify a marker to be removed")
  
  if (is.character(mname)) {
    # converting "character" to "numeric"
    mrk <- which(w$marnames==mname)
    
    # checking if marker exists
    if (length(mrk)==0)
      stop("marker ", mname, " not found in this linkage group")

    # returns the same extracted group, with marker 'mrk' removed
    extracted.group <- list(marnames=w$marnames[-mrk], number=w$number,
                            oriname=w$oriname, name=w$name)
    class(extracted.group) <- "extracted.group"
    cat("  Marker", mname, "removed\n")
    extracted.group
  }
  else stop("'mname' must be of type ", dQuote("character"))
}

