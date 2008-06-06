#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: make.map.R                                                    #
# Contains: make.map                                                  #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 11/07/2007                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

make.map <-
function(w, ordered.group=NULL) {
  # checking for correct objects
  if (any(is.na(match(c("n.mar", "LOD", "max.rf", "recomb", "phases",
                        "analysis", "flags", "arbitr", "segr.type"),
                      names(w))))) 
        stop("this is not an object of class 'rf.2pts'")
  if (is.null(ordered.group))
    stop("you must specify the ordered linkage group")
  
  # starting to make map
  if (is.vector(ordered.group) && is.character(ordered.group)) {
    if (length(ordered.group)<2)
      stop("there must be at least two markers in the linkage group")

    # getting markers names
    markers <- match(ordered.group,colnames(w$recomb))

    # checking if markers really exist
    if (any(is.na(markers)))
      stop("marker ", ordered.group[which(is.na(markers))][1],
           " not found in object 'w'")

    # checking for repeated markers
    if (length(unique(ordered.group)) != length(ordered.group))
      stop("the map given contains repeated markers")

    # making map
    order <- 1:length(markers)
    recomb <- matrix(NA,length(markers),length(markers))
    phases <- matrix(NA,length(markers),length(markers))
    for (i in 1:(length(markers)-1))
      for (j in (i+1):length(markers)) {
        recomb[i,j] <- recomb[j,i] <- w$recomb[markers[i],markers[j]]
        phases[i,j] <- phases[j,i] <- w$phases[markers[i],markers[j]]
      }

    # results
    final.map <- list(order=order, recomb=recomb, marnames=ordered.group,
                      number=-1, name="", LOD=-1, max.rf=-1,
                      phases=phases)
    class(final.map) <- "map"
    final.map
  }
  else stop("'ordered.group' must be a vector of marker names (type ",
            dQuote("character"), ")")
}

