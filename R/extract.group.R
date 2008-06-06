#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: extract.group.R                                               #
# Contains: extract.group, print.extracted.group                      #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 11/07/2007                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

extract.group <-
function(w, i=NULL) {
  # checking for correct objects
  if (any(is.na(match(c("marnames", "n.mar", "LOD", "max.rf", "n.groups",
                        "groups", "name"), names(w))))) 
        stop("this is not an object of class 'group'")
  if (is.null(i))
    stop("you must specify a linkage group")

  # checking if linkage group 'i' exists
  if (i > w$n.groups) stop("linkage group ", i, " not found")

  # "extracting" group
  group <- which(w$groups==i)

  # results
  extracted.group <- list(marnames=w$marnames[group], number=i,
                          oriname=w$name, name=as.character(sys.call())[2])
  class(extracted.group) <- "extracted.group"
  extracted.group
}



# print method for object class 'extracted.group'
print.extracted.group <-
function(x,...) {
  # checking for correct object
  if (any(is.na(match(c("marnames", "number", "oriname", "name"),
                      names(x))))) 
    stop("this is not an object of class 'extracted.group'")
  cat(paste("  Group ", x$number, " extracted from object \"",
            x$name, "\"\n    ",sep=""))
  cat(x$marnames, "\n")
}

