#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: marker.type.R                                                     #
# Contains: marker.type                                                   #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007-9, Gabriel R A Margarido                         #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

marker.type <-
function(w) {
  # checking for correct objects
  if(!any(class(w)=="sequence")) stop(deparse(substitute(w))," is not an object of class 'sequence'")

  for(i in 1:length(w$seq.num))
    # printing marker type
    cat("  Marker", w$seq.num[i], "(", get(w$twopt)$marnames[w$seq.num[i]], ") is of type", get(w$data.name)$segr.type[w$seq.num[i]], "\n")
}

# end of file
