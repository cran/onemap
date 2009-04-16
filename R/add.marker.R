#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: codif.data.R                                                  #
# Contains: codif.data                                                #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2009, Marcelo Mollinari                               #
#                                                                     #
# First version: 02/27/2009                                           #
# Last update: 02/27/2009                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function add markers to a sequence
add.marker<-function(w, mrks)
  {
    if (!any(class(w) == "sequence")) 
      stop(sQuote(deparse(substitute(w))), " is not an object of class 'sequence'")
    seq.num<-c(w$seq.num,mrks)
    return(make.seq(get(w$twopt),seq.num, twopt=w$twopt))
  }
