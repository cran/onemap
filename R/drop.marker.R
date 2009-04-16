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

# This function remove markers from a sequence
drop.marker<-function(w, mrks)
  {
    if (!any(class(w) == "sequence")) 
      stop(deparse(substitute(w)), " is not an object of class 'sequence'")
    seq.num<-w$seq.num;count<-numeric()
    for(i in 1:length(mrks)){
      if(all(seq.num!=mrks[i])) count<-c(count, i)
      seq.num<-seq.num[seq.num!=mrks[i]]
    }
    if(length(count)> 0){
      mrklist<-paste(mrks[count], collapse = ", ")
      msg <- sprintf(ngettext(length(count),
                              "marker %s was not in the sequence",
                              "markers %s were not in the sequence"), mrklist)
      warning(msg, domain=NA)
      mrks<-mrks[-count]
    }
    return(make.seq(get(w$twopt),seq.num,twopt=w$twopt))
  }
