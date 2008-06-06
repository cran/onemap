#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: rcd.R                                                         #
# Contains: rcd                                                       #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 11/07/2007                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

rcd <-
function(r) {
  # checking for correct object
  if (!is.matrix(r))
    stop ("'r' must be a matrix")
  if (dim(r)[1]!=dim(r)[2])
    stop ("'r' must be a square matrix")
  if (!all(is.na(diag(r))))
    stop ("the diagonal of the matrix must be filled with NA's")
  
  markers <- dim(r)[1]
  
  # x defines the non-positioned markers
  x <- 1:markers

  # the group starts with the closest two markers
  i <- which(r==r[which.min(r)])
  i <- i[sample(length(i),1)]

  # 'first' and 'last' are the markers on the edges of the ordered group
  first <- ceiling(i/markers)
  last <- i-((first-1)*markers)
  
  # the two markers are set next to each other  
  order <- c(first,last)

  # markers already ordered are marked as NaN
  x[first] <- NaN
  x[last] <- NaN

  # extending the chain
  while (length(order) < markers) {
    # get the markers closest to the left end of the group
    j <- which(r[first,][x]==r[first,which.min(r[first,][x])])
    # randomly choose one of them
    if (length(j) > 1) j <- j[sample(length(j),1)]

    # get the markers closest to the right end of the group
    k <- which(r[last,][x]==r[last,which.min(r[last,][x])])
    # randomly choose one of them
    if (length(k) > 1) k <- k[sample(length(k),1)]

    if (r[first,j] < r[last,k]) {  # place a marker at the left side
      order <- c(j,order)
      x[j] <- NaN
      first <- j
    }
    else if (r[first,j] > r[last,k]) {  # or at the right side
      order <- c(order,k)
      x[k] <- NaN
      last <- k } 
    else {  # if the distance is the same, randomly choose one side
      rand <- sample(2,1) 
      if (rand==1) { order <- c(j,order); x[j] <- NaN; first <- j }
      else { order <- c(order,k); x[k] <- NaN; last <- k }
    }
  }
  # end of chain

  names(order) <- NULL

  # return 'order' or 'reversed order' according to the following
  # criterion: markers on the start of the input file are placed on
  # "top" of the group 
  corr1 <- cor.test(order,1:markers,method="spearman",
                    alternative="two.sided")$estimate
  corr2 <- cor.test(rev(order),1:markers,method="spearman",
                    alternative="two.sided")$estimate
  if (corr1 > corr2) order
  else if (corr1 < corr2) rev(order)
  else {
    if (which(order==1) < which(rev(order)==1)) order
    else if (which(order==1) > which(rev(order)==1)) rev(order)
    else {
      if (which(order==2) < which(rev(order)==2)) order
      else rev(order)
    }
  }
}

