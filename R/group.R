#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: group.R                                                       #
# Contains: group, print.group                                        #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 11/07/2007                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

group <-
function(w) {
  # checking for correct object
  if (any(is.na(match(c("n.mar", "LOD", "max.rf", "recomb", "phases",
                        "analysis", "flags", "arbitr", "segr.type"),
                      names(w))))) 
    stop("this is not an object of class 'rf.2pts'")
  
  # 'groups' indicates the linkage group to which each marker is associated
  groups <- rep(NA,w$n.mar)
  
  # 'group.count' is the current linkage group
  group.count <- 1
  
  for (m in 1:w$n.mar) {
    if (is.na(groups[m])) {  # marker 'm' is not associated to any group
      
      # check which markers are linked with 'm'
      grouping <- c(m,which(w$phases[m,]=="C/C" | w$phases[m,]=="C/R" |
                            w$phases[m,]=="R/C" | w$phases[m,]=="R/R"))
      
      if (length(grouping) > 1) {  # 'm' is linked with at least one marker
        flag <- 1  # 'flag' indicates if a new marker has been added to the group
        while (flag) {
          flag <- 0
          
          for (i in grouping) {
            # detect all markers linked to those already in group
            group_parc <- which(w$phases[i,]=="C/C" | w$phases[i,]=="C/R" |
                                w$phases[i,]=="R/C" | w$phases[i,]=="R/R")
            
            # check if markers in 'group_parc' are already in the current group
            for (j in group_parc)
              if (all(grouping != j)) {
                grouping <- c(grouping,j)
                flag <- 1
              }
          }
        }
        
        # finishing the current linkage group
        groups[grouping] <- group.count
        group.count <- group.count + 1
      }
    }
  }
  
  # results
  marnames <- colnames(w$recomb)
  final.groups <- list(marnames=marnames, n.mar=w$n.mar, LOD=w$LOD,
                       max.rf=w$max.rf, n.groups=max(groups,na.rm=TRUE),
                       groups=groups, name=as.character(sys.call())[2])
  class(final.groups) <- "group"
  final.groups
}



# print method for object class 'group'
print.group <-
function(x, detailed=TRUE,...) {
  # checking for correct object
  if (any(is.na(match(c("marnames", "n.mar", "LOD", "max.rf",
                        "n.groups", "groups", "name"),
                      names(x))))) 
    stop("this is not an object of class 'group'")
  
  cat("  This is an object of class 'group'\n")
  cat(paste("  It was generated from the object \"", x$name,
            "\"\n\n",sep=""))
  
  # criteria
  cat("  Criteria used to assign markers to groups:\n")
  cat("    LOD =", x$LOD, ", Maximum recombination fraction =",
      x$max.rf, "\n")

  # printing summary
  cat("\n  No. markers:           ", x$n.mar, "\n")
  cat("  No. groups:            ", x$n.groups, "\n")
  cat("  No. linked markers:    ", sum(!is.na(x$groups)), "\n")
  cat("  No. unlinked markers:  ", sum(is.na(x$groups)), "\n")
  
  if (detailed) {
    # printing detailed results (markers in each linkage group)
    cat("\n  Printing groups:")
    for (i in 1:x$n.groups) {
      cat("\n  Group", i, ":\n    ")
      cat(x$marnames[which(x$groups==i)], "\n")
    }
    if (any(is.na(x$groups))) {
      cat("\n  Unlinked markers:\n    ")
      cat(x$marnames[which(is.na(x$groups))], "\n")
    }
  }
}


