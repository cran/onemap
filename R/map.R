#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: map.R                                                         #
# Contains: map, print.map                                            #
#                                                                     #
# Written by Gabriel Rodrigues Alves Margarido                        #
# copyright (c) 2007, Gabriel R A Margarido                           #
#                                                                     #
# First version: 11/07/2007                                           #
# Last update: 11/07/2007                                             #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

map <-
function(w, y, LOD=NULL, max.rf=NULL) {
  # checking for correct objects
  if (any(is.na(match(c("marnames", "number", "oriname", "name"),
                      names(w))))) 
    stop("'w' is not an object of class 'extracted.group'")
  if (any(is.na(match(c("n.mar", "LOD", "max.rf", "recomb", "phases",
                        "analysis", "flags", "arbitr", "segr.type"),
                      names(y))))) 
    stop("'y' is not an object of class 'rf.2pts'")
  if (w$oriname != as.character(sys.call())[3])
    stop("this group ('w' argument) was not originally generated from object \"",
         as.character(sys.call())[3], "\"")
  
  # getting marker names and checking if they exist
  markers <- match(w$marnames,colnames(y$recomb))
  if (any(is.na(markers)))
    stop("marker ", w$marnames[which(is.na(markers))][1],
         " not found in object 'y'")
  
  # update two-point analysis with new criteria
  recomb <- matrix(NA,length(markers),length(markers))
  phases <- matrix(NA,length(markers),length(markers))
  if (is.null(LOD) && is.null(max.rf)) {
    LOD <- y$LOD
    max.rf <- y$max.rf
    # extracting linkage phase and recombination fraction values
    for (i in 1:(length(markers)-1))
      for (j in (i+1):length(markers)) {
        recomb[i,j] <- recomb[j,i] <- y$recomb[markers[i],markers[j]]
        phases[i,j] <- phases[j,i] <- y$phases[markers[i],markers[j]]
      }
  }
  else {
    if(is.null(LOD)) LOD <- y$LOD
    if(is.null(max.rf)) max.rf <- y$max.rf
    n.mar <- length(markers)
    goodness <- character(4)
    
    for (i in 2:n.mar) {
      for (j in 1:(i-1)) {
        phase <- NA
        
        # getting previous results
        current <- y$analysis[acum(markers[i]-2)+markers[j],,]
        
        # choosing most probable assignments
        probab <- which(current[,4]>(max(current[,4]-0.005)) &
                        current[,4]<=max(current[,4]))
        for (a in probab) {
          if (current[a,1] <= max.rf) goodness[a] <- "**"
          else goodness[a] <- "-"
        }
        goodness[-probab] <- "-"
        phase <- which(goodness=="**")
        
        # no assignment meets 'max.rf' criterion
        if (length(phase)==0) {
          phases[i,j] <- phases[j,i] <- "ns"
          recomb[i,j] <- recomb[j,i] <- 0.5
        }
        else {  # one or more assignments meet 'max.rf' criterion
          
          # if more than one probable phase, choose the first one
          if (length(phase)>1) phase <- phase[1]
          
          # checking if phase meets 'LOD' criterion
          if (current[phase,4] >= LOD) {
            recomb[i,j] <- recomb[j,i] <- current[phase,1]
            if (phase==1) phases[i,j] <- phases[j,i] <- "C/C"
            else if (phase==2) phases[i,j] <- phases[j,i] <- "C/R"
            else if (phase==3) phases[i,j] <- phases[j,i] <- "R/C"
            else if (phase==4) phases[i,j] <- phases[j,i] <- "R/R"
          }
          else {  # 'LOD' criterion not met
            phases[i,j] <- phases[j,i] <- "ns"
            recomb[i,j] <- recomb[j,i] <- 0.5
          }
        }
      }
    }
  }
  
  # ordering markers in linkage group with Rapid Chain Delineation
  order <- rcd(recomb)

  # results
  final.map <- list(order=order, recomb=recomb, marnames=w$marnames,
                    number=w$number, name=w$name, LOD=LOD,
                    max.rf=max.rf, phases=phases)
  class(final.map) <- "map"
  final.map
}



# print method for object class 'map'
print.map <-
function(x,cumulative=FALSE,...) {
  # checking for correct object
  if (any(is.na(match(c("order", "recomb", "marnames", "number",
                        "name", "LOD", "max.rf", "phases"),
                      names(x))))) 
        stop("this is not an object of class 'map'")
  
  if (x$number==-1) {  # arbitrary map
    cat("  This is an arbitrary map given by the user.\n\n")
  }
  else {  # "regular" map
    cat(paste("  Group ", x$number, " extracted from object \"",
              x$name, "\"\n",sep=""))
    # criteria
    cat("  Criteria used to order markers in this group:\n")
    cat("    LOD =", x$LOD, ", Maximum recombination fraction =",
        x$max.rf, "\n\n")
  }

  # starting to print map
  cat("  Linkage map:\n\n")

  dist <- numeric(length(x$marnames))
  phases <- character(length(x$marnames))
  # obtaining values
  for (i in 1:(length(dist)-1)) {
    dist[i] <- x$recomb[x$order[i],x$order[i+1]]
    phases[i] <- x$phases[x$order[i],x$order[i+1]]
  }
  
  phases[which(phases=="C/C")]<-"coupling/coupling"
  phases[which(phases=="C/R")]<-"coupling/repulsion"
  phases[which(phases=="R/C")]<-"repulsion/coupling"
  phases[which(phases=="R/R")]<-"repulsion/repulsion"
  phases[which(phases=="ns")]<-"non-significant"
  if (cumulative==TRUE) {
    # cumulative distances (required by other softwares)
    dist <- head(dist,n=-1)
    results <- cbind(x$marnames[x$order],
                     round(c(0,cumsum(kosambi(dist))),3),
                     round(c(0,cumsum(haldane(dist))),3),phases)
  } else {
    # interval distances
    results <- cbind(x$marnames[x$order],
                     round(kosambi(dist),3),
                     round(haldane(dist),3),phases)
    results[length(dist),2:3] <- ""
  }
  
  # displaying results
  results <- as.data.frame(results)
  colnames(results) <- c("Markers", "cM Kosambi", "cM Haldane",
                         "Linkage Phases")
  print(results)
}

