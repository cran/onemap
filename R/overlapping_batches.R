#######################################################################
##                                                                     ##
## Package: Onemap                                                     ##
##                                                                     ##
## File: overlapping_batches.R                                         ##
## Contains: map_overlapping_batches, generate_overlapping_batches,    ##
##          pick_batch_sizes                                           ##
##                                                                     ##
## Adaptated by Cristiane Taniguti from the original in                ##
## BatchMap package                                                    ##
## copyright (c) 2019 Cristiane Taniguti                               ##
##                                                                     ##
##                                                                     ##
## First version: 13/11/2019                                           ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################

##' Function to divide the sequence in batches with user defined size
##' 
##' @param input.seq an object of class \code{sequence}.
##' @param size The center size around which an optimum is to be searched
##' @param overlap The desired overlap between batches
##' 
generate_overlapping_batches <- function(input.seq, size = 50, overlap = 15)
{
  start <- 1
  end <- size
  current <- 1
  res <- list()
  while(end <= length(input.seq$seq.num))
  {
    res[[current]] <- input.seq$seq.num[start:end]
    current <- current + 1
    start <- end - overlap
    if(end == length(input.seq$seq.num)) break
    end <- end + size - overlap - 1
    if(end > length(input.seq$seq.num)) end <- length(input.seq$seq.num)
  }
  sizes <- unlist(lapply(res, length))
  if(length(sizes) < 2)
  {
    warning("You should at least have two overlapping batches.",
            " Reconsider the size parameter.")
  }
  if(any(sizes/size > 1.25))
  {
    warning("One group is 25% bigger than the group size. ",
            "Consider adjusting parameters.")
  }
  return(res)
}

##' Picking optimal batch size values
##'
##' Suggest an optimal batch size value for use in
##' \code{\link[onemap]{map_overlapping_batches}}
##'
##' @param input.seq an object of class \code{sequence}.
##' @param size The center size around which an optimum is to be searched
##' @param overlap The desired overlap between batches
##' @param around The range around the center which is maximally allowed
##' to be searched.
##' @return An integer value for the size which most evenly divides batches. In
##' case of ties, bigger batch sizes are preferred.
##' @author Bastian Schiffthaler, \email{bastian.schiffthaler@umu.se}
##' @seealso \code{\link[onemap]{map_overlapping_batches}}
##'
##' @keywords utilities
##' @examples
##'
##'   LG <- structure(list(seq.num = seq(1,800)), class = "sequence")
##'   batchsize <- pick_batch_sizes(LG, 50, 19)
##' 
##' @export
pick_batch_sizes <- function(input.seq, size = 50, overlap = 15, around = 5)
{
  test.sizes <- c(size, (size - 1):(size - around), (size + 1):(size + around))
  rm.num <- which(test.sizes <= (overlap+1))
  if(length(rm.num) > 0) test.sizes <- test.sizes[-rm.num]
  all.batches <- lapply(test.sizes, function(s){
    generate_overlapping_batches(input.seq, s, overlap)
  })
  x <- unlist(lapply(all.batches, function(f){
    ran <- range(unlist(lapply(f,length)))
    ran[2] - ran[1]
  }))
  x <- which(x == min(x))
  size <- test.sizes[x[length(x)]] #prefer larger maps
  if(length(size) > 0) return(size) else stop("You should at least have two overlapping batches.",
                                              " Reconsider the size parameter.")
}

##' Mapping overlapping batches
##'
##' Apply the batch mapping algorithm using overlapping windows.
##'
##' This algorithm implements the overlapping batch maps for high density
##' marker sets. The mapping problem is reduced to a number of subsets (batches)
##' which carry information forward in order to more accurately estimate
##' recombination fractions and phasing. It is a adapted version of
##' map.overlapping.batches function of BatchMap package. The main differences are
##' that this onemap version do not have the option to reorder the markers 
##' according to ripple algorithm and, if the it finds markers that do not reach the linkage
##' criterias, the algorithm remove the problematic marker and repeat the analysis.
##' Than, the output map can have few markers compared with the input.seq.
##'
##' @param input.seq an object of class \code{sequence}.
##' @param size The center size around which an optimum is to be searched
##' @param overlap The desired overlap between batches
##' @param phase_cores The number of parallel processes to use when estimating
##' the phase of a marker. (Should be no more than 4)
##' @param verbose A logical, if TRUE its output progress status
##' information.
##' @param seeds A vector of phase information used as seeds for the first
##' batch
##' @param rm_unlinked When some pair of markers do not follow the linkage criteria, 
##' if \code{TRUE} one of the markers is removed and map is performed again.
##' @param tol tolerance for the C routine, i.e., the value used to evaluate
##' convergence.
##' @param parallelization.type one of the supported cluster types. This should 
#' be either PSOCK (default) or FORK.
##' @param max.gap the marker will be removed if it have gaps higher than this defined threshold in both sides
##' 
##' @return An object of class \code{sequence}, which is a list containing the
##' following components: \item{seq.num}{a \code{vector} containing the
##' (ordered) indices of markers in the sequence, according to the input file.}
##' \item{seq.phases}{a \code{vector} with the linkage phases between markers
##' in the sequence, in corresponding positions. \code{-1} means that there are
##' no defined linkage phases.} \item{seq.rf}{a \code{vector} with the
##' recombination frequencies between markers in the sequence. \code{-1} means
##' that there are no estimated recombination frequencies.}
##' \item{seq.like}{log-likelihood of the corresponding linkage map.}
##' \item{data.name}{name of the object of class \code{outcross} with the raw
##' data.} \item{twopt}{name of the object of class \code{rf.2pts} with the
##' 2-point analyses.} 
##' @seealso \code{\link[onemap]{pick_batch_sizes}}, \code{\link[onemap]{map}}
##'
##' @keywords utilities
##' @export
map_overlapping_batches <- function(input.seq, size = 50, overlap = 15,
                                    phase_cores = 1, verbose = FALSE, 
                                    seeds = NULL, tol=10E-5, rm_unlinked = TRUE, 
                                    max.gap=FALSE, parallelization.type = "PSOCK")
{
  #TODO: error checks...
  #Create initial set of batches
  batches <- generate_overlapping_batches(input.seq, size, overlap)
  if(verbose)
  {
    message("Have ", length(batches), " batches.")
    message("The number of markers in the final batch is: ",
            length(batches[[length(batches)]]))
    message("Processing batch 1...")
  }
  LGs <- list()
  #The first batch is run in full again to get all necessary data (phases etc.)
  if(is.null(seeds))
  {
    LG <- map(input.seq = make_seq(input.seq$twopt, batches[[1]]), 
              phase_cores = phase_cores, 
              rm_unlinked = rm_unlinked, 
              tol=tol, parallelization.type = parallelization.type)
    if(inherits(LG,"integer")){
      rmed.mks <- batches[[1]][which(!batches[[1]] %in% LG)]
      warning(cat("Markers", rmed.mks,"did not reached the OneMap default criteria. They are probably segregating independently 
                  and new vector of marker numbers was generated without them. Use function map_avoid_unlinked to remove these markers automatically.\n"))
      new.seq <- input.seq$seq.num[-which(input.seq$seq.num == rmed.mks)]
      return(new.seq)
    }
  } else {
    LG <- seeded_map(input.seq = make_seq(input.seq$twopt, batches[[1]],
                                          twopt = input.seq$twopt), 
                     phase_cores = phase_cores,
                     verbose = verbose, seeds = seeds, tol=tol,
                     parallelization.type = parallelization.type)
    if(inherits(LG,"integer")){
      rmed.mks <- batches[[1]][which(!batches[[1]] %in% LG)]
      warning(cat("Markers", rmed.mks,"did not reached the OneMap default criteria. They are probably segregating independently 
                  and new vector of marker numbers was generated without them. Use function map_avoid_unlinked to remove these markers automatically.\n"))
      new.seq <- input.seq$seq.num[-which(input.seq$seq.num == rmed.mks)]
      return(new.seq)
    }
  }
  
  round <- 1
  increment <- 0
  
  LGs[[1]] <- LG
  #Start processing all following batches
  for(i in 2:length(batches))
  {
    if(verbose) #Print previous batch-map segment
    {
      print(LGs[[i - 1]])
      message("Processing batch ",i,"...")
    }
    #Need to use a seeded map in order to not mess with the overlapping area
    #which we trust more from the previous batch (as that had more information)
    seeds <- tail(LGs[[i - 1]]$seq.phases, overlap)
    batches[[i]][1:(overlap+1)] <- tail(LGs[[i - 1]]$seq.num, overlap + 1)
    LG <- seeded_map(input.seq = make_seq(input.seq$twopt,
                                          batches[[i]],
                                          twopt = input.seq$twopt),
                     verbose = verbose,
                     seeds = seeds, rm_unlinked = rm_unlinked,
                     tol=tol, parallelization.type = parallelization.type)
    if(inherits(LG,"integer")){
      rmed.mks <- batches[[1]][which(!batches[[1]] %in% LG)]
      warning(cat("Markers", rmed.mks,"did not reached the OneMap default criteria. They are probably segregating independently 
                  and new vector of marker numbers was generated without them. Use function map_avoid_unlinked to remove these markers automatically.\n"))
      new.seq <- input.seq$seq.num[-which(input.seq$seq.num == rmed.mks)]
      return(new.seq)
    }
    LGs[[i]] <- LG
  }
  #Initialize final order, phases and rfs with the first batch
  final.seq <- LGs[[1]]$seq.num
  final.phase <- LGs[[1]]$seq.phases
  final.rf <- LGs[[1]]$seq.rf
  #Iteratively add data from other batches to the first sequence
  for(i in 2:length(batches))
  {
    start <- length(final.seq) - overlap #start position of overlap with next
    #Add marker order from the next batch to the sequence starting from the
    #start of the overlap
    final.seq[start:length(final.seq)] <- head(LGs[[i]]$seq.num, overlap + 1)
    final.seq <- c(final.seq,
                   LGs[[i]]$seq.num[(overlap + 2):length(LGs[[i]]$seq.num)])
    #Add phases and RFs. We need to shift the indices left by 1
    start <- length(final.phase) - overlap + 1
    final.phase[start:length(final.phase)] <- head(LGs[[i]]$seq.phases, overlap)
    final.phase <- c(final.phase,
                     LGs[[i]]$seq.phases[(overlap + 1):length(LGs[[i]]$seq.phases)])
    final.rf[start:length(final.rf)] <- head(LGs[[i]]$seq.rf, overlap)
    final.rf <- c(final.rf,
                  LGs[[i]]$seq.rf[(overlap + 1):length(LGs[[i]]$seq.rf)])
    
  }
  if(verbose)
  {
    message("Final call to map...")
  }
  #Create final sequence and run
  #final.rf is currently only used for debugging purposes
  s <- make_seq(input.seq$twopt, final.seq, final.phase, input.seq$twopt)
  mp <- map(input.seq = s, rm_unlinked = rm_unlinked, tol=tol, parallelization.type = parallelization.type)
  
  if(max.gap){ # the marker will be removed if it have gaps higher than the threshold in both sides
    idx <- which(get(get(".map.fun", envir=.onemapEnv))(mp$seq.rf) > max.gap)
    rm.seq <- vector()
    for(i in 1:length(idx)){
      if(idx[i] == 1){
        rm.seq <- c(rm.seq, 1)
      } else if(idx[i] == length(mp$seq.rf) -1){
        rm.seq <- c(rm.seq, idx[i] + 1)
        if(idx[i-1] == idx[i] -1)
          rm.seq <- c(rm.seq, idx[i])
      } else if(i-1 != 0){
        if(idx[i-1] == idx[i] -1) {
          rm.seq <- c(rm.seq, idx[i])
        }
      }
    }
    
    if(length(rm.seq)> 0){
      new.seq <- make_seq(mp$twopt, mp$seq.num[-rm.seq])
      # Make exception if number of remaining markers is smaller than the batch size
      if(verbose) cat("Markers", mp$seq.num[rm.seq], "were remove because they cause gaps higher than ", max.gap, " cM with both neighboors markers.\n")
      mp <- map_overlapping_batches(input.seq = new.seq,
                                    size = size, overlap = overlap,
                                    phase_cores = phase_cores, verbose = verbose , 
                                    seeds = seeds, tol=tol, rm_unlinked = rm_unlinked, max.gap=max.gap,
                                    parallelization.type = parallelization.type)
    }
  }
  return(mp)
}
