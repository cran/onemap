#######################################################################
#                                                                     #
# Package: onemap                                                     #
#                                                                     #
# File: onemap_read_vcfR.R                                            #
# Contains: onemap_read_vcfR write_onemap_raw                         #
#                                                                     #
# Written by Cristiane Hayumi Taniguti                                #
#                                                                     #                                            #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################


##' Convert vcf file to onemap object
##'
##' Converts data from a vcf file to onemap initial object, while identify 
##' the appropriate marker segregation patterns.
##'
##' Only biallelic SNPs and indels for diploid variant sites are considered.
##'
##' Genotype information on the parents is required for all cross types. For
##' full-sib progenies, both outbred parents must be genotyped. For backcrosses,
##' F2 intercrosses and recombinant inbred lines, the \emph{original inbred
##' lines} must be genotyped. Particularly for backcross progenies, the
##' \emph{recurrent line must be provided as the first parent} in the function
##' arguments.
##'
##' Marker type is determined based on parental genotypes. Variants for which parent
##' genotypes cannot be determined are discarded.
##'
##' Reference sequence ID and position for each variant site are also stored.
##'
##' @param vcf string defining the path to VCF file;
##' @param vcfR.object object of class vcfR;
##' @param cross type of cross. Must be one of: \code{"outcross"} for full-sibs;
##' \code{"f2 intercross"} for an F2 intercross progeny; \code{"f2 backcross"};
##' \code{"ri self"} for recombinant inbred lines by self-mating; or
##' \code{"ri sib"} for recombinant inbred lines by sib-mating.
##' @param parent1 \code{string} specifying sample ID of the first parent. If f2 backcross population, define here the ID of the backcrossed parent.
##' @param parent2 \code{string} specifying sample ID of the second parent.
##' @param f1 \code{string} if you are working with f2 intercross or backcross populations you may have f1 parents in you vcf, specify its ID here
##' @param only_biallelic if TRUE (default) only biallelic markers are considered, if FALSE multiallelic markers are included.
##' @param output_info_rds define a name for the file with alleles information.
##' @param verbose A logical, if TRUE it output progress status
##' information.
##' 
##' 
##' @return An object of class \code{onemap}, i.e., a list with the following
##' components: \item{geno}{a matrix with integers indicating the genotypes
##' read for each marker. Each column contains data for a marker and each row
##' represents an individual.} \item{n.ind}{number of individuals.}
##' \item{n.mar}{number of markers.} \item{segr.type}{a vector with the
##' segregation type of each marker, as \code{strings}.} \item{segr.type.num}{a
##' vector with the segregation type of each marker, represented in a
##' simplified manner as integers, i.e. 1 corresponds to markers of type
##' \code{"A"}; 2 corresponds to markers of type \code{"B1.5"}; 3 corresponds
##' to markers of type \code{"B2.6"}; 4 corresponds to markers of type
##' \code{"B3.7"}; 5 corresponds to markers of type \code{"C.8"}; 6 corresponds
##' to markers of type \code{"D1"} and 7 corresponds to markers of type
##' \code{"D2"}. Markers for F2 intercrosses are coded as 1; all other crosses
##' are left as \code{NA}.} \item{input}{the name of the input file.}
##' \item{n.phe}{number of phenotypes.} \item{pheno}{a matrix with phenotypic
##' values. Each column contains data for a trait and each row represents an
##' individual.} \item{error}{matrix containing HMM emission probabilities}
##' 
##' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
##' 
##' @seealso \code{read_onemap} for a description of the output object of class onemap.
##' 
##' @importFrom rebus number_range
##' @importFrom vcfR read.vcfR extract.gt masplit
##' 
##' @examples
##' \donttest{
##' data <- onemap_read_vcfR(vcf=system.file("extdata/vcf_example_out.vcf.gz", package = "onemap"),
##'                  cross="outcross",
##'                  parent1=c("P1"),
##'                  parent2=c("P2"))
##' }
##'                 
##'@export                  
onemap_read_vcfR <- function(vcf=NULL,
                             vcfR.object = NULL,
                             cross = NULL,
                             parent1 =NULL,
                             parent2 =NULL,
                             f1=NULL,
                             only_biallelic = TRUE,
                             output_info_rds = NULL,
                             verbose=TRUE){
  
  if (is.null(vcf) & is.null(vcfR.object)) {
    stop("You must specify one vcf file.")
  }
  if (is.null(parent1) || is.null(parent2)) {
    stop("You must specify samples as parents 1 and 2.")
  }
  
  if(is.null(vcfR.object)){
    vcfR.obj <- read.vcfR(vcf, verbose = F)
  } else vcfR.obj <- vcfR.object
  
  if(is.null(cross)) stop("Define a cross type: outcross, f2 intercross, f2 backcross, ri self, ri sib")
  
  n.mk <- dim(vcfR.obj@gt)[1]
  n.ind <- dim(vcfR.obj@gt)[2]-1
  INDS <- dimnames(vcfR.obj@gt)[[2]][-1]
  CHROM <- vcfR.obj@fix[,1]
  POS <- as.numeric(vcfR.obj@fix[,2])
  REF <- vcfR.obj@fix[,4]
  ALT <- vcfR.obj@fix[,5]
  
  if(is.vector(vcfR.obj@gt)){
    jump <- 1
  } else if(dim(vcfR.obj@gt)[1] == 0){
    jump <- 1
  } else jump <- 0
  
  if(jump == 1){
    warning("Input vcfR.objR object do not have markers. An empty object onemap will be generated.")
    onemap.obj <- empty_onemap_obj(vcfR.obj, P1, P2, cross)
    return(onemap.obj)
  }
  
  # Checking marker segregation according with parents
  P1 <- which(dimnames(vcfR.obj@gt)[[2]]==parent1) -1 
  P2 <- which(dimnames(vcfR.obj@gt)[[2]]==parent2) -1
  if(length(P1)==0 | length(P2)==0) stop("One or both parents names could not be found in your data")
  
  MKS <- vcfR.obj@fix[,3]
  if (any(MKS == "." | is.na(MKS))) {
    MKS <- paste0(vcfR.obj@fix[,1],"_", vcfR.obj@fix[,2])
    # Add tag if is duplicated positions (split form of mnps)
    if(any(duplicated(MKS))){
      z <- 1
      for(i in 2:length(MKS)) {
        if(MKS[i] == paste0(strsplit(MKS[i-1], "_")[[1]][1:2], collapse = "_")) {
          z <- z + 1
          MKS[i] <- paste0(MKS[i], "_",z)
        } 
      }
    }
  }
  
  # Geno matrix
  GT_matrix <- extract.gt(vcfR.obj)
  
  # This function do not consider phased genotypes
  GT_matrix[grep("[.]", GT_matrix)] <- "./."
  GT_matrix[is.na(GT_matrix)] <- "./."
  
  GT_names <- names(table(GT_matrix))
  
  phased <- any(grepl("[|]", GT_names))
  if(phased)
    GT_matrix <- gsub("[|]", "/", as.matrix(GT_matrix))
  
  GT_names <- names(table(GT_matrix))
  
  GT_names_up <- strsplit(GT_names, "/")
  max.alleles <- max(as.numeric(do.call(c, GT_names_up[-1])))
  
  if(phased){
    if(length(grep("[.]", GT_names_up)) > 0){
      idx.mis <- grep("[.]", GT_names_up)
      GT_names_up[[idx.mis]] <- 0 # avoiding warning
    } else idx.mis <- "nomis"
    
    GT_names_up <- sapply(GT_names_up, function(x) paste(sort(as.numeric(x)), collapse = "/"))
    
    if(idx.mis != "nomis")
      GT_names_up[idx.mis] <- "./."
    
    only_diff <- which(GT_names_up != GT_names)
    if(length(only_diff) > 0){
      repl <- GT_names_up[only_diff]
      sear <- GT_names[only_diff]
      for(i in 1:length(sear)){
        GT_matrix[which(GT_matrix == sear[i])] <- repl[i]
      }
    }
  }
  
  # keep only biallelic
  if(only_biallelic | cross != "outcross"){
    if(max.alleles > 1){
      rx <- number_range(2, max.alleles)
      rm_multi <- which(apply(GT_matrix, 1, function(x) any(grepl(rx, x))))
      if(length(rm_multi) > 0){
        GT_matrix <- GT_matrix[-rm_multi,]
        CHROM <- CHROM[-rm_multi]
        POS <- POS[-rm_multi]
        MKS <- MKS[-rm_multi]
        REF <- REF[-rm_multi]
        ALT <- ALT[-rm_multi]
      }
    }
  }
  n.mk <- nrow(GT_matrix)
  
  alleles <- strsplit(ALT, ",")
  for(i in 1:length(alleles)) alleles[[i]] <- c(REF[i],alleles[[i]])
  
  mk.type <- mk.type.num <- rep(NA, n.mk)
  if (cross == "outcross"){
    P1_1 <- sapply(strsplit(GT_matrix[,P1], "/"), "[", 1)
    P1_2 <- sapply(strsplit(GT_matrix[,P1], "/"), "[", 2)
    P2_1 <- sapply(strsplit(GT_matrix[,P2], "/"), "[", 1)
    P2_2 <- sapply(strsplit(GT_matrix[,P2], "/"), "[", 2)
    
    # avoid warning
    P1_1_t <- gsub("[.]", NA, P1_1)
    P1_2_t <- gsub("[.]", NA, P1_2)
    P2_1_t <- gsub("[.]", NA, P2_1)
    P2_2_t <- gsub("[.]", NA, P2_2)
    
    P1_1_allele <- unlist(Map("[",alleles,as.numeric(P1_1_t) + 1))
    P1_2_allele <- unlist(Map("[",alleles,as.numeric(P1_2_t) + 1))
    P2_1_allele <- unlist(Map("[",alleles,as.numeric(P2_1_t) + 1))
    P2_2_allele <- unlist(Map("[",alleles,as.numeric(P2_2_t) + 1))
    
    names(P1_1_allele) <- names(P1_2_allele) <- names(P2_1_allele) <- names(P2_2_allele) <- rownames(GT_matrix)
    
    # Marker types
    GT_parents <- cbind(P1_1, P1_2,P2_1, P2_2)
    idx <- which(P1_1 == "." | P2_1 == "." |  P1_2 == "." | P2_2 == ".")
    GT_parents[idx,] <- NA
    
    idx <- apply(GT_parents, 1, function(x) length(x) == length(unique(x)))
    mk.type[idx] <- "A.1"
    mk.type.num[idx] <- 1
    idx <- apply(GT_parents, 1, function(x) (length(x) -1) == length(unique(x)))
    idx.sub <- which(P1_1[idx] != P1_2[idx] & P2_1[idx] != P2_2[idx])
    mk.type[idx][idx.sub] <- "A.2"
    mk.type.num[idx][idx.sub] <- 1
    idx.sub <- which(P1_1[idx] != P1_2[idx] & P2_1[idx] == P2_2[idx])
    mk.type[idx][idx.sub] <- "D1.9"
    mk.type.num[idx][idx.sub] <- 6
    idx.sub <- which(P1_1[idx] == P1_2[idx] & P2_1[idx] != P2_2[idx])
    mk.type[idx][idx.sub] <- "D2.14"
    mk.type.num[idx][idx.sub] <- 7
    idx <- apply(GT_parents, 1, function(x) (length(x) -2) == length(unique(x)))
    idx.sub <- which(P1_1[idx] != P1_2[idx] & P2_1[idx] != P2_2[idx])
    mk.type[idx][idx.sub] <- "B3.7"
    mk.type.num[idx][idx.sub] <- 4
    idx.sub <- which(P1_1[idx] != P1_2[idx] & P2_1[idx] == P2_2[idx])
    mk.type[idx][idx.sub] <- "D1.10"
    mk.type.num[idx][idx.sub] <- 6
    idx.sub <- which(P1_1[idx] == P1_2[idx] & P2_1[idx] != P2_2[idx])
    mk.type[idx][idx.sub] <- "D2.15"
    mk.type.num[idx][idx.sub] <- 7
    
    # It informs to user why markers are being removed
    idx <- which(GT_matrix[,P1] == "./." | GT_matrix[,P2] == "./." |  GT_matrix[,P2] == "." | GT_matrix[,P1] == ".")
    if(verbose){
      if (length(idx) > 0)
        cat(length(MKS[idx]), "Markers were removed of the dataset because one or both of parents have no informed genotypes (are missing data)\n")
    }
    idx <- which(P1_1 == P1_2 & P2_1 == P2_2)
    if(verbose) {
      if (length(idx) > 0)
        cat( length(MKS[idx]), "Markers were removed from the dataset because both of parents are homozygotes, these markers are considered non-informative in outcrossing populations.\n")
    }
    # Excluding non-informative markers
    rm_mk <- which(is.na(mk.type))
    if(length(rm_mk)!=0){
      GT_matrix <- GT_matrix[-rm_mk,]
      P1_1 <- P1_1[-rm_mk]
      P1_2 <- P1_2[-rm_mk]
      P2_1 <- P2_1[-rm_mk]
      P2_2 <- P2_2[-rm_mk]
      MKS <- MKS[-rm_mk]
      n.mk <- n.mk - length(rm_mk)
      CHROM <-CHROM[-rm_mk]
      POS <- POS[-rm_mk]
      REF <- REF[-rm_mk]
      ALT <- ALT[-rm_mk]
      mk.type <- mk.type[-rm_mk]
      mk.type.num <- mk.type.num[-rm_mk]
    } 
    
    if(is.vector(GT_matrix)){
      jump <- 1
    } else if(dim(GT_matrix)[1]==0){
      jump <- 1
    } else jump <- 0
    
    if(jump == 1){
      warning("Input vcfR object do not have markers. An empty object onemap will be generated.")
      
      onemap.obj <- empty_onemap_obj(vcfR.obj, P1, P2, cross)
      return(onemap.obj)
    }
    GT_matrix <- as.matrix(GT_matrix)
    # Codification for OneMap
    idx <- which(mk.type=="A.1" | mk.type=="A.2")
    cat <- paste0(P1_1[idx], "/", P2_1[idx])
    cat.rev <- paste0(P2_1[idx], "/", P1_1[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 1
    cat <- paste0(P1_1[idx], "/", P2_2[idx])
    cat.rev <- paste0(P2_2[idx], "/", P1_1[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 2
    cat <- paste0(P1_2[idx], "/", P2_1[idx])
    cat.rev <- paste0(P2_1[idx], "/", P1_2[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 3
    cat <- paste0(P1_2[idx], "/", P2_2[idx])
    cat.rev <- paste0(P2_2[idx], "/", P1_2[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 4
    
    idx <- which(mk.type=="B3.7")
    cat <- paste0(P1_1[idx], "/", P2_1[idx]) 
    cat.rev <- paste0(P2_1[idx], "/", P1_1[idx]) 
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 1
    cat <- paste0(P1_1[idx], "/", P2_2[idx]) 
    cat.rev <- paste0(P2_2[idx], "/", P1_1[idx]) 
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 2
    cat <- paste0(P1_2[idx], "/", P2_2[idx]) 
    cat.rev <- paste0(P2_2[idx], "/", P1_2[idx]) 
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 3
    
    idx <- which(mk.type=="D1.10")
    idx.sub <- which(P1_1[idx] == P2_1[idx])
    cat <- paste0(P1_1[idx][idx.sub], "/", P2_1[idx][idx.sub])
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_1[idx][idx.sub]) 
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 1
    cat <- paste0(P1_2[idx][idx.sub], "/", P2_1[idx][idx.sub]) 
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_2[idx][idx.sub]) 
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 2
    
    idx.sub <- which(P1_2[idx] == P2_1[idx])
    cat <- paste0(P1_2[idx][idx.sub], "/", P2_1[idx][idx.sub])
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_2[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 1
    cat <- paste0(P1_1[idx][idx.sub], "/", P2_1[idx][idx.sub])
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_1[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 2
    
    idx <- which(mk.type=="D1.9")
    cat <- paste0(P1_1[idx], "/", P2_1[idx])
    cat.rev <- paste0(P2_1[idx], "/", P1_1[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 1
    cat <- paste0(P1_2[idx], "/", P2_1[idx])
    cat.rev <- paste0(P2_1[idx], "/", P1_2[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 2
    
    idx <- which(mk.type=="D2.15" )
    idx.sub <- which(P1_1[idx] == P2_1[idx])
    cat <- paste0(P1_1[idx][idx.sub], "/", P2_1[idx][idx.sub])
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_1[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 1
    cat <- paste0(P1_1[idx][idx.sub], "/", P2_2[idx][idx.sub])
    cat.rev <- paste0(P2_2[idx][idx.sub], "/", P1_2[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 2
    
    idx.sub <- which(P1_1[idx] == P2_2[idx])
    cat <- paste0(P1_2[idx][idx.sub], "/", P2_2[idx][idx.sub])
    cat.rev <- paste0(P2_2[idx][idx.sub], "/", P1_2[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 1
    cat <- paste0(P1_1[idx][idx.sub], "/", P2_1[idx][idx.sub])
    cat.rev <- paste0(P2_1[idx][idx.sub], "/", P1_1[idx][idx.sub])
    GT_matrix[idx[idx.sub],][which(GT_matrix[idx[idx.sub],] == cat | GT_matrix[idx[idx.sub],] == cat.rev)] <- 2
    
    idx <- which(mk.type=="D2.14")
    cat <- paste0(P1_1[idx], "/", P2_1[idx])
    cat.rev <- paste0(P2_1[idx], "/", P1_1[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 1
    cat <- paste0(P1_1[idx], "/", P2_2[idx])
    cat.rev <- paste0(P2_2[idx], "/", P1_1[idx])
    GT_matrix[idx,][which(GT_matrix[idx,] == cat | GT_matrix[idx,] == cat.rev)] <- 2
    
    GT_matrix[grepl("/", GT_matrix)] <- 0
    GT_matrix[grepl("[.]", GT_matrix)] <- 0
    GT_parents <- cbind(P1_1, P1_2, P2_1, P2_2)
  } else if(cross== "f2 intercross"){
    # Marker type
    mk.type[which(GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "1/1")] <- "A.H.B.1"
    mk.type[which(GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "0/0")] <- "A.H.B.2"
    GT_parents <- GT_matrix[,c(P1,P2)]
    
    P1_1_allele <- GT_parents[,1]
    P1_1_allele[which(GT_parents[,1] == "0/0")] <- REF[which(GT_parents[,1] == "0/0")]
    P1_1_allele[which(GT_parents[,1] == "1/1")] <- ALT[which(GT_parents[,1] == "1/1")]
    P2_1_allele <- P1_1_allele
    
    P2_2_allele <- GT_parents[,2]
    P2_2_allele[which(GT_parents[,2] == "0/0")] <- REF[which(GT_parents[,2] == "0/0")]
    P2_2_allele[which(GT_parents[,2] == "1/1")] <- ALT[which(GT_parents[,2] == "1/1")]
    P1_2_allele <- P2_2_allele
    
    # Informs to user why markers are being removed
    idx <- which(GT_matrix[,P1] == "./." | GT_matrix[,P2] == "./.")
    if(verbose){
      if (length(idx) > 0)
        cat(length(MKS[idx]), "Markers were removed of the dataset because one or both of parents have no informed genotypes (are missing data)\n")
      
      idx <- which((GT_matrix[,P1] == "0/1" | GT_matrix[,P2] == "0/1"))
      if (length(idx) > 0)
        cat(length(MKS[idx]), "Markers were removed from the dataset because one or both of the parents are heterozygotes, we do not expect heterozygotes parents in F2 populations.\n") 
      
      idx <- which((GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "0/0") | (GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "1/1"))
      if (length(idx) > 0)
        cat(length(MKS[idx]), "Markers were removed from the dataset because they are monomorphic for the parents, these markers are not informative for the genetic map.\n") 
    }
    # Excluding non-informative markers
    rm_mk <- which(is.na(mk.type))
    if(length(rm_mk)!=0){
      GT_matrix <- GT_matrix[-rm_mk,]
      MKS <- MKS[-rm_mk]
      n.mk <- n.mk - length(rm_mk)
      CHROM <- CHROM[-rm_mk]
      POS <- POS[-rm_mk]
      REF <- REF[-rm_mk]
      ALT <- ALT[-rm_mk]
      mk.type <- mk.type[-rm_mk]
      mk.type.num <- mk.type.num[-rm_mk]
      GT_parents <- GT_parents[-rm_mk,]
    } 
    
    if(is.vector(GT_matrix)){
      jump <- 1
    } else if(dim(GT_matrix)[1]==0){
      jump <- 1
    } else jump <- 0
    
    if(jump == 1){
      warning("Input vcfR object do not have markers. An empty object onemap will be generated.")
      
      onemap.obj <- empty_onemap_obj(vcfR.obj, P1, P2, cross)
      return(onemap.obj)
    }
    
    # Codification for OneMap
    GT_matrix[-which(GT_matrix == "1/1" | GT_matrix == "0/0" | GT_matrix == "0/1")] <- 0
    GT_matrix[which(GT_matrix == "0/1")] <- 2
    
    idx <- which(mk.type=="A.H.B.1")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 1
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 3
    
    idx <- which(mk.type=="A.H.B.2")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 3
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 1
    
    mk.type <- mk.type.num <- rep("A.H.B", n.mk)
    mk.type.num[mk.type=="A.H.B"] <- 4
    
  } else if(cross=="f2 backcross"){
    mk.type[which(GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "1/1")] <- "A.H.1"
    mk.type[which(GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "0/0")] <- "A.H.2"
    GT_parents <- GT_matrix[,c(P1,P2)]
    
    P1_1_allele <- GT_parents[,1]
    P1_1_allele[which(GT_parents[,1] == "0/0")] <- REF[which(GT_parents[,1] == "0/0")]
    P1_1_allele[which(GT_parents[,1] == "1/1")] <- ALT[which(GT_parents[,1] == "1/1")]
    P2_1_allele <- P1_1_allele
    
    P2_2_allele <- GT_parents[,2]
    P2_2_allele[which(GT_parents[,2] == "0/0")] <- REF[which(GT_parents[,2] == "0/0")]
    P2_2_allele[which(GT_parents[,2] == "1/1")] <- ALT[which(GT_parents[,2] == "1/1")]
    P1_2_allele <- P2_2_allele
    
    # Informs to user why markers are being removed
    if(verbose) {
      idx <- which(GT_matrix[,P1] == "./." | GT_matrix[,P2] == "./.")
      if (length(idx) > 0)
        cat(length(MKS[idx]), "Markers were removed of the dataset because one or both of parents have no informed genotypes (are missing data)\n")
      
      idx <- which((GT_matrix[,P1] == "0/1" | GT_matrix[,P2] == "0/1"))
      if (length(idx) > 0)
        cat(length(MKS[idx]), "Markers were removed from the dataset because one or both of the parents are heterozygotes, we do not expect heterozygotes parents in F2 populations.\n") 
      
      idx <- which((GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "0/0") | (GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "1/1"))
      if (length(idx) > 0)
        cat(length(MKS[idx]), "Markers were removed from the dataset because they are monomorphic for the parents, these markers are not informative for the genetic map.\n") 
    }
    # Excluding non-informative markers
    rm_mk <- which(is.na(mk.type))
    if(length(rm_mk)!=0){
      GT_matrix <- GT_matrix[-rm_mk,]
      MKS <- MKS[-rm_mk]
      n.mk <- n.mk - length(rm_mk)
      CHROM <- CHROM[-rm_mk]
      POS <- POS[-rm_mk]
      REF <- REF[-rm_mk]
      ALT <- ALT[-rm_mk]
      mk.type <- mk.type[-rm_mk]
      mk.type.num <- mk.type.num[-rm_mk]
      GT_parents <- GT_parents[-rm_mk,]
    } 
    
    if(is.vector(GT_matrix)){
      jump <- 1
    } else if(dim(GT_matrix)[1]==0){
      jump <- 1
    } else jump <- 0
    
    if(jump == 1){
      warning("Input vcfR object do not have markers. An empty object onemap will be generated.")
      
      onemap.obj <- empty_onemap_obj(vcfR.obj, P1, P2, cross)
      return(onemap.obj)
    }
    
    GT_matrix[-which(GT_matrix == "1/1" | GT_matrix == "0/0" | GT_matrix == "0/1")] <- 0
    GT_matrix[which(GT_matrix == "0/1")] <- 2
    
    idx <- which(mk.type=="A.H.1")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 1
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 0
    
    idx <- which(mk.type=="A.H.2")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 0
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 1
    
    
    mk.type <- mk.type.num <- rep("A.H", n.mk)
    mk.type.num[mk.type=="A.H"] <- 8
    
  } else if(cross=="ri self" || cross=="ri sib"){
    # Marker type
    mk.type[which(GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "1/1")] <- "A.B.1"
    mk.type[which(GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "0/0")] <- "A.B.2"
    GT_parents <- GT_matrix[,c(P1,P2)]
    
    P1_1_allele <- GT_parents[,1]
    P1_2_allele <- GT_parents[,1]
    P1_1_allele[which(GT_parents[,1] == "0/0")] <- REF[which(GT_parents[,1] == "0/0")]
    P1_1_allele[which(GT_parents[,1] == "1/1")] <- ALT[which(GT_parents[,1] == "1/1")]
    P1_2_allele <- P1_1_allele
    P2_1_allele <- GT_parents[,2]
    P2_2_allele <- GT_parents[,2]
    P2_1_allele[which(GT_parents[,2] == "0/0")] <- REF[which(GT_parents[,2] == "0/0")]
    P2_1_allele[which(GT_parents[,2] == "1/1")] <- ALT[which(GT_parents[,2] == "1/1")]
    P2_2_allele <- P2_1_allele
    
    # Informs to user why markers are being removed
    if(verbose) {
      idx <- which(GT_matrix[,P1] == "./." | GT_matrix[,P2] == "./.")
      if (length(idx) > 0)
        cat(length(MKS[idx]), "Markers were removed of the dataset because one or both of parents have no informed genotypes (are missing data)\n")
      
      idx <- which((GT_matrix[,P1] == "0/1" | GT_matrix[,P2] == "0/1"))
      if (length(idx) > 0)
        cat(length(MKS[idx]), "Markers were removed from the dataset because one or both of the parents are heterozygotes, we do not expect heterozygotes parents in RILs populations.\n") 
      
      idx <- which((GT_matrix[,P1] == "0/0" & GT_matrix[,P2] == "0/0") | (GT_matrix[,P1] == "1/1" & GT_matrix[,P2] == "1/1"))
      if (length(idx) > 0)
        cat(length(MKS[idx]), "Markers were removed from the dataset because they are monomorphic for the parents, these markers are not informative for the genetic map.\n") 
    }
    # Excluding non-informative markers
    rm_mk <- which(is.na(mk.type))
    if(length(rm_mk)!=0){
      GT_matrix <- GT_matrix[-rm_mk,]
      MKS <- MKS[-rm_mk]
      n.mk <- n.mk - length(rm_mk)
      CHROM <- CHROM[-rm_mk]
      POS <- POS[-rm_mk]
      REF <- REF[-rm_mk]
      ALT <- ALT[-rm_mk]
      mk.type <- mk.type[-rm_mk]
      mk.type.num <- mk.type.num[-rm_mk]
      GT_parents <- GT_parents[-rm_mk,]
    }
    
    if(is.vector(GT_matrix)){
      jump <- 1
    } else if(dim(GT_matrix)[1]==0){
      jump <- 1
    } else jump <- 0
    
    if(jump == 1){
      warning("Input vcfR object do not have markers. An empty object onemap will be generated.")
      
      onemap.obj <- empty_onemap_obj(vcfR.obj, P1, P2, cross)
      return(onemap.obj)
    }
    
    # Onemap codification
    GT_matrix[-which(GT_matrix == "1/1" | GT_matrix == "0/0" | GT_matrix == "0/1")] <- 0
    GT_matrix[which(GT_matrix == "0/1")] <- 0
    
    idx <- which(mk.type=="A.B.1")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 1
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 3
    
    idx <- which(mk.type=="A.B.2")
    GT_matrix[idx,][which(GT_matrix[idx,] == "0/0")] <- 3
    GT_matrix[idx,][which(GT_matrix[idx,] == "1/1")] <- 1
    
    mk.type <- mk.type.num <- rep("A.B", n.mk)
    mk.type.num[mk.type=="A.B"] <- 9
    
  }
  GT_matrix[is.na(GT_matrix)] <- 0
  
  if(is.vector(GT_matrix)){
    jump <- 1
  } else if(dim(GT_matrix)[1]==0){
    jump <- 1
  } else jump <- 0
  
  if(jump == 1){
    warning("Input vcfR.objR object do not have markers. An empty object onemap will be generated.")
    
    onemap.obj <- empty_onemap_obj(vcfR.obj, P1, P2, cross)
    return(onemap.obj)
  }
  
  # Removing parents
  if(is.null(f1)){
    GT_matrix <- apply(GT_matrix[,-c(P1,P2), drop=F],2,as.numeric)
    if(!inherits(GT_matrix, "matrix")) GT_matrix <- t(as.matrix(GT_matrix)) # If there is only one marker
    colnames(GT_matrix)  <-  INDS[-c(P1,P2)] 
  } else{
    F1 <- which(dimnames(vcfR.obj@gt)[[2]]==f1) - 1
    GT_matrix <- apply(GT_matrix[,-c(P1,P2,F1), drop=F],2,as.numeric)
    if(!inherits(GT_matrix, "matrix")) GT_matrix <- t(as.matrix(GT_matrix))
    colnames(GT_matrix)  <-  INDS[-c(P1,P2,F1)] 
  }
  rownames(GT_matrix)  <- MKS
  
  ref_alt_alleles <- data.frame(P1_1_allele = P1_1_allele[match(MKS, names(P1_1_allele))], # smaller number in the VCF codification Ex: 0 if 0/1; 2 if 2/4
                                P1_2_allele = P1_2_allele[match(MKS, names(P1_2_allele))], # larger number in the VCF codification Ex: 1 if 0/1; 4 if 2/4
                                P2_1_allele = P2_1_allele[match(MKS, names(P2_1_allele))],
                                P2_2_allele = P2_2_allele[match(MKS, names(P2_2_allele))])
  
  legacy_crosses <- setNames(c("outcross", "f2", "backcross", "riself", "risib"), 
                             c("outcross", "f2 intercross", "f2 backcross", "ri self", "ri sib"))
  
  onemap.obj <- structure(list(geno= t(GT_matrix),
                               n.ind = if(!inherits(GT_matrix, "matrix")) length(GT_matrix) else dim(GT_matrix)[2],
                               n.mar = n.mk,
                               segr.type = mk.type,
                               segr.type.num = as.numeric(mk.type.num),
                               n.phe = 0,
                               pheno = NULL,
                               CHROM = CHROM,
                               POS = POS,
                               ref_alt_alleles = ref_alt_alleles,
                               input = "vcf"),
                          class=c("onemap",legacy_crosses[cross]))
  
  onemap.obj  <- rm_dupli_mks(onemap.obj)
  new.onemap.obj <- create_probs(onemap.obj, global_error = 10^-5)
  
  if(!is.null(output_info_rds)){
    info <- data.frame(CHROM = onemap.obj$CHROM, 
                       POS = onemap.obj$POS, 
                       ID = colnames(onemap.obj$geno), 
                       REF = REF, 
                       ALT = ALT)
    saveRDS(info, file = output_info_rds)
  }
  
  return(new.onemap.obj)
}


##' Convert onemap object to onemap raw file
##' 
##' Converts onemap R object to onemap input file. The input file brings information about the mapping population:
##' First line: cross type, it can be "outcrossing", "f2 intercross", "f2 backcross", "ri self" or "ri sib".
##' Second line:  number of individuals, number of markers, presence (1) or absence (0) of chromossome and position of the markers, and number of phenotypes mesured.
##' Third line: Individuals/sample names; 
##' Followed lines: marker name, marker type and genotypes. One line for each marker.
##' Final lines: chromossome, position and phenotypes informations. 
##' See more about input file format at vignettes.
##' 
##' @param onemap.obj object of class `onemap`
##' 
##' @param file.name a character for the onemap raw file name. Default is "out.raw"
##' 
##' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
##' @seealso \code{read_onemap} for a description of the output object of class onemap.
##' 
##' @return a onemap input file
##' 
##' @examples
##' 
##' data(onemap_example_out)
##' write_onemap_raw(onemap_example_out, file.name = paste0(tempfile(), ".raw"))
##' 
##'@export                  
write_onemap_raw <- function(onemap.obj=NULL, 
                             file.name = NULL){
  
  if(inherits(onemap.obj, "outcross")){
    cross <- "outcross"
  } else if(inherits(onemap.obj, "f2")){
    cross <- "f2 intercross"
  } else if(inherits(onemap.obj, "backcross")){
    cross <- "f2 backcross"
  } else if(inherits(onemap.obj, "riself")){
    cross <- "ri self"
  } else if(inherits(onemap.obj, "risib")){
    cross <- "ri sib"
  }
  
  if(is.null(file.name)) file.name = paste0(tempfile(), ".raw")
  
  fileConn<-file(file.name, "w")
  head1 <- paste("data type", cross)
  head2 <- paste(onemap.obj$n.ind,
                 onemap.obj$n.mar,
                 as.numeric(!is.null(onemap.obj$CHROM)),
                 as.numeric(!is.null(onemap.obj$POS)), 
                 onemap.obj$n.phe)
  ind.names <- rownames(onemap.obj$geno)
  if(is.null(ind.names))
    ind.names <- paste0("ID", 1:onemap.obj$n.ind)
  
  geno.mat <- onemap.obj$geno
  if(is.vector(geno.mat)) geno.mat <- matrix(geno.mat)
  
  if(inherits(onemap.obj, "outcross")){
    
    geno.mat[which(geno.mat == 0)] <- "-"
    
    idx <- which(onemap.obj$segr.type == "A.1")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "ac"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "ad"
    geno.mat[,idx][which(geno.mat[,idx]== 3)] <- "bc"
    geno.mat[,idx][which(geno.mat[,idx]== 4)] <- "bd"
    
    idx <- which(onemap.obj$segr.type == "A.2")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "ac"
    geno.mat[,idx][which(geno.mat[,idx]== 3)] <- "ba"
    geno.mat[,idx][which(geno.mat[,idx]== 4)] <- "bc"
    
    idx <- which(onemap.obj$segr.type == "A.3")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "ac"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 3)] <- "bc"
    geno.mat[,idx][which(geno.mat[,idx]== 4)] <- "b"
    
    idx <- which(onemap.obj$segr.type == "A.4")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "ab"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 3)] <- "b"
    geno.mat[,idx][which(geno.mat[,idx]== 4)] <- "o"
    
    idx <- which(onemap.obj$segr.type == "B1.5" | onemap.obj$segr.type == "B2.6" | onemap.obj$segr.type == "B3.7")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "ab"
    geno.mat[,idx][which(geno.mat[,idx]== 3)] <- "b"
    
    idx <- which(onemap.obj$segr.type == "D1.9" | onemap.obj$segr.type == "D2.14")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "ac"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "bc"
    
    idx <- which(onemap.obj$segr.type == "D1.10" | onemap.obj$segr.type == "D2.15")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "ab"
    
    idx <- which(onemap.obj$segr.type == "D1.11"  | onemap.obj$segr.type == "D2.16")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "b"
    
    idx <- which(onemap.obj$segr.type == "D1.12" | onemap.obj$segr.type == "D2.17")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "ab"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "a"
    
    idx <- which(onemap.obj$segr.type == "D1.13" | onemap.obj$segr.type == "D2.18" | onemap.obj$segr.type == "C.8")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "o"
  }
  if(inherits(onemap.obj, c("f2","backcross"))){
    
    geno.mat[which(geno.mat == 0)] <- "-"
    
    idx <- which(onemap.obj$segr.type == "A.H.B" | onemap.obj$segr.type == "A.H")
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "ab"
    geno.mat[,idx][which(geno.mat[,idx]== 3)] <- "b"
    
    idx <- which(onemap.obj$segr.type == "D.B")
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "b"
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "d"
    
    idx <- which(onemap.obj$segr.type == "C.A")
    geno.mat[,idx][which(geno.mat[,idx]== 2)] <- "a"
    geno.mat[,idx][which(geno.mat[,idx]== 1)] <- "c"
  }
  if(inherits(onemap.obj, c("riself", "risib"))){
    
    geno.mat[which(geno.mat == 0)] <- "-"
    geno.mat[which(geno.mat == 1)] <- "a"
    geno.mat[which(geno.mat == 3)] <- "b"
  }
  
  mat <- data.frame(paste0(rep("*", onemap.obj$n.mar),colnames(onemap.obj$geno)), 
                    onemap.obj$segr.type, t(geno.mat))
  colnames(mat) <- rownames(mat) <- NULL
  mat <- apply(mat, 1, function(x) paste(x, collapse = " "))
  
  writeLines(c(head1, head2, paste(ind.names, collapse = " "), mat),con =  fileConn)
  
  if(onemap.obj$n.phe > 0){
    onemap.obj$pheno[which(is.na(onemap.obj$pheno))] <- "-"
    fen <- paste0("*", paste(colnames(onemap.obj$pheno), apply(onemap.obj$pheno,2, function(x) paste(x, collapse = " "))))
    writeLines( fen, con = fileConn)
  }
  
  if(length(onemap.obj$CHROM)>0){
    onemap.obj$pheno[which(is.na(onemap.obj$pheno))] <- "-"
    chrom <- paste(paste0("*", "CHROM"), paste(onemap.obj$CHROM, collapse = " "))
    writeLines(chrom, con = fileConn)
  }
  
  if(length(onemap.obj$POS)>0){
    onemap.obj$pheno[which(is.na(onemap.obj$pheno))] <- "-"
    pos <- paste(paste0("*", "POS"), paste(onemap.obj$POS, collapse = " "))
    writeLines(pos, con = fileConn)
  }
  close(fileConn)
}
