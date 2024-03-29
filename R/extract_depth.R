globalVariables(c("gt.onemap.alt.ref", "gt.vcf.alt.ref", "parents"))

#' Extract allele counts of progeny and parents of vcf file
#' 
#' Uses vcfR package and onemap object to generates list of vectors with
#' reference allele count and total counts for each marker and genotypes 
#' included in onemap object (only available for biallelic sites)
#' 
#' @param vcfR.object object output from vcfR package
#' @param onemap.object onemap object output from read_onemap, read_mapmaker or onemap_read_vcf function
#' @param vcf.par vcf format field that contain allele counts informations, the implemented are: AD, DPR, GQ, PL, GL. 
#' AD and DPR return a list with allele depth information. GQ returns a matrix with error probability for each genotype. 
#' PL return a data.frame with genotypes probabilities for every genotype. 
#' @param parent1 parent 1 identification in vcfR object
#' @param parent2 parent 2 identification in vcfR object
#' @param f1 if your cross type is f2, you must define the F1 individual
#' @param recovering TRUE/FALSE, if TRUE evaluate all markers from vcf file, if FALSE evaluate only markers in onemap object
#' 
#' @return list containing the following components: \item{palt}{a \code{matrix} with parent 1 and 2 
#' alternative allele counts.} \item{pref}{a \code{matrix} with parent 1 and 2 
#' reference allele counts.} \item{psize}{a \code{matrix} with parent 1 and 2 
#' total allele counts.}\item{oalt}{a \code{matrix} with progeny 
#' alternative allele counts.}\item{oref}{a \code{matrix} with progeny 
#' reference allele counts.}\item{osize}{a \code{matrix} with progeny 
#' total allele counts.}\item{n.mks}{total number of markers.} 
#' 
#' \item{n.ind}{total number of individuals in progeny.} \item{inds}{progeny individuals identification.}
#' \item{mks}{markers identification.} \item{onemap.object}{same onemap.object inputted}
#' 
#' @import dplyr
#' @importFrom vcfR extract.gt masplit
#' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu} 
#' @export
extract_depth <- function(vcfR.object=NULL,
                          onemap.object= NULL,
                          vcf.par = c("GQ","AD", "DPR, PL", "GL"),
                          parent1="P1",
                          parent2="P2",
                          f1="F1",
                          recovering = FALSE){
  if(is.null(vcfR.object))
    stop("You must specify one vcfR object.")
  
  if(!inherits(vcfR.object,"vcfR"))
    stop("You must specify one vcfR object.")
  
  if(is.null(onemap.object))
    stop("You must specify one onemap object.")
  
  if(!inherits(onemap.object,"onemap"))
    stop("You must specify one onemap object.")
  
  # Infos                                                                                                        
  ind <- rownames(onemap.object$geno)
  IND <- colnames(vcfR.object@gt)[-1]
  mks <- colnames(onemap.object$geno)
  MKS <- vcfR.object@fix[,3]
  n.mks <- length(mks)
  n.ind <- length(ind)
  N.MKs <- dim(vcfR.object@gt)[1]
  N.IND <- dim(vcfR.object@gt)[2]-1
  chr.pos.vcf <- paste0(vcfR.object@fix[,1], "_",as.numeric(vcfR.object@fix[,2]))
  chr.pos.onemap <- paste0(onemap.object$CHROM,"_", onemap.object$POS)
  
  # If there are no marker names                                                                                 
  if(anyNA(MKS))
    MKS <- chr.pos.vcf
  
  parents <- c(which(IND == parent1),which(IND == parent2))
  
  if(recovering==FALSE){
    rm.mks <- which(!(chr.pos.vcf %in% chr.pos.onemap))
    rm.ind <- which(IND[-parents] %in% ind==FALSE)                                                                             
    CHROM <- onemap.object$CHROM
    POS <- onemap.object$POS
  } else {
    CHROM <- vcfR.object@fix[,1]
    POS <- vcfR.object@fix[,2]
    rm.mks <- NULL
    rm.ind <- NULL
  }
  
  par_matrix <- extract.gt(vcfR.object, element = vcf.par)
  
  # Replacing missing data with compatible format                                                                
  if(length(which(par_matrix == ".")) > 0 | length(which(is.na(par_matrix))) > 0 ){
    if(vcf.par=="GQ") {
      par_matrix[which(par_matrix == ".")] <- NA
    } else if (vcf.par == "PL" | vcf.par == "GL") { 
      par_matrix[which(par_matrix == "." | is.na(par_matrix))] <- "0,0,0"
    } else {
      par_matrix[which(par_matrix == ".")] <- "0,0"
    }
  }
  if(length(rm.mks)>0 & length(rm.ind)>0){
    par_matrix <- par_matrix[-rm.mks, -rm.ind]
    IND <- IND[-rm.ind]
    MKS <- MKS[-rm.mks]
  } else if(length(rm.mks)>0){
    par_matrix <- par_matrix[-rm.mks,]
    MKS <- MKS[-rm.mks]
  } else if(length(rm.ind)>0){
    par_matrix <- par_matrix[,-rm.ind]
    IND <- IND[-rm.ind]
  }
  
  n.ind <- N.IND - length(rm.ind)
  n.mks <- N.MKs - length(rm.mks)
  # The probabilities must be calculated if AD or DPR parameters were chosen                                    
  if(vcf.par=="AD" | vcf.par=="DPR"){
    ref_matrix <- masplit(par_matrix, record=1, sort=0)
    alt_matrix <- masplit(par_matrix, record=2, sort=0)
  } else if(vcf.par=="GQ"){
    error_matrix <- 10^(-apply(par_matrix,1,as.numeric)/10)
    idx <- which(IND %in% c(parent1, parent2, f1))
    error_matrix <- error_matrix[-idx,]
    rownames(error_matrix) <- ind
    return(error_matrix)
  } else if (vcf.par == "PL" | vcf.par == "GL") {
    idx <- which(IND %in% c(parent1, parent2, f1))
    if(vcf.par == "PL"){
      probs <- par_matrix  %>% .[,-idx] %>% strsplit(., ",") %>% 
        do.call(rbind, .) %>% apply(., 2,as.numeric) %>% 
        apply(., 2, function(x) 10^(-x/10)) 
    } else {
      probs <- par_matrix %>% .[,-idx] %>% strsplit(., ",") %>% 
        do.call(rbind, .) %>% apply(., 2,as.numeric) %>% 
        apply(., 2, function(x) 10^(x)) 
    }
    
    probs[apply(probs, 1, function(x) all(x == 1)),] <- NA
    
    sums <- apply(probs, 1, sum)
    comb <- expand.grid(MKS, IND[-idx])
    
    # Probs sum 1
    probs <- probs/sums 
    probs[apply(probs, 1, function(x) all(is.na(x))),] <- 1
    
    rownames(probs) <- paste0(comb$Var1, "_", comb$Var2)
    return(probs)
  } 
  
  if(vcf.par!="GQ" & vcf.par!="PL"){
    if(vcf.par=="DPR"){
      size_matrix <- ref_matrix
      ref_matrix <- size_matrix - alt_matrix
      alt_matrix <- size_matrix - ref_matrix
    } else if(vcf.par=="AD"){
      size_matrix <- ref_matrix + alt_matrix
    }
    
    # Splitting offspring and parents                                                                               
    idx <- parents
    palt <- alt_matrix[,idx]
    pref <- ref_matrix[,idx]
    psize <- size_matrix[,idx]
    
    if(recovering==TRUE){
      IND <- IND[-c(idx)]
      oalt <- alt_matrix[,-idx]
      oref <- ref_matrix[,-idx]
      osize <- size_matrix[,-idx]
    } else {
      IND <- IND[-parents]
      oalt <- alt_matrix[,-parents]
      oref <- ref_matrix[,-parents]
      osize <- size_matrix[,-parents]
    }
    
    n.ind <- dim(oref)[2]
    n.mks <- dim(oref)[1]
    
    structure(list(palt=palt,
                   pref=pref,
                   psize=psize,
                   oalt=oalt,
                   oref=oref,
                   osize=osize,
                   n.mks=n.mks,
                   n.ind=n.ind,
                   inds = IND,
                   mks = MKS,
                   CHROM = CHROM,
                   POS = POS,
                   onemap.object = onemap.object))
  }
}

