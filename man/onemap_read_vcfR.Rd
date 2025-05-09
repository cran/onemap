% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onemap_read_vcfR.R
\name{onemap_read_vcfR}
\alias{onemap_read_vcfR}
\title{Convert vcf file to onemap object}
\usage{
onemap_read_vcfR(
  vcf = NULL,
  vcfR.object = NULL,
  cross = NULL,
  parent1 = NULL,
  parent2 = NULL,
  f1 = NULL,
  only_biallelic = TRUE,
  output_info_rds = NULL,
  verbose = TRUE
)
}
\arguments{
\item{vcf}{string defining the path to VCF file;}

\item{vcfR.object}{object of class vcfR;}

\item{cross}{type of cross. Must be one of: \code{"outcross"} for full-sibs;
\code{"f2 intercross"} for an F2 intercross progeny; \code{"f2 backcross"};
\code{"ri self"} for recombinant inbred lines by self-mating; or
\code{"ri sib"} for recombinant inbred lines by sib-mating.}

\item{parent1}{\code{string} specifying sample ID of the first parent. If f2 backcross population, define here the ID of the backcrossed parent.}

\item{parent2}{\code{string} specifying sample ID of the second parent.}

\item{f1}{\code{string} if you are working with f2 intercross or backcross populations you may have f1 parents in you vcf, specify its ID here}

\item{only_biallelic}{if TRUE (default) only biallelic markers are considered, if FALSE multiallelic markers are included.}

\item{output_info_rds}{define a name for the file with alleles information.}

\item{verbose}{A logical, if TRUE it output progress status
information.}
}
\value{
An object of class \code{onemap}, i.e., a list with the following
components: \item{geno}{a matrix with integers indicating the genotypes
read for each marker. Each column contains data for a marker and each row
represents an individual.} \item{n.ind}{number of individuals.}
\item{n.mar}{number of markers.} \item{segr.type}{a vector with the
segregation type of each marker, as \code{strings}.} \item{segr.type.num}{a
vector with the segregation type of each marker, represented in a
simplified manner as integers, i.e. 1 corresponds to markers of type
\code{"A"}; 2 corresponds to markers of type \code{"B1.5"}; 3 corresponds
to markers of type \code{"B2.6"}; 4 corresponds to markers of type
\code{"B3.7"}; 5 corresponds to markers of type \code{"C.8"}; 6 corresponds
to markers of type \code{"D1"} and 7 corresponds to markers of type
\code{"D2"}. Markers for F2 intercrosses are coded as 1; all other crosses
are left as \code{NA}.} \item{input}{the name of the input file.}
\item{n.phe}{number of phenotypes.} \item{pheno}{a matrix with phenotypic
values. Each column contains data for a trait and each row represents an
individual.} \item{error}{matrix containing HMM emission probabilities}
}
\description{
Converts data from a vcf file to onemap initial object, while identify 
the appropriate marker segregation patterns.
}
\details{
Only biallelic SNPs and indels for diploid variant sites are considered.

Genotype information on the parents is required for all cross types. For
full-sib progenies, both outbred parents must be genotyped. For backcrosses,
F2 intercrosses and recombinant inbred lines, the \emph{original inbred
lines} must be genotyped. Particularly for backcross progenies, the
\emph{recurrent line must be provided as the first parent} in the function
arguments.

Marker type is determined based on parental genotypes. Variants for which parent
genotypes cannot be determined are discarded.

Reference sequence ID and position for each variant site are also stored.
}
\examples{
\donttest{
data <- onemap_read_vcfR(vcf=system.file("extdata/vcf_example_out.vcf.gz", package = "onemap"),
                 cross="outcross",
                 parent1=c("P1"),
                 parent2=c("P2"))
}
                
}
\seealso{
\code{read_onemap} for a description of the output object of class onemap.
}
\author{
Cristiane Taniguti, \email{chtaniguti@tamu.edu}
}
