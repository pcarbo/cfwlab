\name{cfw.map}
\docType{data}
\alias{cfw.map}
\title{Information on 79,824 SNPs genotyped in CFW mice}
\description{
  Data frame containing information about 79,824 single nucleotide
  polymorphisms (SNPs) on chromosomes 1-19 genotyped in the CFW outbred
  mouse cohort. All selected SNPs have >2\% minor allele frequency and
  "imputation quality" >95\%. All genomic positions are based on Mouse
  Genome Assembly 38 from the NCBI database (mm10, December 2011).
}
\usage{data(cfw.map)}
\format{

The columns of \code{cfw.map} are as follows:
  
  \describe{
  
    \item{id}{Unique label assigned to each SNP. This is either the
      dbSNP database id (release 137), or if the SNP is not registered
      in dbSNP, then it is an internal id of the form "cfw-N-XXXXXX",
      where N is the chromsome number, and XXXXXX is the base-pair
      position of the SNP on the chromosome.}

    \item{chr}{Chromosome number.}

    \item{pos}{Base-pair position on chromosome.}

    \item{ref}{Reference allele. (This always matches the dbSNP reference
      allele.)}
    
    \item{alt}{Alternative allele. This usually matches one of the
      alternative alleles in the corresponding dbSNP entry.}

    \item{quality}{Imputation quality of each SNP. Here we define
      "imputation quality" as the proportion of samples in which the
      probability of the maximum-probability genotype is >0.5.}
  }
}
\details{
  These data were compiled using script \code{compile.gwas.data.R} in
  the \code{inst} package directory.
}
\source{
  Data were retrieved and compiled from
  \url{http://dx.doi.org/10.5061/dryad.2rs41}.
}
\references{
  See \link{cfwlab-package} for a list of references.
}
\seealso{
  \code{\link{cfw.geno}}
}
\examples{
# View the first few rows of the marker data table.
data(cfw.map)
head(cfw.map)
}
\keyword{datasets}
