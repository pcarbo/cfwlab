\name{cfw.geno}
\docType{data}
\alias{cfw.geno}
\title{TITLE GOES HERE}
\description{
  DESCRIPTION GOES HERE.
}
\usage{data(cfw.geno)}
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
\source{
  Data were retrieved and prepared from
  \url{http://dx.doi.org/10.5061/dryad.2rs41}.
}
\references{
  See \link{cfwlab-package} for a list of references.
}
\seealso{
  \code{\link{cfw.map}}, \code{\link{cfw.pheno}}
}
\examples{
# View the first several rows (samples) and columns (SNPs) of the
# genotype matrix.
data(cfw.geno)
geno[1:8,1:5]
}
\keyword{datasets}
