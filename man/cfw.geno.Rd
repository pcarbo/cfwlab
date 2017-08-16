\name{cfw.geno}
\docType{data}
\alias{cfw.geno}
\title{Genotype data for 1,038 mice at 79,824 SNPs}
\description{
  Genotype data for CFW mice obtained by genotype-by-sequencing (GBS)
  technology.
}
\usage{data(cfw.geno)}
\format{
  The genotype data are stored in an 1,038 x 79,824 matrix
  \code{cfw.geno}. Rows are labeled by sample id, and columns are
  labeled by marker (SNP) id. 
  
  Genotypes are encoded as "dosages"; specifically, the expected
  number of times the alternative allele is observed in the
  genotype. This will either be an integer (0, 1 or 2), or a real number
  between 0 and 2 when there is some uncertainty in the estimate of the
  genotype.
}
\details{
  These data were compiled using script \code{compile.gwas.data.R} in
  the \code{inst} package directory.
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
cfw.geno[1:8,1:5]

# All genotype dosages should be between 0 and 2.
range(cfw.geno)
}
\keyword{datasets}
