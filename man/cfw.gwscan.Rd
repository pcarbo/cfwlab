\name{cfw.gwscan}
\docType{data}
\alias{cfw.gwscan}
\title{Selected QTL mapping results from CFW study}
\description{
  Results of genome-wide scans for 8 physiological phenotypes from CFW
  study. Association p-values were computed using the program GEMMA.
}
\usage{data(cfw.gwscan)}
\format{
  Individual columns of the table \code{cfw.gwscan} are as follows:
  
  \describe{
    \item{id}{Unique label assigned to each SNP.}

    \item{chr}{Chromosome number.}

    \item{pos}{Base-pair position on chromosome.}

    \item{TA}{p-values for association with tibialis anterior (TA)
      muscle weight, controlling for tibia length.}

    \item{EDL}{p-values for association with extensor digitorum longus
      (EDL) muscle weight, controlling for tibia length.}

    \item{soleus}{p-values for association with soleus muscle weight,
      controlling for tibia length.}

    \item{plant}{p-values for association with plantaris muscle weight,
      controlling for tibia length.}

    \item{gastroc}{p-values for association with gastrocnemius
      ("gastroc") muscle weight, controlling for tibia length.}
    
    \item{tibia}{p-values for association with tibia length, controlling
      for body weight.}

    \item{abBMD}{p-values for association with "abnormally high"
      bone-mass density.}

    \item{testis}{p-values for association with testis weight,
      controlling for body weight.}
  }
}
\details{
  These QTL mapping results were generated using \code{gemma} version
  0.96 and script \code{compile.gwscans.R} in the \code{inst} package
  directory.
}
\references{
  See \link{cfwlab-package} for a list of references.
}
\seealso{
  \code{\link{cfw.map}}, \code{\link{cfw.geno}}, \code{\link{cfw.pheno}}
}
\examples{
# Get the strongest association (minimum p-value) for each phenotype
# analyzed.
data(cfw.gwscan)
10^(-sapply(cfw.gwscan[-(1:3)],max))
}
\keyword{datasets}
