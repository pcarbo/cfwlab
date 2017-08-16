\name{cfw.gwscan}
\docType{data}
\alias{cfw.gwscan}
\title{Selected QTL mapping results from CFW study}
\description{
  Results of genome-wide scans for 8 physiological phenotypes from CFW
  study. Association p-values were computed using the program GEMMA.
}
\usage{data(cfw.pheno)}
\format{
  Individual columns of the table \code{cfw.gwscan} are as follows:
  
  \describe{
    \item{id}{Unique label assigned to each SNP.}

    \item{chr}{Chromosome number.}

    \item{pos}{Base-pair position on chromosome.}

    \item{TA}{}

    \item{EDL}{}

    \item{soleus}{}

    \item{plant}{}

    \item{gastroc}{}
    
    \item{tibia}{}

    \item{abBMD}{}

    \item{testis}{}
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
