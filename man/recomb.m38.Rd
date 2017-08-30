\name{recomb.m38}
\docType{data}
\alias{recomb.m38}
\title{Estimated fine-scale recombination rates in the mouse genome}
\description{
  Data frame containing recombination rate estimates at 240,532
  positions on chromosomes 1--19 in the mouse genome. All genomic
  positions are based on Mouse Genome Assembly 38 from the NCBI database
  (mm10, December 2011).
}
\usage{data(recomb.m38)}
\format{
  The columns of \code{recomb.m38} are as follows:
  
  \describe{
    \item{chr}{Chromosome number.}

    \item{pos}{Base-pair position on chromosome.}

    \item{recomb.rate}{Estimated fine-scale recombination rate in units
      4*Ne*r/kb, where e*Ne*r is the sex-averaged, population-scaled
      recombination rate. See the Brunschwig \emph{et al} paper for
      details.}
}
\details{
  These data were compiled using script \code{compile.recomb.data.R} in
  the \code{inst} package directory.

  The recombination rate data are from Table S1 of Brunschwig \emph{et
  al}. Since these data were based on Mouse Genome Assembly 37 (mm9), the
  LiftOver online utility from the UCSC Genome Browser website was
  applied to these data to convert to Genome Assembly 38 (mm10)
  co-coordinates.
}
\references{
  Brunschwig, H., Levi, L., Ben-David, E., Williams, R. W., Yakir, B. and
  Shifman, S. (2012). Fine-scale maps of recombination rates and hotspots
  in the mouse genome. \emph{Genetics} \bold{191}, 757--764.
  \url{http://doi.org/10.1534/genetics.112.141036}
}
\seealso{
  \code{\link{cfw.map}}, \code{\link{genes.m38}}
}
\keyword{datasets}
