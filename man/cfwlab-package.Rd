\name{cfwlab-package}
\alias{cfwlab-package}
\alias{cfwlab}
\docType{package}
\title{CFW genome-wide association study data}
\description{
This R package contains data from Parker \emph{et al}, 2016. These data
were collected as part of a large study to assess the viability of using
Carworth Farms White (CFW) mice for mapping genes and genetic loci
underlying complex traits relevant to the study of human disease and
psychology.

For code implementing QTL mapping of physiological, behavioral and
gene expression phenotypes, and other analyses of these data, see
\url{https://github.com/pcarbo/cfw}.

To use these data for your research, please cite the references below.
} \references{
  Parker, C. C., Gopalakrishnan, S., Carbonetto, P., Gonzales, N.M.,
  Leung, E., Park, Y.J., Aryee, E., Davis, J., Blizard, D.A.,
  Ackert-Bicknell, C.L., Lionikas, A., Pritchard, J.K. and Palmer,
  A.A. (2016) Genome-wide association study of behavioral, physiological
  and gene expression traits in outbred CFW mice. \emph{Nature Genetics}
  \bold{48}, 919--926. \url{http://dx.doi.org/10.1038/ng.3609}

  Parker, C. C., Gopalakrishnan, S., Carbonetto, P., Gonzales, N.M.,
  Leung, E., Park, Y.J., Aryee, E., Davis, J., Blizard, D.A.,
  Ackert-Bicknell, C.L., Lionikas, A., Pritchard, J.K. and Palmer,
  A.A. (2016) Data from: Genome-wide association study of behavioral,
  physiological and gene expression traits in outbred CFW
  mice. \emph{Dryad Digital Repository}. 
  \url{http://dx.doi.org/10.5061/dryad.2rs41}
}
\author{Peter Carbonetto and John Novembre}
\seealso{
  \url{https://github.com/pcarbo/cfwlab}
}
\examples{
  # Load all the data from this package.
  data(cfw.pheno)
  data(hmdp.pheno)
  data(cfw.map)
  data(cfw.geno)
  data(cfw.gwscan)
  data(genes.m38)
  data(recomb.m38)
}