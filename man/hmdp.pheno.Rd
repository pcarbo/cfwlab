\name{hmdp.pheno}
\docType{data}
\alias{hmdp.pheno}
\title{Bone-mass density data from HMDP study.}
\description{
  Data frame containing bone-mass density (BMD) data from Table S1 of
  Farber \emph{et al}, 2011. These data are part of a QTL mapping study
  for BMD traits using mice from the Hybrid Mouse Diversity Panel
  (HMDP). The BMD measurements were taken in 878 mice from 97 HMDP
  strains.
}
\usage{data(hmdp.pheno)}
\format{
  The three BMD measurements are total body areal BMD
  (\code{"totalbody"}), femur areal BMD (\code{"femur"}), and lumbar
  spine areal BMD (\code{"spine"}). Units are mg/cm^2.
}
\details{
  These data were compiled using script \code{compile.hmdp.data.R} in
  the \code{inst} package directory.
}
\references{
  Farber, C. R., Bennett, B. J., Orozco, L., Zou, W., Lira, A.,
  Kostem, E., Kang, H. M., Furlotte, N., Berberyan, A., Ghazalpour, A.,
  Suwanwela, J., Drake, T. A., Eskin, E., Wang, Q. T., Teitelbaum,
  S. L. and Lusis, A. J. (2011) Mouse genome-wide association and systems
  genetics identify Asxl2 as a regulator of bone mineral density and
  osteoclastogenesis. \emph{PLoS Genetics} \bold{7}, e1002038.
  \url{https://doi.org/10.1371/journal.pgen.1002038}
}
\examples{
# View the first few rows of the phenotype data table.
data(hmdp.pheno)
head(hmdp.pheno)
}
\keyword{datasets}
