\name{genes.m38}
\docType{data}
\alias{genes.m38}
\title{RefSeq gene annotations for mouse genome assembly 38 (mm10)}
\description{
  Data frame containing 23,452 protein-coding gene annotations for
  release 38 of the Mouse Genome Assembly (mm10). This includes genes on
  autosomal chromosomes (1--19), X, Y and the mitochondrial DNA ("MT").
}
\usage{data(genes.m38)}
\format{
  The columns of \code{genes.m38} are as follows:
  
  \describe{
    \item{gene.symbol}{Gene symbol.}

    \item{gene.id}{NCBI GeneID.}
    
    \item{chr}{Chromosome.}

    \item{start}{Start position of transcript.}

    \item{stop}{Stop position of transcript.}
  }
}
\details{
  These data were compiled using script \code{compile.refseq.data.R} in
  the \code{inst} package directory.
}
\source{
  Data were retrieved and compiled from the RefSeq Genome Assemblies FTP site,
  \url{ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq}.
}
\seealso{
  \code{\link{cfw.map}}
}
\examples{
# View the first few genes annotated to chromosome 1.
data(genes.m38)
head(genes.m38)
}
\keyword{datasets}
