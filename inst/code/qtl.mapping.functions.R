# This function removes outliers from the phenotype data for the given
# phenotype, optionally conditioning on covariates. If covariates are
# specified, outliers are determined according to the residual of the
# phenotype after regressing on the covariates.
remove.outliers <- function (pheno, phenotype, covariates, outliers) {

  # Specify the function for removing the outliers.
  is.outlier <- function (x) {
    y           <- outliers(x)
    y[is.na(y)] <- FALSE
    return(y)
  }
  
  # If we are conditioning on one or more covariates, get the
  # residuals of the phenotype conditioned on these covariates.
  if (length(covariates) > 0) {
    f <- formula(paste(phenotype,"~",paste(covariates,collapse="+")))
    r <- resid(lm(f,pheno,na.action = na.exclude))
  } else
    r <- pheno[[phenotype]]

  # Remove the outlying data points.
  pheno[is.outlier(r),phenotype] <- NA

  # Return the updated phenotype data table.
  return(pheno)
}

# This function writes the phenotype data to a file in the format used
# by GEMMA. Each line of the file contains one phenotype observation.
write.gemma.pheno <- function (file, phenotype, pheno) {
  y <- pheno[phenotype]
  if (is.numeric(y))
    y <- round(y,digits = 6)
  write.table(y,file,quote = FALSE,row.names = FALSE,col.names = FALSE)
}

# This function writes the covariate data to a file in the format used
# by GEMMA. Each line corresponds to a sample. We must include an
# additional covariate for the intercept.
write.gemma.covariates <- function (file, covariates, pheno) {
  if (is.null(covariates)) {
    write.table(data.frame(rep(1,nrow(pheno))),file,sep = " ",quote = FALSE,
                row.names = FALSE,col.names = FALSE)
  } else {
    write.table(cbind(1,data.frame(lapply(pheno[covariates],function (x) {
      if (is.numeric(x))
        round(x,digits=6)
      else
        x
    }))),file,sep = " ",quote = FALSE,row.names = FALSE,col.names = FALSE)
  }
}

# This function writes the SNP information to a space-delimited text
# file in the format used by GEMMA. The file contains one line per
# SNP, with three columns: id, base-pair position and chromosome.
write.gemma.map <- function (file, map)
  write.table(map[c("id","pos","chr")],file,sep = " ",quote = FALSE,
              row.names = FALSE,col.names = FALSE)

# Store the genotype dosages as a space-delimited text file in the
# format used by GEMMA, in which we have one row per SNP, and one
# column per sample. The first three columns give the SNP id and
# the two SNP alleles.
write.gemma.geno <- function (file, geno, map) {
  geno <- t(geno)
  geno <- as.data.frame(geno,check.names = FALSE)
  geno <- round(geno,digits = 3)
  geno <- cbind(map[c("id","ref","alt")],geno)
  write.table(geno,file,sep = " ",quote = FALSE,row.names = FALSE,
              col.names = FALSE)
}

# This function reads in the GEMMA association results from GEMMA, and
# returns a data frame containing 4 columns: chromosome number
# ("chr"); base-pair position ("pos"); SNP id ("id"); and the negative
# base-10 logarithm of the p-value ("log10p").
read.gemma.results <- function (file) {
  out <- read.table(file,sep = "\t",header = TRUE,check.names = FALSE,
                       quote = "",stringsAsFactors = FALSE)
  out <- out[c("chr","ps","rs","p_lrt")]
  out <- transform(out,p_lrt = -log10(p_lrt))
  colnames(out) <- c("chr","pos","id","log10p")
  return(out)
}
