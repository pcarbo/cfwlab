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
