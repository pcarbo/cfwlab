# TO DO: Explain here what this script does, and how to use it.
#
# NOTES:
#
#  * First run compile.gwas.data.R before running this script.
#
source("qtl.mapping.functions.R")

# This data structure provides information about all the QTL
# analyses.
analyses <- list(
  TA      = list(pheno="TA",cov=c("SW16","tibia"),
                 outliers=function (x) x < (-18)),
  EDL     = list(pheno="EDL",cov=c("SW16","tibia"),
                 outliers=function (x) x < (-5) | x > 4),
  soleus  = list(pheno="soleus",cov=c("SW16","tibia"),
                 outliers=function (x) x < (-4) | x > 4),
  plant   = list(pheno="plantaris",cov=c("SW16","tibia"),
                 outliers=function (x) x < (-9) | x > 8),
  gastroc = list(pheno="gastroc",cov=c("SW16","tibia"),
                 outliers=function (x) x < (-40) | x > 50),
  tibia   = list(pheno="tibia",cov=c("SW6","SW16","sacweight"),
                 outliers=function (x) x < (-1.5)),
  BMD     = list(pheno="BMD",cov="SW16",outliers=function (x) x > 0.14),
  abBMD   = list(pheno="abBMD",cov="SW16",outliers=NULL),
  testis  = list(pheno="testisweight",cov="sacweight",
                 outliers=function (x) x < (-0.075)))

# Initialize the table containing the QTL mapping results.
cfw.qtls <- NULL

# Repeat for each QTL analysis.
for (which.analysis in names(analyses)) {
  cat("Running",which.analysis,"analysis:\n")

  # Get the phenotype and covariates used in the QTL mapping.
  analysis   <- analyses[[which.analysis]]
  phenotype  <- analysis$pheno
  covariates <- analysis$cov
  outliers   <- analysis$outliers

  # Load the phenotype data.
  cat(" - Loading and preparing phenotype data.\n")
  load("../../data/cfw.pheno.RData")

  # Create some binary covariates.
  cfw.pheno <-
    cbind(cfw.pheno,
          data.frame(SW6 = factor(as.integer(cfw.pheno$round == "SW6")),
                     SW16 = factor(as.integer(cfw.pheno$round == "SW16"))))

  # Retain only the table columns needed for the data analysis.
  cfw.pheno <- cfw.pheno[c(phenotype,covariates)]
  
  # Discard outliers.
  if (!is.null(outliers))
    cfw.pheno <- remove.outliers(cfw.pheno,phenotype,covariates,outliers)

  # Load the genotype data.
  cat(" - Loading genotype data.\n")
  load("../../data/cfw.map.RData")
  load("../../data/cfw.geno.RData")
  
  # Only analyze samples (i.e. rows of the genotype and phenotype
  # matrices) for which the phenotype and all the covariates are
  # observed.
  rows      <- which(!apply(is.na(cfw.pheno),1,any))
  cfw.pheno <- cfw.pheno[rows,]
  cfw.geno  <- cfw.geno[rows,]

  # Compute association p-values using GEMMA.
  gwscan.gemma <- run.gemma(phenotype,covariates,pheno,X,map,
                            gemmadir,gemma.exe,chromosomes)

  # Add a new row to the QTL mapping results.
  
  stop()
}

# Save results to file.
# TO DO.
