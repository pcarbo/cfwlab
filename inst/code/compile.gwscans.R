# Compile genome-wide scans for several physiological traits.
# Association p-values are computed for all SNPs using GEMMA. To run
# this script, you must first generate the phenotype and genotype data
# files using script compile.gwas.data.R, and store these data files
# in the "data" folder of this repository.
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
  abBMD   = list(pheno="abBMD",cov="SW16",outliers=NULL),
  testis  = list(pheno="testisweight",cov="sacweight",
                 outliers=function (x) x < (-0.075)))

# Load the SNP data.
load("../../data/cfw.map.RData")

# Initialize the table containing the QTL mapping results.
cfw.gwscan        <- as.data.frame(matrix(0,nrow(cfw.map),length(analyses)))
names(cfw.gwscan) <- names(analyses)
cfw.gwscan        <- cbind(cfw.map[c("id","chr","pos")],cfw.gwscan)
    
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
  
  # Only analyze samples (i.e., rows) for which the phenotype and all
  # the covariates are observed.
  rows      <- which(!apply(is.na(cfw.pheno),1,any))
  cfw.pheno <- cfw.pheno[rows,]
  cfw.geno  <- cfw.geno[rows,]

  # Give summary of QTL mapping analysis.
  cat(" - Mapping QTLs for",phenotype,"in",nrow(cfw.pheno),"mice, ")
  if (!is.null(covariates)) {
    cat("controlling for ",paste(covariates,collapse=" + "),".\n",sep="")
  } else {
    cat("with no covariates included.\n")
  }

  # Write the phenotype data to a text file.
  cat(" - Writing phenotype data file.\n")
  write.gemma.pheno("pheno.txt",phenotype,cfw.pheno)

  # Write the covariate data to a text file.
  cat(" - Writing covariate data file.\n")
  write.gemma.covariates("cov.txt",covariates,cfw.pheno)

  # Write out the genotypes and SNP information.
  cat(" - Writing genotype data to file.\n");
  write.gemma.map("map.txt",cfw.map)
  write.gemma.geno("geno.txt",cfw.geno,cfw.map)

  # Run GEMMA.
  # chromosome using the kinship matrix computed using all the
  # markers *not* on the chromosome.
  cat(" - Computing p-values for",nrow(cfw.map),"candidate SNPs.\n")
  system("./gemma -g geno.txt -a map.txt -p pheno.txt -c cov.txt -lm 2",
         ignore.stdout = TRUE)

  # Load the results of the GEMMA association analysis.
  cat(" - Reading GEMMA results from file.\n")
  cfw.gwscan[[which.analysis]] <-
      read.gemma.results("output/result.assoc.txt")$log10p
}

# Save results to file.
save(list = "cfw.gwscan",file = "cfw.gwscan.RData")
