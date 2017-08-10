# Compile the phenotype and genotype data from files stored in theq
# Data Dryad repository.
library(curl)
suppressMessages(library(data.table))

# Read in the phenotype data from the Data Dryad repository.
cat("Reading phenotype data.\n")
pheno <- read.csv(curl(paste("https://datadryad.org/bitstream/handle",
                             "10255/dryad.117919/pheno.csv",sep = "/")),
                  quote = "",header = TRUE,check.names = FALSE,
                  stringsAsFactors = FALSE,comment.char = "#")

# Convert some of the columns to factors.
cat("Preparing phenotype data.\n")
pheno <- transform(pheno,
                   id            = as.character(id),
                   round         = factor(round,paste0("SW",1:25)),
                   FCbox         = factor(FCbox),
                   PPIbox        = factor(PPIbox),
                   methcage      = factor(methcage),
                   methcycle     = factor(methcycle),
                   discard       = factor(discard),
                   mixup         = factor(mixup),
                   earpunch      = factor(earpunch),
                   abnormalbone  = factor(abnormalbone),
                   experimenters = factor(experimenters))

# Convert some of the columns to double precision.
pheno <- transform(pheno,fastglucose = as.double(fastglucose))

# Remove rows marked as "discard" and as possible sample mixups.
pheno <- subset(pheno,discard == "no" & mixup == "no")

# Create a new binary trait that is equal to 1 when the mouse has an
# "abnormal" bone-mineral density, and 0 otherwise.
pheno <- cbind(pheno,data.frame(abBMD = factor(as.numeric(pheno$BMD > 90))))

# Create a binary trait from the p120b1 trait that indicates that
# the mouse does not respond to the pulse, possibly indicating
# deafness. I choose mice with p120b1 values on the "long tail" of
# the distribution for this phenotype.
pheno <- cbind(pheno,data.frame(deaf = factor(as.numeric(pheno$p120b1 < 50))))

# Remove most of the behavioural phenotype data.
pheno <- pheno[c("id","round","cageid","FCbox","PPIbox","methcage",
                 "methcycle","earpunch","glucoseage","methage",
                 "FCage","PPIage","sacage","bw0","bw1","bw2","bw3",
                 "PPIweight","sacweight","BMD","TA","EDL","gastroc",
                 "plantaris","soleus","tibia","abnormalbone",
                 "experimenters","testisweight","taillength",
                 "fastglucose","PreTrainD1","AvToneD1","AvContextD2",
                 "AvAltContextD3","AvToneD3","abBMD","deaf",
                 "pp3PPIavg","pp6PPIavg","pp12PPIavg")]

# Read the marker (i.e., SNP) data. 
cat("Reading marker data.\n")
bases   <- c("A","T","G","C")
map <- read.table(curl(paste("https://datadryad.org/bitstream/handle",
                             "10255/dryad.117921/map.txt",sep = "/")),
                  sep = " ",header = TRUE,stringsAsFactors = FALSE)
map <- transform(map,
                 chr = factor(chr,1:19),
                 ref = factor(ref,bases),
                 alt = factor(alt,bases))

# Read the genotype data. To avoid downloading a large file every time
# this script is run, you will need to download and uncompress
# geno.txt.gz manually from the Data Dryad webpage:
#
#   http://dx.doi.org/10.5061/dryad.2rs41
#
cat("Reading genotype data.\n")
n    <- nrow(map)
geno <- fread("geno.txt",sep = " ",header = TRUE,showProgress = FALSE,
              colClasses = c("character","character",rep("double",n)))
class(geno)    <- "data.frame"
rownames(geno) <- geno$id
discard        <- factor(geno$discard)
geno           <- as.matrix(geno[-(1:2)])

# Discard genotype samples from mislabeled flowcell samples.
cat("Preparing & filtering genotype data.\n")
geno <- geno[discard == "no",]

# Align the phenotypes and genotypes
ids   <- intersect(pheno$id,rownames(geno))
pheno <- pheno[match(ids,pheno$id),]
geno  <- geno[match(ids,rownames(geno)),]

# Returns the minor allele frequency given a vector of genotypes
# encoded as allele counts.
compute.maf <- function (geno) {
  f <- mean(geno,na.rm = TRUE)/2
  return(min(f,1-f))
}

# Discard SNPs with low "imputation quality" assessed by inspecting
# the genotype probabilities. Retain SNPs for which: (1) at least 95%
# of the samples have a maximum probability genotype greater than than
# 0.5; (2) the minor allele frequency is greater than 2%.
f       <- apply(geno,2,compute.maf)
markers <- which(map$quality > 0.95 & f > 0.02)
map     <- map[markers,]
geno    <- geno[,markers]

# Save the data to .RData files.
cat("Saving data to .RData files.\n")
cfw.pheno <- pheno
cfw.map   <- map
cfw.geno  <- geno
save(list = "cfw.pheno",file = "cfw.pheno.RData")
save(list = "cfw.map",  file = "cfw.map.RData")
save(list = "cfw.geno", file = "cfw.geno.RData")
