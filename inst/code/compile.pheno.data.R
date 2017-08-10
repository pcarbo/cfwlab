# Compile the phenotype data from the CSV file stored in the Data
# Dryad repository.
library(curl)

# Read in the phenotype data from the Data Dryad repository.
cfw.pheno <- read.csv(curl(paste("https://datadryad.org/bitstream/handle",
                                 "10255/dryad.117919/pheno.csv",sep = "/")),
                      quote = "",header = TRUE,check.names = FALSE,
                      stringsAsFactors = FALSE,comment.char = "#")


# Convert some of the columns to factors.
cfw.pheno <- transform(cfw.pheno,
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
cfw.pheno <- transform(cfw.pheno,fastglucose = as.double(fastglucose))

# Remove rows marked as "discard" and as possible sample mixups.
cfw.pheno <- subset(cfw.pheno,discard == "no" & mixup == "no")

# Create a new binary trait that is equal to 1 when the mouse has an
# "abnormal" bone-mineral density, and 0 otherwise.
cfw.pheno <- cbind(cfw.pheno,
                   data.frame(abBMD = factor(as.numeric(cfw.pheno$BMD > 90))))

# Create a binary trait from the p120b1 trait that indicates that
# the mouse does not respond to the pulse, possibly indicating
# deafness. I choose mice with p120b1 values on the "long tail" of
# the distribution for this phenotype.
cfw.pheno <-
  cbind(cfw.pheno,
        data.frame(deaf = factor(as.numeric(cfw.pheno$p120b1 < 50))))

# Remove most of the behavioural phenotype data.
cfw.pheno <- cfw.pheno[c("id","round","cageid","FCbox","PPIbox","methcage",
                         "methcycle","earpunch","glucoseage","methage",
                         "FCage","PPIage","sacage","bw0","bw1","bw2","bw3",
                         "PPIweight","sacweight","BMD","TA","EDL","gastroc",
                         "plantaris","soleus","tibia","abnormalbone",
                         "experimenters","testisweight","taillength",
                         "fastglucose","PreTrainD1","AvToneD1","AvContextD2",
                         "AvAltContextD3","AvToneD3","abBMD","deaf",
                         "pp3PPIavg","pp6PPIavg","pp12PPIavg")]

# Same the phenotype data table to an .RData file.
save(list = "cfw.pheno",file = "cfw.pheno.RData")
