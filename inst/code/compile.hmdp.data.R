# Compile the HMDP phenotype data from the CSV file.

# Load the phenotype data from the HMDP study.
hmdp.pheno <- read.csv("../data/farber2011.TableS1.csv",header = TRUE,
                       stringsAsFactors = FALSE,comment.char = "#")

# Adjust the BMD measurement units, and convert "sex" into a factor.
hmdp.pheno <- transform(hmdp.pheno,
                        sex       = factor(sex,c("M","F")),
                        totalbody = 1000 * totalbody,
                        femur     = 1000 * femur,
                        spine     = 1000 * spine)

# Same the HMDP phenotype data table to an .RData file.
save(list = "hmdp.pheno",file = "hmdp.pheno.RData")
