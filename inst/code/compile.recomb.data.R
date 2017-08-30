# Compile the recombination rate data from Brunschwig et al (2012).

# Load the recombination rate data from Brunschwig et al (2012),
# "lifted over" to the Mouse Genome Assembly 38 (mm10).
recomb.m38 <- read.table("../data/brunschwig_tableS1_mm10.txt.gz",sep = " ",
                         stringsAsFactors = TRUE,header = TRUE)
recomb.m38 <- transform(recomb.m38,chr = factor(chr))

# Save the data to an .RData file.
save(list = "recomb.m38", file = "recomb.m38.RData")
