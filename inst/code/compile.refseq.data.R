# Download and compile all (protein coding) gene annotations for NCBI
# release 38 of the Mouse Genome Assembly.
library(data.table)
library(curl)

# Download the feature_table.tar.gz file from the RefSeq database, and
# the accompanying README. This file contains annotation data for
# genes, RNAs and another genomic elements. 
cat("Downloading and processing RefSeq annotations.\n")
curl_download(paste("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq",
                    "vertebrate_mammalian/Mus_musculus",
                    "latest_assembly_versions/GCF_000001635.25_GRCm38.p5",
                    "README.txt",sep = "/"),"README.txt")
curl_download(paste("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq",
                    "vertebrate_mammalian/Mus_musculus",
                    "latest_assembly_versions/GCF_000001635.25_GRCm38.p5",
                    "GCF_000001635.25_GRCm38.p5_feature_table.txt.gz",sep="/"),
              "feature_table.txt.gz")
system(paste("gzcat feature_table.txt.gz | grep gene |",
             "grep protein_coding > gene_table.txt"))

# Read in the gene annotations from the tab-delimited text file. See
# the accompanying README.txt file for more information about this file.
cat("Reading RefSeq annotations.\n")
genes.m38        <- fread("gene_table.txt",sep = "\t",stringsAsFactors = FALSE)
class(genes.m38) <- "data.frame"
genes.m38        <- genes.m38[c(15,16,6,8,9)]
names(genes.m38) <- c("gene.symbol","gene.id","chr","start","end")
genes.m38        <- transform(genes.m38,chr = factor(chr,c(1:19,"X","Y","MT")))
                       
# Save the data to an .RData file.
save(list = "genes.m38",file = "genes.m38.RData")
