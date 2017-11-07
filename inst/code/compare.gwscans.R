# A short script to compare results from two different runs of the R
# script, compile.gwscan.R. This could be used, for example, to
# compare the results of GEMMA using two different versions of the GNU
# Scientific Library (GSL).

# LOAD RESULTS
# ------------
# Load the results from running compile.gwscans.R in two different
# settings (e.g., using two different versions of the GEMMA software).
cat("Loading results from running compile.gwscan.R.\n")
load("cfw.gwscan1.RData")
gwscan1 <- cfw.gwscan
load("cfw.gwscan2.RData")
gwscan2 <- cfw.gwscan
rm(cfw.gwscan)

# COMPARE RESULTS
# ---------------
# Check that all the SNPs are the same.
cat("Comparing results.\n")
if (all(gwscan1$id  == gwscan2$id) &
    all(gwscan1$chr == gwscan2$chr) &
    all(gwscan1$pos == gwscan2$pos)) {
  cat("All SNPs match.\n")
} else {
  stop("SNPs do not match")
}

# Compare the (log10) p-values.
cat("Largest difference in log10 p-values for each phenotype:\n")
print(t(t(apply(abs(gwscan1[-(1:3)] - gwscan2[-(1:3)]),2,max))))
