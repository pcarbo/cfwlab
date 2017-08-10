# cfwlab

R package containing data from the
[CFW genome-wide association study](https://dx.doi.org/10.1038/ng.3609),
designed for use in courses. These data were compiled from the
[Data Dryad repository](http://dx.doi.org/10.5061/dryad.2rs41).

## Quick Start

1. Install the cfwlab package in R or RStudio. The easiest way to
   install the package is with devtools:

   ```R
   install.packages("devtools")
   library(devtools)
   install_github("pcarbo/cfwlab")
   ```

2. Load the package into your R environment:

   ```R
   library(cfwlab)
   ```

3. Load and inspect the data sets, e.g.,

   ```R
   data(cfw.pheno)
   head(cfw.pheno)
   ```

4. Read the package documentation:

   ```R
   help(package = cfwlab)
   ```

## Credits

The cfwlab package was developed by:

[Peter Carbonetto](http://pcarbo.github.io) &
[John Novembre](http://jnpopgen.org)<br>
Department of Human Genetics<br>
University of Chicago
