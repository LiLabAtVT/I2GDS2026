#!/bin/bash
# Check which R packages needed for the downstream analysis are available on ARC
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2 >/dev/null 2>&1
module list 2>&1 | tail -3
Rscript -e '
pk <- c("phyloseq","vegan","microbiome","DirichletMultinomial","randomForest","microeco",
        "file2meco","SpiecEasi","igraph","ggplot2","pheatmap","agricolae","reshape2",
        "magrittr","dplyr","tidyr","ape","ggrepel","patchwork","lefser","Maaslin2")
ok <- sapply(pk, requireNamespace, quietly = TRUE)
print(ok); cat("R", R.version.string, "\nlib:", .libPaths(), "\n")'
