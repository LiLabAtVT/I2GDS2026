# Step 0 (R): install packages not provided by the ARC R-bundle-Bioconductor/3.20 module
# Run with: module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2; Rscript 00_install_R_packages.R
lib <- Sys.getenv("R_LIBS_USER"); dir.create(lib, recursive=TRUE, showWarnings=FALSE); .libPaths(c(lib,.libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=8)
need <- function(p) !requireNamespace(p, quietly=TRUE)
for (p in c("phyloseq","microbiome")) if (need(p)) BiocManager::install(p, lib=lib, update=FALSE, ask=FALSE)
for (p in c("microeco","file2meco","remotes","ggalluvial","ggpubr")) if (need(p)) install.packages(p, lib=lib)
if (need("SpiecEasi")) remotes::install_github("zdk123/SpiecEasi@5f396da85baa114b31c13d9744c05387a1b04c23", lib=lib, upgrade="never")  # Oct-2023 commit (v1.1.3); newer versions need R>=4.5
pk <- c("phyloseq","vegan","microbiome","DirichletMultinomial","randomForest","microeco","file2meco","SpiecEasi","igraph")
print(sapply(pk, requireNamespace, quietly=TRUE))
