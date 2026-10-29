#!/bin/bash
#SBATCH --job-name=R_install
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code/logs/%x_%j.out
# Install R packages used by the paper that are not in the ARC R-bundle modules
set -euo pipefail
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
Rscript - <<'EOF'
lib <- .libPaths()[1]; dir.create(lib, recursive = TRUE, showWarnings = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib)
need <- function(p) !requireNamespace(p, quietly = TRUE)
if (need("microbiome")) BiocManager::install("microbiome", lib = lib, update = FALSE, ask = FALSE)
for (p in c("microeco", "file2meco", "remotes")) if (need(p)) install.packages(p, lib = lib)
if (need("SpiecEasi")) remotes::install_github("zdk123/SpiecEasi@v1.1.1", lib = lib, upgrade = "never")  # latest needs R>=4.5
pk <- c("microbiome", "microeco", "file2meco", "SpiecEasi")
print(sapply(pk, requireNamespace, quietly = TRUE))
for (p in pk) cat(p, as.character(packageVersion(p)), "\n")
EOF
