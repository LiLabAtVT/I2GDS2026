#!/bin/bash
#SBATCH -J r_stats -A introtogds -p normal_q -t 12:00:00 -c 16 --mem=64G -o logs/04_05_R_%j.log
# Steps 4-5: R statistics and SparCC networks
set -euo pipefail
cd /projects/intro2gds/I2GDS2026/Team1/Kamal/Claude_Code
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
Rscript scripts/04_stats.R
Rscript scripts/05_network.R
