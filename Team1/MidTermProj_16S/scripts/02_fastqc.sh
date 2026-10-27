#!/bin/bash
#SBATCH -J fastqc -A introtogds -p normal_q -t 2:00:00 -c 8 --mem=16G -o logs/02_fastqc_%j.log
# Step 2: read quality check (paper: FastQC, Andrews 2010)
set -euo pipefail
cd /projects/intro2gds/I2GDS2026/Team1/Kamal/Claude_Code
module load FastQC/0.12.1-Java-11
mkdir -p results/fastqc
fastqc -t 8 -o results/fastqc raw_data/*.fastq.gz
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
Rscript scripts/02b_fastqc_summary.R
