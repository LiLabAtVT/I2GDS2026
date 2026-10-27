#!/bin/bash
# Reproduce everything: download on the login node, then chain the Slurm jobs.
set -euo pipefail
cd /projects/intro2gds/I2GDS2026/Team1/Kamal/Claude_Code
bash scripts/01_download.sh
sbatch scripts/02_fastqc.sh
j3=$(sbatch --parsable scripts/03_qiime2.sh)
sbatch --dependency=afterok:$j3 scripts/04_05_R.sh
