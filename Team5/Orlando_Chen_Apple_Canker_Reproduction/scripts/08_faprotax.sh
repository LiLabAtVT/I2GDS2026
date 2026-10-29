#!/bin/bash
#SBATCH --job-name=faprotax
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code/logs/08_faprotax_%j.log
# Step 8: FAPROTAX functional prediction (R env canker_r, microeco 1.1.0)
set -euo pipefail
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate canker_r
Rscript /projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code/scripts/08_faprotax.R
