#!/bin/bash
#SBATCH --job-name=comp_ml
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code/logs/07_composition_ml_%j.log
# Step 7: composition, Random Forest, LEfSe, DMM, Erwinia ASVs (R env canker_r: R 4.3.3, randomForest 4.7-1.2,
# DirichletMultinomial 1.44.0, microbiomeMarker 1.8.0, phyloseq 1.46.0)
set -euo pipefail
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate canker_r
Rscript /projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code/scripts/07_composition_ml.R
