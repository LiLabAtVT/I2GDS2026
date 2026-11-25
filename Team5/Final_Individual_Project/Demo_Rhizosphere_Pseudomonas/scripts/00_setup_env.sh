#!/bin/bash
#SBATCH --job-name=demo_env
#SBATCH --time=01:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=demo_env_%j.out
# Creates the conda environment ./env (inside this folder) from environment.yml. Run once:
#   sbatch --account=<your_allocation> scripts/00_setup_env.sh      (on ARC; or plain `bash scripts/00_setup_env.sh` elsewhere)
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
source "$(conda info --base)/etc/profile.d/conda.sh"
command -v mamba >/dev/null && SOLVER=mamba || SOLVER=conda
$SOLVER env create -p ./env -f environment.yml -y
conda activate ./env
Rscript -e 'for (p in c("ANCOMBC","DESeq2","CVXR","vegan","ggplot2")) cat(p, as.character(packageVersion(p)), "\n")'
blastn -version | head -1
