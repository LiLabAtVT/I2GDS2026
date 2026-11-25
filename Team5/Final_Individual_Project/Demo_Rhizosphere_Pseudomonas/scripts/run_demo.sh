#!/bin/bash
#SBATCH --job-name=soy_demo
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=run_demo_%j.out
# Runs Part 1 (differential abundance) then Part 2 (Pseudomonas variants). From this folder:
#   sbatch --account=<your_allocation> scripts/run_demo.sh      (on ARC; or plain `bash scripts/run_demo.sh` elsewhere)
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ./env
NC=${SLURM_CPUS_PER_TASK:-4}
echo "== Part 1: differential abundance ($(date))"
Rscript scripts/01_differential_abundance.R input_dataset output $NC
echo "== Part 2: Pseudomonas variants ($(date))"
Rscript scripts/02_pseudomonas_variants.R input_dataset output $NC
echo "== done ($(date)); compare output/ with expected_output/"
