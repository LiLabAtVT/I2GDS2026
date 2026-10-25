#!/bin/bash
#SBATCH --job-name=fb_dmm
# Slurm account is not hard-coded: run  export SBATCH_ACCOUNT=<your_account>  once before sbatch
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/10_erwinia_dmm_%j.log

# Step 10: E. amylovora ASVs + DMM (paper Fig. 5B-D)
#   10a ASV numbering + BLAST (fast)  ->  10b Ea per sample + DMM k = 1..7 (7 fits in parallel)
# Submit from the project root directory:  sbatch scripts/10_erwinia_dmm.sh
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
[ -s results/10_erwinia_dmm/blast_top5.tsv ] || bash scripts/10a_asv_ids_blast.sh
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate fb_r
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=$CONDA_PREFIX/lib/R/library
Rscript scripts/10b_erwinia_dmm.R
echo "Step 10 done"
