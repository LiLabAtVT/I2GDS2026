#!/bin/bash
#SBATCH --job-name=fb_dada2_sum
# Slurm account is not hard-coded: run  export SBATCH_ACCOUNT=<your_account>  once before sbatch
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/03b_dada2_summary_%j.log

# Step 3c: summaries + exports + plots of the DADA2 output
# Called at the end of 03_dada2_denoise.sh, or alone:  sbatch scripts/03b_dada2_summary.sh
set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
OUT=results/03_dada2
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate qiime2-2024.10
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
export R_LIBS_USER=$CONDA_PREFIX/lib/R/library

# ---- 3c summaries ----
qiime metadata tabulate --m-input-file $OUT/denoising-stats.qza --o-visualization $OUT/denoising-stats.qzv
qiime feature-table summarize --i-table $OUT/table.qza --m-sample-metadata-file metadata/metadata.tsv --o-visualization $OUT/table.qzv
qiime feature-table tabulate-seqs --i-data $OUT/rep-seqs.qza --o-visualization $OUT/rep-seqs.qzv

qiime tools export --input-path $OUT/denoising-stats.qza --output-path $OUT/export
qiime tools export --input-path $OUT/table.qza           --output-path $OUT/export
qiime tools export --input-path $OUT/rep-seqs.qza        --output-path $OUT/export
biom convert -i $OUT/export/feature-table.biom -o $OUT/export/feature-table.tsv --to-tsv

Rscript scripts/03_dada2_plots.R
echo "Step 3c done"
