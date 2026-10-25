#!/bin/bash
#SBATCH --job-name=fb_dada2
# Slurm account is not hard-coded: run  export SBATCH_ACCOUNT=<your_account>  once before sbatch
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/03_dada2_%j.log

# ---------------------------------------------------------------
# Step 3: QIIME2 (v2024.10) + DADA2 CCS mode (PacBio HiFi full-length 16S)
#   paper: "DADA2 was used within QIIME2 to remove quality control and
#           sequence noise and to detect and remove chimeras"
#   Parameters (not given in the paper) = standard PacBio HiFi settings
#   (Callahan et al. 2019 NAR; QIIME2 denoise-ccs):
#     primers 27F / 1492R, auto-orientation, length 1000-1600 bp,
#     maxEE 2, truncQ 2, consensus chimera removal, min fold parent 3.5
# Submit from the project root directory:  sbatch scripts/03_dada2_denoise.sh
# ---------------------------------------------------------------
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
OUT=results/03_dada2
mkdir -p $OUT results/pdf

module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate qiime2-2024.10
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
export R_LIBS_USER=$CONDA_PREFIX/lib/R/library
export TMPDIR=$PWD/tmp_q2; mkdir -p $TMPDIR
qiime --version | head -1

# ---- 3a import ----
[ -s $OUT/demux.qza ] || qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path metadata/manifest.tsv \
  --output-path $OUT/demux.qza

# ---- 3b DADA2 CCS (skipped if already done) ----
[ -s $OUT/table.qza ] || qiime dada2 denoise-ccs \
  --i-demultiplexed-seqs $OUT/demux.qza \
  --p-front AGRGTTYGATYMTGGCTCAG \
  --p-adapter RGYTACCTTGTTACGACTT \
  --p-min-len 1000 --p-max-len 1600 \
  --p-max-ee 2 --p-trunc-q 2 \
  --p-chimera-method consensus \
  --p-min-fold-parent-over-abundance 3.5 \
  --p-n-threads 32 \
  --o-table $OUT/table.qza \
  --o-representative-sequences $OUT/rep-seqs.qza \
  --o-denoising-stats $OUT/denoising-stats.qza \
  --verbose

# ---- 3c summaries, exports, plots (separate script so it can be re-run alone) ----
bash scripts/03b_dada2_summary.sh
rm -rf $TMPDIR
echo "Step 3 done"
