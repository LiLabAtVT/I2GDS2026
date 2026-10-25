#!/bin/bash
#SBATCH --job-name=fb_alpha
# Slurm account is not hard-coded: run  export SBATCH_ACCOUNT=<your_account>  once before sbatch
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/05_alpha_%j.log

# ---------------------------------------------------------------
# Step 5: rarefaction + alpha diversity in QIIME2 (paper Fig. 2B)
#   paper: "ASV abundances were rarefied to the same sequencing depth";
#          "Shannon, Simpson and Chao1 were calculated using QIIME2";
#          ANOVA in R (F test, P < 0.05)
#   Two tables (user decision after Step 4):
#     paper    = all ASVs (incl. chloroplast/mitochondria), depth 14,881 (= min sample, keeps 48)
#     filtered = bacteria only,                              depth 505    (keeps 32 cankers, 0 healthy)
#   core-metrics-phylogenetic also writes the Bray-Curtis PCoA used in Step 6.
# Submit from the project root directory:  sbatch scripts/05_alpha_diversity.sh
# ---------------------------------------------------------------
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
T3=results/03_dada2
T4=results/04_taxonomy
OUT=results/05_alpha
mkdir -p $OUT results/pdf

module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate qiime2-2024.10
export TMPDIR=$PWD/tmp_q2; mkdir -p $TMPDIR

# tree with ALL ASVs (needed for Faith PD on the paper-style table)
[ -s $T4/rooted-tree-all.qza ] || qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $T3/rep-seqs.qza --p-n-threads 16 \
  --o-alignment $T4/aligned-rep-seqs-all.qza --o-masked-alignment $T4/masked-aligned-rep-seqs-all.qza \
  --o-tree $T4/unrooted-tree-all.qza --o-rooted-tree $T4/rooted-tree-all.qza

run_set () {   # $1 name  $2 table  $3 tree  $4 depth
  local D=$OUT/$1
  rm -rf $D
  qiime diversity core-metrics-phylogenetic \
    --i-table $2 --i-phylogeny $3 --p-sampling-depth $4 \
    --m-metadata-file metadata/metadata.tsv --p-n-jobs-or-threads 16 \
    --output-dir $D
  for m in chao1 simpson; do
    qiime diversity alpha --i-table $D/rarefied_table.qza --p-metric $m --o-alpha-diversity $D/${m}_vector.qza
  done
  mkdir -p $D/export
  for v in observed_features shannon chao1 simpson faith_pd evenness; do
    qiime tools export --input-path $D/${v}_vector.qza --output-path $D/export/$v
    qiime diversity alpha-group-significance --i-alpha-diversity $D/${v}_vector.qza \
      --m-metadata-file metadata/metadata.tsv --o-visualization $D/${v}_group_significance.qzv
  done
  qiime diversity alpha-rarefaction --i-table $2 --i-phylogeny $3 --p-max-depth $5 \
    --p-metrics observed_features --p-metrics shannon \
    --m-metadata-file metadata/metadata.tsv --o-visualization $D/alpha_rarefaction.qzv
  qiime tools export --input-path $D/alpha_rarefaction.qzv --output-path $D/export/alpha_rarefaction
  qiime tools export --input-path $D/rarefied_table.qza --output-path $D/export
  biom convert -i $D/export/feature-table.biom -o $D/export/rarefied_table.tsv --to-tsv
}
run_set paper    $T3/table.qza          $T4/rooted-tree-all.qza 14881 14881
run_set filtered $T4/table-filtered.qza $T4/rooted-tree.qza     505   5000
conda deactivate

conda activate fb_r
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=$CONDA_PREFIX/lib/R/library
Rscript scripts/05_alpha_plots.R
rm -rf $TMPDIR
echo "Step 5 done"
