#!/bin/bash
#SBATCH --job-name=fb_taxonomy
# Slurm account is not hard-coded: run  export SBATCH_ACCOUNT=<your_account>  once before sbatch
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=logs/04_taxonomy_%j.log

# ---------------------------------------------------------------
# Step 4: taxonomy (paper: "Sequences were assigned to species using the
#         Silva 138 Reference Taxonomy Database"), host/organelle filtering,
#         phylogenetic tree (for Faith PD), class-level composition (Fig. 2A)
#   classifier: ref/silva-138-99-nb-classifier.qza (QIIME2 data resources,
#               sklearn 1.4.2 build for QIIME2 2024.10)
# Submit from the project root directory:  sbatch scripts/04_taxonomy_tree.sh
# ---------------------------------------------------------------
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
IN=results/03_dada2
OUT=results/04_taxonomy
mkdir -p $OUT/export results/pdf

module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate qiime2-2024.10
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
export R_LIBS_USER=$CONDA_PREFIX/lib/R/library
export TMPDIR=$PWD/tmp_q2; mkdir -p $TMPDIR

# ---- 4a classify (Naive Bayes, SILVA 138 99%) ----
[ -s $OUT/taxonomy.qza ] || qiime feature-classifier classify-sklearn \
  --i-classifier ref/silva-138-99-nb-classifier.qza \
  --i-reads $IN/rep-seqs.qza \
  --p-n-jobs 4 \
  --o-classification $OUT/taxonomy.qza
qiime metadata tabulate --m-input-file $OUT/taxonomy.qza --o-visualization $OUT/taxonomy.qzv

# ---- 4b remove mitochondria, chloroplast, unassigned / no phylum ----
qiime taxa filter-table \
  --i-table $IN/table.qza --i-taxonomy $OUT/taxonomy.qza \
  --p-include p__ --p-exclude mitochondria,chloroplast \
  --o-filtered-table $OUT/table-filtered.qza
qiime taxa filter-seqs \
  --i-sequences $IN/rep-seqs.qza --i-taxonomy $OUT/taxonomy.qza \
  --p-include p__ --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences $OUT/rep-seqs-filtered.qza
qiime feature-table summarize --i-table $OUT/table-filtered.qza \
  --m-sample-metadata-file metadata/metadata.tsv --o-visualization $OUT/table-filtered.qzv

# ---- 4c phylogenetic tree (MAFFT + FastTree, midpoint rooted) ----
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $OUT/rep-seqs-filtered.qza --p-n-threads 16 \
  --o-alignment $OUT/aligned-rep-seqs.qza --o-masked-alignment $OUT/masked-aligned-rep-seqs.qza \
  --o-tree $OUT/unrooted-tree.qza --o-rooted-tree $OUT/rooted-tree.qza

# ---- 4d barplot + class-level table ----
qiime taxa barplot --i-table $OUT/table-filtered.qza --i-taxonomy $OUT/taxonomy.qza \
  --m-metadata-file metadata/metadata.tsv --o-visualization $OUT/taxa-barplot.qzv

# ---- exports for R ----
qiime tools export --input-path $OUT/taxonomy.qza          --output-path $OUT/export
qiime tools export --input-path $OUT/table-filtered.qza    --output-path $OUT/export
qiime tools export --input-path $OUT/rep-seqs-filtered.qza --output-path $OUT/export
qiime tools export --input-path $OUT/rooted-tree.qza       --output-path $OUT/export
biom convert -i $OUT/export/feature-table.biom -o $OUT/export/feature-table.tsv --to-tsv
# unfiltered table too, to report what was removed
qiime tools export --input-path $IN/table.qza --output-path $OUT/export/unfiltered
biom convert -i $OUT/export/unfiltered/feature-table.biom -o $OUT/export/unfiltered/feature-table.tsv --to-tsv

Rscript scripts/04_taxonomy_plots.R
rm -rf $TMPDIR
echo "Step 4 done"
