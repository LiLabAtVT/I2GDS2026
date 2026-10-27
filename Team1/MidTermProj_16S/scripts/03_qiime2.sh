#!/bin/bash
#SBATCH -J qiime2 -A introtogds -p normal_q -t 24:00:00 -c 32 --mem=128G -o logs/03_qiime2_%j.log
# Step 3: QIIME2 2024.10 — import, DADA2 (CCS) denoise + chimera removal, SILVA 138 taxonomy,
# tree, rarefied core diversity (Chao1, observed, Shannon, Simpson, Faith PD, Bray-Curtis PCoA).
set -euo pipefail
cd /projects/intro2gds/I2GDS2026/Team1/Kamal/Claude_Code
module load QIIME2/2024.10-amplicon
export TMPDIR=${TMPDIR:-/localscratch/$USER}
Q=results/qiime2; mkdir -p $Q refs
M=metadata/sample-metadata.tsv

[ -s $Q/demux.qza ] || qiime tools import --type 'SampleData[SequencesWithQuality]' --input-format SingleEndFastqManifestPhred33V2 \
  --input-path metadata/manifest.tsv --output-path $Q/demux.qza

# Paper gives no DADA2 params; these are the standard full-length PacBio settings (Callahan et al. 2019).
# 27F front, reverse-complement of 1492R as 3' adapter.
qiime dada2 denoise-ccs --i-demultiplexed-seqs $Q/demux.qza \
  --p-front AGRGTTYGATYMTGGCTCAG --p-adapter AAGTCGTAACAAGGTARCY \
  --p-min-len 1000 --p-max-len 1600 --p-n-threads 32 --verbose \
  --o-table $Q/table.qza --o-representative-sequences $Q/rep-seqs.qza --o-denoising-stats $Q/dada2-stats.qza

# SILVA 138 (99%, full-length) Naive Bayes classifier for the container's sklearn 1.4.2
[ -s refs/silva-138-99-nb-classifier.qza ] || curl -sSL -o refs/silva-138-99-nb-classifier.qza \
  https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza
qiime feature-classifier classify-sklearn --i-classifier refs/silva-138-99-nb-classifier.qza \
  --i-reads $Q/rep-seqs.qza --p-n-jobs 8 --o-classification $Q/taxonomy.qza

qiime phylogeny align-to-tree-mafft-fasttree --i-sequences $Q/rep-seqs.qza --p-n-threads 32 \
  --o-alignment $Q/aligned.qza --o-masked-alignment $Q/masked.qza --o-tree $Q/unrooted.qza --o-rooted-tree $Q/rooted-tree.qza

# Rarefy to the shallowest sample so all 48 samples are kept (paper does not state the depth).
qiime tools export --input-path $Q/table.qza --output-path $Q/export
apptainer exec /apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif \
  biom summarize-table -i $Q/export/feature-table.biom -o $Q/table-summary.txt
DEPTH=$(awk '/^ Min:/{print int($2)}' $Q/table-summary.txt); echo "rarefaction depth: $DEPTH"
rm -rf $Q/core
qiime diversity core-metrics-phylogenetic --i-table $Q/table.qza --i-phylogeny $Q/rooted-tree.qza \
  --p-sampling-depth $DEPTH --m-metadata-file $M --output-dir $Q/core
for m in chao1 simpson; do
  qiime diversity alpha --i-table $Q/core/rarefied_table.qza --p-metric $m --o-alpha-diversity $Q/core/${m}_vector.qza
done

# Export for R
for a in dada2-stats taxonomy rep-seqs rooted-tree core/rarefied_table; do
  qiime tools export --input-path $Q/$a.qza --output-path $Q/export/$(basename $a)
done
for v in observed_features shannon simpson chao1 faith_pd evenness; do
  qiime tools export --input-path $Q/core/${v}_vector.qza --output-path $Q/export/alpha_$v
done
apptainer exec /apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif \
  biom convert -i $Q/export/rarefied_table/feature-table.biom -o $Q/export/rarefied_table.tsv --to-tsv
apptainer exec /apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif \
  biom convert -i $Q/export/feature-table.biom -o $Q/export/feature-table.tsv --to-tsv
echo done
