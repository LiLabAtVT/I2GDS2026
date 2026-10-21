#!/bin/bash
# Step 4: QIIME2 2024.10 -- import, DADA2 (CCS mode) denoising + chimera removal,
#         SILVA 138 taxonomy, phylogeny (for Faith PD), and export of tables for R.
set -euo pipefail
source $HOME/fireblight_canker_reanalysis/scripts/qiime_env.sh
Q=$P/qiime; cd $Q

# 4a. manifest (sample-id -> fastq) from metadata
awk -F'\t' -v d=$P/data/raw 'BEGIN{OFS="\t";print "sample-id","absolute-filepath"} NR>1{print $1, d"/"$2".fastq.gz"}' \
  $P/metadata/sample_metadata.tsv > manifest.tsv
[ -f demux.qza ] || qiime tools import --type 'SampleData[SequencesWithQuality]' \
  --input-format SingleEndFastqManifestPhred33V2 --input-path manifest.tsv --output-path demux.qza
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv

# 4b. DADA2 denoise-ccs. Primers from the paper (27F / 1492R).
#     Length/EE/quality cutoffs are not reported in the paper -> DADA2 PacBio HiFi defaults
#     (Callahan et al. 2019): keep 1000-1600 bp, maxEE=2, minQ=0 (trunc-q 0).
[ -f table.qza ] || qiime dada2 denoise-ccs --i-demultiplexed-seqs demux.qza \
  --p-front AGRGTTYGATYMTGGCTCAG --p-adapter RGYTACCTTGTTACGACTT \
  --p-min-len 1000 --p-max-len 1600 --p-max-ee 2 --p-trunc-q 0 \
  --p-n-threads 8 \
  --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats dada2-stats.qza
qiime metadata tabulate --m-input-file dada2-stats.qza --o-visualization dada2-stats.qzv

# 4c. Taxonomy: SILVA 138 99% full-length Naive Bayes classifier (paper: Silva 138)
CLF=$P/refs/silva-138-99-nb-classifier.qza
[ -f $CLF ] || curl -sL -o $CLF https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza
[ -f taxonomy.qza ] || qiime feature-classifier classify-sklearn --i-classifier $CLF \
  --i-reads rep-seqs.qza --p-n-jobs 2 --p-reads-per-batch 2000 --o-classification taxonomy.qza
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza \
  --m-metadata-file $P/metadata/sample_metadata.tsv --o-visualization taxa-barplot.qzv

# 4d. Phylogeny (needed for Faith's phylogenetic diversity, reported in the paper)
[ -f rooted-tree.qza ] || qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza \
  --p-n-threads 8 --o-alignment aligned.qza --o-masked-alignment masked.qza \
  --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

# 4e. Export for downstream R analyses
rm -rf export && mkdir export
for a in table taxonomy rep-seqs rooted-tree dada2-stats; do qiime tools export --input-path $a.qza --output-path export/$a; done
q2 biom convert -i export/table/feature-table.biom -o export/feature-table.tsv --to-tsv
echo "STEP4 DONE"
