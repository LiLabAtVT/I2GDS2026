#!/bin/bash
#SBATCH --job-name=diversity
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code/logs/06_diversity_%j.log
# Step 6: rarefaction, alpha/beta diversity and PERMANOVA (paper: rarefied in QIIME2; Chao1, observed,
# Shannon, Simpson, phylogenetic diversity; Bray-Curtis PCoA; PERMANOVA health status x cultivar).
#   A = paper-like: unfiltered table (chloroplast kept), all 48 samples, depth = smallest sample (14,881)
#   B = sensitivity: bacteria only (no chloroplast/mitochondria), 32 canker samples, depth = smallest (505)
set -euo pipefail
B=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code
META=$B/raw_data/sample_metadata.tsv
D4=$B/results/04_dada2; D5=$B/results/05_taxonomy; OUT=$B/results/06_diversity
mkdir -p $OUT
module load apptainer/1.4.0
SIF=/apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif
export TMPDIR=${TMPDIR:-/tmp}; export MPLCONFIGDIR=$TMPDIR/mpl
mkdir -p $HOME/qiime2/data_root
Q() { apptainer exec --bind $HOME/qiime2:/home/qiime2 --bind $HOME/qiime2/data_root:/data --bind /projects --bind $TMPDIR $SIF "$@"; }
N=$SLURM_CPUS_PER_TASK

# Phylogeny for Faith's PD (all ASVs; subsets of the table use the same tree)
Q qiime phylogeny align-to-tree-mafft-fasttree --i-sequences $D4/rep-seqs.qza --p-n-threads $N \
  --o-alignment $OUT/aligned-rep-seqs.qza --o-masked-alignment $OUT/masked-aligned-rep-seqs.qza \
  --o-tree $OUT/unrooted-tree.qza --o-rooted-tree $OUT/rooted-tree.qza

# Bacteria-only canker table for analysis B
Q qiime feature-table filter-samples --i-table $D5/table-no-host.qza --m-metadata-file $META \
  --p-where "[health_status]='Canker'" --o-filtered-table $OUT/table-B-canker-bacteria.qza

run_set() {  # $1 = label, $2 = table, $3 = depth, $4 = PERMANOVA formula
  local L=$1 T=$2 DEPTH=$3 F=$4 C=$OUT/core-$1
  rm -rf $C
  Q qiime diversity core-metrics-phylogenetic --i-table $T --i-phylogeny $OUT/rooted-tree.qza \
    --p-sampling-depth $DEPTH --m-metadata-file $META --p-n-jobs-or-threads $N --output-dir $C
  for m in chao1 simpson; do
    Q qiime diversity-lib alpha-passthrough --i-table $C/rarefied_table.qza --p-metric $m --o-vector $C/${m}_vector.qza
  done
  mkdir -p $C/export
  for v in observed_features shannon evenness faith_pd chao1 simpson; do
    Q qiime diversity alpha-group-significance --i-alpha-diversity $C/${v}_vector.qza --m-metadata-file $META \
      --o-visualization $C/${v}-group-significance.qzv
    Q qiime tools export --input-path $C/${v}_vector.qza --output-path $C/export/$v
  done
  Q qiime tools export --input-path $C/bray_curtis_pcoa_results.qza --output-path $C/export/bray_curtis_pcoa
  Q qiime tools export --input-path $C/bray_curtis_distance_matrix.qza --output-path $C/export/bray_curtis_dm
  Q qiime diversity adonis --i-distance-matrix $C/bray_curtis_distance_matrix.qza --m-metadata-file $META \
    --p-formula "$F" --p-permutations 999 --p-n-jobs $N --o-visualization $C/adonis_bray_curtis.qzv
  Q qiime tools export --input-path $C/adonis_bray_curtis.qzv --output-path $C/export/adonis
  for col in health_status cultivar; do
    [ "$L" = B ] && [ $col = health_status ] && continue
    Q qiime diversity beta-group-significance --i-distance-matrix $C/bray_curtis_distance_matrix.qza \
      --m-metadata-file $META --m-metadata-column $col --p-pairwise --p-permutations 999 \
      --o-visualization $C/bray_curtis_${col}_significance.qzv
    Q qiime tools export --input-path $C/bray_curtis_${col}_significance.qzv --output-path $C/export/permanova_$col
  done
}
run_set A $D4/table.qza 14881 "health_status*cultivar"
run_set B $OUT/table-B-canker-bacteria.qza 505 "cultivar"

Q python $B/scripts/06_diversity_report.py
echo "Step 6 done"
