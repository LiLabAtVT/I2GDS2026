#!/bin/bash
#SBATCH --job-name=q2_dada2
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code/logs/%x_%j.out
# Step 5: QIIME2 v2024.10 (same version as the paper)
#   import -> DADA2 denoise-ccs (primer removal, denoising, chimera removal)
#   -> SILVA 138 taxonomy -> remove host organelles -> phylogeny -> export for R
set -euo pipefail
BASE=/projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code
cd "$BASE"
module load QIIME2/2024.10-amplicon
export TMPDIR=${TMPDIR:-/tmp}
# on ARC 'qiime' is a function running an Apptainer image; run biom from the same image
biom() { apptainer exec --bind /projects --bind "$TMPDIR" \
  /apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif biom "$@"; }
NT=$SLURM_CPUS_PER_TASK
Q=results/02_qiime2
mkdir -p $Q raw_data/reference

# --- manifest ---
{
  printf "sample-id\tabsolute-filepath\n"
  tail -n +3 metadata/sample_metadata.tsv | cut -f1 | while read -r s; do
    printf "%s\t%s/raw_data/fastq/%s.fastq.gz\n" "$s" "$BASE" "$s"
  done
} > $Q/manifest.tsv

# --- import ---
[[ -s $Q/demux.qza ]] || qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path $Q/manifest.tsv --output-path $Q/demux.qza
qiime demux summarize --i-data $Q/demux.qza --o-visualization $Q/demux.qzv

# --- DADA2 for PacBio CCS/HiFi; primers 27F / 1492R from the paper ---
[[ -s $Q/table_raw.qza ]] || qiime dada2 denoise-ccs \
  --i-demultiplexed-seqs $Q/demux.qza \
  --p-front AGRGTTYGATYMTGGCTCAG \
  --p-adapter RGYTACCTTGTTACGACTT \
  --p-min-len 1000 --p-max-len 1600 \
  --p-trunc-len 0 --p-max-ee 2 \
  --p-chimera-method consensus \
  --p-n-threads $NT \
  --o-table $Q/table_raw.qza \
  --o-representative-sequences $Q/rep_seqs_raw.qza \
  --o-denoising-stats $Q/dada2_stats.qza
qiime metadata tabulate --m-input-file $Q/dada2_stats.qza --o-visualization $Q/dada2_stats.qzv

# --- SILVA 138 full-length 16S Naive Bayes classifier (QIIME2 2024.10 = scikit-learn 1.4.2) ---
CLS=raw_data/reference/silva-138-99-nb-classifier.qza
[[ -s $CLS ]] || wget -q -O $CLS https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza \
  || { rm -f $CLS; echo "classifier download failed"; exit 1; }
[[ -s $Q/taxonomy.qza ]] || qiime feature-classifier classify-sklearn \
  --i-classifier $CLS --i-reads $Q/rep_seqs_raw.qza \
  --p-n-jobs 8 --p-reads-per-batch 2000 \
  --o-classification $Q/taxonomy.qza

# --- remove host-derived chloroplast / mitochondria and unassigned-at-phylum reads ---
qiime taxa filter-table --i-table $Q/table_raw.qza --i-taxonomy $Q/taxonomy.qza \
  --p-include p__ --p-exclude mitochondria,chloroplast \
  --o-filtered-table $Q/table.qza
qiime feature-table filter-seqs --i-data $Q/rep_seqs_raw.qza --i-table $Q/table.qza \
  --o-filtered-data $Q/rep_seqs.qza
qiime feature-table summarize --i-table $Q/table.qza \
  --m-sample-metadata-file metadata/sample_metadata.tsv --o-visualization $Q/table.qzv

# --- phylogeny (for Faith's phylogenetic diversity) ---
[[ -s $Q/rooted_tree.qza ]] || qiime phylogeny align-to-tree-mafft-fasttree --i-sequences $Q/rep_seqs.qza \
  --p-n-threads $NT \
  --o-alignment $Q/aligned.qza --o-masked-alignment $Q/masked.qza \
  --o-tree $Q/unrooted_tree.qza --o-rooted-tree $Q/rooted_tree.qza

# --- export plain-text files for R ---
E=$Q/export; rm -rf $E; mkdir -p $E
for a in table_raw table taxonomy rep_seqs rooted_tree dada2_stats; do
  qiime tools export --input-path $Q/$a.qza --output-path $E/$a
done
biom convert -i $E/table/feature-table.biom -o $E/feature_table.tsv --to-tsv
biom convert -i $E/table_raw/feature-table.biom -o $E/feature_table_raw.tsv --to-tsv

# per-sample read counts after organelle filtering (used to pick the rarefaction depth)
qiime tools export --input-path $Q/table.qzv --output-path $E/table_summary
biom summarize-table -i $E/table/feature-table.biom -o $E/table_summary/summary.txt
echo "== Reads per sample after filtering:"; cat $E/table_summary/summary.txt
