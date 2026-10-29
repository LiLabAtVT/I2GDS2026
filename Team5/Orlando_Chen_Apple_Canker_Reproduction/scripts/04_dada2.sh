#!/bin/bash
#SBATCH --job-name=dada2
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code/logs/04_dada2_%j.log
# Step 4: Import PacBio HiFi reads into QIIME2 2024.10 and denoise with DADA2 (denoise-ccs).
# Primers 27F/1492R (paper); length 1000-1600 bp, maxEE 2, consensus chimera removal (DADA2 CCS defaults).
set -euo pipefail
B=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code
OUT=$B/results/04_dada2; mkdir -p $OUT
module load apptainer/1.4.0
SIF=/apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif
export TMPDIR=${TMPDIR:-/tmp}; export MPLCONFIGDIR=$TMPDIR/mpl
mkdir -p $HOME/qiime2/data_root
Q() { apptainer exec --bind $HOME/qiime2:/home/qiime2 --bind $HOME/qiime2/data_root:/data --bind /projects --bind $TMPDIR $SIF "$@"; }
Q qiime --version

# 1. Manifest + import
{ echo -e "sample-id\tabsolute-filepath"
  for f in $B/raw_data/fastq/*.fastq.gz; do echo -e "$(basename $f .fastq.gz)\t$f"; done; } > $OUT/manifest.tsv
Q qiime tools import --type 'SampleData[SequencesWithQuality]' \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path $OUT/manifest.tsv --output-path $OUT/demux.qza
Q qiime demux summarize --i-data $OUT/demux.qza --o-visualization $OUT/demux.qzv

# 2. DADA2 CCS
Q qiime dada2 denoise-ccs --i-demultiplexed-seqs $OUT/demux.qza \
  --p-front AGRGTTYGATYMTGGCTCAG --p-adapter RGYTACCTTGTTACGACTT \
  --p-min-len 1000 --p-max-len 1600 --p-max-ee 2 \
  --p-chimera-method consensus --p-n-threads $SLURM_CPUS_PER_TASK \
  --o-table $OUT/table.qza --o-representative-sequences $OUT/rep-seqs.qza \
  --o-denoising-stats $OUT/denoising-stats.qza --verbose

# 3. Summaries + exports
Q qiime feature-table summarize --i-table $OUT/table.qza --m-sample-metadata-file $B/raw_data/sample_metadata.tsv --o-visualization $OUT/table.qzv
Q qiime feature-table tabulate-seqs --i-data $OUT/rep-seqs.qza --o-visualization $OUT/rep-seqs.qzv
Q qiime metadata tabulate --m-input-file $OUT/denoising-stats.qza --o-visualization $OUT/denoising-stats.qzv
mkdir -p $OUT/export
Q qiime tools export --input-path $OUT/denoising-stats.qza --output-path $OUT/export
Q qiime tools export --input-path $OUT/table.qza --output-path $OUT/export
Q qiime tools export --input-path $OUT/rep-seqs.qza --output-path $OUT/export
Q biom convert -i $OUT/export/feature-table.biom -o $OUT/export/feature-table.tsv --to-tsv

# 4. Report PDF
Q python $B/scripts/04_dada2_report.py
echo "Step 4 done"
