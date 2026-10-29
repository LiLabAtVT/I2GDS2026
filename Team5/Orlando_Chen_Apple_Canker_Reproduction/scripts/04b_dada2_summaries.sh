#!/bin/bash
#SBATCH --job-name=dada2_sum
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code/logs/04b_dada2_summaries_%j.log
# Step 4b: summaries, exports and report PDF from the finished DADA2 run (job 7773425).
# Same as section 3-4 of 04_dada2.sh; rerun after fixing the metadata ID header (sample_id -> sample-id).
set -euo pipefail
B=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code
OUT=$B/results/04_dada2
module load apptainer/1.4.0
SIF=/apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif
export TMPDIR=${TMPDIR:-/tmp}; export MPLCONFIGDIR=$TMPDIR/mpl
mkdir -p $HOME/qiime2/data_root
Q() { apptainer exec --bind $HOME/qiime2:/home/qiime2 --bind $HOME/qiime2/data_root:/data --bind /projects --bind $TMPDIR $SIF "$@"; }

Q qiime feature-table summarize --i-table $OUT/table.qza --m-sample-metadata-file $B/raw_data/sample_metadata.tsv --o-visualization $OUT/table.qzv
Q qiime feature-table tabulate-seqs --i-data $OUT/rep-seqs.qza --o-visualization $OUT/rep-seqs.qzv
Q qiime metadata tabulate --m-input-file $OUT/denoising-stats.qza --o-visualization $OUT/denoising-stats.qzv
mkdir -p $OUT/export
Q qiime tools export --input-path $OUT/denoising-stats.qza --output-path $OUT/export
Q qiime tools export --input-path $OUT/table.qza --output-path $OUT/export
Q qiime tools export --input-path $OUT/rep-seqs.qza --output-path $OUT/export
Q biom convert -i $OUT/export/feature-table.biom -o $OUT/export/feature-table.tsv --to-tsv

Q python $B/scripts/04_dada2_report.py
echo "Step 4b done"
