#!/bin/bash
#SBATCH --job-name=taxonomy
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code/logs/05_taxonomy_%j.log
# Step 5: Taxonomy with the QIIME2 pre-trained naive Bayes SILVA 138 (99%, full-length) classifier,
# as in the paper ("Silva 138 Reference Taxonomy Database", docs.qiime2.org data-resources).
# The paper describes no host filtering, so the unfiltered table stays primary; a table without
# mitochondria/chloroplast/unassigned ASVs is written alongside for comparison.
set -euo pipefail
B=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code
IN=$B/results/04_dada2; OUT=$B/results/05_taxonomy; REF=$B/reference
mkdir -p $OUT $REF
module load apptainer/1.4.0
SIF=/apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif
export TMPDIR=${TMPDIR:-/tmp}; export MPLCONFIGDIR=$TMPDIR/mpl
mkdir -p $HOME/qiime2/data_root
Q() { apptainer exec --bind $HOME/qiime2:/home/qiime2 --bind $HOME/qiime2/data_root:/data --bind /projects --bind $TMPDIR $SIF "$@"; }

# 1. Classifier (must match the container's scikit-learn version)
SKL=$(Q python -c "import sklearn; print(sklearn.__version__)"); echo "scikit-learn in container: $SKL"
CLF=$REF/silva-138-99-nb-classifier_sklearn-$SKL.qza
[ -s $CLF ] || wget -q -O $CLF https://data.qiime2.org/classifiers/sklearn-$SKL/silva/silva-138-99-nb-classifier.qza
md5sum $CLF

# 2. Classify ASVs
Q qiime feature-classifier classify-sklearn --i-classifier $CLF --i-reads $IN/rep-seqs.qza \
  --p-n-jobs 4 --o-classification $OUT/taxonomy.qza --verbose
Q qiime metadata tabulate --m-input-file $OUT/taxonomy.qza --o-visualization $OUT/taxonomy.qzv
Q qiime taxa barplot --i-table $IN/table.qza --i-taxonomy $OUT/taxonomy.qza \
  --m-metadata-file $B/raw_data/sample_metadata.tsv --o-visualization $OUT/taxa-barplot.qzv

# 3. Comparison table without host organelles / unassigned
Q qiime taxa filter-table --i-table $IN/table.qza --i-taxonomy $OUT/taxonomy.qza \
  --p-include d__Bacteria --p-exclude mitochondria,chloroplast --o-filtered-table $OUT/table-no-host.qza
Q qiime feature-table filter-seqs --i-data $IN/rep-seqs.qza --i-table $OUT/table-no-host.qza \
  --o-filtered-data $OUT/rep-seqs-no-host.qza
Q qiime tools export --input-path $OUT/taxonomy.qza --output-path $OUT/export

# 4. Report PDF
Q python $B/scripts/05_taxonomy_report.py
echo "Step 5 done"
