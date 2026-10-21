#!/bin/bash
# Reproduce the whole reanalysis (run on a compute node, e.g. inside an ARC interactive job):
#   bash ~/fireblight_canker_reanalysis/scripts/run_all.sh
set -euo pipefail
P=$HOME/fireblight_canker_reanalysis; S=$P/scripts; L=$P/logs
bash $S/01_get_metadata.sh                > $L/01_metadata.log 2>&1
bash $S/02_download_raw.sh                > $L/02_download.log 2>&1
bash $S/03_fastqc.sh                      > $L/03_fastqc.log 2>&1
bash $S/04_qiime2_dada2_taxonomy.sh       > $L/04_qiime2.log 2>&1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1   # multi-threaded BLAS on a 128-core node with 8 allocated CPUs makes SparCC ~300x slower
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
Rscript $S/00_install_R_packages.R        > $L/00_install_R.log 2>&1
for s in 05_build_phyloseq_qc 06_composition_alpha_beta 07_taxa_RF_LEfSe_DMM 08_FAPROTAX 09_SparCC_networks; do
  Rscript $S/$s.R > $L/$s.log 2>&1
done
bash $S/10_finalize_README.sh
