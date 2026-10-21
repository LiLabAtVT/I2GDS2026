#!/bin/bash
# Run step 9 (SparCC networks) + README tables detached from the login session; safe across VPN/ssh disconnects.
P=$HOME/fireblight_canker_reanalysis
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1   # multi-threaded BLAS on a 128-core node with 8 allocated CPUs makes SparCC ~300x slower
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2 >/dev/null 2>&1
Rscript $P/scripts/09_SparCC_networks.R > $P/logs/09_SparCC_networks.log 2>&1 && echo "STEP9 OK $(date)" >> $P/logs/09_SparCC_networks.log || echo "STEP9 FAILED $(date)" >> $P/logs/09_SparCC_networks.log
bash $P/scripts/10_finalize_README.sh >> $P/logs/09_SparCC_networks.log 2>&1
