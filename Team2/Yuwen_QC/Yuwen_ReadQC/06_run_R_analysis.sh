#!/bin/bash
#SBATCH --job-name=R_analysis
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code/logs/%x_%j.out
# Step 6: all downstream statistics and figures in R (PDFs -> results/pdf)
#   usage: sbatch [--export=ALL,STEPS=06c:06d] scripts/06_run_R_analysis.sh   (steps separated by ':')
set -euo pipefail
cd /projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
# one BLAS thread per process: forked workers (SpiecEasi::sparccboot) deadlock with threaded BLAS
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
STEPS=${STEPS:-"06a:06b:06c:06d:06e:06f:06g"}
STEPS=${STEPS//:/ }
for s in $STEPS; do
  f=$(ls scripts/${s}_*.R)
  echo "================ $f  ($(date))"
  Rscript "$f"
done
echo "PDFs:"; ls -1 results/pdf
