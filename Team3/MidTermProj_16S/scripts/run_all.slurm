#!/bin/bash
# Slurm batch version of run_all.sh -- reproduces the whole reanalysis as one job.
# Submit from anywhere:   sbatch ~/fireblight_canker_reanalysis/scripts/run_all.slurm
# Check progress:         squeue -u $USER ; tail -f ~/fireblight_canker_reanalysis/logs/slurm_run_all_<jobid>.out
#SBATCH --job-name=fireblight_all
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8            # DADA2 / MAFFT / FastTree use --p-n-threads 8
#SBATCH --mem=62G
#SBATCH --time=24:00:00              # interactive run took ~8 h (step 9 SparCC ~1.5 h)
#SBATCH --output=/home/%u/fireblight_canker_reanalysis/logs/slurm_run_all_%j.out
#SBATCH --error=/home/%u/fireblight_canker_reanalysis/logs/slurm_run_all_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=elowenc@vt.edu

set -eo pipefail                     # no -u: Lmod's module function references unset variables
type module >/dev/null 2>&1 || source /etc/profile.d/modules.sh

P=$HOME/fireblight_canker_reanalysis; S=$P/scripts; L=$P/logs
mkdir -p $L
step() { echo "[$(date '+%F %T')] $1"; }

step "job $SLURM_JOB_ID on $(hostname), $SLURM_CPUS_PER_TASK CPUs"
step "01 metadata";  bash $S/01_get_metadata.sh          > $L/01_metadata.log 2>&1
step "02 download";  bash $S/02_download_raw.sh          > $L/02_download.log 2>&1
step "03 fastqc";    bash $S/03_fastqc.sh                > $L/03_fastqc.log 2>&1
step "04 qiime2";    bash $S/04_qiime2_dada2_taxonomy.sh > $L/04_qiime2.log 2>&1

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1   # multi-threaded BLAS on a 128-core node with 8 allocated CPUs makes SparCC ~300x slower
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
step "00 install R packages"; Rscript $S/00_install_R_packages.R > $L/00_install_R.log 2>&1
for s in 05_build_phyloseq_qc 06_composition_alpha_beta 07_taxa_RF_LEfSe_DMM 08_FAPROTAX 09_SparCC_networks; do
  step "$s"; Rscript $S/$s.R > $L/$s.log 2>&1
done
step "10 finalize README"; bash $S/10_finalize_README.sh
step "ALL DONE"
