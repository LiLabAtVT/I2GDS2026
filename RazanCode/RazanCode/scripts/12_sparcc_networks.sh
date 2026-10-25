#!/bin/bash
#SBATCH --job-name=fb_sparcc
# Slurm account is not hard-coded: run  export SBATCH_ACCOUNT=<your_account>  once before sbatch
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/12_sparcc_%j.log

# ---------------------------------------------------------------
# Step 12: SparCC co-occurrence networks (paper Figs. 7-8, Table 1)
#   paper used the MENAP/iNAP web pipeline (SparCC); here FastSpar 1.0 (c++ SparCC) with 1,000 bootstraps for
#   pseudo P-values, then igraph for topology (12b_network_analysis.R). Edges: |r| >= 0.2 and P < 0.01.
# Submit from the project root directory:  sbatch scripts/12_sparcc_networks.sh
# ---------------------------------------------------------------
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh

conda activate fb_r
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=$CONDA_PREFIX/lib/R/library
Rscript scripts/12a_network_inputs.R
conda deactivate

conda activate fastspar
NBOOT=1000
for net in canker_T5 canker_T0 healthy_T1; do
  D=results/12_networks/$net
  echo "== $net"
  if [ -s $D/pvalues.tsv ]; then echo "already done"; continue; fi
  rm -rf $D/boot $D/boot_cor; mkdir -p $D/boot $D/boot_cor
  fastspar --otu_table $D/input.tsv --correlation $D/correlation.tsv --covariance $D/covariance.tsv \
           --iterations 50 --threads 16 --seed 2026 --yes > $D/fastspar.log 2>&1 \
    || { echo "fastspar failed for $net (see $D/fastspar.log)"; continue; }
  fastspar_bootstrap --otu_table $D/input.tsv --number $NBOOT --prefix $D/boot/data --threads 16 --seed 2026 >> $D/fastspar.log 2>&1
  for f in $D/boot/data_*.tsv; do
    b=$(basename $f .tsv)
    fastspar --otu_table $f --correlation $D/boot_cor/cor_$b.tsv --covariance $D/boot_cor/cov_$b.tsv \
             --iterations 5 --threads 16 --yes >> $D/fastspar_boot.log 2>&1
  done
  fastspar_pvalues --otu_table $D/input.tsv --correlation $D/correlation.tsv --prefix $D/boot_cor/cor_data_ \
                   --permutations $NBOOT --outfile $D/pvalues.tsv --threads 16 >> $D/fastspar.log 2>&1
  rm -rf $D/boot $D/boot_cor $D/fastspar_boot.log     # bootstrap files are large; P-values kept
  echo "done $net"
done
conda deactivate

conda activate fb_r
Rscript scripts/12b_network_analysis.R
echo "Step 12 done"
