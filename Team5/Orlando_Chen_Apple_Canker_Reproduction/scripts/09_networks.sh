#!/bin/bash
#SBATCH --job-name=sparcc
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code/logs/09_networks_%j.log
# Step 9: SparCC networks with FastSpar (SparCC re-implementation; 1000 bootstraps for pseudo P values, as MENAP's SparCC)
set -euo pipefail
B=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code
O=$B/results/09_networks
N=$SLURM_CPUS_PER_TASK; NBOOT=1000
source $HOME/miniconda3/etc/profile.d/conda.sh
# FastSpar in its own env (installed here, on the compute node)
[ -d $HOME/miniconda3/envs/fastspar ] || mamba create -y -n fastspar -c conda-forge -c bioconda fastspar
conda activate fastspar; fastspar --version || true; conda deactivate

conda activate canker_r
Rscript $B/scripts/09_networks.R prep
conda deactivate

conda activate fastspar
for net in canker_T5 canker_T0 healthy_T1; do
  d=$O/$net; rm -rf $d/boot_counts $d/boot_cor; mkdir -p $d/boot_counts $d/boot_cor
  echo "== $net"; date
  fastspar --otu_table $d/counts.tsv --correlation $d/median_correlation.tsv --covariance $d/median_covariance.tsv \
    --iterations 50 --exclude_iterations 10 --threshold 0.1 --threads $N --seed 2026 > $d/fastspar.log
  fastspar_bootstrap --otu_table $d/counts.tsv --number $NBOOT --prefix $d/boot_counts/data --threads $N --seed 2026
  ls $d/boot_counts/*.tsv | xargs -P $N -I{} sh -c \
    'b=$(basename {}); fastspar --otu_table {} --correlation '"$d"'/boot_cor/cor_$b --covariance '"$d"'/boot_cor/cov_$b --iterations 5 --threads 1 --seed 2026 >/dev/null'
  fastspar_pvalues --otu_table $d/counts.tsv --correlation $d/median_correlation.tsv --prefix $d/boot_cor/cor_data_ \
    --permutations $NBOOT --outfile $d/pvalues.tsv --threads $N
  rm -rf $d/boot_counts $d/boot_cor/cov_*   # keep bootstrap correlations only
  tar -czf $d/boot_cor.tar.gz -C $d boot_cor && rm -rf $d/boot_cor
done
conda deactivate

conda activate canker_r
Rscript $B/scripts/09_networks.R analyze
echo "Step 9 done"
