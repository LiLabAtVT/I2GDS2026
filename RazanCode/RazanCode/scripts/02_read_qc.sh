#!/bin/bash
#SBATCH --job-name=fb_readqc
# Slurm account is not hard-coded: run  export SBATCH_ACCOUNT=<your_account>  once before sbatch
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/02_read_qc_%j.log

# ---------------------------------------------------------------
# Step 2: raw read QC (paper: "assessed the quality ... using FastQC
#         and determined a cut-off score for low quality")
#   2a FastQC on every run  + MultiQC summary (PDF plots exported)
#   2b seqkit stats (per sample) + per-read length / mean Q table
#   2c primer / orientation check: 27F and 1492R in the first 40 bp
#   2d R plots -> results/pdf/02_*.pdf
# Submit from the project root directory:  sbatch scripts/02_read_qc.sh
# ---------------------------------------------------------------
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
OUT=results/02_read_qc
mkdir -p $OUT/fastqc $OUT/multiqc results/pdf

module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh

# ---- 2a FastQC + MultiQC ----
module load FastQC/0.12.1-Java-11
fastqc -t 16 --memory 2048 -o $OUT/fastqc raw_data/fastq/*.fastq.gz
module unload FastQC/0.12.1-Java-11
module load MultiQC/1.22.3-foss-2023b
multiqc -f --export -o $OUT/multiqc $OUT/fastqc
cp $OUT/multiqc/multiqc_plots/pdf/*.pdf results/pdf/ 2>/dev/null || true
for f in results/pdf/fastqc_*.pdf; do mv "$f" "results/pdf/02_multiqc_$(basename "$f")"; done 2>/dev/null || true
module unload MultiQC/1.22.3-foss-2023b

# ---- 2b seqkit stats ----
conda activate sra
seqkit stats -a -T -j 16 raw_data/fastq/*.fastq.gz > $OUT/seqkit_stats.tsv
# per-read length + mean quality (id, length, avgQ)
for f in raw_data/fastq/*.fastq.gz; do
  s=$(basename $f .fastq.gz)
  seqkit fx2tab -n -i -l -q $f | awk -v s=$s 'BEGIN{OFS="\t"}{print s,$2,$3}'
done | gzip > $OUT/per_read_len_qual.tsv.gz

# ---- 2c primer / orientation check (IUPAC-aware, first 40 bp) ----
F27=AGRGTTYGATYMTGGCTCAG      # 27F
R1492=RGYTACCTTGTTACGACTT     # 1492R (as written 5'->3')
echo -e "sample\ttotal\tfwd_27F_at_start\trev_1492R_at_start" > $OUT/primer_check.tsv
for f in raw_data/fastq/*.fastq.gz; do
  s=$(basename $f .fastq.gz)
  tot=$(seqkit stats -T $f | awk 'NR==2{print $4}')
  nf=$(seqkit subseq -r 1:40 $f | seqkit grep -s -d -j 16 -p $F27 | seqkit stats -T | awk 'NR==2{print $4}')
  nr=$(seqkit subseq -r 1:40 $f | seqkit grep -s -d -j 16 -p $R1492 | seqkit stats -T | awk 'NR==2{print $4}')
  echo -e "$s\t$tot\t$nf\t$nr" >> $OUT/primer_check.tsv
done
conda deactivate

# ---- 2d plots ----
conda activate qiime2-2024.10
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
export R_LIBS_USER=$CONDA_PREFIX/lib/R/library
Rscript scripts/02_read_qc_plots.R
echo "Step 2 done"
