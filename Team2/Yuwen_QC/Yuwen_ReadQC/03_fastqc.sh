#!/bin/bash
#SBATCH --job-name=read_qc
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code/logs/%x_%j.out
# Step 3: raw read QC (paper: "assessed the quality ... using FastQC (Andrews 2010)
#         and determined a cut-off score for low quality")
#   3a FastQC on every run + MultiQC summary (PDF plots exported)
#   3b seqkit stats (per sample) + per-read length / mean Q table
#   3c primer / orientation check: 27F and 1492R in the first 40 bp (IUPAC-aware)
#   3d R plots -> results/pdf/03_*.pdf
#   usage: sbatch scripts/03_fastqc.sh
set -euo pipefail
BASE=/projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code
cd "$BASE"
T=${SLURM_CPUS_PER_TASK:-4}
FQ=results/01_fastqc        # FastQC zips (also read by scripts/06a_qc_summary.R)
OUT=results/01_read_qc
mkdir -p "$FQ" "$OUT/multiqc" results/pdf

# ---- 3a FastQC + MultiQC ----
module load FastQC/0.12.1-Java-11
fastqc -t "$T" --memory 2048 -o "$FQ" raw_data/fastq/*.fastq.gz
echo "FastQC reports: $(ls $FQ/*_fastqc.zip | wc -l)"
module reset
module load MultiQC/1.22.3-foss-2023b
multiqc -f --export -o "$OUT/multiqc" "$FQ"
for f in "$OUT"/multiqc/multiqc_plots/pdf/*.pdf; do
  [ -e "$f" ] && cp "$f" "results/pdf/03_multiqc_$(basename "$f")"
done
module reset

# ---- 3b seqkit stats ----
# seqkit is not an ARC module; keep a small conda env inside the project
module load Miniconda3
set +u; source "$(conda info --base)/etc/profile.d/conda.sh"
ENV="$BASE/envs/seqkit"
[ -x "$ENV/bin/seqkit" ] || conda create -y -p "$ENV" --override-channels -c conda-forge -c bioconda seqkit
conda activate "$ENV"; set -u
seqkit stats -a -T -j "$T" raw_data/fastq/*.fastq.gz > "$OUT/seqkit_stats.tsv"
# per-read length + mean quality (sample, length, avgQ)
for f in raw_data/fastq/*.fastq.gz; do
  s=$(basename "$f" .fastq.gz)
  seqkit fx2tab -n -i -l -q "$f" | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$2,$3}'
done | gzip > "$OUT/per_read_len_qual.tsv.gz"

# ---- 3c primer / orientation check (IUPAC-aware, first 40 bp, + strand only) ----
# HiFi CCS reads come in both orientations: 27F at the start = forward read,
# 1492R at the start = reverse read. DADA2 denoise-ccs reorients and trims them.
F27=AGRGTTYGATYMTGGCTCAG      # 27F
R1492=RGYTACCTTGTTACGACTT     # 1492R (5'->3')
echo -e "sample\ttotal\tfwd_27F_at_start\trev_1492R_at_start" > "$OUT/primer_check.tsv"
for f in raw_data/fastq/*.fastq.gz; do
  s=$(basename "$f" .fastq.gz)
  tot=$(seqkit stats -T "$f" | awk 'NR==2{print $4}')
  nf=$(seqkit subseq -r 1:40 "$f" | seqkit grep -s -d -P -j "$T" -p "$F27"   | seqkit stats -T | awk 'NR==2{print $4}')
  nr=$(seqkit subseq -r 1:40 "$f" | seqkit grep -s -d -P -j "$T" -p "$R1492" | seqkit stats -T | awk 'NR==2{print $4}')
  echo -e "$s\t$tot\t$nf\t$nr" >> "$OUT/primer_check.tsv"
done
set +u; conda deactivate; set -u
module reset

# ---- 3d plots ----
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
Rscript scripts/03b_read_qc_plots.R
echo "Step 3 done"
