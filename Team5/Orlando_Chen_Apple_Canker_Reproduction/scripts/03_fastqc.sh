#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code/logs/03_fastqc_%j.log
# Step 3: FastQC on raw HiFi reads + per-sample read summary and QC PDF.
set -euo pipefail
B=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code
OUT=$B/results/03_fastqc; mkdir -p $OUT $B/pdf
module load FastQC/0.12.1-Java-11 apptainer/1.4.0
fastqc --version
fastqc -t $SLURM_CPUS_PER_TASK -o $OUT $B/raw_data/fastq/*.fastq.gz
SIF=/apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif
apptainer exec --bind /projects $SIF python $B/scripts/03_read_summary.py
echo "Step 3 done"
