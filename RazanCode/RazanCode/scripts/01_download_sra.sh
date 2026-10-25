#!/bin/bash
#SBATCH --job-name=fb_download
# Slurm account is not hard-coded: run  export SBATCH_ACCOUNT=<your_account>  once before sbatch
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=06:00:00
#SBATCH --output=logs/01_download_%j.log

# ---------------------------------------------------------------
# Step 1: download raw PacBio HiFi 16S reads for PRJNA1348523
# Paper: BakhshiGanje et al. 2026, Phytobiomes J. (PBIOMES-10-25-0086-R)
# Input : metadata/SRR_Acc_List.txt (48 runs)
# Output: raw_data/fastq/SRRxxxxxxxx.fastq.gz
# Submit from the project root directory:  sbatch scripts/01_download_sra.sh
# ---------------------------------------------------------------

module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate sra
set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p raw_data/fastq raw_data/sra_tmp

while read -r SRR; do
  if [ -s "raw_data/fastq/${SRR}.fastq.gz" ]; then
    echo "SKIP $SRR (already downloaded)"; continue
  fi
  echo "Downloading $SRR"
  prefetch "$SRR" --max-size 100G -O raw_data/sra_tmp
  fasterq-dump "raw_data/sra_tmp/$SRR" --threads 8 -t raw_data/sra_tmp -O raw_data/fastq/
  pigz -p 8 raw_data/fastq/${SRR}*.fastq
  rm -rf "raw_data/sra_tmp/$SRR"
done < metadata/SRR_Acc_List.txt

rmdir raw_data/sra_tmp 2>/dev/null || true

# checksum + simple count check
n_exp=$(wc -l < metadata/SRR_Acc_List.txt)
n_got=$(ls raw_data/fastq/*.fastq.gz | wc -l)
echo "Expected $n_exp files, found $n_got"
(cd raw_data/fastq && md5sum *.fastq.gz > ../md5sums.txt)
du -sh raw_data/fastq
