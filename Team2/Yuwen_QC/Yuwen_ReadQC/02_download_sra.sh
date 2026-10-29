#!/bin/bash
#SBATCH --job-name=dl_sra
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code/logs/%x_%j.out
# Step 3: download the 48 PacBio HiFi 16S runs of PRJNA1348523 and build sample metadata
set -euo pipefail
BASE=/projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code
cd "$BASE"
module load SRA-Toolkit/3.2.0-gompi-2024a pigz/2.8-GCCcore-13.3.0

# --- 1. De-duplicated run table -> sample metadata (QIIME2 format) ---
# runinfo columns: 1 Run, 4 spots, 26 BioSample, 30 SampleName (e.g. "Aztec Fuji_1_Canker_HiFi.fastq")
awk -F, 'NR>1 && $1 ~ /^SRR/ && !seen[$1]++ {print $1"\t"$4"\t"$26"\t"$30}' metadata/sra_runinfo.csv \
  | sort -k4,4 > metadata/runs_unique.tsv

{
  printf "sample-id\trun\tbiosample\tcultivar\treplicate\thealth\tcultivar_health\tspots\n"
  printf "#q2:types\tcategorical\tcategorical\tcategorical\tcategorical\tcategorical\tcategorical\tnumeric\n"
  while IFS=$'\t' read -r run spots bios name; do
    name=${name%_HiFi.fastq}
    health=${name##*_}; tmp=${name%_*}; rep=${tmp##*_}; cult=${tmp%_*}
    cult=${cult/cultiovar/cultivar}                 # fix typo in SRA submission
    cid=$(echo "$cult" | tr -d " '")                # e.g. AztecFuji, HewesCrab
    printf "%s_%s_%s\t%s\t%s\t%s\t%s\t%s\t%s_%s\t%s\n" \
      "$cid" "$rep" "$health" "$run" "$bios" "$cult" "$rep" "$health" "$cid" "$health" "$spots"
  done < metadata/runs_unique.tsv
} > metadata/sample_metadata.tsv
echo "Samples in metadata: $(($(wc -l < metadata/sample_metadata.tsv) - 2))"

# --- 2. Download reads ---
mkdir -p raw_data/sra_cache raw_data/fastq
tail -n +3 metadata/sample_metadata.tsv | while IFS=$'\t' read -r sid run rest; do
  out=raw_data/fastq/${sid}.fastq.gz
  [[ -s $out ]] && { echo "skip $sid"; continue; }
  prefetch "$run" -O raw_data/sra_cache --max-size 5G
  fasterq-dump raw_data/sra_cache/"$run" -O raw_data/fastq -e "$SLURM_CPUS_PER_TASK" -t "${TMPDIR:-/tmp}"
  mv raw_data/fastq/"$run".fastq raw_data/fastq/"$sid".fastq
  pigz -p "$SLURM_CPUS_PER_TASK" raw_data/fastq/"$sid".fastq
  echo "done $sid ($run)"
done

# --- 3. Verify read counts against SRA ---
tail -n +3 metadata/sample_metadata.tsv | while IFS=$'\t' read -r sid run bios cult rep health ch spots; do
  n=$(( $(pigz -dc raw_data/fastq/${sid}.fastq.gz | wc -l) / 4 ))
  printf "%s\t%s\t%s\t%s\n" "$sid" "$spots" "$n" "$([[ $n == "$spots" ]] && echo OK || echo MISMATCH)"
done > metadata/download_check.tsv
grep -c OK metadata/download_check.tsv | xargs echo "Runs verified OK:"
rm -rf raw_data/sra_cache
