#!/bin/bash
#SBATCH --job-name=dl_reads
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code/logs/02_download_%j.log
# Step 2: Download PacBio HiFi 16S FASTQs for PRJNA1348523 from the ENA mirror, rename to sample IDs, verify MD5.
set -euo pipefail
B=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code
R=$B/raw_data; mkdir -p $R/fastq; cd $R

curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1348523&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,fastq_md5&format=tsv" > ena_filereport.tsv

# run -> sample_id map from Step 1 metadata (skip header + #q2:types line)
awk -F'\t' 'NR>2{print $2"\t"$1}' sample_metadata.tsv > run2sample.tsv

: > md5_check.txt
tail -n +2 ena_filereport.tsv | while IFS=$'\t' read -r run url bytes md5; do
  sid=$(awk -F'\t' -v r="$run" '$1==r{print $2}' run2sample.tsv)
  f=fastq/${sid}.fastq.gz
  for try in 1 2 3; do
    wget -q -O "$f" "https://$url" && [ "$(md5sum < "$f" | cut -d' ' -f1)" == "$md5" ] && break
    echo "retry $try: $run"; sleep 10
  done
  got=$(md5sum < "$f" | cut -d' ' -f1)
  [ "$got" == "$md5" ] && st=OK || st=FAIL
  echo -e "$run\t$sid\t$md5\t$got\t$st" | tee -a md5_check.txt
done
echo "OK: $(grep -c OK md5_check.txt) / $(tail -n +2 ena_filereport.tsv | wc -l)"
