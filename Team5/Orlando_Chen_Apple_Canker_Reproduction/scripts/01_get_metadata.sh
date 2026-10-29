#!/bin/bash
# Step 1: Retrieve SRA run metadata for BioProject PRJNA1348523 (BakhshiGanje et al. 2026, Phytobiomes J.)
# Uses NCBI E-utilities via curl (read-only). Output: raw_data/SraRunInfo.csv
set -euo pipefail
B=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code
OUT=$B/raw_data
EU=https://eutils.ncbi.nlm.nih.gov/entrez/eutils

# Run table (one row per SRR)
curl -s "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA1348523" -o $OUT/SraRunInfo.csv
# Fallback via esearch/efetch if the above returned nothing
if [ ! -s $OUT/SraRunInfo.csv ]; then
  WEB=$(curl -s "$EU/esearch.fcgi?db=sra&term=PRJNA1348523&usehistory=y" | grep -oP '(?<=<WebEnv>).*(?=</WebEnv>)')
  KEY=$(curl -s "$EU/esearch.fcgi?db=sra&term=PRJNA1348523&usehistory=y" | grep -oP '(?<=<QueryKey>).*(?=</QueryKey>)')
  curl -s "$EU/efetch.fcgi?db=sra&query_key=$KEY&WebEnv=$WEB&rettype=runinfo&retmode=text" -o $OUT/SraRunInfo.csv
fi

echo "Runs: $(($(wc -l < $OUT/SraRunInfo.csv)-1))"
