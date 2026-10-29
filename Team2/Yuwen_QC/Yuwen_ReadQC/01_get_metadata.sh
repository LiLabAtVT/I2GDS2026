#!/bin/bash
# Step 2: retrieve SRA run/sample metadata for BioProject PRJNA1348523
# (BakhshiGanje et al. 2026, Phytobiomes J. 10:277-292)
set -euo pipefail
BASE=/projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code
cd "$BASE/metadata"

EUTILS=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
BP_UID=1348523   # NCBI uid of PRJNA1348523

# SRA uids linked to the BioProject (first <Id> is the BioProject itself)
ids=$(curl -s "$EUTILS/elink.fcgi?dbfrom=bioproject&db=sra&id=$BP_UID" \
      | grep -oP '(?<=<Id>)\d+' | tail -n +2 | paste -sd,)
echo "SRA entries linked: $(echo "$ids" | tr ',' '\n' | wc -l)"

# Run table
curl -s "$EUTILS/efetch.fcgi?db=sra&id=$ids&rettype=runinfo&retmode=text" -o sra_runinfo.csv
# Full experiment XML (contains sample attributes: cultivar, health status, ...)
curl -s "$EUTILS/efetch.fcgi?db=sra&id=$ids&rettype=full&retmode=xml" -o sra_experiments.xml

echo "Runs in runinfo: $(($(grep -c . sra_runinfo.csv) - 1))"
