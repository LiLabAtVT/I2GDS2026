#!/bin/bash
# Step 1: download PRJNA1348523 PacBio HiFi reads (via ENA mirror) and build QIIME2 metadata/manifest.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p raw_data metadata
# SRA toolkit can't resolve accessions on this cluster; ENA mirrors the same runs as fastq.gz with md5
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1348523&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv" -o metadata/ena_filereport.tsv
tail -n +2 metadata/ena_filereport.tsv | while IFS=$'\t' read -r run url md5 bytes; do
  out=raw_data/$run.fastq.gz
  [ -s "$out" ] && echo "$md5  $out" | md5sum -c --quiet && continue
  curl -sS --retry 3 -o "$out" "https://$url"
  echo "$md5  $out" | md5sum -c --quiet || { echo "md5 FAIL $run"; exit 1; }
done
# metadata: sample-id = SRR, cultivar / replicate / status parsed from SampleName (col 30)
python3 - <<'PY'
import csv
rows = list(csv.DictReader(open("metadata/runinfo.csv")))
with open("metadata/sample-metadata.tsv", "w") as m, open("metadata/manifest.tsv", "w") as f:
    m.write("sample-id\tsample_name\tcultivar\treplicate\tstatus\n")
    f.write("sample-id\tabsolute-filepath\n")
    import os
    for r in rows:
        cv, rep, st, _ = r["SampleName"].split("_")
        cv = cv.replace("cultiovar", "cultivar").replace("'", "")
        m.write(f"{r['Run']}\t{cv}_{rep}_{st}\t{cv}\t{rep}\t{st}\n")
        f.write(f"{r['Run']}\t{os.path.abspath('raw_data/'+r['Run']+'.fastq.gz')}\n")
PY
echo "done: $(ls raw_data/*.fastq.gz | wc -l) files"
