#!/usr/bin/env python3
"""Step 1a: build QIIME2 sample metadata + manifest from the SRA run table.

Input : metadata/runinfo.csv  (NCBI SRA RunInfo for PRJNA1348523)
Output: metadata/metadata.tsv  (QIIME2 metadata)
        metadata/manifest.tsv  (QIIME2 SingleEndFastqManifestPhred33V2)
Run from the project root directory:  python3 scripts/00_make_metadata.py
"""
import csv
import os

base = os.getcwd()
rows = list(csv.DictReader(open("metadata/runinfo.csv")))

# NOTE: "sample_name" is a reserved QIIME2 metadata column, so the short name is called "label"
# SRA SampleName looks like "Aztec Fuji_1_Canker_HiFi.fastq"
fix = {"Unknown cultiovar": "Unknown cultivar"}  # typo in the SRA record
short = {"Aztec Fuji": "Fuji", "Cripps Pink": "CP", "Gala": "Gala",
         "Granniwinkle": "Nic", "Hewe's Crab": "HW", "Unknown cultivar": "UN"}

out = []
for r in rows:
    cultivar, rep, status = r["SampleName"].split("_")[:3]
    cultivar = fix.get(cultivar, cultivar)
    out.append({
        "sample-id": r["Run"],
        "label": f"{short[cultivar]}{rep}_{status}",
        "cultivar": cultivar,
        "health_status": status,
        "replicate": rep,
        "group": f"{cultivar}_{status}",
        "instrument": r["Model"],
        "biosample": r["BioSample"],
    })
out.sort(key=lambda d: (d["health_status"], d["cultivar"], int(d["replicate"])))

cols = list(out[0].keys())
with open("metadata/metadata.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    f.write("#q2:types\t" + "\t".join("categorical" for _ in cols[1:]) + "\n")
    for d in out:
        f.write("\t".join(d[c] for c in cols) + "\n")

with open("metadata/manifest.tsv", "w") as f:
    f.write("sample-id\tabsolute-filepath\n")
    for d in out:
        f.write(f"{d['sample-id']}\t{base}/raw_data/fastq/{d['sample-id']}.fastq.gz\n")

n_c = sum(d["health_status"] == "Canker" for d in out)
print(f"{len(out)} samples: {n_c} Canker, {len(out) - n_c} Healthy")
