#!/usr/bin/env python3
"""Step 1b: Build QIIME2 sample metadata + a run->sample map from SraRunInfo.csv.
Fixes the SRA typo 'Unknown cultiovar' -> 'Unknown cultivar'; sample IDs have no spaces/apostrophes."""
import csv
B = "/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code/raw_data"
short = {"Aztec Fuji":"Fuji","Cripps Pink":"CP","Gala":"Gala","Granniwinkle":"Nic",
         "Hewe's Crab":"HW","Unknown cultivar":"UN"}   # abbreviations used in the paper's Fig. 5C
rows = [r for r in csv.DictReader(open(f"{B}/SraRunInfo.csv")) if r["Run"].startswith("SRR")]
out = []
for r in rows:
    name = r["SampleName"].replace("_HiFi.fastq", "").replace("cultiovar", "cultivar")
    cultivar, rep, status = name.rsplit("_", 2)
    sid = f"{short[cultivar]}{rep}_{status}"
    out.append(dict(sample_id=sid, run=r["Run"], biosample=r["BioSample"], cultivar=cultivar,
                    replicate=rep, health_status=status, platform=r["Model"], raw_reads=r["spots"]))
out.sort(key=lambda d: (d["health_status"], d["cultivar"], int(d["replicate"])))
cols = list(out[0])
with open(f"{B}/sample_metadata.tsv", "w") as f:
    f.write("\t".join(["sample-id"] + cols[1:]) + "\n")  # QIIME2 requires "sample-id" as ID header
    f.write("#q2:types\t" + "\t".join("numeric" if c == "raw_reads" else "categorical" for c in cols[1:]) + "\n")
    for d in out: f.write("\t".join(d[c] for c in cols) + "\n")
with open(f"{B}/SRR_Acc_List.txt", "w") as f:
    f.write("\n".join(d["run"] for d in out) + "\n")
from collections import Counter
print(Counter((d["health_status"], d["cultivar"]) for d in out))
