#!/usr/bin/env python3
"""Step 2a: turn the step 1 sample sheet into the standard sheet used by step 2.

Why: different step 1 scripts name their columns differently (for example
"status" instead of "health_status", or the SRR run number stored in
"sample-id"). This script finds the columns it needs, checks that every
sample has a FASTQ file, and writes one clean table with fixed column names:

    sample-id   run   cultivar   health_status

The rest of step 2 only reads that clean table, so it does not care which
step 1 script was used.

How each column is found (first match wins; set the environment variable to
choose a column yourself):
    run            RUN_COL,      else "run",      else the column holding SRR/ERR/DRR numbers
    sample-id      SAMPLE_COL,   else "sample-id" (or "label" if "sample-id" is the run number)
    cultivar       CULTIVAR_COL, else "cultivar"
    health_status  STATUS_COL,   else "health_status", else "status"
Rows starting with "#" (e.g. QIIME 2's "#q2:types" line) are skipped.

Usage: 02a_prepare_samples.py <step1_sheet.tsv> <fastq_dir> <out.tsv>
"""
import csv
import os
import re
import sys

sheet, fq_dir, out = sys.argv[1:4]


def fail(msg):
    # Print a clear message to the log and stop the job
    sys.exit(f"ERROR (02a_prepare_samples): {msg}")


# Read the sheet, skipping comment/directive rows that start with "#"
with open(sheet, newline="") as fh:
    rows = [r for r in csv.reader(fh, delimiter="\t") if r and not r[0].startswith("#")]
if len(rows) < 2:
    fail(f"{sheet} has no samples")
header, data = rows[0], rows[1:]
cols = {name: [r[i] if i < len(r) else "" for r in data] for i, name in enumerate(header)}


def pick(env_var, candidates, what):
    # Use the column the user named, otherwise the first candidate that exists
    if os.environ.get(env_var):
        if os.environ[env_var] not in cols:
            fail(f"{env_var}={os.environ[env_var]} is not a column. Columns are: {', '.join(header)}")
        return os.environ[env_var]
    for c in candidates:
        if c in cols:
            return c
    fail(f"no {what} column found (looked for {', '.join(candidates)}). "
         f"Columns are: {', '.join(header)}. Set {env_var}=<column name>.")


# Run numbers look like SRR36091169 (SRA), ERR... (ENA) or DRR... (DDBJ)
is_run = re.compile(r"^[SED]RR\d+$")
run_like = [c for c in header if all(is_run.match(v) for v in cols[c])]
run_col = pick("RUN_COL", ["run"] + run_like, "run (SRR number)")
# If sample-id just repeats the run number, a "label" column has the readable name
sample_default = ["label", "sample-id"] if cols.get("sample-id") == cols[run_col] else ["sample-id", "label"]
sample_col = pick("SAMPLE_COL", sample_default + [run_col], "sample name")
cultivar_col = pick("CULTIVAR_COL", ["cultivar"], "cultivar")
status_col = pick("STATUS_COL", ["health_status", "status"], "health status")

# Every sample needs a FASTQ file named <run>.fastq.gz
missing = [r for r in cols[run_col] if not os.path.isfile(os.path.join(fq_dir, f"{r}.fastq.gz"))]
if missing:
    fail(f"{len(missing)} of {len(data)} samples have no FASTQ in {fq_dir} "
         f"(expected <run>.fastq.gz), e.g. {', '.join(missing[:3])}")

with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["sample-id", "run", "cultivar", "health_status"])
    for i in range(len(data)):
        w.writerow([cols[sample_col][i], cols[run_col][i], cols[cultivar_col][i], cols[status_col][i]])

print(f"samples: {len(data)}  columns used: sample-id={sample_col}, run={run_col}, "
      f"cultivar={cultivar_col}, health_status={status_col}")
print(f"health status groups: {sorted(set(cols[status_col]))}")
