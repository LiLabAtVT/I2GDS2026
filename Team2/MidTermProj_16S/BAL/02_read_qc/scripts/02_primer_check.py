#!/usr/bin/env python3
"""Step 2b: count, per sample, whether reads start with the 27F or the 1492R primer.

Why: PacBio sequences each DNA molecule from whichever end it happens to
start at, so a read can run 27F -> 1492R ("forward") or 1492R -> 27F
("reverse"). Knowing the mix tells us DADA2 must flip reads into one
direction (its denoise-ccs mode does this). Reads with no primer at the
start are likely to be dropped in step 3.

PacBio full-length 16S primers (degenerate):
  27F   AGRGTTYGATYMTGGCTCAG
  1492R RGYTACCTTGTTACGACTT
A read is "forward" if 27F occurs in its first 60 bp, "reverse" if 1492R does,
"none" otherwise (primer absent or too many errors).

This is a check only; it does not change or remove any reads.

Usage: 02_primer_check.py <samples.tsv> <fastq_dir> <out.tsv>
(samples.tsv is the clean sheet written by 02a_prepare_samples.py)
"""
import csv
import gzip
import re
import sys

# "Degenerate" primers contain letters that stand for more than one base
# (IUPAC codes), because the primer is a mix designed to match many bacteria.
# e.g. R = A or G, Y = C or T. This table turns each letter into a pattern
# that matches any of its bases: R -> [AG].
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "M": "[AC]",
         "K": "[GT]", "S": "[CG]", "W": "[AT]", "N": "[ACGT]"}
# Build a search pattern for each primer, one letter at a time
F27 = re.compile("".join(IUPAC[b] for b in "AGRGTTYGATYMTGGCTCAG"))
R1492 = re.compile("".join(IUPAC[b] for b in "RGYTACCTTGTTACGACTT"))
# Only look at the start of each read, where a primer should be
WINDOW = 60

# Read the three file paths given on the command line
meta, fq_dir, out = sys.argv[1:4]
# Load the sample sheet (one row per sample)
samples = list(csv.DictReader(open(meta), delimiter="\t"))

with open(out, "w") as fh:
    # Header of the output table
    fh.write("sample-id\trun\tforward\treverse\tnone\ttotal\n")
    for s in samples:
        # Running counts for this sample
        n = {"forward": 0, "reverse": 0, "none": 0}
        # Open the compressed FASTQ directly ("rt" = read as text)
        with gzip.open(f"{fq_dir}/{s['run']}.fastq.gz", "rt") as fq:
            for i, line in enumerate(fq):
                # FASTQ has 4 lines per read (ID, sequence, +, quality);
                # the sequence is line 2 of each group, so skip the rest
                if i % 4 != 1:
                    continue
                # First 60 bases of the read
                head = line[:WINDOW]
                if F27.search(head):
                    n["forward"] += 1
                elif R1492.search(head):
                    n["reverse"] += 1
                else:
                    n["none"] += 1
        tot = sum(n.values())
        # One row per sample in the output table
        fh.write(f"{s['sample-id']}\t{s['run']}\t{n['forward']}\t{n['reverse']}\t{n['none']}\t{tot}\n")
        # Progress message for the log
        print(s["sample-id"], n, flush=True)
