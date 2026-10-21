#!/bin/bash
# Step 3: read quality check with FastQC (paper: Andrews 2010) + per-sample read-length summary
set -euo pipefail
P=$HOME/fireblight_canker_reanalysis
module load FastQC/0.12.1-Java-11
mkdir -p $P/results/fastqc
fastqc -q -t 8 -o $P/results/fastqc $P/data/raw/*.fastq.gz
