#!/bin/bash
# Step 0: list the ARC software modules relevant to this re-analysis
module -t spider 2>&1 | grep -iE '^(SRA-Toolkit|sra-tools|R|R-bundle|Miniforge|Miniconda|Anaconda|QIIME2|FastQC|pigz|LEfSe)' | sort -u
