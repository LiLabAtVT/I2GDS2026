#!/bin/bash
# Show per-sample read depth of the organelle-filtered ASV table (to choose a rarefaction depth)
cd /projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code/results/02_qiime2/export
ls table_summary | head -30
f=feature_table.tsv
echo "ASVs (filtered): $(($(wc -l < $f) - 2))   ASVs (raw): $(($(wc -l < feature_table_raw.tsv) - 2))"
awk -F'\t' 'NR==2{for(i=2;i<=NF;i++) n[i]=$i} NR>2{for(i=2;i<=NF;i++) s[i]+=$i}
  END{for(i in n) print s[i]"\t"n[i]}' $f | sort -n
echo "== raw (before organelle removal) total reads per sample, lowest 5:"
awk -F'\t' 'NR==2{for(i=2;i<=NF;i++) n[i]=$i} NR>2{for(i=2;i<=NF;i++) s[i]+=$i}
  END{for(i in n) print s[i]"\t"n[i]}' feature_table_raw.tsv | sort -n | head -5
echo "== top taxa removed as organelles:"
grep -iE 'chloroplast|mitochondria' taxonomy/taxonomy.tsv | wc -l
