#!/bin/bash
# Step 11a: BLAST every bacterial ASV vs NCBI 16S type strains -> species label where best hit >= 99% identity
#   (SILVA 138 NB gives no species names; FAPROTAX assigns e.g. plant_pathogen at species level: Erwinia amylovora, ...)
# Run from the project root directory:  bash scripts/11a_blast_all_asvs.sh
set -eo pipefail
OUT=results/11_faprotax; mkdir -p $OUT
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate sra
blastn -query results/04_taxonomy/export/dna-sequences.fasta -db ref/ncbi_16S/16S_ribosomal_RNA \
  -max_target_seqs 5 -num_threads 8 -outfmt "6 qseqid sseqid pident length qlen slen mismatch evalue bitscore stitle" \
  > $OUT/blast_all_asvs_top5.tsv
# best hit per ASV (highest bitscore, then identity); species = first two words of the title
sort -t$'\t' -k1,1 -k9,9gr -k3,3gr $OUT/blast_all_asvs_top5.tsv | awk -F'\t' 'BEGIN{OFS="\t"; print "feature_id","pident","aln_len","qlen","species"}
  !seen[$1]++ {split($10,w," "); print $1,$3,$4,$5,w[1]" "w[2]}' > $OUT/blast_all_asvs_best.tsv
echo "ASVs with a hit: $(($(wc -l < $OUT/blast_all_asvs_best.tsv)-1));  >=99% identity: $(awk -F'\t' 'NR>1 && $2>=99' $OUT/blast_all_asvs_best.tsv | wc -l)"
