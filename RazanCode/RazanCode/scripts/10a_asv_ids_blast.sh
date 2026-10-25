#!/bin/bash
# Step 10a: readable ASV numbers + BLAST species calls (paper: E. amylovora ASVs and DMM drivers identified by NCBI BLAST)
#   - ASVs (bacteria-only table) numbered by total abundance: ASV0001 = most abundant. Our IDs cannot be mapped to the
#     paper's numbers (ASV1032 etc.), which are not published with sequences.
#   - BLAST: every Erwinia ASV + the 100 most abundant ASVs vs NCBI 16S_ribosomal_RNA (type strains), top 5 hits each.
# Run from the project root directory:  bash scripts/10a_asv_ids_blast.sh
set -eo pipefail
OUT=results/10_erwinia_dmm
mkdir -p $OUT
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh

conda activate fb_r
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=$CONDA_PREFIX/lib/R/library
Rscript - <<'EOF'
tab <- as.matrix(read.delim("results/04_taxonomy/export/feature-table.tsv", skip = 1, check.names = FALSE, row.names = 1))
tax <- read.delim("results/04_taxonomy/export/taxonomy_ranks.tsv")
fa  <- readLines("results/04_taxonomy/export/dna-sequences.fasta")
seqs <- setNames(fa[c(FALSE, TRUE)], sub("^>", "", fa[c(TRUE, FALSE)]))
tot <- sort(rowSums(tab), decreasing = TRUE)
map <- data.frame(asv_num = sprintf("ASV%04d", seq_along(tot)), feature_id = names(tot), total_reads = tot,
                  genus = sub("^g__", "", tax$Genus[match(names(tot), tax$asv)]),
                  taxonomy = tax$taxon[match(names(tot), tax$asv)], length = nchar(seqs[names(tot)]))
write.table(map, "results/10_erwinia_dmm/asv_id_map.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
q <- unique(c(map$feature_id[grepl("Erwinia", map$genus)], head(map$feature_id, 100)))
writeLines(paste0(">", map$asv_num[match(q, map$feature_id)], "\n", seqs[q]), "results/10_erwinia_dmm/blast_query.fasta")
cat(length(q), "query ASVs written;", sum(grepl("Erwinia", map$genus)), "Erwinia ASVs\n")
EOF
conda deactivate

conda activate sra
blastn -query $OUT/blast_query.fasta -db ref/ncbi_16S/16S_ribosomal_RNA -max_target_seqs 5 -num_threads 4 \
  -outfmt "6 qseqid sseqid pident length qlen slen mismatch evalue bitscore stitle" > $OUT/blast_top5.tsv
echo "BLAST hits: $(wc -l < $OUT/blast_top5.tsv)"
