# Part 2. Pseudomonas is the strongest shared recruit (Part 1). Which Pseudomonas variants win, and does it depend on site or genotype?
#   2a  name the variants: BLAST against NCBI 16S type strains, keep ties
#   2b  how much Pseudomonas is there in each compartment
#   2c  the winning variant (most abundant ASV) in every rhizosphere plot
#   2d  were the winners already in the bulk soil? which variant dominates bulk soil?
#   2e  does genotype decide the Pseudomonas variant composition? (PERMANOVA, permutations within block)
#   2f  figure: Pseudomonas variants in every rhizosphere plot
# Usage: Rscript scripts/02_pseudomonas_variants.R <input_dataset dir> <output dir> <threads>
a <- commandArgs(TRUE); IN <- a[1]; OUT <- a[2]; NT <- as.integer(ifelse(length(a) > 2, a[3], 4))
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(vegan)})
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) read.delim(file.path(IN, f), check.names = FALSE)
md <- rd("sample_metadata.tsv")
pc <- rd("pseudomonas_asv_counts_rarefied80k.tsv"); PM <- t(as.matrix(pc[, -1])); colnames(PM) <- pc$asv; PM <- PM[md$sample_id, ]  # samples x ASVs
info <- rd("pseudomonas_asv_info.tsv")
DEPTH <- 80000                                     # every sample was rarefied to 80,000 reads

## ---- 2a BLAST naming (parameters as in the study) ----
QCOV <- 0.95; TIE <- 0.25                          # >= 95% query coverage; hits within 0.25% identity of the best are ties
bdir <- file.path(OUT, "blastdb"); dir.create(bdir, showWarnings = FALSE)
ref <- file.path(IN, "reference/ncbi16S_typestrains_subset.fasta")
system2("makeblastdb", c("-in", ref, "-dbtype", "nucl", "-out", file.path(bdir, "typestrains")), stdout = FALSE)
q <- file.path(OUT, "query.fasta")
fa <- readLines(file.path(IN, "pseudomonas_asv_sequences.fasta")); hdr <- sub("^>", "", fa[c(TRUE, FALSE)])
keep <- hdr %in% info$asv[info$blasted_in_study]  # the variants named in the study (the 69 most abundant / enriched)
writeLines(c(rbind(fa[c(TRUE, FALSE)][keep], fa[c(FALSE, TRUE)][keep])), q)
bl_out <- file.path(OUT, "2a_blast_hits.tsv")
system2("blastn", c("-query", q, "-db", file.path(bdir, "typestrains"), "-outfmt", shQuote("6 qseqid sseqid pident length qlen mismatch gapopen evalue bitscore stitle"),
                    "-max_target_seqs", "100", "-num_threads", NT, "-out", bl_out))
bl <- read.delim(bl_out, header = FALSE, col.names = c("qseqid", "sseqid", "pident", "length", "qlen", "mismatch", "gapopen", "evalue", "bitscore", "stitle"))
bl$stitle <- sub("^[A-Z]{2}_[0-9]+\\.[0-9]+\\s+", "", bl$stitle)   # drop the leading accession from the FASTA title
ty <- bl %>% mutate(cov = length / qlen, genus = sub("^\\s*(\\S+)\\s+(\\S+).*$", "\\1", stitle), epithet = sub("^\\s*(\\S+)\\s+(\\S+).*$", "\\2", stitle)) %>%
  filter(cov >= QCOV) %>% group_by(qseqid) %>% filter(pident >= max(pident) - TIE) %>%
  summarise(identity = max(pident), n_tied = n_distinct(paste(genus, epithet)),
            closest = paste(head(unique(paste0(substr(genus, 1, 1), ". ", epithet)), 4), collapse = ", "),
            non_pseudomonas_tie = any(genus != "Pseudomonas"), .groups = "drop") %>%
  mutate(name = paste0(sub(",.*", "", closest), " (", qseqid, ")")) %>% rename(asv = qseqid) %>%
  arrange(as.integer(sub("ASV", "", asv)))
write.table(ty, file.path(OUT, "2a_variant_names.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
nm <- setNames(ty$name, ty$asv)

## ---- 2b Pseudomonas share per compartment ----
md$pse_pct <- 100 * rowSums(PM) / DEPTH
share <- md %>% group_by(site, compartment) %>% summarise(n = n(), mean_pct = round(mean(pse_pct), 2), median_pct = round(median(pse_pct), 2),
                                                          min_pct = round(min(pse_pct), 2), max_pct = round(max(pse_pct), 2), .groups = "drop")
write.table(share, file.path(OUT, "2b_pseudomonas_share.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 2c winner in every rhizosphere plot ----
win <- md %>% filter(compartment == "Rhizosphere") %>% rowwise() %>%
  mutate(winner = colnames(PM)[which.max(PM[sample_id, ])], winner_share_of_pseudomonas = round(100 * max(PM[sample_id, ]) / sum(PM[sample_id, ]), 1)) %>%
  ungroup() %>% mutate(winner_name = ifelse(winner %in% names(nm), nm[winner], winner))
bulk <- md %>% filter(compartment == "Bulk") %>% select(site_plot, bulk_id = sample_id)
win <- win %>% left_join(bulk, by = "site_plot") %>% rowwise() %>%
  mutate(winner_reads_in_own_bulk = ifelse(is.na(bulk_id), NA, PM[bulk_id, winner])) %>% ungroup()
write.table(win %>% select(site, plot, block, genotype_label, genotype, pse_pct, winner, winner_name, winner_share_of_pseudomonas, winner_reads_in_own_bulk),
            file.path(OUT, "2c_winner_per_plot.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ws <- win %>% group_by(site) %>% summarise(plots = n(), distinct_winners = n_distinct(winner), top_winner = names(sort(table(winner), TRUE))[1],
  top_winner_plots = max(table(winner)), median_winner_share = median(winner_share_of_pseudomonas),
  winner_found_in_own_bulk = sum(winner_reads_in_own_bulk > 0, na.rm = TRUE), plots_with_bulk = sum(!is.na(winner_reads_in_own_bulk)), .groups = "drop")

## ---- 2d dominant Pseudomonas variant in bulk soil ----
bd <- md %>% filter(compartment == "Bulk") %>% rowwise() %>%
  mutate(pse_reads = sum(PM[sample_id, ]), dominant = colnames(PM)[which.max(PM[sample_id, ])],
         dominant_share = round(100 * max(PM[sample_id, ]) / pse_reads, 1), ASV3_reads = PM[sample_id, "ASV3"]) %>% ungroup()
write.table(bd %>% select(site, plot, block, pse_reads, dominant, dominant_share, ASV3_reads), file.path(OUT, "2d_bulk_dominant_per_plot.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
bs <- bd %>% group_by(site) %>% summarise(bulk_plots = n(), median_pse_reads = median(pse_reads),
  dominant_counts = paste(names(sort(table(dominant), TRUE)), sort(table(dominant), TRUE), sep = ":", collapse = ", "),
  ASV3_detected = sum(ASV3_reads > 0), .groups = "drop")

## ---- 2e genotype effect on variant composition ----
MINREADS <- 200; NPERM <- 999
set.seed(1)
pg <- bind_rows(lapply(c("Kentland", "Warsaw"), function(s) {
  d <- md %>% filter(site == s, compartment == "Rhizosphere"); x <- PM[d$sample_id, ]; ok <- rowSums(x) >= MINREADS; d <- d[ok, ]; x <- x[ok, ]
  x <- x[, colSums(x) > 0]; dd <- vegdist(x / rowSums(x), "bray")
  ad <- adonis2(dd ~ genotype, data = d, permutations = how(blocks = factor(d$block), nperm = NPERM))
  data.frame(site = s, samples = nrow(d), R2 = round(ad$R2[1], 3), p = ad$`Pr(>F)`[1]) }))
write.table(pg, file.path(OUT, "2e_genotype_permanova.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 2f figure ----
GL <- md %>% filter(compartment == "Rhizosphere")
top8 <- names(sort(colSums(PM[GL$sample_id, ]), decreasing = TRUE))[1:8]
pp <- GL %>% select(sample_id, site, block, genotype_label) %>% bind_cols(as.data.frame(PM[GL$sample_id, top8]) %>% setNames(nm[top8])) %>%
  mutate(Others = rowSums(PM[GL$sample_id, setdiff(colnames(PM), top8), drop = FALSE])) %>%
  pivot_longer(c(unname(nm[top8]), "Others"), names_to = "variant", values_to = "reads") %>%
  mutate(pct = 100 * reads / DEPTH, variant = factor(variant, c(unname(nm[top8]), "Others")), block = factor(block, 1:4))
PAL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948", "#C8C8C8")
g <- ggplot(pp, aes(block, pct, fill = variant)) + geom_col(width = .85, colour = "white", linewidth = .2) + facet_grid(site ~ genotype_label) +
  scale_fill_manual(values = PAL, name = "Pseudomonas ASV\n(closest type strain,\nputative)") + scale_x_discrete(drop = FALSE) +
  labs(x = "Block (1-4) within each genotype (A-C high root:shoot, D-F low)", y = "Pseudomonas, % of all reads",
       title = "Pseudomonas variants in each rhizosphere plot, by genotype and block") +
  theme_bw(12) + theme(legend.text = element_text(face = "italic"), strip.text = element_text(face = "bold"))
ggsave(file.path(OUT, "2_fig_pseudomonas_variants_per_plot.png"), g, width = 13, height = 6.5, dpi = 200)

txt <- c("PART 2 - PSEUDOMONAS VARIANTS", "",
         sprintf("2a BLAST: %d variants named (>= %.0f%% query coverage; ties within %.2f%% identity); %d tie with a non-Pseudomonas type strain",
                 nrow(ty), 100 * QCOV, TIE, sum(ty$non_pseudomonas_tie)),
         capture.output(print(as.data.frame(ty %>% filter(asv %in% top8) %>% select(asv, identity, n_tied, closest)), row.names = FALSE)),
         "", "2b Pseudomonas, % of all reads:", capture.output(print(as.data.frame(share), row.names = FALSE)),
         "", "2c Winner (most abundant Pseudomonas ASV) per rhizosphere plot:", capture.output(print(as.data.frame(ws), row.names = FALSE)),
         "", "2d Dominant Pseudomonas ASV in bulk soil (Pseudomonas is rare in bulk soil: see median reads):", capture.output(print(as.data.frame(bs), row.names = FALSE)),
         "", sprintf("2e PERMANOVA, Pseudomonas variant composition ~ genotype (Bray-Curtis, samples with >= %d Pseudomonas reads, %d permutations within block):", MINREADS, NPERM),
         capture.output(print(pg, row.names = FALSE)))
writeLines(txt, file.path(OUT, "2_summary.txt")); cat(txt, sep = "\n")
cat("\nPart 2 done. Outputs in", OUT, "\n")
