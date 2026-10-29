# Step 6b: build phyloseq object from QIIME2 exports, rename ASVs, rarefy
source("scripts/R_common.R")
library(ape)

# ---- metadata ----
meta <- read.delim("metadata/sample_metadata.tsv", check.names = FALSE)
meta <- meta[!grepl("^#", meta[[1]]), ]
names(meta)[1] <- "sample_id"            # phyloseq would mangle "sample-id"
rownames(meta) <- meta$sample_id
meta$cultivar <- factor(meta$cultivar, CULTIVARS)
meta$health   <- factor(meta$health, c("Canker", "Healthy"))
meta$group    <- factor(paste(meta$cultivar, meta$health, sep = " | "))

# ---- ASV table (organelle-filtered) ----
otu <- read.delim(file.path(EXP, "feature_table.tsv"), skip = 1, row.names = 1, check.names = FALSE)
otu <- as.matrix(otu)

# ---- taxonomy (SILVA 138) ----
tx <- read.delim(file.path(EXP, "taxonomy/taxonomy.tsv"), check.names = FALSE)
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tt <- t(sapply(strsplit(tx$Taxon, ";\\s*"), function(v) { length(v) <- 7; v }))
tt <- sub("^[a-z]__", "", tt); tt[tt == "" | grepl("uncultured|unidentified|metagenome", tt, TRUE)] <- NA
dimnames(tt) <- list(tx[["Feature ID"]], ranks)

# ---- tree ----
tree <- read_tree(file.path(EXP, "rooted_tree/tree.nwk"))

ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tt[rownames(otu), ]),
               sample_data(meta[colnames(otu), ]), phy_tree(tree))

# ---- rename hash IDs to ASV0001.. (ordered by total abundance) ----
ord <- names(sort(taxa_sums(ps), decreasing = TRUE))
ps <- prune_taxa(ord, ps)
ps <- phyloseq(otu_table(ps)[ord, ], tax_table(ps)[ord, ], sample_data(ps), phy_tree(ps))
new_id <- sprintf("ASV%04d", seq_along(ord))
map <- data.frame(ASV = new_id, feature_id = ord)
taxa_names(ps) <- new_id

seqs <- Biostrings::readDNAStringSet(file.path(EXP, "rep_seqs/dna-sequences.fasta"))
seqs <- seqs[map$feature_id]; names(seqs) <- map$ASV
Biostrings::writeXStringSet(seqs, tab_out("ASV_sequences.fasta"))
map$sequence <- as.character(seqs)
write.csv(cbind(map[, 1:2], as.data.frame(tax_table(ps))), tab_out("ASV_taxonomy.csv"), row.names = FALSE)

# ---- summary: reads & ASVs per sample (Supplementary Fig. S1) ----
ss <- as(sample_data(ps), "data.frame") %>%
  mutate(reads = sample_sums(ps), n_ASV = colSums(otu_table(ps) > 0))
write.csv(ss, tab_out("reads_ASVs_per_sample.csv"), row.names = FALSE)
cat(sprintf("Filtered table: %d ASVs, %s reads in %d samples\n",
            ntaxa(ps), format(sum(ss$reads), big.mark = ","), nsamples(ps)))
print(ss %>% group_by(health) %>%
        summarise(n = n(), ASV_min = min(n_ASV), ASV_mean = round(mean(n_ASV)), ASV_max = max(n_ASV),
                  reads_min = min(reads), reads_median = median(reads)))

# ---- rarefaction (depth: env RARE_DEPTH, default = smallest sample) ----
# Healthy wood is almost entirely host chloroplast/mitochondria: 8 healthy samples have 0 bacterial
# reads and the rest 7-177, versus 505-76,056 in cankers. Split approach:
#   * canker-only analyses (cultivar effects): rarefy cankers to the smallest canker depth (505)
#   * canker vs healthy: unrarefied counts / relative abundance of all samples with reads
empty <- setdiff(meta$sample_id, sample_names(ps))
write.csv(meta[empty, ], tab_out("samples_without_bacterial_reads.csv"), row.names = FALSE)
cat("Samples with 0 bacterial reads after organelle removal:", length(empty), "-", paste(empty, collapse = ", "), "\n")

ps_can <- prune_samples(sample_data(ps)$health == "Canker", ps)
depth <- as.integer(Sys.getenv("RARE_DEPTH", min(sample_sums(ps_can))))
ps_rare <- rarefy_even_depth(ps_can, sample.size = depth, rngseed = SEED, replace = FALSE, verbose = FALSE)
cat(sprintf("Cankers rarefied to %d reads/sample: %d samples, %d ASVs retained\n",
            depth, nsamples(ps_rare), ntaxa(ps_rare)))

# rarefaction curves
rc <- vegan::rarecurve(t(as(otu_table(ps), "matrix")), step = 20, tidy = TRUE)
names(rc) <- c("sample", "sampled", "ASVs_observed")
rc <- left_join(rc, ss, by = c(sample = "sample_id"))
p_rc <- ggplot(rc, aes(sampled, ASVs_observed, group = sample, colour = health)) +
  geom_line(alpha = .7) + geom_vline(xintercept = depth, linetype = 2) +
  scale_colour_manual(values = HEALTH_COL) + scale_x_log10() +
  labs(title = "Rarefaction curves (dashed = canker rarefaction depth)", x = "Reads sampled (log10)")
p_s1 <- ggplot(ss, aes(cultivar, n_ASV, fill = health)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(.8)) +
  geom_point(position = position_jitterdodge(jitter.width = .15, dodge.width = .8), size = 1) +
  scale_fill_manual(values = HEALTH_COL) + scale_y_log10() +
  labs(x = NULL, y = "Unique ASVs per sample (log10)", title = "Supplementary Fig. S1: ASVs per sample") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
pdf(pdf_out("S1_ASVs_per_sample_rarefaction.pdf"), width = 11, height = 5)
print(p_s1 | p_rc)
dev.off()

saveRDS(ps, file.path(RDS, "ps_filtered.rds"))
saveRDS(ps_rare, file.path(RDS, "ps_canker_rarefied.rds"))
