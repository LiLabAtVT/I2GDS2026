# Step 6e: Fig. 5B/D (Dirichlet-multinomial mixture clustering of the 32 cankers, k = 1..7,
#          k = 4 reported) and Fig. 5C (Erwinia ASVs per cultivar)
source("scripts/R_common.R")
library(DirichletMultinomial); library(pheatmap)
ps      <- readRDS(file.path(RDS, "ps_filtered.rds"))
ps_rare <- readRDS(file.path(RDS, "ps_canker_rarefied.rds"))   # cankers, rarefied to 505
asv_tax <- read.csv(tab_out("ASV_taxonomy.csv"), row.names = 1)
asv_lab <- function(ids) {
  g <- asv_tax[ids, "Genus"]; s <- asv_tax[ids, "Species"]
  paste0(ids, " (", ifelse(is.na(g), "unclassified", g), ifelse(is.na(s), "", paste0(" ", s)), ")")
}

# ================= Fig. 5B/D: DMM on canker ASV counts =================
can <- prune_samples(sample_data(ps)$health == "Canker", ps)
can <- prune_taxa(taxa_sums(can) > 0, can)
cnt <- t(as(otu_table(can), "matrix"))                          # samples x ASVs
cnt <- cnt[, colSums(cnt > 0) >= 2]                             # drop singletons-by-sample
cat(sprintf("DMM input: %d canker samples x %d ASVs\n", nrow(cnt), ncol(cnt)))
set.seed(SEED)
fits <- lapply(1:7, function(k) dmn(cnt, k, verbose = FALSE, seed = SEED))
saveRDS(fits, file.path(RDS, "dmm_fits_k1-7.rds"))
fitstat <- data.frame(k = 1:7, Laplace = sapply(fits, laplace), AIC = sapply(fits, AIC), BIC = sapply(fits, BIC))
write.csv(fitstat, tab_out("Fig5B_DMM_model_fit.csv"), row.names = FALSE)
print(fitstat)
k_best <- which.min(fitstat$Laplace); K <- 4                     # paper reports k = 4
fit <- fits[[K]]
cl  <- apply(mixture(fit), 1, which.max)
md  <- as(sample_data(can), "data.frame")
memb <- data.frame(sample = names(cl), cluster = cl, cultivar = md[names(cl), "cultivar"])
write.csv(memb, tab_out("Fig5B_DMM_k4_cluster_membership.csv"), row.names = FALSE)
print(table(cluster = memb$cluster, cultivar = memb$cultivar))

# drivers: ASVs with the largest fitted weights in each component
fv <- fitted(fit, scale = TRUE)
drv <- bind_rows(lapply(seq_len(K), function(j) {
  o <- order(fv[, j], decreasing = TRUE)[1:10]
  data.frame(community_type = j, ASV = rownames(fv)[o], weight = fv[o, j],
             label = asv_lab(rownames(fv)[o]))
}))
write.csv(drv, tab_out("Fig5D_DMM_top10_drivers.csv"), row.names = FALSE)

p_fit <- ggplot(fitstat, aes(k, Laplace)) + geom_line() + geom_point(size = 2) +
  geom_point(data = fitstat[K, ], colour = "red", size = 3) +
  labs(title = "DMM model fit (lower = better)", subtitle = sprintf("minimum Laplace at k = %d; k = 4 shown as in the paper", k_best))
top_asv <- unique(drv$ASV)
rel <- sweep(cnt, 1, rowSums(cnt), "/")[, top_asv]
ord_s <- order(cl)
ann <- data.frame(cluster = factor(cl), cultivar = md[names(cl), "cultivar"], row.names = names(cl))
pdf(pdf_out("Fig5B_D_DMM_clustering.pdf"), width = 11, height = 8)
print(p_fit)
draw_pheatmap(t(sqrt(rel[ord_s, ])), cluster_cols = FALSE, cluster_rows = TRUE,
         annotation_col = ann, annotation_colors = list(cultivar = CULT_COL[levels(droplevels(ann$cultivar))]),
         labels_row = asv_lab(top_asv), fontsize_row = 6, show_colnames = TRUE, fontsize_col = 6,
         gaps_col = cumsum(table(cl[ord_s]))[-K],
         color = colorRampPalette(c("white", "#fdae61", "#d7191c", "#67000d"))(100),
         main = "Fig. 5B: DMM community types (k = 4) - sqrt relative abundance of driver ASVs")
for (j in seq_len(K)) {
  dj <- drv[drv$community_type == j, ]
  print(ggplot(dj, aes(weight, factor(label, rev(label)))) + geom_col(fill = "#2c3e50") +
          labs(x = "Fitted weight", y = NULL, title = sprintf("Fig. 5D: top 10 drivers, community type %d", j)))
}
dev.off()

# ================= Fig. 5C: Erwinia ASVs per cultivar (cankers) =================
can_r <- prune_samples(sample_data(ps_rare)$health == "Canker", ps_rare)
g <- as.character(tax_table(can_r)[, "Genus"])
erw <- prune_taxa(taxa_names(can_r)[!is.na(g) & g == "Erwinia"], can_r)
erw <- prune_taxa(taxa_sums(erw) > 0, erw)
cat("Erwinia ASVs in cankers:", ntaxa(erw), "\n")
writeLines(capture.output(print(table(as(tax_table(erw), "matrix")[, "Species"], useNA = "ifany"))),
           tab_out("Fig5C_Erwinia_species_labels.txt"))
E <- t(as(otu_table(erw), "matrix"))
emd <- as(sample_data(erw), "data.frame")
es <- data.frame(emd[, c("cultivar", "replicate")], Erwinia_reads = rowSums(E),
                 Erwinia_ASVs = rowSums(E > 0), Erwinia_Shannon = vegan::diversity(E, "shannon"))
write.csv(es, tab_out("Fig5C_Erwinia_per_sample.csv"))
esum <- es %>% group_by(cultivar) %>%
  summarise(across(c(Erwinia_reads, Erwinia_ASVs, Erwinia_Shannon), list(mean = mean, sd = sd)), .groups = "drop")
write.csv(esum, tab_out("Fig5C_Erwinia_per_cultivar.csv"), row.names = FALSE)
print(esum)
el <- es %>% pivot_longer(c(Erwinia_reads, Erwinia_ASVs, Erwinia_Shannon))
p5c <- ggplot(el, aes(cultivar, value, fill = cultivar)) +
  stat_summary(fun = mean, geom = "col", alpha = .8) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .3) +
  geom_jitter(width = .15, size = 1) + facet_wrap(~ name, scales = "free_y") +
  scale_fill_manual(values = CULT_COL) +
  labs(x = NULL, y = NULL, title = "Fig. 5C: Erwinia ASVs in cankers per cultivar (mean +/- SE; rarefied)") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
top_e <- names(sort(taxa_sums(erw), decreasing = TRUE))[1:min(10, ntaxa(erw))]
ed <- psmelt(prune_taxa(top_e, erw))
p5c2 <- ggplot(ed, aes(Sample, Abundance, fill = OTU)) + geom_col() +
  facet_grid(~ cultivar, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Reads (rarefied)", fill = "Erwinia ASV", title = "Top Erwinia ASVs per canker sample") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size = 6), strip.text = element_text(size = 7))
pdf(pdf_out("Fig5C_Erwinia_ASVs_by_cultivar.pdf"), width = 12, height = 5)
print(p5c); print(p5c2)
dev.off()

# FASTA of Erwinia + DMM driver ASVs for BLAST confirmation (the paper used NCBI BLAST)
seqs <- Biostrings::readDNAStringSet(tab_out("ASV_sequences.fasta"))
Biostrings::writeXStringSet(seqs[unique(c(top_e, drv$ASV))], tab_out("ASVs_for_BLAST.fasta"))
