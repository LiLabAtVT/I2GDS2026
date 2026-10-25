# Step 10b: E. amylovora (Ea) ASVs per canker sample (paper Fig. 5C) + Dirichlet-multinomial mixtures (Fig. 5B, 5D)
#   paper: DMM (Holmes et al. 2012), count data fitted for k = 1-7, results for k = 4 shown (clusters of 4, 5, 11, 12 samples,
#          n = 32); 4 Ea ASVs are the main community drivers, followed by Curtobacterium flaccumfaciens, Novosphingobium fluoreni,
#          Methylobacterium, Microbacterium (species by BLAST)
#   Table (user decision): relaxed = bacteria only; DMM on the 32 cankers, ASVs with >= 50 reads in cankers.
#   Needs results of scripts/10a_asv_ids_blast.sh (ASV numbers + BLAST).
# Run via:  sbatch scripts/10_erwinia_dmm.sh
suppressPackageStartupMessages({ library(DirichletMultinomial); library(ggplot2); library(dplyr); library(tidyr); library(parallel) })
set.seed(2026)

OUT  <- "results/10_erwinia_dmm"
meta <- read.delim("metadata/metadata.tsv", comment.char = "#"); names(meta)[1] <- "sample"
map  <- read.delim(file.path(OUT, "asv_id_map.tsv"))
tab  <- as.matrix(read.delim("results/04_taxonomy/export/feature-table.tsv", skip = 1, check.names = FALSE, row.names = 1))
rownames(tab) <- map$asv_num[match(rownames(tab), map$feature_id)]
cult_order <- c("Cripps Pink", "Aztec Fuji", "Gala", "Hewe's Crab", "Granniwinkle", "Unknown cultivar")
cult_col <- c("Cripps Pink" = "#1B9E77", "Aztec Fuji" = "#D95F02", "Gala" = "#7570B3",
              "Hewe's Crab" = "#E7298A", "Granniwinkle" = "#66A61E", "Unknown cultivar" = "#E6AB02")

# ---- BLAST best hit per ASV -> species label ----
bl <- read.delim(file.path(OUT, "blast_top5.tsv"), header = FALSE,
                 col.names = c("asv", "sseqid", "pident", "length", "qlen", "slen", "mismatch", "evalue", "bitscore", "stitle"))
best <- bl |> group_by(asv) |> arrange(desc(bitscore), desc(pident), .by_group = TRUE) |> slice(1) |> ungroup() |>
  mutate(species = sapply(strsplit(stitle, " "), function(w) paste(w[1:2], collapse = " ")),
         coverage = round(100 * length / qlen, 1))
map$blast_species <- best$species[match(map$asv_num, best$asv)]
map$blast_pident  <- best$pident[match(map$asv_num, best$asv)]
map$label <- ifelse(!is.na(map$blast_species), paste0(map$blast_species, " (", map$blast_pident, "%)"), map$genus)
write.table(map, file.path(OUT, "asv_id_map_blast.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ck  <- meta[meta$health_status == "Canker", ]
ck  <- ck[order(match(ck$cultivar, cult_order), as.numeric(ck$replicate)), ]
cnt <- tab[, ck$sample]
rel <- 100 * sweep(cnt, 2, colSums(cnt), "/")

# ---- Fig. 5C: Ea ASVs per canker sample ----
ea <- map$asv_num[grepl("Erwinia", map$genus) & grepl("amylovora", map$blast_species)]
cat("E. amylovora ASVs (BLAST):", paste(ea, collapse = ", "), "\n")
top4 <- head(ea, 4)
ea_df <- data.frame(sample = ck$sample, label = sub("_Canker", "", ck$label), cultivar = ck$cultivar)
for (a in top4) ea_df[[a]] <- rel[a, ]
ea_df$`Other Ea ASVs` <- colSums(rel[setdiff(ea, top4), , drop = FALSE])
ea_df$total_Ea <- rowSums(ea_df[, c(top4, "Other Ea ASVs")])
ea_df$n_Ea_ASVs <- colSums(cnt[ea, ] > 0)
ea_df$shannon_Ea <- apply(cnt[ea, ], 2, function(x) { p <- x[x > 0] / sum(x); if (sum(x) == 0) 0 else -sum(p * log(p)) })
write.table(ea_df, file.path(OUT, "ea_asvs_per_canker.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
long <- ea_df |> pivot_longer(c(all_of(top4), "Other Ea ASVs"), names_to = "ASV", values_to = "pct") |>
  mutate(label = factor(label, levels = ea_df$label[order(ea_df$total_Ea)]),
         ASV = factor(ASV, levels = rev(c(top4, "Other Ea ASVs"))))
pdf("results/pdf/10_erwinia_amylovora_asvs_per_canker.pdf", width = 11, height = 5)
print(ggplot(long, aes(label, pct, fill = ASV)) + geom_col() +
  scale_fill_manual(values = setNames(c("#1F4E79", "#E07B39", "#1B7837", "#2E9FDF", "grey70"), c(top4, "Other Ea ASVs"))) +
  labs(x = NULL, y = "% of bacterial reads", fill = NULL, title = "Erwinia amylovora ASVs per canker sample (BLAST-confirmed)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "top"))
cs <- ea_df |> group_by(cultivar) |> summarise(mean_Ea_pct = mean(total_Ea), mean_n_ASVs = mean(n_Ea_ASVs),
                                               mean_shannon = mean(shannon_Ea), .groups = "drop") |>
  mutate(cultivar = factor(cultivar, levels = cult_order))
print(ggplot(ea_df, aes(factor(cultivar, levels = cult_order), total_Ea, colour = cultivar)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 2) + scale_colour_manual(values = cult_col) +
  labs(x = NULL, y = "E. amylovora, % of bacterial reads", title = "E. amylovora abundance by cultivar (cankers)") +
  theme_bw() + theme(legend.position = "none"))
dev.off()
kw <- kruskal.test(total_Ea ~ cultivar, data = ea_df)

# ---- DMM ----
keep <- rowSums(cnt) >= 50
X <- t(cnt[keep, ])                                             # samples x ASVs
cat("DMM input:", nrow(X), "samples x", ncol(X), "ASVs (>= 50 reads in cankers)\n")
fits <- mclapply(1:7, function(k) dmn(X, k, verbose = FALSE, seed = 2026), mc.cores = 7)
fit_tab <- data.frame(k = 1:7, Laplace = sapply(fits, laplace), AIC = sapply(fits, AIC), BIC = sapply(fits, BIC))
write.table(fit_tab, file.path(OUT, "dmm_model_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
kbest <- which.min(fit_tab$Laplace)
saveRDS(fits, file.path(OUT, "dmm_fits.rds"))

f4 <- fits[[4]]
cl <- apply(mixture(f4), 1, which.max)
ck$dmm_k4 <- cl[ck$sample]; ck$dmm_best <- apply(mixture(fits[[kbest]]), 1, which.max)[ck$sample]
write.table(ck[, c("sample", "label", "cultivar", "dmm_k4", "dmm_best")], file.path(OUT, "dmm_clusters.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p0 <- fitted(fits[[1]], scale = TRUE)
p4 <- fitted(f4, scale = TRUE)
colnames(p4) <- paste0("type", 1:4)
drv <- data.frame(ASV = rownames(p4), p4, overall_diff = rowSums(abs(p4 - as.vector(p0))))
drv$label <- map$label[match(drv$ASV, map$asv_num)]
drv <- drv |> arrange(desc(overall_diff))
write.table(drv, file.path(OUT, "dmm_k4_drivers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# community type containing the most Erwinia-dominated samples = the paper's "community type 1"
ea_rank <- sapply(1:4, function(i) mean(ea_df$total_Ea[match(names(cl)[cl == i], ea_df$sample)]))
ct <- which.max(ea_rank)
d1 <- data.frame(ASV = rownames(p4), value = p4[, ct] - p0[, 1]) |> arrange(desc(value)) |> head(10)
d1$label <- map$label[match(d1$ASV, map$asv_num)]
d1$ASV <- factor(d1$ASV, levels = rev(d1$ASV))

pdf("results/pdf/10_dmm_clustering.pdf", width = 10, height = 6)
print(ggplot(fit_tab, aes(k, Laplace)) + geom_line() + geom_point(size = 2) +
  geom_point(data = fit_tab[kbest, ], colour = "red", size = 3.5) +
  labs(title = sprintf("DMM model fit (lowest Laplace: k = %d)", kbest), x = "Number of Dirichlet components (k)") + theme_bw())
heatmapdmn(X, fits[[1]], f4, ntaxa = 30)
title(sprintf("DMM k = 4 - cluster sizes: %s", paste(table(factor(cl, levels = 1:4)), collapse = ", ")), line = 3)
print(ggplot(d1, aes(value, ASV)) + geom_col(fill = "grey30") +
  geom_text(aes(x = 0, label = paste0("  ", label)), hjust = 0, colour = "white", size = 3) +
  labs(x = "Difference in fitted proportion vs 1-component model", y = NULL,
       title = sprintf("Top 10 drivers: community type %d (Erwinia-dominated, k = 4)", ct)) + theme_bw())
print(ggplot(ck, aes(factor(dmm_k4), fill = cultivar)) + geom_bar() + scale_fill_manual(values = cult_col) +
  labs(x = "DMM community type (k = 4)", y = "Canker samples", title = "Cultivar composition of DMM clusters") + theme_bw())
dev.off()

# ---- summary ----
cat("\nEa ASVs total % per canker: min", round(min(ea_df$total_Ea), 1), " median", round(median(ea_df$total_Ea), 1),
    " max", round(max(ea_df$total_Ea), 1), "\n")
print(cs); cat("Kruskal-Wallis Ea % ~ cultivar: P =", signif(kw$p.value, 3), "\n")
cat("\nDMM fit:\n"); print(fit_tab)
cat("k = 4 cluster sizes:", table(cl), " (paper: 4, 5, 11, 12)\n")
print(table(ck$dmm_k4, ck$cultivar))
cat("\nTop 10 drivers, community type", ct, ":\n"); print(d1[, c("ASV", "value", "label")], row.names = FALSE)
cat("\nTop 10 overall drivers (k = 4):\n"); print(head(drv[, c("ASV", "overall_diff", "label")], 10), row.names = FALSE)
