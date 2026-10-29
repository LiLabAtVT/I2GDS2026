#!/usr/bin/env Rscript
# Step 7: community composition and machine-learning analyses on the bacteria-only table
# (chloroplast, mitochondria and non-bacterial ASVs removed; not rarefied; relative abundances).
# Reproduces paper Fig. 2A (classes), 4A (top-20 genus heatmap, canker), 4B (Random Forest by cultivar),
# 5A (LEfSe canker vs healthy), 5B/D (Dirichlet-multinomial k = 1-7, k = 4 shown), 5C (Erwinia ASVs).
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(pheatmap); library(randomForest)
  library(DirichletMultinomial); library(phyloseq); library(microbiomeMarker); library(parallel)
})
set.seed(2026)
B <- "/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code"
O <- file.path(B, "results/07_composition_ml"); dir.create(O, showWarnings = FALSE, recursive = TRUE)
NCORE <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
CULT <- c("Cripps Pink", "Aztec Fuji", "Gala", "Hewe's Crab", "Granniwinkle", "Unknown cultivar")
CCOL <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"), CULT)

# ---- Inputs -----------------------------------------------------------------
meta <- read.delim(file.path(B, "raw_data/sample_metadata.tsv"), check.names = FALSE)
meta <- meta[meta[["sample-id"]] != "#q2:types", ]
rownames(meta) <- meta[["sample-id"]]
meta$cultivar <- factor(meta$cultivar, CULT); meta$health_status <- factor(meta$health_status, c("Canker", "Healthy"))
ft <- as.matrix(read.delim(file.path(B, "results/04_dada2/export/feature-table.tsv"), skip = 1,
                           row.names = 1, check.names = FALSE, comment.char = ""))
tax <- read.delim(file.path(B, "results/05_taxonomy/taxonomy_by_rank.tsv"), row.names = 1, check.names = FALSE)
tax <- tax[rownames(ft), ]
keep <- tax$Domain == "Bacteria" & !(tax$Order %in% "Chloroplast") & !(tax$Family %in% "Mitochondria")
ft <- ft[keep, ]; tax <- tax[keep, ]
for (r in c("Phylum", "Class", "Order", "Family", "Genus")) tax[[r]][is.na(tax[[r]]) | tax[[r]] == ""] <- "Unassigned"

# Rename ASVs ASV0001... by total abundance; keep the lookup
ord <- order(-rowSums(ft)); ft <- ft[ord, ]; tax <- tax[ord, ]
ids <- sprintf("ASV%04d", seq_len(nrow(ft)))
write.table(data.frame(ASV = ids, feature_id = rownames(ft), tax[, c("Phylum", "Class", "Order", "Family", "Genus", "Confidence")],
                       total_reads = rowSums(ft)), file.path(O, "asv_lookup.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
rownames(ft) <- ids; rownames(tax) <- ids

# Drop samples with no bacterial reads; order samples as in the paper
dropped <- colnames(ft)[colSums(ft) == 0]
ft <- ft[, colSums(ft) > 0]
meta <- meta[colnames(ft), ]
so <- rownames(meta)[order(meta$health_status, meta$cultivar, as.integer(meta$replicate))]
ft <- ft[, so]; meta <- meta[so, ]
rel <- sweep(ft, 2, colSums(ft), "/") * 100
cat("Bacterial ASVs:", nrow(ft), " samples kept:", ncol(ft), " dropped (0 bacterial reads):", paste(dropped, collapse = ", "), "\n")
isC <- meta$health_status == "Canker"

collapse <- function(m, rank) rowsum(m, tax[rownames(m), rank])
top_other <- function(g, n) { k <- names(sort(rowMeans(g), decreasing = TRUE))[seq_len(min(n, nrow(g)))]
  rbind(g[k, , drop = FALSE], Other = colSums(g[setdiff(rownames(g), k), , drop = FALSE])) }

pdf(file.path(B, "pdf/07_composition_ml.pdf"), width = 14, height = 8)
grid::grid.newpage()
grid::grid.text(paste0("Step 7 - composition and machine learning (bacteria only, not rarefied)\n\n",
  "Bacterial ASVs: ", nrow(ft), "   Samples: ", sum(isC), " canker + ", sum(!isC), " healthy\n",
  "Healthy samples dropped (0 bacterial reads after chloroplast removal):\n", paste(dropped, collapse = ", "),
  "\n\nASV IDs are ours (ranked by abundance); see results/07_composition_ml/asv_lookup.tsv"), gp = grid::gpar(fontsize = 13))

# ---- 1. Fig. 2A: class-level relative abundance ----------------------------
cls <- top_other(collapse(rel, "Class"), 20)
write.table(round(cls, 4), file.path(O, "class_relative_abundance.tsv"), sep = "\t", quote = FALSE)
d <- as.data.frame(cls) %>% mutate(Class = rownames(cls)) %>% pivot_longer(-Class, names_to = "sample", values_to = "pct") %>%
  left_join(meta %>% mutate(sample = rownames(meta)) %>% select(sample, cultivar, health_status), by = "sample") %>%
  mutate(sample = factor(sample, so), Class = factor(Class, rev(rownames(cls))))
pal <- c(setNames(grDevices::hcl.colors(nrow(cls) - 1, "Dark 3"), rownames(cls)[-nrow(cls)]), Other = "grey80")
print(ggplot(d, aes(sample, pct, fill = Class)) + geom_col(width = 0.9) + scale_fill_manual(values = pal) +
  facet_grid(~ health_status + cultivar, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Relative abundance (%)", title = "Fig. 2A reproduction: class-level composition (bacteria only)") +
  theme_bw(10) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
                       strip.text = element_text(size = 7), legend.text = element_text(size = 8)))

# ---- 2. Fig. 4A: top-20 genus heatmap, canker samples ----------------------
genC <- collapse(rel[, isC], "Genus"); genC <- genC[rownames(genC) != "Unassigned", ]
top20 <- names(sort(rowMeans(genC), decreasing = TRUE))[1:20]
write.table(round(genC[top20, ], 4), file.path(O, "canker_top20_genera.tsv"), sep = "\t", quote = FALSE)
mC <- meta[isC, ]
ph <- pheatmap(log10(genC[top20, ] + 0.01), cluster_rows = FALSE, cluster_cols = FALSE,
  gaps_col = cumsum(table(droplevels(mC$cultivar))), annotation_col = data.frame(Cultivar = mC$cultivar, row.names = rownames(mC)),
  annotation_colors = list(Cultivar = CCOL[levels(droplevels(mC$cultivar))]),
  color = colorRampPalette(c("navy", "lightyellow", "firebrick"))(100), breaks = seq(-2, 2, length.out = 101),
  legend_breaks = c(-2, -1, 0, 1, 2), legend_labels = c("0.01", "0.1", "1", "10", "100"),
  main = "Fig. 4A reproduction: top 20 genera in canker samples (% relative abundance, log scale)", fontsize_col = 7, silent = TRUE)
grid::grid.newpage(); grid::grid.draw(ph$gtable)

# ---- 3. Fig. 4B: Random Forest, canker genus abundances ~ cultivar ---------
X <- t(genC); X <- X[, colSums(X > 0) >= 2, drop = FALSE]; colnames(X) <- make.names(colnames(X))
y <- droplevels(mC$cultivar)
rf <- randomForest(x = X, y = y, ntree = 1000, importance = TRUE)
imp <- data.frame(genus = colnames(X), MeanDecreaseGini = rf$importance[, "MeanDecreaseGini"]) %>% arrange(desc(MeanDecreaseGini))
write.table(imp, file.path(O, "random_forest_importance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
oob <- tail(rf$err.rate[, "OOB"], 1)
cat("Random Forest OOB error:", round(oob, 3), "\n"); print(rf$confusion)
top8 <- head(imp, 8)
ab <- as.data.frame(X[, top8$genus, drop = FALSE]) %>% mutate(cultivar = y) %>% pivot_longer(-cultivar, names_to = "genus") %>%
  group_by(genus, cultivar) %>% summarise(mean = mean(value), se = sd(value) / sqrt(n()), .groups = "drop") %>%
  mutate(genus = factor(genus, rev(top8$genus)))
p1 <- ggplot(top8, aes(MeanDecreaseGini, factor(genus, rev(genus)))) + geom_col(fill = "grey40") +
  labs(y = NULL, title = sprintf("Random Forest (cultivar, canker; OOB error %.0f%%)", 100 * oob)) + theme_bw(10)
p2 <- ggplot(ab, aes(mean, genus, fill = cultivar)) + geom_col(position = position_dodge(0.8), width = 0.8) +
  geom_errorbar(aes(xmin = pmax(mean - se, 0), xmax = mean + se), position = position_dodge(0.8), width = 0.3) +
  scale_fill_manual(values = CCOL) + labs(x = "Mean relative abundance (%)", y = NULL, title = "Relative abundance by cultivar") + theme_bw(10)
grid::grid.newpage(); vp <- grid::viewport(layout = grid::grid.layout(1, 2))
grid::pushViewport(vp); print(p1, vp = grid::viewport(layout.pos.col = 1)); print(p2, vp = grid::viewport(layout.pos.col = 2)); grid::popViewport()

# ---- 4. Fig. 5A: LEfSe, canker vs healthy (genus level) --------------------
ps <- phyloseq(otu_table(ft, taxa_are_rows = TRUE), tax_table(as.matrix(tax[, c("Phylum", "Class", "Order", "Family", "Genus")])),
               sample_data(meta[, c("cultivar", "health_status")]))
# LEfSe bootstraps the LDA step; with only 8 healthy samples some random draws return no markers,
# so the seed is fixed here and stability is checked over 10 seeds.
# marker_table is an S4 data.frame subclass; convert to a plain data.frame so write.table/ggplot work
plain_df <- function(m) { d <- setNames(as.data.frame(m@.Data, stringsAsFactors = FALSE), m@names); d$ef_lda <- as.numeric(d$ef_lda); d }
lefse_run <- function(seed) { set.seed(seed)
  tryCatch(plain_df(marker_table(run_lefse(ps, group = "health_status", taxa_rank = "Genus", norm = "CPM",
                                                 kw_cutoff = 0.05, lda_cutoff = 2))), error = function(e) NULL) }
mt <- lefse_run(2026); if (is.null(mt)) mt <- data.frame(feature = character(), enrich_group = character(), ef_lda = numeric())
stab <- bind_rows(lapply(1:10, function(s) { m <- lefse_run(s); if (is.null(m)) data.frame(seed = s, feature = NA, enrich_group = NA)
                                             else data.frame(seed = s, feature = m$feature, enrich_group = m$enrich_group) }))
stab_sum <- stab %>% filter(!is.na(feature)) %>% count(feature, enrich_group, name = "seeds_detected_of_10") %>% arrange(desc(seeds_detected_of_10))
write.table(stab_sum, file.path(O, "lefse_genus_stability_10seeds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("LEfSe markers per seed:", paste(sapply(1:10, function(s) sum(stab$seed == s & !is.na(stab$feature))), collapse = " "), "\n")
write.table(mt, file.path(O, "lefse_genus.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("LEfSe markers:", nrow(mt), "\n"); print(head(mt[order(-mt$ef_lda), ], 30))
if (nrow(mt) > 0) {
  mt2 <- head(mt[order(-mt$ef_lda), ], 35)
  print(ggplot(mt2, aes(ef_lda, reorder(feature, ef_lda), colour = enrich_group, size = ef_lda)) + geom_point() +
    scale_colour_manual(values = c(Canker = "#d6604d", Healthy = "#4daf4a")) + scale_size(range = c(2, 6)) +
    labs(x = "LDA score (log10)", y = NULL, colour = "Enriched in", size = "LDA",
         title = sprintf("Fig. 5A reproduction: LEfSe genus level (KW P<0.05, LDA>2; %d canker vs %d healthy)", sum(isC), sum(!isC))) + theme_bw(10))
}

# ---- 5. Fig. 5B/D: Dirichlet-multinomial mixtures, canker samples (ASV level) ----
cnt <- t(ft[, isC]); cnt <- cnt[, colSums(cnt > 0) >= 4, drop = FALSE]   # ASVs in >= 4 of 32 canker samples
cat("DMM input:", nrow(cnt), "samples x", ncol(cnt), "ASVs\n")
fits <- mclapply(1:7, function(k) { set.seed(2026); dmn(cnt, k, verbose = FALSE) }, mc.cores = min(7, NCORE))
fitstat <- data.frame(k = 1:7, Laplace = sapply(fits, laplace), AIC = sapply(fits, AIC), BIC = sapply(fits, BIC))
write.table(fitstat, file.path(O, "dmm_model_fit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(fitstat)
print(ggplot(pivot_longer(fitstat, -k), aes(k, value)) + geom_line() + geom_point() + facet_wrap(~ name, scales = "free_y") +
  labs(title = sprintf("DMM model fit, k = 1-7 (best by Laplace: k = %d; paper presents k = 4)", fitstat$k[which.min(fitstat$Laplace)])) + theme_bw(11))
k4 <- fits[[4]]
cl <- apply(mixture(k4), 1, which.max)
cls_tab <- data.frame(sample = rownames(cnt), cluster = cl, cultivar = meta[rownames(cnt), "cultivar"])
write.table(cls_tab, file.path(O, "dmm_k4_clusters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("DMM k=4 cluster sizes (paper: 4, 5, 11, 12):", paste(sort(table(cl)), collapse = ", "), "\n"); print(table(cls_tab$cluster, cls_tab$cultivar))
# Drivers: ASVs whose fitted Dirichlet weight in each community type differs most from the k = 1 fit
p0 <- fitted(fits[[1]], scale = TRUE); p4 <- fitted(k4, scale = TRUE)
drivers <- bind_rows(lapply(seq_len(ncol(p4)), function(j) {
  dd <- p4[, j] - p0[, 1]; top <- head(order(-abs(dd)), 10)
  data.frame(community_type = j, n_samples = sum(cl == j), ASV = rownames(p4)[top], difference = dd[top],
             genus = tax[rownames(p4)[top], "Genus"], family = tax[rownames(p4)[top], "Family"])
}))
write.table(drivers, file.path(O, "dmm_k4_top10_drivers.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(drivers)
print(ggplot(drivers, aes(abs(difference), reorder(paste(ASV, genus), abs(difference)))) + geom_col(fill = "grey35") +
  facet_wrap(~ paste0("Community type ", community_type, " (n = ", n_samples, ")"), scales = "free") +
  labs(x = "|difference from k = 1 fit|", y = NULL, title = "Fig. 5D reproduction: top 10 drivers per DMM community type (k = 4)") + theme_bw(9))
topasv <- names(sort(colSums(cnt), decreasing = TRUE))[1:100]
hm <- pheatmap(t(sqrt(cnt[order(cl), topasv])), cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = FALSE,
  annotation_col = data.frame(DMM_type = factor(cl[order(cl)]), Cultivar = meta[rownames(cnt)[order(cl)], "cultivar"], row.names = rownames(cnt)[order(cl)]),
  annotation_colors = list(Cultivar = CCOL), gaps_col = cumsum(table(cl)), fontsize_col = 7,
  main = "Fig. 5B reproduction: top 100 ASVs (sqrt counts) in canker samples ordered by DMM community type (k = 4)", silent = TRUE)
grid::grid.newpage(); grid::grid.draw(hm$gtable)

# ---- 6. Fig. 5C: Erwinia ASVs per canker sample ---------------------------
erw <- rownames(tax)[tax$Genus == "Erwinia"]
erwC <- rel[erw, isC, drop = FALSE]
top4 <- names(sort(rowSums(ft[erw, isC, drop = FALSE]), decreasing = TRUE))[1:4]
e <- as.data.frame(t(erwC[top4, ])) %>% mutate(sample = rownames(.)) %>% pivot_longer(-sample, names_to = "ASV", values_to = "pct")
e$sample <- factor(e$sample, names(sort(colSums(erwC[top4, ]))))
print(ggplot(e, aes(sample, pct, fill = ASV)) + geom_col() + theme_bw(10) +
  labs(x = NULL, y = "% of bacterial reads", title = "Fig. 5C reproduction: the four most abundant Erwinia ASVs per canker sample") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)))
shan <- function(v) { p <- v[v > 0] / sum(v); if (length(p) == 0) 0 else -sum(p * log(p)) }
erwsum <- data.frame(sample = colnames(erwC), cultivar = meta[colnames(erwC), "cultivar"],
                     n_Erwinia_ASVs = colSums(ft[erw, isC, drop = FALSE] > 0), Erwinia_pct = colSums(erwC),
                     Erwinia_Shannon = apply(ft[erw, isC, drop = FALSE], 2, shan))
bycult <- erwsum %>% group_by(cultivar) %>% summarise(samples = n(), total_Erwinia_ASVs = length(unique(unlist(lapply(sample, function(s) erw[ft[erw, s] > 0])))),
  mean_ASVs_per_sample = mean(n_Erwinia_ASVs), mean_Erwinia_pct = mean(Erwinia_pct), mean_Shannon = mean(Erwinia_Shannon)) %>% arrange(desc(total_Erwinia_ASVs))
write.table(erwsum, file.path(O, "erwinia_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bycult, file.path(O, "erwinia_by_cultivar.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Erwinia ASVs total:", length(erw), "\n"); print(as.data.frame(bycult))
grid::grid.newpage()
grid::grid.text(paste0("Erwinia ASVs by cultivar (paper ranking: Cripps Pink > Gala > Aztec Fuji > Granniwinkle > Hewe's Crab > Unknown)\n\n",
  paste(capture.output(print(as.data.frame(bycult), digits = 3, row.names = FALSE)), collapse = "\n")),
  gp = grid::gpar(fontsize = 11, fontfamily = "mono"))
dev.off()
writeLines(capture.output(sessionInfo()), file.path(O, "sessionInfo.txt"))
cat("Step 7 done\n")
