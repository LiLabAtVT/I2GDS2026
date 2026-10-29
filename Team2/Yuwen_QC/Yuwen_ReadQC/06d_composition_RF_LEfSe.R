# Step 6d: Fig. 4A (top-20 genus heatmap in cankers), Fig. 4B (Random Forest by cultivar),
#          Fig. 5A (LEfSe canker vs healthy)
source("scripts/R_common.R")
library(pheatmap); library(randomForest); library(microeco); library(file2meco)
ps      <- readRDS(file.path(RDS, "ps_filtered.rds"))
ps_rare <- readRDS(file.path(RDS, "ps_canker_rarefied.rds"))   # cankers, rarefied to 505

# ---- genus-level relative abundance, canker samples ----
can <- prune_samples(sample_data(ps_rare)$health == "Canker", ps_rare)
can <- prune_taxa(taxa_sums(can) > 0, can)
gen <- tax_glom(can, "Genus", NArm = FALSE)
taxa_names(gen) <- make.unique(tax_label(tax_table(gen), "Genus"))
gen_rel <- transform_sample_counts(gen, function(x) 100 * x / sum(x))
G <- t(as(otu_table(gen_rel), "matrix"))                     # samples x genera (%)
md <- as(sample_data(gen_rel), "data.frame")

# ================= Fig. 4A: heatmap of 20 most abundant genera =================
top20 <- names(sort(colMeans(G), decreasing = TRUE))[1:20]
H <- sapply(CULTIVARS, function(cv) colMeans(G[md$cultivar == cv, top20, drop = FALSE]))
write.csv(H, tab_out("Fig4A_top20_genera_mean_relab_canker.csv"))
pdf(pdf_out("Fig4A_heatmap_top20_genera_canker.pdf"), width = 7, height = 7)
draw_pheatmap(log10(H + 0.01), cluster_cols = FALSE, cluster_rows = FALSE,
         color = colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(100),
         display_numbers = round(H, 2), number_color = "black", fontsize_number = 6,
         main = "Fig. 4A: top 20 genera in cankers\n(mean % relative abundance; colour = log10 %)")
dev.off()

# ================= Fig. 4B: Random Forest (canker, cultivar as class) =================
keep <- colnames(G)[colSums(G > 0) >= 2]                        # genera seen in >= 2 cankers
X <- G[, keep]; colnames(X) <- make.names(colnames(X), unique = TRUE)
y <- droplevels(md$cultivar)
set.seed(SEED)
rf <- randomForest(x = X, y = y, ntree = 1000, importance = TRUE)
print(rf)
imp <- data.frame(genus = keep, MeanDecreaseGini = rf$importance[, "MeanDecreaseGini"]) %>%
  arrange(-MeanDecreaseGini)
write.csv(imp, tab_out("Fig4B_randomForest_importance.csv"), row.names = FALSE)
write.csv(rf$confusion, tab_out("Fig4B_randomForest_confusion.csv"))
top8 <- head(imp, 8)
ra <- data.frame(G[, top8$genus, drop = FALSE] / 100, cultivar = y, check.names = FALSE) %>%
  pivot_longer(-cultivar, names_to = "genus") %>%
  group_by(genus, cultivar) %>% summarise(relab = mean(value), .groups = "drop") %>%
  mutate(genus = factor(genus, rev(top8$genus)))
p4b1 <- ggplot(top8, aes(MeanDecreaseGini, factor(genus, rev(genus)))) +
  geom_point(size = 3, colour = "#2c3e50") + labs(y = NULL, x = "Mean Decrease Gini")
p4b2 <- ggplot(ra, aes(relab, genus, fill = cultivar)) + geom_col() +
  scale_fill_manual(values = CULT_COL) + labs(y = NULL, x = "Relative abundance") +
  theme(axis.text.y = element_blank())
pdf(pdf_out("Fig4B_randomForest_canker_by_cultivar.pdf"), width = 10, height = 4.5)
print((p4b1 | p4b2) + plot_annotation(
  title = "Fig. 4B: Random Forest - genera that discriminate cankers of each cultivar",
  subtitle = sprintf("ntree = 1000; OOB error = %.1f%%", 100 * tail(rf$err.rate[, "OOB"], 1))))
dev.off()

# ================= Fig. 5A: LEfSe (canker vs healthy) =================
mt <- to_meco(ps)
saveRDS(mt, file.path(RDS, "microtable_filtered.rds"))

set.seed(SEED)
lf <- trans_diff$new(dataset = mt, method = "lefse", group = "health",
                     alpha = 0.05, lefse_subgroup = NULL, taxa_level = "all", p_adjust_method = "none")
res <- lf$res_diff
write.csv(res, tab_out("Fig5A_LEfSe_canker_vs_healthy.csv"), row.names = FALSE)
cat("LEfSe features (P < 0.05):", nrow(res), "; LDA > 4:", sum(res$LDA > 4), "\n")
pdf(pdf_out("Fig5A_LEfSe_canker_vs_healthy.pdf"), width = 10, height = 9)
# Fig. 5A: LEfSe, canker vs healthy (LDA score log10 > 4); microeco returns composite plots
print(lf$plot_diff_bar(threshold = 4, use_number = 1:40, color_values = unname(HEALTH_COL)))
print(lf$plot_diff_abund(use_number = 1:20, color_values = unname(HEALTH_COL)))
dev.off()
