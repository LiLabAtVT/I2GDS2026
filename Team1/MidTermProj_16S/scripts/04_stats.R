# Step 4: downstream statistics mirroring the paper (alpha ANOVA, Bray-Curtis PCoA + PERMANOVA,
# genus heatmap, Random Forest, LEfSe, FAPROTAX, Dirichlet-multinomial). All figures -> results/pdf/.
suppressPackageStartupMessages({
  library(ggplot2); library(reshape2); library(dplyr); library(magrittr); library(vegan)
  library(randomForest); library(microeco); library(DirichletMultinomial); library(pheatmap)
})
setwd("/projects/intro2gds/I2GDS2026/Team1/Kamal/Claude_Code")
set.seed(123)
ex <- "results/qiime2/export"; out <- "results/pdf"; tab <- "results/tables"
dir.create(out, TRUE, FALSE); dir.create(tab, TRUE, FALSE)

meta <- read.delim("metadata/sample-metadata.tsv", row.names = 1)
meta$status <- factor(meta$status, c("Healthy", "Canker"))
meta$cultivar <- factor(meta$cultivar)
otu <- as.matrix(read.delim(file.path(ex, "rarefied_table.tsv"), skip = 1, row.names = 1, check.names = FALSE))
meta <- meta[colnames(otu), ]
tax <- read.delim(file.path(ex, "taxonomy/taxonomy.tsv"), row.names = 1)
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxm <- t(sapply(strsplit(tax$Taxon, "; "), function(x) { x <- sub("^[a-z]__", "", x); length(x) <- 7; x }))
dimnames(taxm) <- list(rownames(tax), ranks); taxm[is.na(taxm) | taxm == ""] <- "Unassigned"
taxm <- as.data.frame(taxm)[rownames(otu), ]
cat("Samples:", ncol(otu), " ASVs (rarefied):", nrow(otu), " depth:", colSums(otu)[1], "\n")

## 1. Alpha diversity (QIIME2 metrics) + ANOVA ---------------------------------------------
alpha <- do.call(cbind, lapply(c("observed_features", "chao1", "shannon", "simpson", "faith_pd", "evenness"), function(m)
  read.delim(file.path(ex, paste0("alpha_", m), "alpha-diversity.tsv"), row.names = 1)))
alpha <- cbind(meta, alpha[rownames(meta), ])
write.csv(alpha, file.path(tab, "alpha_diversity.csv"))
metrics <- c("observed_features", "chao1", "shannon", "simpson", "faith_pd", "pielou_evenness")
aov_res <- do.call(rbind, lapply(metrics, function(m) {
  st <- summary(aov(as.formula(paste(m, "~ status")), alpha))[[1]]
  cvC <- summary(aov(as.formula(paste(m, "~ cultivar")), subset(alpha, status == "Canker")))[[1]]
  cvH <- summary(aov(as.formula(paste(m, "~ cultivar")), subset(alpha, status == "Healthy")))[[1]]
  data.frame(metric = m, F_status = st$`F value`[1], P_status = st$`Pr(>F)`[1],
             F_cultivar_canker = cvC$`F value`[1], P_cultivar_canker = cvC$`Pr(>F)`[1],
             F_cultivar_healthy = cvH$`F value`[1], P_cultivar_healthy = cvH$`Pr(>F)`[1])
}))
write.csv(aov_res, file.path(tab, "alpha_anova.csv"), row.names = FALSE); print(aov_res)
al <- melt(alpha[, c("cultivar", "status", metrics)], id.vars = c("cultivar", "status"))
pdf(file.path(out, "Fig2B_alpha_diversity.pdf"), 12, 7)
print(ggplot(al, aes(status, value, fill = status)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = 1) +
  facet_wrap(~variable, scales = "free_y", nrow = 2) + theme_bw() + labs(title = "Alpha diversity: canker vs healthy (rarefied)"))
print(ggplot(al, aes(cultivar, value, fill = status)) + geom_boxplot(outlier.shape = NA) +
  facet_wrap(~variable, scales = "free_y", nrow = 2) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Alpha diversity by cultivar"))
dev.off()

## 2. Composition barplots (Fig 2A style) ---------------------------------------------------
rel <- sweep(otu, 2, colSums(otu), "/")
agg <- function(rank) rowsum(rel, taxm[[rank]])
barplot_rank <- function(rank, n = 15) {
  a <- agg(rank); top <- names(sort(rowMeans(a), TRUE))[1:min(n, nrow(a))]
  a <- rbind(a[top, , drop = FALSE], Other = colSums(a[!rownames(a) %in% top, , drop = FALSE]))
  d <- melt(a, varnames = c("taxon", "sample")); d <- cbind(d, meta[as.character(d$sample), c("cultivar", "status")])
  d$taxon <- factor(d$taxon, rownames(a))
  ggplot(d, aes(sample, value, fill = taxon)) + geom_col() + facet_grid(~status + cultivar, scales = "free_x", space = "free_x") +
    theme_bw() + theme(axis.text.x = element_blank(), strip.text = element_text(size = 6)) +
    labs(y = "Relative abundance", title = paste("Top", n, rank))
}
pdf(file.path(out, "Fig2A_relative_abundance.pdf"), 16, 7)
for (r in c("Phylum", "Class", "Genus")) print(barplot_rank(r))
dev.off()

## 3. Bray-Curtis PCoA + PERMANOVA (Fig 3) --------------------------------------------------
bc <- vegdist(t(otu), "bray")
perm <- adonis2(bc ~ status * cultivar, data = meta, permutations = 999, by = "terms")
print(perm); write.csv(as.data.frame(perm), file.path(tab, "permanova_bray.csv"))
pc <- cmdscale(bc, k = 2, eig = TRUE); ve <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)
pd <- cbind(meta, PC1 = pc$points[, 1], PC2 = pc$points[, 2])
pdf(file.path(out, "Fig3_braycurtis_PCoA.pdf"), 8, 6)
print(ggplot(pd, aes(PC1, PC2, colour = status, shape = cultivar)) + geom_point(size = 3) + stat_ellipse(aes(group = status)) +
  theme_bw() + labs(x = paste0("PCoA1 (", ve[1], "%)"), y = paste0("PCoA2 (", ve[2], "%)"),
  title = sprintf("Bray-Curtis PCoA | PERMANOVA status R2=%.3f P=%.3f", perm$R2[1], perm$`Pr(>F)`[1])))
dev.off()

## 4. Genus heatmap of canker samples + Random Forest by cultivar (Fig 4) -------------------
can <- rownames(meta)[meta$status == "Canker"]
gen <- agg("Genus")[, can]
topg <- names(sort(rowMeans(gen), TRUE))[1:30]
pdf(file.path(out, "Fig4A_canker_genus_heatmap.pdf"), 10, 9)
pheatmap(log10(gen[topg, ] + 1e-4), annotation_col = meta[can, "cultivar", drop = FALSE],
         labels_col = meta[can, "sample_name"], main = "Top 30 genera in canker tissue (log10 rel. abundance)")
dev.off()
gk <- gen[rowSums(gen > 0) >= 2 & rownames(gen) != "Unassigned", ]
rf <- randomForest(x = t(gk), y = droplevels(meta[can, "cultivar"]), ntree = 1000, importance = TRUE)
print(rf); imp <- importance(rf); imp <- imp[order(-imp[, "MeanDecreaseGini"]), ]
write.csv(imp, file.path(tab, "randomforest_canker_cultivar_importance.csv"))
top <- rownames(imp)[1:20]
pdf(file.path(out, "Fig4B_randomforest_canker_cultivar.pdf"), 11, 7)
varImpPlot(rf, n.var = 20, main = sprintf("Random Forest (canker, cultivar) OOB error = %.1f%%", 100 * tail(rf$err.rate[, 1], 1)))
mg <- melt(gk[top, ], varnames = c("genus", "sample")); mg$cultivar <- meta[as.character(mg$sample), "cultivar"]
mg <- aggregate(value ~ genus + cultivar, mg, mean); mg$genus <- factor(mg$genus, rev(top))
print(ggplot(mg, aes(cultivar, genus, fill = value)) + geom_tile() + scale_fill_viridis_c(trans = "sqrt") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(fill = "mean rel. ab.", title = "Top-20 RF genera by cultivar"))
dev.off()

## 5. LEfSe canker vs healthy (Fig 5) -------------------------------------------------------
mt <- microtable$new(otu_table = as.data.frame(otu), tax_table = taxm, sample_table = meta)
mt$tidy_dataset()
ld <- trans_diff$new(dataset = mt, method = "lefse", group = "status", alpha = 0.05, lefse_subgroup = NULL)
write.csv(ld$res_diff, file.path(tab, "lefse_status.csv"), row.names = FALSE)
pdf(file.path(out, "Fig5_LEfSe_canker_vs_healthy.pdf"), 10, 12)
print(ld$plot_diff_bar(threshold = 2, use_number = 1:40))
print(ld$plot_diff_cladogram(use_taxa_num = 100, use_feature_num = 40, clade_label_level = 5))
dev.off()

## 6. FAPROTAX functional prediction (microeco) ---------------------------------------------
fn <- trans_func$new(mt); fn$cal_spe_func(prok_database = "FAPROTAX"); fn$cal_spe_func_perc(abundance_weighted = TRUE)
fp <- fn$res_spe_func_perc; fp <- fp[, colSums(fp) > 0, drop = FALSE]
write.csv(fp, file.path(tab, "faprotax_function_percent.csv"))
fpk <- fp[, order(-colMeans(fp))[1:min(30, ncol(fp))], drop = FALSE]
pdf(file.path(out, "Fig6_FAPROTAX_functions.pdf"), 12, 9)
pheatmap(t(fpk), annotation_col = meta[, c("status", "cultivar")], labels_col = meta[rownames(fpk), "sample_name"],
         fontsize_col = 6, main = "FAPROTAX predicted functions (% abundance-weighted)")
fl <- melt(as.matrix(fpk), varnames = c("sample", "func")); fl$status <- meta[as.character(fl$sample), "status"]
print(ggplot(fl, aes(func, value, fill = status)) + geom_boxplot(outlier.size = .5) + coord_flip() + theme_bw() +
  labs(y = "% of community", x = NULL))
dev.off()

## 7. Dirichlet-multinomial clustering, k = 1..7, report k = 4 (Fig 7) ----------------------
gcount <- t(rowsum(otu, taxm$Genus)); gcount <- gcount[, colSums(gcount) > 0]
fits <- lapply(1:7, function(k) dmn(gcount, k, verbose = FALSE))
lap <- sapply(fits, laplace); aic <- sapply(fits, AIC); bic <- sapply(fits, BIC)
write.csv(data.frame(k = 1:7, Laplace = lap, AIC = aic, BIC = bic), file.path(tab, "dmm_fit.csv"), row.names = FALSE)
f4 <- fits[[4]]; cl <- apply(mixture(f4), 1, which.max)
dmm_tab <- table(cluster = cl, group = paste(meta[names(cl), "status"], meta[names(cl), "cultivar"]))
write.csv(as.data.frame.matrix(dmm_tab), file.path(tab, "dmm_k4_assignments.csv"))
fd <- fitted(f4, scale = TRUE); drv <- unique(unlist(lapply(1:4, function(i) names(sort(fd[, i], TRUE))[1:8])))
pdf(file.path(out, "Fig7_DirichletMultinomial.pdf"), 9, 8)
plot(1:7, lap, type = "b", xlab = "Number of Dirichlet components (k)", ylab = "Model fit",
     main = "DMM model fit (Laplace; lower = better)"); lines(1:7, aic, type = "b", col = 2); lines(1:7, bic, type = "b", col = 4)
legend("topright", c("Laplace", "AIC", "BIC"), col = c(1, 2, 4), lty = 1)
pheatmap(fd[drv, ], cluster_cols = FALSE, main = "DMM k = 4: driver genera per component")
pheatmap(as.matrix(as.data.frame.matrix(dmm_tab)), cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
         number_format = "%d", main = "Samples per DMM component (k = 4)")
dev.off()
cat("04_stats.R done\n")
