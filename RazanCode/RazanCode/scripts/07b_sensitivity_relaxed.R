# Step 7b (sensitivity): can the paper's healthy-sample results be explained by
#   "chloroplast/mitochondria removed, every sample with >= 1 bacterial read kept, NOT rarefied"?
#   Alpha from raw counts (vegan: observed, Chao1, Shannon ln, Simpson), Bray-Curtis on relative abundance,
#   PCoA + PERMANOVA status x cultivar -> compare with paper Fig. 2B / Fig. 3 / PERMANOVA.
# Run from the project root directory (env fb_r):
#   R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=$CONDA_PREFIX/lib/R/library Rscript scripts/07b_sensitivity_relaxed.R
suppressPackageStartupMessages({ library(vegan); library(ggplot2); library(dplyr); library(tidyr) })
set.seed(2026)

OUT <- "results/07b_sensitivity_relaxed"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
meta <- read.delim("metadata/metadata.tsv", comment.char = "#"); names(meta)[1] <- "sample"
meta$health_status <- factor(meta$health_status, levels = c("Canker", "Healthy"))
cult_col <- c("Aztec Fuji" = "#F8766D", "Cripps Pink" = "#B79F00", "Gala" = "#00BA38",
              "Granniwinkle" = "#00BFC4", "Hewe's Crab" = "#619CFF", "Unknown cultivar" = "#F564E3")

tab <- as.matrix(read.delim("results/04_taxonomy/export/feature-table.tsv", skip = 1, check.names = FALSE, row.names = 1))
x   <- t(tab); x <- x[rowSums(x) > 0, ]
md  <- meta[match(rownames(x), meta$sample), ]
cat("Samples kept:", nrow(x), "(", sum(md$health_status == "Canker"), "canker,", sum(md$health_status == "Healthy"), "healthy )\n")
cat("Healthy kept:", paste(md$label[md$health_status == "Healthy"], rowSums(x)[md$health_status == "Healthy"], sep = "=", collapse = ", "), "\n")

# ---- alpha (raw counts) ----
er <- estimateR(x)
al <- data.frame(sample = rownames(x), reads = rowSums(x), Observed = er["S.obs", ], Chao1 = er["S.chao1", ],
                 Shannon = diversity(x, "shannon"), Simpson = diversity(x, "simpson")) |> left_join(md, by = "sample")
write.table(al, file.path(OUT, "alpha_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
alpha_sum <- al |> pivot_longer(c(Observed, Chao1, Shannon, Simpson), names_to = "metric") |>
  group_by(metric, health_status) |>
  summarise(mean = round(mean(value), 2), min = round(min(value), 2), max = round(max(value), 2), .groups = "drop")
alpha_p <- sapply(c("Observed", "Chao1", "Shannon", "Simpson"),
                  function(m) signif(wilcox.test(al[[m]] ~ al$health_status)$p.value, 3))

# ---- beta (relative abundance) ----
rel <- x / rowSums(x)
bc  <- vegdist(rel, "bray")
pc  <- cmdscale(bc, k = 2, eig = TRUE); pv <- 100 * pc$eig / sum(pc$eig[pc$eig > 0])
sc  <- data.frame(sample = rownames(x), PC1 = pc$points[, 1], PC2 = pc$points[, 2]) |> left_join(md, by = "sample")
a1 <- adonis2(bc ~ health_status * cultivar, data = md, permutations = 999, by = "terms")
a2 <- adonis2(bc ~ cultivar, data = md, permutations = 999)
bd <- permutest(betadisper(bc, md$health_status), permutations = 999)
perm <- rbind(data.frame(model = "status * cultivar", term = rownames(a1), R2 = round(a1$R2, 4), F = round(a1$F, 3), P = a1$`Pr(>F)`),
              data.frame(model = "cultivar alone", term = rownames(a2), R2 = round(a2$R2, 4), F = round(a2$F, 3), P = a2$`Pr(>F)`))
write.table(perm, file.path(OUT, "permanova.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf("results/pdf/07b_sensitivity_relaxed.pdf", width = 9, height = 6)
print(ggplot(sc, aes(PC1, PC2)) + geom_point(aes(colour = cultivar, shape = health_status), size = 2.6) +
  stat_ellipse(aes(group = health_status), linetype = 2) +
  scale_colour_manual(values = cult_col) + scale_shape_manual(values = c(Canker = 16, Healthy = 17)) +
  labs(x = sprintf("[%.1f%%]", pv[1]), y = sprintf("[%.1f%%]", pv[2]), colour = "Cultivar", shape = "Health status",
       title = "Bray-Curtis PCoA - bacteria only, all samples with >=1 read, relative abundance (not rarefied)") +
  theme_bw() + coord_equal())
print(al |> pivot_longer(c(Observed, Chao1, Shannon, Simpson), names_to = "metric") |>
  ggplot(aes(cultivar, value, colour = cultivar)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15) +
  facet_grid(metric ~ health_status, scales = "free") + scale_colour_manual(values = cult_col) +
  labs(x = NULL, y = NULL, title = "Alpha diversity (raw counts) - bacteria only, not rarefied") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none"))
dev.off()

cat("\nAlpha by status:\n"); print(alpha_sum, n = 20)
cat("Wilcoxon canker vs healthy P:\n"); print(alpha_p)
cat("\nPCoA axes:", round(pv[1], 1), "/", round(pv[2], 1), " (paper 23.7 / 6.9)\n")
cat("\nPERMANOVA:\n"); print(perm, row.names = FALSE)
cat("\nbetadisper status: F =", round(bd$tab$F[1], 2), " P =", bd$tab$`Pr(>F)`[1], "\n")
