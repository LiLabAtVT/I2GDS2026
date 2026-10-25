# Step 7: top-20 genus heatmap in canker samples (paper Fig. 4A)
#   paper: "Abundances of the top 20 most abundant bacterial genus levels" per canker sample, grouped by
#          cultivar, log colour scale (0.01-10 %); top 5 = Erwinia, Sphingomonas, Massilia, Pantoea, Novosphingobium
#   Relative abundance from the unrarefied tables (paper does not say rarefied); two versions:
#     paper    = all ASVs (chloroplast/mitochondria kept)
#     filtered = bacteria only
# Run from the project root directory (env fb_r):
#   R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=$CONDA_PREFIX/lib/R/library Rscript scripts/07_genus_heatmap.R
suppressPackageStartupMessages({ library(pheatmap); library(dplyr) })

OUT <- "results/07_composition"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
meta <- read.delim("metadata/metadata.tsv", comment.char = "#"); names(meta)[1] <- "sample"
tax  <- read.delim("results/04_taxonomy/export/taxonomy_ranks.tsv")
paper_top5 <- c("Erwinia", "Sphingomonas", "Massilia", "Pantoea", "Novosphingobium")
cult_order <- c("Cripps Pink", "Aztec Fuji", "Gala", "Hewe's Crab", "Granniwinkle", "Unknown cultivar")  # paper order
cult_col <- c("Cripps Pink" = "#1B9E77", "Aztec Fuji" = "#D95F02", "Gala" = "#7570B3",
              "Hewe's Crab" = "#E7298A", "Granniwinkle" = "#66A61E", "Unknown cultivar" = "#E6AB02")

genus_rel <- function(tab) {
  g <- sub("^g__", "", tax$Genus[match(rownames(tab), tax$asv)])
  g[is.na(g) | g == ""] <- "Unclassified"
  gt <- rowsum(tab, g)
  100 * sweep(gt, 2, pmax(colSums(gt), 1), "/")
}

make <- function(set, file) {
  tab <- as.matrix(read.delim(file, skip = 1, check.names = FALSE, row.names = 1))
  md  <- meta[meta$health_status == "Canker", ]
  md  <- md[order(match(md$cultivar, cult_order), as.numeric(md$replicate)), ]
  rel <- genus_rel(tab)[, md$sample]
  write.table(round(rel, 5), file.path(OUT, sprintf("genus_relative_abundance_cankers_%s.tsv", set)),
              sep = "\t", quote = FALSE)
  keep <- setdiff(rownames(rel), "Unclassified")
  top  <- names(sort(rowMeans(rel[keep, ]), decreasing = TRUE))[1:20]
  m    <- rel[top, ]
  colnames(m) <- md$label
  ann  <- data.frame(Cultivar = md$cultivar, row.names = md$label)
  gaps <- cumsum(table(factor(md$cultivar, levels = cult_order)))
  lm   <- log10(pmax(m, 0.01))                                   # floor at 0.01 % like the paper's scale
  pdf(sprintf("results/pdf/07_genus_heatmap_cankers_%s.pdf", set), width = 11, height = 6.5)
  pheatmap(lm, cluster_rows = FALSE, cluster_cols = FALSE, gaps_col = gaps[-length(gaps)],
           annotation_col = ann, annotation_colors = list(Cultivar = cult_col),
           color = colorRampPalette(c("#08306B", "#4292C6", "#FFFFBF", "#F46D43", "#A50026"))(100),
           breaks = seq(-2, 2, length.out = 101),
           legend_breaks = -2:2, legend_labels = c("0.01", "0.1", "1", "10", "100"),
           border_color = "grey60", fontsize_col = 7, fontsize_row = 9,
           main = sprintf("Top 20 genera in canker samples (%% relative abundance, log scale) - %s table", set))
  dev.off()
  mean_top <- data.frame(set = set, rank = 1:20, genus = top, mean_pct = round(rowMeans(rel[top, ]), 2),
                         in_paper_top5 = top %in% paper_top5)
  per_cult <- sapply(cult_order, function(cv) round(rowMeans(rel[top, md$sample[md$cultivar == cv], drop = FALSE]), 2))
  write.table(cbind(mean_top, per_cult), file.path(OUT, sprintf("top20_genera_cankers_%s.tsv", set)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  mean_top
}

res <- bind_rows(make("paper", "results/03_dada2/export/feature-table.tsv"),
                 make("filtered", "results/04_taxonomy/export/feature-table.tsv"))
print(res, row.names = FALSE)
