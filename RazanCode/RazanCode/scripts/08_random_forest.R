# Step 8: Random Forest - genera that discriminate apple cultivars in canker tissue (paper Fig. 4B)
#   paper: randomForest (Liaw & Wiener 2002) on relative abundance of canker microbiota; top 8 by Mean Decrease Gini:
#          Friedmanniella, Aurantisolimonas, Deinococcus (Cripps Pink); Rhizobium group (Aztec Fuji);
#          Cutibacterium, Gibbsiella (Gala); Nocardioides, Shimwellia (Hewe's Crab)
#   Table (user decision): relaxed = bacteria only, not rarefied; here the 32 canker samples, genus relative abundance.
# Run from the project root directory (env fb_r):
#   R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=$CONDA_PREFIX/lib/R/library Rscript scripts/08_random_forest.R
suppressPackageStartupMessages({ library(randomForest); library(ggplot2); library(dplyr); library(tidyr); library(patchwork) })
set.seed(2026)

OUT <- "results/08_random_forest"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
meta <- read.delim("metadata/metadata.tsv", comment.char = "#"); names(meta)[1] <- "sample"
tax  <- read.delim("results/04_taxonomy/export/taxonomy_ranks.tsv")
tab  <- as.matrix(read.delim("results/04_taxonomy/export/feature-table.tsv", skip = 1, check.names = FALSE, row.names = 1))
cult_order <- c("Cripps Pink", "Aztec Fuji", "Gala", "Hewe's Crab", "Granniwinkle", "Unknown cultivar")
cult_col <- c("Cripps Pink" = "#1B9E77", "Aztec Fuji" = "#D95F02", "Gala" = "#7570B3",
              "Hewe's Crab" = "#E7298A", "Granniwinkle" = "#66A61E", "Unknown cultivar" = "#E6AB02")
paper_rf <- c("Friedmanniella", "Aurantisolimonas", "Deinococcus", "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
              "Cutibacterium", "Gibbsiella", "Nocardioides", "Shimwellia")

md  <- meta[meta$health_status == "Canker", ]
g   <- sub("^g__", "", tax$Genus[match(rownames(tab), tax$asv)]); g[is.na(g) | g == ""] <- "Unclassified"
gt  <- rowsum(tab[, md$sample], g)
rel <- sweep(gt, 2, colSums(gt), "/")                         # genus relative abundance (0-1)
rel <- rel[rowSums(rel > 0) > 0 & rownames(rel) != "Unclassified", ]
X   <- as.data.frame(t(rel)); colnames(X) <- make.names(colnames(X), unique = TRUE)
key <- setNames(rownames(rel), colnames(X))
y   <- factor(md$cultivar, levels = cult_order)
cat("Samples:", nrow(X), " genera:", ncol(X), "\n"); print(table(y))

rf <- randomForest(x = X, y = y, ntree = 1000, importance = TRUE)
print(rf)
imp <- as.data.frame(importance(rf)); imp$genus <- key[rownames(imp)]
imp <- imp |> arrange(desc(MeanDecreaseGini))
# cultivar in which each genus is most abundant (colour, as in the paper)
mean_by_cult <- sapply(cult_order, function(cv) rowMeans(rel[, md$cultivar == cv, drop = FALSE]))
imp$top_cultivar <- cult_order[apply(mean_by_cult[imp$genus, , drop = FALSE], 1, which.max)]
imp$in_paper_top8 <- imp$genus %in% paper_rf
write.table(imp |> select(genus, top_cultivar, MeanDecreaseGini, MeanDecreaseAccuracy, all_of(cult_order), in_paper_top8),
            file.path(OUT, "rf_importance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(rf$confusion), file.path(OUT, "rf_confusion.tsv"), sep = "\t", quote = FALSE)

top8 <- head(imp, 8)
top8$genus_lab <- factor(paste0("g__", top8$genus), levels = rev(paste0("g__", top8$genus)))
ab <- lapply(top8$genus, function(gn) data.frame(genus = gn, cultivar = md$cultivar, rel = rel[gn, ])) |> bind_rows() |>
  group_by(genus, cultivar) |> summarise(mean = mean(rel), se = sd(rel) / sqrt(n()), .groups = "drop") |>
  mutate(genus_lab = factor(paste0("g__", genus), levels = levels(top8$genus_lab)),
         cultivar = factor(cultivar, levels = rev(cult_order)))

p1 <- ggplot(top8, aes(MeanDecreaseGini, genus_lab, fill = top_cultivar)) + geom_col(colour = "black") +
  scale_fill_manual(values = cult_col, name = "Most abundant in") + labs(x = "Mean Decrease Gini", y = NULL) + theme_bw()
p2 <- ggplot(ab, aes(mean, genus_lab, fill = cultivar)) +
  geom_col(position = position_dodge(0.8), width = 0.8) +
  geom_errorbar(aes(xmin = mean, xmax = mean + se), position = position_dodge(0.8), width = 0.3) +
  scale_fill_manual(values = cult_col, name = "Group") + labs(x = "Relative abundance", y = NULL) +
  theme_bw() + theme(axis.text.y = element_blank())
pdf("results/pdf/08_random_forest_cultivar.pdf", width = 11, height = 5.5)
print(p1 + p2 + plot_layout(widths = c(1, 1.2)) +
      plot_annotation(title = sprintf("Random Forest (1,000 trees): canker genera predicting cultivar - OOB error %.1f%%",
                                      100 * tail(rf$err.rate[, "OOB"], 1))))
top20 <- head(imp, 20); top20$genus <- factor(top20$genus, levels = rev(top20$genus))
print(ggplot(top20, aes(MeanDecreaseGini, genus, fill = top_cultivar)) + geom_col(colour = "black") +
  scale_fill_manual(values = cult_col, name = "Most abundant in") +
  labs(x = "Mean Decrease Gini", y = NULL, title = "Top 20 genera by Random Forest importance") + theme_bw())
dev.off()

cat("\nTop 15 by Mean Decrease Gini:\n")
print(head(imp |> transmute(genus, top_cultivar, Gini = round(MeanDecreaseGini, 3), in_paper_top8), 15), row.names = FALSE)
cat("\nPaper's 8 genera - rank here (of", nrow(imp), "), mean % in cankers:\n")
for (gn in paper_rf) cat(sprintf("  %-50s rank %s  mean %.3f%%\n", gn,
    ifelse(gn %in% imp$genus, which(imp$genus == gn), "absent"), ifelse(gn %in% rownames(rel), 100 * mean(rel[gn, ]), 0)))
