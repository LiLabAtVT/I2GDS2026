# Step 9: LEfSe, canker vs healthy bark (paper Fig. 5A)
#   paper: LEfSe (Segata et al. 2011); Erwinia (LDA > 5) enriched in cankers; Bacillus (~5.7), Methylobacterium-Methylorubrum,
#          Sphingomonas, Hymenobacter, Pseudomonas, Comamonadaceae enriched in healthy bark; LDA >= ~3.5 shown
#   Table (user decision): relaxed = bacteria only, 40 samples with >= 1 read (32 canker, 8 healthy), not rarefied.
#   microeco trans_diff(method = "lefse"), genus level, Kruskal-Wallis alpha 0.05, LDA > 2.
# Run from the project root directory (env fb_r):
#   R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=$CONDA_PREFIX/lib/R/library Rscript scripts/09_lefse.R
suppressPackageStartupMessages({ library(microeco); library(ggplot2); library(dplyr) })
set.seed(2026)

OUT <- "results/09_lefse"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
meta <- read.delim("metadata/metadata.tsv", comment.char = "#"); names(meta)[1] <- "sample"
tab  <- read.delim("results/04_taxonomy/export/feature-table.tsv", skip = 1, check.names = FALSE, row.names = 1)
tab  <- tab[, colSums(tab) > 0]                                    # relaxed: every sample with >= 1 bacterial read
tax  <- read.delim("results/04_taxonomy/export/taxonomy_ranks.tsv")
tt   <- tax[match(rownames(tab), tax$asv), c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
names(tt)[1] <- "Kingdom"; rownames(tt) <- rownames(tab)
tt[is.na(tt)] <- ""
st <- meta[match(colnames(tab), meta$sample), ]; rownames(st) <- st$sample

mt <- microtable$new(otu_table = tab, tax_table = tt, sample_table = st)
mt$tidy_dataset()
mt$tax_table <- tidy_taxonomy(mt$tax_table)
mt$cal_abund()                                        # relative abundance per rank (also used by LEfSe)
cat("Samples:", ncol(mt$otu_table), "ASVs:", nrow(mt$otu_table), "\n"); print(table(mt$sample_table$health_status))

ld <- trans_diff$new(dataset = mt, method = "lefse", group = "health_status", taxa_level = "Genus",
                     alpha = 0.05, lefse_subgroup = NULL, p_adjust_method = "none")
res <- ld$res_diff
# microeco assigns "Group" by the higher MEDIAN; the original LEfSe (Segata et al.) assigns the class by the higher
# MEAN. With 8 sparse healthy samples these differ, so both are kept; plots use the original (mean) convention.
ga <- mt$taxa_abund$Genus
ga <- ga[match(res$Taxa, rownames(ga)), , drop = FALSE]
res$Group_microeco_median <- res$Group
res$mean_Canker  <- rowMeans(ga[, mt$sample_table$health_status == "Canker", drop = FALSE])
res$mean_Healthy <- rowMeans(ga[, mt$sample_table$health_status == "Healthy", drop = FALSE])
res$Group <- ifelse(res$mean_Canker >= res$mean_Healthy, "Canker", "Healthy")
ld$res_diff <- res
write.table(res, file.path(OUT, "lefse_genus.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf("results/pdf/09_lefse_canker_vs_healthy.pdf", width = 9, height = max(5, 0.22 * sum(res$LDA > 2) + 2))
print(ld$plot_diff_bar(threshold = 2, color_values = c(Canker = "#E07B79", Healthy = "#4FA37F")) +
        labs(title = "LEfSe: genera enriched in canker vs healthy bark (LDA > 2, KW P < 0.05)"))
top <- res |> filter(LDA > 2) |> arrange(desc(LDA)) |> head(30)
top$Taxa <- factor(sub(".*\\|", "", top$Taxa), levels = rev(unique(sub(".*\\|", "", top$Taxa))))
print(ggplot(top, aes(LDA, Taxa, colour = Group, size = LDA)) + geom_point() +
        scale_colour_manual(values = c(Canker = "#E07B79", Healthy = "#4FA37F"), name = "Enriched in") +
        scale_size_continuous(name = "LDA magnitude") + geom_vline(xintercept = c(2, 4), linetype = 2, colour = "grey50") +
        labs(x = "LDA score (log10)", y = "Bacterial genera", title = "LEfSe (paper Fig. 5A style)") + theme_bw())
dev.off()

cat("\nSignificant genera (LDA > 2):", sum(res$LDA > 2), " canker:", sum(res$LDA > 2 & res$Group == "Canker"),
    " healthy:", sum(res$LDA > 2 & res$Group == "Healthy"), "\n")
print(res |> filter(LDA > 2) |> arrange(Group, desc(LDA)) |>
        transmute(Group, Genus = sub(".*\\|", "", Taxa), LDA = round(LDA, 2), P = signif(P.unadj, 3)), row.names = FALSE)
paper <- c("Bacillus", "Erwinia", "Methylobacterium-Methylorubrum", "Sphingomonas", "Hymenobacter", "Pseudomonas",
           "Pantoea", "Massilia", "Kosakonia", "Novosphingobium", "Curtobacterium")
cat("\nPaper genera in our result:\n")
for (g in paper) { r <- res[grepl(paste0("g__", g, "$"), res$Taxa), ]
  cat(sprintf("  %-32s %s\n", g, if (nrow(r)) sprintf("%s  LDA %.2f  P %.3g", r$Group, r$LDA, r$P.unadj) else "not significant")) }
