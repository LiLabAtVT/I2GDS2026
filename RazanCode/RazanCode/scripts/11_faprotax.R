# Step 11: FAPROTAX functional prediction (paper Fig. 6)
#   paper: FAPROTAX (Louca et al. 2016) via microeco in R; heatmap of functions per sample, canker vs healthy x cultivar,
#          grouped as Energy source / C-cycle / N-cycle / Other; chemoheterotrophy, fermentation, ureolysis, plant pathogen
#          in cankers; methanol oxidation, methylotrophy, methanotrophy, hydrocarbon degradation in healthy bark
#   Table (user decision): relaxed = bacteria only, 40 samples with >= 1 read, not rarefied.
#   Two scores per sample: % of ASVs with the function (unweighted, "number of ASVs" as in the paper) and
#   abundance-weighted % (share of reads).
# Run from the project root directory (env fb_r):
#   R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=$CONDA_PREFIX/lib/R/library Rscript scripts/11_faprotax.R
suppressPackageStartupMessages({ library(microeco); library(pheatmap); library(dplyr) })

OUT <- "results/11_faprotax"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
meta <- read.delim("metadata/metadata.tsv", comment.char = "#"); names(meta)[1] <- "sample"
tab  <- read.delim("results/04_taxonomy/export/feature-table.tsv", skip = 1, check.names = FALSE, row.names = 1)
tab  <- tab[, colSums(tab) > 0]
tax  <- read.delim("results/04_taxonomy/export/taxonomy_ranks.tsv")
tt   <- tax[match(rownames(tab), tax$asv), c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
names(tt)[1] <- "Kingdom"; rownames(tt) <- rownames(tab); tt[is.na(tt)] <- ""
# species from BLAST (scripts/11a_blast_all_asvs.sh): best hit >= 99% identity AND same genus as SILVA
# (SILVA gives no species; FAPROTAX assigns plant_pathogen etc. at species level, e.g. Erwinia amylovora)
bb <- read.delim(file.path(OUT, "blast_all_asvs_best.tsv"))
bb <- bb[match(rownames(tt), bb$feature_id), ]
bgen <- sub(" .*", "", bb$species)
ok <- !is.na(bb$pident) & bb$pident >= 99 & mapply(function(b, s) nzchar(b) && grepl(b, s, fixed = TRUE), bgen, tt$Genus)
tt$Species[ok] <- paste0("s__", bb$species[ok])
cat("Species added from BLAST:", sum(ok), "ASVs (", round(100 * sum(tab[ok, ]) / sum(tab), 1), "% of reads )\n")
write.table(data.frame(feature_id = rownames(tt), Genus = tt$Genus, Species = tt$Species), file.path(OUT, "species_labels_used.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
st <- meta[match(colnames(tab), meta$sample), ]; rownames(st) <- st$sample
mt <- microtable$new(otu_table = tab, tax_table = tt, sample_table = st)
mt$tidy_dataset(); mt$tax_table <- tidy_taxonomy(mt$tax_table)

tf <- trans_func$new(mt)
tf$cal_spe_func(prok_database = "FAPROTAX")
cat("ASVs with >= 1 FAPROTAX function:", sum(rowSums(tf$res_spe_func) > 0), "of", nrow(tf$res_spe_func), "\n")
groups <- tf$func_group_list$FAPROTAX

cult_order <- c("Cripps Pink", "Aztec Fuji", "Gala", "Hewe's Crab", "Granniwinkle", "Unknown cultivar")
cult_col <- c("Cripps Pink" = "#1B9E77", "Aztec Fuji" = "#D95F02", "Gala" = "#7570B3",
              "Hewe's Crab" = "#E7298A", "Granniwinkle" = "#66A61E", "Unknown cultivar" = "#E6AB02")
ord <- st[order(st$health_status, match(st$cultivar, cult_order), as.numeric(st$replicate)), ]

run <- function(weighted) {
  tag <- if (weighted) "abundance_weighted" else "pct_ASVs"
  tf$cal_spe_func_perc(abundance_weighted = weighted)
  P <- 100 * as.matrix(tf$res_spe_func_perc)                         # samples x functions (proportion -> %)
  P <- P[ord$sample, colSums(P) > 0, drop = FALSE]
  write.table(round(P, 4), file.path(OUT, sprintf("faprotax_%s.tsv", tag)), sep = "\t", quote = FALSE)
  # canker vs healthy per function
  w <- lapply(colnames(P), function(f) {
    x <- P[, f]; g <- ord$health_status
    data.frame(func = f, mean_Canker = mean(x[g == "Canker"]), mean_Healthy = mean(x[g == "Healthy"]),
               P = suppressWarnings(wilcox.test(x ~ g)$p.value))
  }) |> bind_rows()
  w$P_BH <- p.adjust(w$P, "BH")
  w$group <- sapply(w$func, function(f) { k <- names(groups)[sapply(groups, function(v) f %in% v)]; if (length(k)) k[1] else "Other" })
  w <- w |> arrange(group, P)
  write.table(w, file.path(OUT, sprintf("faprotax_canker_vs_healthy_%s.tsv", tag)), sep = "\t", quote = FALSE, row.names = FALSE)

  fo <- w |> mutate(group = factor(group, levels = unique(c("Energy source", "C-cycle", "N-cycle", "S-cycle", names(groups), "Other")))) |>
    arrange(group, func)
  M <- t(P[, fo$func, drop = FALSE]); colnames(M) <- ord$label
  rownames(M) <- ifelse(fo$P_BH < 0.05, paste0(fo$func, " *"), fo$func)
  gaps_r <- cumsum(table(droplevels(fo$group))); gaps_c <- cumsum(table(paste(ord$health_status)))
  pdf(sprintf("results/pdf/11_faprotax_heatmap_%s.pdf", tag), width = 13, height = max(7, 0.2 * nrow(M) + 3))
  pheatmap(M, cluster_rows = FALSE, cluster_cols = FALSE, gaps_row = gaps_r[-length(gaps_r)], gaps_col = gaps_c[-length(gaps_c)],
           annotation_col = data.frame(Status = ord$health_status, Cultivar = ord$cultivar, row.names = ord$label),
           annotation_row = data.frame(Group = as.character(fo$group), row.names = rownames(M)),
           annotation_colors = list(Cultivar = cult_col, Status = c(Canker = "#B2182B", Healthy = "#2166AC")),
           color = colorRampPalette(c("white", "#F4C2CB", "#B03060", "#7A0A2E"))(100),
           border_color = "grey80", fontsize_row = 7, fontsize_col = 6,
           main = sprintf("FAPROTAX functions (%s) - * = canker vs healthy Wilcoxon BH P < 0.05",
                          if (weighted) "% of reads" else "% of ASVs"))
  dev.off()
  w
}
w_u <- run(FALSE)
w_w <- run(TRUE)

paper_funcs <- c("aerobic_chemoheterotrophy", "anaerobic_chemoheterotrophy", "fermentation", "ureolysis", "plant_pathogen",
                 "methanol_oxidation", "methylotrophy", "methanotrophy", "hydrocarbon_degradation", "phototrophy", "nitrogen_fixation")
for (nm in c("pct_ASVs", "abundance_weighted")) {
  w <- if (nm == "pct_ASVs") w_u else w_w
  cat("\n====", nm, "====  functions detected:", nrow(w), " BH-significant:", sum(w$P_BH < 0.05), "\n")
  print(w |> filter(func %in% paper_funcs) |>
          transmute(func, group, Canker = round(mean_Canker, 2), Healthy = round(mean_Healthy, 2),
                    higher = ifelse(mean_Canker > mean_Healthy, "Canker", "Healthy"), P = signif(P, 2), P_BH = signif(P_BH, 2)),
        row.names = FALSE)
}
