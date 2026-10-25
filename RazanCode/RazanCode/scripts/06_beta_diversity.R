# Step 6: beta diversity (paper Fig. 3) - Bray-Curtis PCoA + PERMANOVA
#   paper: Bray-Curtis PCoA of healthy vs canker; PERMANOVA: health status R2 = 17.2%, F = 9.35, P = 0.001;
#          cultivar alone P = 0.223; status x cultivar R2 = 12.4%, F = 1.43, P = 0.016; pairwise P = 0.001
#   Input: rarefied tables from Step 5 (paper = all ASVs @14,881; filtered = bacteria only @505, cankers only)
# Run from the project root directory (env fb_r):
#   R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=$CONDA_PREFIX/lib/R/library Rscript scripts/06_beta_diversity.R
suppressPackageStartupMessages({ library(vegan); library(ggplot2); library(dplyr) })
set.seed(2026)

OUT <- "results/06_beta"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
meta <- read.delim("metadata/metadata.tsv", comment.char = "#")
names(meta)[1] <- "sample"
meta$health_status <- factor(meta$health_status, levels = c("Canker", "Healthy"))
cult_col <- c("Aztec Fuji" = "#F8766D", "Cripps Pink" = "#B79F00", "Gala" = "#00BA38",
              "Granniwinkle" = "#00BFC4", "Hewe's Crab" = "#619CFF", "Unknown cultivar" = "#F564E3")

tidy_adonis <- function(a, set, model) {
  data.frame(set = set, model = model, term = rownames(a), Df = a$Df, SumOfSqs = round(a$SumOfSqs, 4),
             R2 = round(a$R2, 4), F = round(a$F, 3), P = a$`Pr(>F)`, row.names = NULL)
}

run_beta <- function(set, depth) {
  tab <- as.matrix(read.delim(sprintf("results/05_alpha/%s/export/rarefied_table.tsv", set),
                              skip = 1, check.names = FALSE, row.names = 1))
  x  <- t(tab)                                   # samples x ASVs
  md <- meta[match(rownames(x), meta$sample), ]
  bc <- vegdist(x, method = "bray")
  pc <- cmdscale(bc, k = 3, eig = TRUE)
  pv <- 100 * pc$eig / sum(pc$eig[pc$eig > 0])
  sc <- data.frame(sample = rownames(x), PC1 = pc$points[, 1], PC2 = pc$points[, 2]) |> left_join(md, by = "sample")
  write.table(sc, file.path(OUT, sprintf("pcoa_scores_%s.tsv", set)), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(as.matrix(bc), file.path(OUT, sprintf("bray_curtis_%s.tsv", set)), sep = "\t", quote = FALSE)

  res <- list(); disp <- list()
  two <- nlevels(droplevels(md$health_status)) == 2
  if (two) {
    res[[1]] <- tidy_adonis(adonis2(bc ~ health_status * cultivar, data = md, permutations = 999, by = "terms"),
                            set, "status * cultivar (sequential)")
    res[[2]] <- tidy_adonis(adonis2(bc ~ cultivar, data = md, permutations = 999), set, "cultivar alone")
    res[[3]] <- tidy_adonis(adonis2(bc ~ health_status, data = md, permutations = 999), set, "status alone")
    bd <- betadisper(bc, md$health_status)
    pt <- permutest(bd, permutations = 999)
    disp[[1]] <- data.frame(set = set, grouping = "health_status", F = round(pt$tab$F[1], 3), P = pt$tab$`Pr(>F)`[1],
                            mean_dist = paste(names(bd$group.distances), round(bd$group.distances, 3), collapse = "; "))
  }
  ck <- md$health_status == "Canker"
  bck <- as.dist(as.matrix(bc)[ck, ck])
  res[[length(res) + 1]] <- tidy_adonis(adonis2(bck ~ cultivar, data = md[ck, ], permutations = 999), set, "cultivar within cankers")
  # instrument check: Hewe's Crab (Sequel IIe) is the only cultivar on the second instrument
  res[[length(res) + 1]] <- tidy_adonis(adonis2(bck ~ instrument, data = md[ck, ], permutations = 999), set,
                                        "instrument within cankers (= Hewe's Crab vs rest)")
  bd2 <- betadisper(bck, md$cultivar[ck]); pt2 <- permutest(bd2, permutations = 999)
  disp[[length(disp) + 1]] <- data.frame(set = set, grouping = "cultivar (cankers)", F = round(pt2$tab$F[1], 3),
                                         P = pt2$tab$`Pr(>F)`[1],
                                         mean_dist = paste(names(bd2$group.distances), round(bd2$group.distances, 3), collapse = "; "))
  # pairwise cultivar PERMANOVA within cankers (BH-adjusted)
  cv <- sort(unique(md$cultivar[ck])); pw <- list()
  for (i in 1:(length(cv) - 1)) for (j in (i + 1):length(cv)) {
    k <- ck & md$cultivar %in% cv[c(i, j)]
    a <- adonis2(as.dist(as.matrix(bc)[k, k]) ~ cultivar, data = md[k, ], permutations = 999)
    pw[[length(pw) + 1]] <- data.frame(set = set, a = cv[i], b = cv[j], R2 = round(a$R2[1], 3), F = round(a$F[1], 3), P = a$`Pr(>F)`[1])
  }
  pw <- bind_rows(pw); pw$P_BH <- round(p.adjust(pw$P, "BH"), 4)

  p <- ggplot(sc, aes(PC1, PC2)) +
    geom_point(aes(colour = cultivar, shape = health_status), size = 2.6) +
    scale_colour_manual(values = cult_col) + scale_shape_manual(values = c(Canker = 16, Healthy = 17), drop = FALSE) +
    labs(x = sprintf("[%.1f%%]", pv[1]), y = sprintf("[%.1f%%]", pv[2]), colour = "Cultivar", shape = "Health status",
         title = sprintf("Bray-Curtis PCoA - %s table (rarefied to %s reads, n = %d)", set, format(depth, big.mark = ","), nrow(sc))) +
    theme_bw() + coord_equal()
  if (two) p <- p + stat_ellipse(aes(group = health_status), linetype = 2, level = 0.95)
  else p <- p + stat_ellipse(aes(group = cultivar, colour = cultivar), linetype = 3, level = 0.95)
  pdf(sprintf("results/pdf/06_bray_curtis_pcoa_%s.pdf", set), width = 8, height = 6.5)
  print(p)
  if (two) print(p + facet_wrap(~ health_status) + guides(shape = "none"))
  dev.off()

  list(permanova = bind_rows(res), dispersion = bind_rows(disp), pairwise = pw,
       axes = data.frame(set = set, PC1 = round(pv[1], 1), PC2 = round(pv[2], 1)))
}

r <- list(paper = run_beta("paper", 14881), filtered = run_beta("filtered", 505))
for (k in c("permanova", "dispersion", "pairwise", "axes")) {
  z <- bind_rows(lapply(r, `[[`, k))
  write.table(z, file.path(OUT, paste0(k, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n==== ", k, " ====\n"); print(z, row.names = FALSE)
}
