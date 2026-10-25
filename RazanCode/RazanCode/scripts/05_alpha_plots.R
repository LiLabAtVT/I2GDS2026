# Step 5 (R): ANOVA + Tukey letters on QIIME2 alpha metrics, Fig. 2B-style boxplots,
#             rarefaction curves, and paper-style class barplot (Fig. 2A, all ASVs)
# Run from the project root directory (called by scripts/05_alpha_diversity.sh, env fb_r)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(agricolae) })

meta <- read.delim("metadata/metadata.tsv", comment.char = "#")
names(meta)[1] <- "sample"
meta$health_status <- factor(meta$health_status, levels = c("Canker", "Healthy"))
metrics <- c(simpson = "Simpson", shannon = "Shannon (log2)", observed_features = "Observed ASVs",
             chao1 = "Chao1", faith_pd = "Faith PD", evenness = "Pielou evenness")
cult_col <- c("Aztec Fuji" = "#1B9E77", "Cripps Pink" = "#D95F02", "Gala" = "#7570B3",
              "Granniwinkle" = "#E7298A", "Hewe's Crab" = "#66A61E", "Unknown cultivar" = "#E6AB02")

analyse <- function(set, depth) {
  D <- file.path("results/05_alpha", set)
  al <- lapply(names(metrics), function(m) {
    x <- read.delim(file.path(D, "export", m, "alpha-diversity.tsv"))
    names(x) <- c("sample", "value"); x$metric <- m; x
  }) |> bind_rows() |> left_join(meta, by = "sample")
  al$group <- paste(al$cultivar, al$health_status, sep = " | ")
  write.table(al |> select(sample, label, cultivar, health_status, metric, value) |>
                pivot_wider(names_from = metric, values_from = value),
              file.path(D, "alpha_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  stats <- list(); letters_df <- list()
  two_status <- length(unique(al$health_status)) == 2
  for (m in names(metrics)) {
    x <- al[al$metric == m, ]
    fit1 <- aov(value ~ group, data = x)                     # 12 (or 6) groups, one-way (Fig. 2B letters)
    hsd  <- HSD.test(fit1, "group", group = TRUE)$groups
    letters_df[[m]] <- data.frame(group = rownames(hsd), letter = hsd$groups, metric = m)
    a1 <- summary(fit1)[[1]]
    row <- data.frame(set = set, metric = m, test = "one-way ANOVA: group",
                      F = a1$`F value`[1], P = a1$`Pr(>F)`[1])
    if (two_status) {
      a2 <- summary(aov(value ~ health_status * cultivar, data = x))[[1]]
      rn <- trimws(rownames(a2))
      row <- rbind(row, data.frame(set = set, metric = m, test = paste("two-way ANOVA:", rn[1:3]),
                                   F = a2$`F value`[1:3], P = a2$`Pr(>F)`[1:3]))
      w <- wilcox.test(value ~ health_status, data = x)
      row <- rbind(row, data.frame(set = set, metric = m, test = "Wilcoxon: Canker vs Healthy",
                                   F = NA, P = w$p.value))
    } else {
      a2 <- summary(aov(value ~ cultivar, data = x))[[1]]
      row <- rbind(row, data.frame(set = set, metric = m, test = "one-way ANOVA: cultivar (cankers)",
                                   F = a2$`F value`[1], P = a2$`Pr(>F)`[1]))
    }
    stats[[m]] <- row
  }
  stats <- bind_rows(stats); stats$F <- round(stats$F, 3); stats$P <- signif(stats$P, 3)
  write.table(stats, file.path(D, "alpha_anova.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  letters_df <- bind_rows(letters_df)
  write.table(letters_df, file.path(D, "alpha_tukey_letters.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # Fig. 2B style: cultivars ordered by mean within each health status, Tukey letters on top
  pdf(sprintf("results/pdf/05_alpha_diversity_%s.pdf", set), width = 9, height = 5.5)
  for (m in names(metrics)) {
    x <- al[al$metric == m, ]
    ord <- x |> group_by(health_status, cultivar, group) |> summarise(mu = mean(value), .groups = "drop") |>
      arrange(health_status, desc(mu))
    x$group <- factor(x$group, levels = ord$group)
    lab <- letters_df[letters_df$metric == m, ] |>
      left_join(x |> group_by(group, health_status) |> summarise(top = max(value), .groups = "drop"), by = "group")
    lab$group <- factor(lab$group, levels = ord$group)
    ymax <- max(x$value)
    p <- ggplot(x, aes(group, value, colour = cultivar)) +
      geom_boxplot(outlier.shape = NA, fill = NA) +
      geom_jitter(width = 0.15, size = 1.3) +
      geom_text(data = lab, aes(group, top + 0.06 * ymax, label = letter), colour = "black", size = 3.5) +
      scale_colour_manual(values = cult_col) +
      scale_x_discrete(labels = function(z) sub(" \\|.*", "", z)) +
      facet_grid(~ health_status, scales = "free_x", space = "free_x") +
      labs(x = NULL, y = metrics[[m]], colour = "Cultivar",
           title = sprintf("%s - %s table, rarefied to %s reads", metrics[[m]], set, format(depth, big.mark = ",")),
           subtitle = "Letters: Tukey HSD after one-way ANOVA across cultivar x health-status groups") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    print(p)
  }
  dev.off()

  # rarefaction curves
  rc <- lapply(c("observed_features", "shannon"), function(m) {
    f <- file.path(D, "export/alpha_rarefaction", paste0(m, ".csv"))
    if (!file.exists(f)) return(NULL)
    r <- read.csv(f, check.names = FALSE)
    dc <- grep("^depth-", names(r), value = TRUE)
    r[, c("sample-id", dc)] |> pivot_longer(all_of(dc), names_to = "k", values_to = "value") |>
      mutate(depth = as.numeric(sub("depth-([0-9]+)_iter.*", "\\1", k)), metric = m) |>
      group_by(`sample-id`, metric, depth) |> summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  }) |> bind_rows()
  if (nrow(rc)) {
    names(rc)[1] <- "sample"
    rc <- left_join(rc, meta, by = "sample")
    pdf(sprintf("results/pdf/05_rarefaction_curves_%s.pdf", set), width = 9, height = 4.5)
    print(ggplot(rc, aes(depth, value, group = sample, colour = health_status)) +
      geom_line(alpha = 0.7) + facet_wrap(~ metric, scales = "free_y") +
      scale_colour_manual(values = c(Canker = "#B2182B", Healthy = "#2166AC")) +
      geom_vline(xintercept = depth, linetype = 2) +
      labs(x = "Sequencing depth", y = NULL, colour = "Health status",
           title = sprintf("Alpha rarefaction (%s table; dashed = rarefaction depth)", set)) + theme_bw())
    dev.off()
  }

  cat("\n====", set, "table, depth", depth, "— samples kept:", length(unique(al$sample)), "====\n")
  print(al |> group_by(metric, health_status) |>
          summarise(mean = round(mean(value), 2), min = round(min(value), 2), max = round(max(value), 2), .groups = "drop"))
  print(stats, row.names = FALSE)
}

analyse("paper", 14881)
analyse("filtered", 505)

# ---- Fig. 2A paper-style: class-level composition on ALL ASVs (incl. chloroplast) ----
tab <- as.matrix(read.delim("results/03_dada2/export/feature-table.tsv", skip = 1, check.names = FALSE, row.names = 1))
tax <- read.delim("results/04_taxonomy/export/taxonomy_ranks.tsv")
cls <- tax$Class[match(rownames(tab), tax$asv)]
cls[is.na(cls) | cls %in% c("c__", "")] <- "Unclassified"
rel <- sweep(rowsum(tab, cls), 2, colSums(tab), "/")
top <- names(sort(rowMeans(rel), decreasing = TRUE))[1:min(21, nrow(rel))]
long <- as.data.frame(rel) |> mutate(Class = rownames(rel)) |>
  pivot_longer(-Class, names_to = "sample", values_to = "rel") |>
  mutate(Class = ifelse(Class %in% top, Class, "Other")) |>
  group_by(sample, Class) |> summarise(rel = sum(rel), .groups = "drop") |> left_join(meta, by = "sample")
long$Class <- factor(long$Class, levels = c("Other", rev(top)))
ord <- meta |> arrange(desc(health_status == "Healthy"), cultivar, as.numeric(replicate))
long$label <- factor(long$label, levels = ord$label)
long$health_status <- factor(long$health_status, levels = c("Healthy", "Canker"))
pal <- c("#1F4E79", "#1B7837", "#E07B39", "#8E44AD", "#2E86C1", "#C0392B", "#F1C40F", "#16A085",
         "#7F8C8D", "#D35400", "#27AE60", "#5B2C6F", "#A93226", "#117864", "#B7950B", "#2874A6",
         "#935116", "#1D8348", "#6C3483", "#B03A2E", "#0E6655", "grey80")
names(pal) <- c(top, "Other")[seq_along(pal)]
pdf("results/pdf/05_class_relative_abundance_all_ASVs.pdf", width = 13, height = 6.5)
print(ggplot(long, aes(label, 100 * rel, fill = Class)) + geom_col(width = 0.9) +
  scale_fill_manual(values = pal) +
  facet_grid(~ health_status + cultivar, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Relative abundance (%) (Class level)", fill = "Class",
       title = "Class-level composition, all ASVs incl. chloroplast/mitochondria (paper-style Fig. 2A)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                     strip.text = element_text(size = 7), legend.text = element_text(size = 7)))
dev.off()
