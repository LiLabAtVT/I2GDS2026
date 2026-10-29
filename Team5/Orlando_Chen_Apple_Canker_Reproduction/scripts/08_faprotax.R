#!/usr/bin/env Rscript
# Step 8: FAPROTAX functional prediction with microeco (paper: "FAPROTAX analysis was performed using the
# microeco package"; Fig. 6). Bacteria-only table (chloroplast/mitochondria removed), not rarefied.
# Output: % of reads assigned to each function per sample (abundance-weighted), heatmap grouped as Fig. 6,
# and tests: canker vs healthy (Wilcoxon) and among cultivars within canker (Kruskal-Wallis).
suppressPackageStartupMessages({ library(microeco); library(dplyr); library(tidyr); library(ggplot2) })
B <- "/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code"
# Optional args: [1] species assignment TSV from Step 8b (ASV hash -> species), [2] output tag (e.g. "08b_faprotax_species")
args <- commandArgs(trailingOnly = TRUE)
SPF <- if (length(args) >= 1) args[1] else NA; TAG <- if (length(args) >= 2) args[2] else "08_faprotax"
O <- file.path(B, "results", TAG); dir.create(O, showWarnings = FALSE, recursive = TRUE)
CULT <- c("Cripps Pink", "Aztec Fuji", "Gala", "Hewe's Crab", "Granniwinkle", "Unknown cultivar")

meta <- read.delim(file.path(B, "raw_data/sample_metadata.tsv"), check.names = FALSE)
meta <- meta[meta[["sample-id"]] != "#q2:types", ]; rownames(meta) <- meta[["sample-id"]]
ft <- as.matrix(read.delim(file.path(B, "results/04_dada2/export/feature-table.tsv"), skip = 1,
                           row.names = 1, check.names = FALSE, comment.char = ""))
tax <- read.delim(file.path(B, "results/05_taxonomy/taxonomy_by_rank.tsv"), row.names = 1, check.names = FALSE)[rownames(ft), ]
keep <- tax$Domain == "Bacteria" & !(tax$Order %in% "Chloroplast") & !(tax$Family %in% "Mitochondria")
ft <- ft[keep, ]; tax <- tax[keep, ]
ft <- ft[, colSums(ft) > 0]
if (!is.na(SPF)) {   # add BLAST species (only unambiguous calls, same genus as SILVA)
  sp <- read.delim(SPF, check.names = FALSE); sp <- sp[sp$assigned == "yes", ]
  i <- match(sp$feature_id, rownames(tax)); tax$Species[i[!is.na(i)]] <- sp$species[!is.na(i)]
  cat("Species added from BLAST:", sum(!is.na(i)), "ASVs\n"); print(table(tax$Species[!is.na(tax$Species)]))
}
meta <- meta[colnames(ft), c("cultivar", "health_status", "replicate")]
meta$cultivar <- factor(meta$cultivar, CULT); meta$health_status <- factor(meta$health_status, c("Canker", "Healthy"))

# microeco expects prefixed ranks (k__, p__, ... g__, s__)
pre <- c(Kingdom = "k__", Phylum = "p__", Class = "c__", Order = "o__", Family = "f__", Genus = "g__", Species = "s__")
tt <- data.frame(Kingdom = tax$Domain, Phylum = tax$Phylum, Class = tax$Class, Order = tax$Order,
                 Family = tax$Family, Genus = tax$Genus, Species = tax$Species, row.names = rownames(tax))
for (r in names(pre)) tt[[r]] <- paste0(pre[[r]], ifelse(is.na(tt[[r]]), "", tt[[r]]))
mt <- microtable$new(sample_table = meta, otu_table = as.data.frame(ft), tax_table = tt)
mt$tidy_dataset()
print(mt)

t1 <- trans_func$new(mt)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)
perc <- t1$res_spe_func_perc                                    # samples x functions, % of reads
write.table(perc, file.path(O, "faprotax_function_percent_per_sample.tsv"), sep = "\t", quote = FALSE)
asv_func <- t1$res_spe_func
cat("ASVs with >=1 FAPROTAX function:", sum(rowSums(asv_func) > 0), "of", nrow(asv_func),
    sprintf("(%.1f%% of bacterial reads)\n", 100 * sum(ft[rownames(asv_func)[rowSums(asv_func) > 0], ]) / sum(ft)))

# Function groups as in Fig. 6 (Energy source / C-cycle / N-cycle / Other)
grp <- t1$func_group_list
grp <- if (!is.null(grp$FAPROTAX)) grp$FAPROTAX else grp
fg <- bind_rows(lapply(names(grp), function(g) data.frame(func = grp[[g]], group = g))) %>% distinct(func, .keep_all = TRUE)
long <- as.data.frame(perc) %>% mutate(sample = rownames(perc)) %>% pivot_longer(-sample, names_to = "func", values_to = "pct") %>%
  left_join(fg, by = "func") %>% mutate(group = ifelse(is.na(group), "Other", group)) %>%
  left_join(meta %>% mutate(sample = rownames(meta)), by = "sample")
present <- long %>% group_by(func) %>% summarise(mx = max(pct)) %>% filter(mx > 0) %>% pull(func)
long <- filter(long, func %in% present)

# Tests
wil <- long %>% group_by(group, func) %>%
  summarise(mean_canker = mean(pct[health_status == "Canker"]), mean_healthy = mean(pct[health_status == "Healthy"]),
            p_wilcox = tryCatch(wilcox.test(pct ~ health_status, exact = FALSE)$p.value, error = function(e) NA), .groups = "drop") %>%
  mutate(p_BH = p.adjust(p_wilcox, "BH")) %>% arrange(p_wilcox)
kw <- long %>% filter(health_status == "Canker") %>% group_by(group, func) %>%
  summarise(p_kruskal_cultivar = tryCatch(kruskal.test(pct ~ cultivar)$p.value, error = function(e) NA), .groups = "drop") %>%
  mutate(p_BH = p.adjust(p_kruskal_cultivar, "BH")) %>% arrange(p_kruskal_cultivar)
bycult <- long %>% group_by(health_status, cultivar, func) %>% summarise(mean_pct = mean(pct), .groups = "drop") %>%
  pivot_wider(names_from = c(health_status, cultivar), values_from = mean_pct)
write.table(wil, file.path(O, "faprotax_canker_vs_healthy_wilcoxon.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(kw, file.path(O, "faprotax_canker_cultivar_kruskal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bycult, file.path(O, "faprotax_mean_by_status_cultivar.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 200); print(as.data.frame(wil), digits = 3); print(as.data.frame(head(kw, 15)), digits = 3)
focus <- c("plant_pathogen", "ureolysis", "fermentation", "methanol_oxidation", "methylotrophy", "hydrocarbon_degradation",
           "aerobic_chemoheterotrophy", "anaerobic_chemoheterotrophy")
print(as.data.frame(filter(bycult, func %in% focus)), digits = 3)
if ("plant_pathogen" %in% colnames(perc)) {
  pp <- long %>% filter(func == "plant_pathogen", health_status == "Canker")
  cat("plant_pathogen among canker cultivars, Kruskal P =", signif(kruskal.test(pct ~ cultivar, data = pp)$p.value, 3), "\n")
  print(as.data.frame(pp %>% group_by(cultivar) %>% summarise(mean = mean(pct), sd = sd(pct), n = n())), digits = 3)
}

# Heatmap as Fig. 6
so <- rownames(meta)[order(meta$health_status, meta$cultivar, as.integer(meta$replicate))]
long$sample <- factor(long$sample, so)
GORD <- c("Energy source", "C-cycle", "N-cycle", "S-cycle", "Other")
long$group <- factor(long$group, c(intersect(GORD, unique(long$group)), setdiff(unique(long$group), GORD)))
pdf(file.path(B, "pdf", paste0(TAG, ".pdf")), width = 16, height = 11)
print(ggplot(long, aes(sample, func, fill = pct)) + geom_tile(colour = "grey85") +
  scale_fill_gradient(low = "white", high = "#8b0a50", name = "Percentage (%)") +
  facet_grid(group ~ health_status + cultivar, scales = "free", space = "free") +
  labs(x = NULL, y = NULL, title = paste0("Fig. 6 reproduction: FAPROTAX functions (% of bacterial reads), bacteria-only table",
                                         if (is.na(SPF)) " (genus-level taxonomy)" else " (+ BLAST species, Step 8b)"),
       subtitle = paste("Healthy samples with no bacterial reads are excluded; n =", sum(meta$health_status == "Canker"), "canker,",
                        sum(meta$health_status == "Healthy"), "healthy")) +
  theme_bw(9) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), strip.text = element_text(size = 7),
                      panel.spacing = unit(0.1, "lines")))
top <- wil %>% filter(!is.na(p_wilcox)) %>% slice_head(n = 15) %>% pull(func)
print(ggplot(filter(long, func %in% top), aes(health_status, pct, fill = health_status)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8) + facet_wrap(~ func, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c(Canker = "#d6604d", Healthy = "#4daf4a")) +
  labs(x = NULL, y = "% of bacterial reads", title = "Functions differing most between canker and healthy (Wilcoxon; see TSV for BH-adjusted P)") +
  theme_bw(9) + theme(legend.position = "none"))
dev.off()
cat("Step 8 done\n")
