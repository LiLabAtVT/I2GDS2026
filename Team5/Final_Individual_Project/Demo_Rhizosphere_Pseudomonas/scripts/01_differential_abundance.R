# Part 1. Which bacteria does the soybean root enrich or deplete, and do three methods agree?
# Genus-level rhizosphere vs bulk soil at each site (each rhizosphere sample paired with the bulk soil of its own plot).
#   Method A  ANCOM-BC2  (bias-corrected, compositional; the method used in the study)
#   Method B  DESeq2     (negative-binomial count model; common alternative)
#   Method C  paired Wilcoxon test on CLR-transformed abundances (simple compositional check)
# Usage: Rscript scripts/01_differential_abundance.R <input_dataset dir> <output dir> <cores>
a <- commandArgs(TRUE); IN <- a[1]; OUT <- a[2]; NC <- as.integer(ifelse(length(a) > 2, a[3], 4))
suppressPackageStartupMessages({library(ANCOMBC); library(DESeq2); library(dplyr); library(tidyr); library(ggplot2)})
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) read.delim(file.path(IN, f), check.names = FALSE)

## ---- input ----
md <- rd("sample_metadata.tsv"); rownames(md) <- md$sample_id
cnt <- rd("genus_counts.tsv"); G <- as.matrix(cnt[, -1]); rownames(G) <- cnt$genus; G <- G[, md$sample_id]
tax <- rd("genus_taxonomy.tsv"); rownames(tax) <- tax$genus
cat(sprintf("Input: %d samples, %d genera\n", ncol(G), nrow(G)))

## ---- parameters (identical to the study) ----
PRV <- 0.10      # a genus must be present in >= 10% of the samples being tested
ALPHA <- 0.05    # Benjamini-Hochberg adjusted p (q) threshold

paired_site <- function(s) {          # one site, plots that have BOTH a bulk and a rhizosphere sample
  d <- md %>% filter(site == s) %>% group_by(site_plot) %>% filter(n() == 2) %>% ungroup() %>% as.data.frame()
  rownames(d) <- d$sample_id; d$compartment <- relevel(factor(d$compartment), "Bulk"); d$site_plot <- factor(d$site_plot); d }

res_all <- list()
for (s in c("Kentland", "Warsaw")) {
  d <- paired_site(s); x <- G[, d$sample_id]; x <- x[rowSums(x) > 0, ]
  cat(sprintf("\n== %s: %d paired plots (%d samples), %d genera with reads\n", s, nrow(d) / 2, nrow(d), nrow(x)))

  ## Method A: ANCOM-BC2 (same call as the study)
  set.seed(123)
  an <- ancombc2(data = x, taxa_are_rows = TRUE, meta_data = d, tax_level = NULL, fix_formula = "compartment + site_plot", rand_formula = NULL,
                 p_adj_method = "BH", pseudo_sens = TRUE, prv_cut = PRV, lib_cut = 1000, s0_perc = 0.05,
                 group = "compartment", struc_zero = TRUE, neg_lb = TRUE, alpha = ALPHA, n_cl = NC, verbose = FALSE,
                 global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                 lme_control = lme4::lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore", check.nlev.gtr.1 = "ignore"))
  ra <- an$res %>% transmute(genus = taxon, ancombc2_lfc = lfc_compartmentRhizosphere, ancombc2_q = q_compartmentRhizosphere,
                             ancombc2_sig = diff_compartmentRhizosphere & passed_ss_compartmentRhizosphere)
  tested <- ra$genus

  ## Method B: DESeq2 on the same genera (counts, plot as a blocking term, poscounts size factors for sparse data)
  dds <- DESeqDataSetFromMatrix(x[tested, ], colData = d, design = ~ site_plot + compartment)
  dds <- DESeq(estimateSizeFactors(dds, type = "poscounts"), quiet = TRUE)
  rb <- as.data.frame(results(dds, contrast = c("compartment", "Rhizosphere", "Bulk"), alpha = ALPHA)) %>%
    transmute(genus = rownames(.), deseq2_log2fc = log2FoldChange, deseq2_q = padj, deseq2_sig = !is.na(padj) & padj < ALPHA)

  ## Method C: paired Wilcoxon on CLR (CLR over all genera and all 94 samples, pseudo-count 0.5)
  L <- log(G + 0.5); C <- sweep(L, 2, colMeans(L), "-")
  b <- d$sample_id[d$compartment == "Bulk"][order(d$site_plot[d$compartment == "Bulk"])]
  r <- d$sample_id[d$compartment == "Rhizosphere"][order(d$site_plot[d$compartment == "Rhizosphere"])]
  rc <- bind_rows(lapply(tested, function(t) { dd <- C[t, r] - C[t, b]
    data.frame(genus = t, clr_diff = mean(dd), p = suppressWarnings(wilcox.test(dd, exact = FALSE)$p.value)) })) %>%
    mutate(clr_q = p.adjust(p, "BH"), clr_sig = clr_q < ALPHA) %>% select(-p)

  m <- ra %>% left_join(rb, by = "genus") %>% left_join(rc, by = "genus") %>%
    mutate(site = s, phylum = tax[genus, "Phylum"], direction = ifelse(ancombc2_lfc > 0, "enriched", "depleted"))
  write.table(m, file.path(OUT, paste0("1_DA_", s, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  res_all[[s]] <- m
}
R <- bind_rows(res_all)

## ---- summaries ----
cnt_tab <- R %>% group_by(site) %>% summarise(genera_tested = n(),
  ancombc2_enriched = sum(ancombc2_sig & ancombc2_lfc > 0, na.rm = TRUE), ancombc2_depleted = sum(ancombc2_sig & ancombc2_lfc < 0, na.rm = TRUE),
  deseq2_enriched = sum(deseq2_sig & deseq2_log2fc > 0, na.rm = TRUE), deseq2_depleted = sum(deseq2_sig & deseq2_log2fc < 0, na.rm = TRUE),
  clr_enriched = sum(clr_sig & clr_diff > 0, na.rm = TRUE), clr_depleted = sum(clr_sig & clr_diff < 0, na.rm = TRUE),
  all3_enriched = sum(ancombc2_sig & deseq2_sig & clr_sig & ancombc2_lfc > 0 & deseq2_log2fc > 0 & clr_diff > 0, na.rm = TRUE),
  r_lfc_ancombc2_vs_deseq2 = round(cor(ancombc2_lfc, deseq2_log2fc, use = "complete.obs"), 2),
  r_lfc_ancombc2_vs_clr = round(cor(ancombc2_lfc, clr_diff, use = "complete.obs"), 2), .groups = "drop")
write.table(cnt_tab, file.path(OUT, "1_method_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

W <- R %>% select(site, genus, phylum, ancombc2_lfc, ancombc2_sig) %>% pivot_wider(names_from = site, values_from = c(ancombc2_lfc, ancombc2_sig)) %>%
  filter(!is.na(ancombc2_lfc_Kentland), !is.na(ancombc2_lfc_Warsaw)) %>%
  mutate(category = case_when(
    ancombc2_sig_Kentland & ancombc2_sig_Warsaw & ancombc2_lfc_Kentland > 0 & ancombc2_lfc_Warsaw > 0 ~ "Enriched at both",
    ancombc2_sig_Kentland & ancombc2_sig_Warsaw & ancombc2_lfc_Kentland < 0 & ancombc2_lfc_Warsaw < 0 ~ "Depleted at both",
    ancombc2_sig_Kentland & ancombc2_sig_Warsaw ~ "Opposite directions",
    ancombc2_sig_Kentland | ancombc2_sig_Warsaw ~ "Significant at one site", TRUE ~ "Not significant"))
write.table(W, file.path(OUT, "1_rhizosphere_effect_both_sites.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
both_up <- W %>% filter(category == "Enriched at both") %>% arrange(desc(ancombc2_lfc_Kentland + ancombc2_lfc_Warsaw))

pse <- R %>% filter(genus == "Pseudomonas") %>% select(site, ancombc2_lfc, ancombc2_q, deseq2_log2fc, deseq2_q, clr_diff, clr_q)
rank_pse <- R %>% filter(ancombc2_sig, ancombc2_lfc > 0, !grepl("unclassified|Unassigned", genus)) %>% group_by(site) %>%
  arrange(desc(ancombc2_lfc), .by_group = TRUE) %>% summarise(pseudomonas_rank_among_named_enriched = which(genus == "Pseudomonas")[1], .groups = "drop")
txt <- c("PART 1 - DIFFERENTIAL ABUNDANCE, RHIZOSPHERE vs OWN-PLOT BULK SOIL (genus level)",
         "", "Significant genera per method (q < 0.05; ANCOM-BC2 also requires the pseudo-count sensitivity test):",
         capture.output(print(as.data.frame(cnt_tab), row.names = FALSE)),
         "", "Genera tested at both sites, by ANCOM-BC2 category:", capture.output(print(table(W$category))),
         "", "Enriched at BOTH sites (ANCOM-BC2):", paste(" ", both_up$genus),
         "", "Pseudomonas in each method (ANCOM-BC2 lfc = natural log; DESeq2 = log2; CLR = mean paired difference):",
         capture.output(print(as.data.frame(pse), digits = 3, row.names = FALSE)),
         capture.output(print(as.data.frame(rank_pse), row.names = FALSE)))
writeLines(txt, file.path(OUT, "1_summary.txt")); cat(txt, sep = "\n")

## ---- figures ----
lab <- W %>% filter(genus %in% c("Pseudomonas", "Candidatus_Udaeobacter"))
COL <- c("Enriched at both" = "#1F77B4", "Depleted at both" = "#17566B", "Opposite directions" = "#E5751F", "Significant at one site" = "grey70", "Not significant" = "grey88")
g1 <- ggplot(W, aes(ancombc2_lfc_Kentland, ancombc2_lfc_Warsaw, colour = category)) + geom_hline(yintercept = 0, colour = "grey50") +
  geom_vline(xintercept = 0, colour = "grey50") + geom_point(size = 1.8, alpha = .8) +
  geom_text(data = lab, aes(label = genus), colour = "black", fontface = "italic", hjust = 1.08, vjust = -.7, size = 3.6) +
  scale_colour_manual(values = COL, name = NULL) + guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Rhizosphere effect at Kentland (ANCOM-BC2 log fold change)", y = "Rhizosphere effect at Warsaw (log fold change)",
       title = "Rhizosphere effect on each genus at two sites (ANCOM-BC2)",
       subtitle = sprintf("%d genera enriched at both sites; %d depleted at both; %d in opposite directions",
                          sum(W$category == "Enriched at both"), sum(W$category == "Depleted at both"), sum(W$category == "Opposite directions"))) +
  theme_bw(12) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "1_fig_quadrant_kentland_vs_warsaw.png"), g1, width = 8, height = 7.5, dpi = 200)
g2 <- ggplot(R, aes(ancombc2_lfc, deseq2_log2fc, colour = ancombc2_sig & deseq2_sig)) + geom_abline(slope = 1 / log(2), linetype = 2, colour = "grey50") +
  geom_point(size = 1.3, alpha = .7) + facet_wrap(~ site) + scale_colour_manual(values = c(`TRUE` = "#E5751F", `FALSE` = "grey70"), name = "Significant in both") +
  labs(x = "ANCOM-BC2 log fold change (natural log)", y = "DESeq2 log2 fold change", title = "Effect sizes from two methods on the same genera",
       subtitle = "Dashed line = equal effect after converting natural log to log2") + theme_bw(12) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "1_fig_ancombc2_vs_deseq2.png"), g2, width = 10, height = 5.5, dpi = 200)
cat("\nPart 1 done. Outputs in", OUT, "\n")
