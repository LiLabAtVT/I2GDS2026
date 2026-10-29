# Step 6c: Fig. 2A (class composition), Fig. 2B (alpha diversity + ANOVA),
#          Fig. 3 (Bray-Curtis PCoA + PERMANOVA)
# Split approach (healthy wood has only 7-177 bacterial reads):
#   canker vs healthy -> unrarefied counts / relative abundance of all samples with reads
#   cultivar effects within cankers -> cankers rarefied to 505 reads
source("scripts/R_common.R")
library(vegan); library(agricolae)
ps     <- readRDS(file.path(RDS, "ps_filtered.rds"))          # 40 samples, unrarefied
ps_can <- readRDS(file.path(RDS, "ps_canker_rarefied.rds"))   # 32 cankers, rarefied

# ================= Fig. 2A: class-level relative abundance =================
cl <- tax_glom(ps, "Class", NArm = FALSE)
taxa_names(cl) <- make.unique(tax_label(tax_table(cl), "Class"))
cl_s <- transform_sample_counts(cl, function(x) 100 * x / sum(x))      # per-sample %
ds <- psmelt(cl_s)
d <- ds %>% group_by(health, OTU) %>% summarise(Abundance = mean(Abundance), .groups = "drop")
top <- d %>% group_by(OTU) %>% summarise(m = max(Abundance)) %>% arrange(-m) %>% slice(1:21) %>% pull(OTU)
lev <- c(rev(top), "Other")
pal <- c(setNames(colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(top)), top), Other = "grey70")
d$Class2  <- factor(ifelse(d$OTU %in% top, d$OTU, "Other"), lev)
ds$Class2 <- factor(ifelse(ds$OTU %in% top, ds$OTU, "Other"), lev)
p2a <- ggplot(d, aes(health, Abundance, fill = Class2)) + geom_col(width = .7) +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "Mean relative abundance (%) (class level)", fill = "Class",
       title = "Fig. 2A: class-level composition, healthy vs canker",
       subtitle = sprintf("mean of per-sample %%; healthy n = %d, canker n = %d (samples with bacterial reads)",
                          sum(sample_data(ps)$health == "Healthy"), sum(sample_data(ps)$health == "Canker")))
p2a_s <- ggplot(ds, aes(Sample, Abundance, fill = Class2)) + geom_col(width = .9) +
  facet_grid(~ health, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = "Relative abundance (%)", fill = "Class", title = "Class composition per sample") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 6))
write.csv(d %>% select(health, Class = OTU, Abundance), tab_out("Fig2A_class_relab.csv"), row.names = FALSE)
pdf(pdf_out("Fig2A_class_composition.pdf"), width = 12, height = 7)
print(p2a); print(p2a_s)
dev.off()

# ================= Fig. 2B: alpha diversity =================
alpha_tab <- function(p) {
  al <- estimate_richness(p, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson"))
  rownames(al) <- sample_names(p)
  if (requireNamespace("picante", quietly = TRUE)) {
    al$PD <- picante::pd(t(as(otu_table(p), "matrix")), phy_tree(p), include.root = FALSE)$PD
  }
  cbind(as(sample_data(p), "data.frame"), reads = sample_sums(p), al)
}
al_all <- alpha_tab(ps)        # unrarefied, canker vs healthy
al_can <- alpha_tab(ps_can)    # rarefied cankers, cultivar effects
metrics <- intersect(c("Simpson", "Shannon", "Observed", "Chao1", "ACE", "PD"), names(al_all))
write.csv(al_all, tab_out("Fig2B_alpha_all_samples_unrarefied.csv"), row.names = FALSE)
write.csv(al_can, tab_out("Fig2B_alpha_cankers_rarefied505.csv"), row.names = FALSE)

# (i) canker vs healthy (unrarefied): Wilcoxon + two-way ANOVA health x cultivar (as in paper)
wil <- bind_rows(lapply(metrics, function(m) {
  w <- wilcox.test(al_all[[m]] ~ al_all$health)
  data.frame(metric = m, median_canker = median(al_all[al_all$health == "Canker", m]),
             median_healthy = median(al_all[al_all$health == "Healthy", m]), W = w$statistic, p = w$p.value)
}))
write.csv(wil, tab_out("Fig2B_alpha_canker_vs_healthy_wilcoxon.csv"), row.names = FALSE)
print(wil)
aov2 <- bind_rows(lapply(metrics, function(m) {
  a <- summary(aov(as.formula(paste(m, "~ health * cultivar")), data = al_all))[[1]]
  data.frame(metric = m, term = trimws(rownames(a)), a, check.names = FALSE)
}))
write.csv(aov2, tab_out("Fig2B_alpha_ANOVA_health_x_cultivar.csv"), row.names = FALSE)

# (ii) cultivar effect within cankers (rarefied) and within healthy (unrarefied): ANOVA + Tukey HSD letters
cult_test <- function(al, label) {
  bind_rows(lapply(metrics, function(m) {
    al2 <- al %>% group_by(cultivar) %>% filter(n() >= 2) %>% ungroup() %>% droplevels()
    fit <- aov(as.formula(paste(m, "~ cultivar")), data = al2)
    a <- summary(fit)[[1]]
    hsd <- HSD.test(fit, "cultivar", group = TRUE)$groups
    data.frame(set = label, metric = m, cultivar = rownames(hsd), letter = hsd$groups,
               F = a[1, "F value"], p_anova = a[1, "Pr(>F)"])
  }))
}
let <- bind_rows(cult_test(al_can, "Canker"),
                 cult_test(al_all[al_all$health == "Healthy", ], "Healthy"))
write.csv(let, tab_out("Fig2B_alpha_cultivar_ANOVA_TukeyHSD.csv"), row.names = FALSE)
print(let %>% distinct(set, metric, F, p_anova))

al_l <- bind_rows(al_can %>% mutate(set = "Canker"), al_all %>% filter(health == "Healthy") %>% mutate(set = "Healthy")) %>%
  select(set, cultivar, all_of(metrics)) %>% pivot_longer(all_of(metrics), names_to = "metric") %>%
  mutate(metric = factor(metric, metrics))
lab_y <- al_l %>% group_by(set, metric, cultivar) %>% summarise(y = max(value), .groups = "drop") %>%
  group_by(set, metric) %>% mutate(y = y + 0.1 * max(y)) %>%
  left_join(let %>% mutate(metric = factor(metric, metrics)), by = c("set", "metric", "cultivar"))
p2b <- ggplot(al_l, aes(cultivar, value, fill = cultivar)) +
  geom_boxplot(outlier.shape = NA, alpha = .8) + geom_jitter(width = .15, size = .8) +
  geom_text(data = lab_y, aes(y = y, label = letter), size = 3, na.rm = TRUE) +
  facet_grid(metric ~ set, scales = "free_y") + scale_fill_manual(values = CULT_COL) +
  labs(x = NULL, y = NULL, title = "Fig. 2B: alpha diversity by cultivar",
       subtitle = "Canker: rarefied to 505 reads; Healthy: unrarefied (7-177 reads). Letters = Tukey HSD (P < 0.05)") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
p2b_h <- al_all %>% select(health, all_of(metrics)) %>% pivot_longer(-health, names_to = "metric") %>%
  mutate(metric = factor(metric, metrics)) %>%
  ggplot(aes(health, value, fill = health)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = .8) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) + scale_fill_manual(values = HEALTH_COL) +
  labs(x = NULL, y = NULL, title = "Alpha diversity: canker vs healthy (unrarefied)",
       subtitle = paste(sprintf("%s P = %.2g", wil$metric, wil$p), collapse = ";  ")) +
  theme(legend.position = "none")
pdf(pdf_out("Fig2B_alpha_diversity.pdf"), width = 10, height = 12)
print(p2b)
dev.off()
pdf(pdf_out("Fig2B_alpha_canker_vs_healthy.pdf"), width = 12, height = 4)
print(p2b_h)
dev.off()

# ================= Fig. 3: Bray-Curtis PCoA + PERMANOVA =================
set.seed(SEED)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
bc <- phyloseq::distance(ps_rel, method = "bray")
md <- as(sample_data(ps_rel), "data.frame")
perm <- adonis2(bc ~ health * cultivar, data = md, permutations = 999, by = "terms")
print(perm)
write.csv(data.frame(term = rownames(perm), as.data.frame(perm), check.names = FALSE),
          tab_out("Fig3_PERMANOVA_health_x_cultivar.csv"), row.names = FALSE)
write.csv(as.data.frame(anova(betadisper(bc, md$health))), tab_out("Fig3_betadisper_health.csv"))
# cultivar effect among cankers only (rarefied)
bc_c <- phyloseq::distance(ps_can, method = "bray")
perm_c <- adonis2(bc_c ~ cultivar, data = as(sample_data(ps_can), "data.frame"), permutations = 999)
print(perm_c)
write.csv(data.frame(term = rownames(perm_c), as.data.frame(perm_c), check.names = FALSE),
          tab_out("Fig3_PERMANOVA_cultivar_cankers_only.csv"), row.names = FALSE)

pcoa_plot <- function(p, dist, title, sub) {
  ord <- ordinate(p, "PCoA", distance = dist)
  ev <- 100 * ord$values$Relative_eig[1:2]
  plot_ordination(p, ord, color = "cultivar", shape = "health") + geom_point(size = 3) +
    scale_colour_manual(values = CULT_COL) +
    labs(x = sprintf("PCoA1 (%.1f%%)", ev[1]), y = sprintf("PCoA2 (%.1f%%)", ev[2]), title = title, subtitle = sub)
}
sub1 <- sprintf("PERMANOVA  health: R2 = %.3f, F = %.2f, P = %.3f;  cultivar: P = %.3f;  health x cultivar: R2 = %.3f, P = %.3f",
                perm["health", "R2"], perm["health", "F"], perm["health", "Pr(>F)"], perm["cultivar", "Pr(>F)"],
                perm["health:cultivar", "R2"], perm["health:cultivar", "Pr(>F)"])
sub2 <- sprintf("PERMANOVA cultivar (cankers, rarefied 505): R2 = %.3f, F = %.2f, P = %.3f",
                perm_c[1, "R2"], perm_c[1, "F"], perm_c[1, "Pr(>F)"])
pdf(pdf_out("Fig3_PCoA_BrayCurtis.pdf"), width = 9, height = 6.5)
print(pcoa_plot(ps_rel, bc, "Fig. 3: PCoA of Bray-Curtis dissimilarity (relative abundance, all samples)", sub1) +
        stat_ellipse(aes(group = health, linetype = health), colour = "grey40"))
print(pcoa_plot(ps_can, bc_c, "PCoA of Bray-Curtis dissimilarity, cankers only (rarefied)", sub2))
dev.off()
