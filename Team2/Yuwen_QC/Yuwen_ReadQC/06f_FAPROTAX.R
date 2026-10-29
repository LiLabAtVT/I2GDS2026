# Step 6f: Fig. 6 - FAPROTAX functional prediction with microeco (as in the paper)
source("scripts/R_common.R")
library(microeco)
# unrarefied, all samples with bacterial reads (canker vs healthy comparison)
ps <- readRDS(file.path(RDS, "ps_filtered.rds"))
mt <- to_meco(ps)

tf <- trans_func$new(mt)
tf$cal_spe_func(prok_database = "FAPROTAX")
# % of ASVs in each sample assigned to each function (the paper reports numbers of ASVs per function)
tf$cal_spe_func_perc(abundance_weighted = FALSE)
fp <- tf$res_spe_func_perc                                        # samples x functions
fp <- fp[, colSums(fp) > 0, drop = FALSE]
write.csv(fp, tab_out("Fig6_FAPROTAX_percent_ASVs_per_sample.csv"))
write.csv(tf$res_spe_func, tab_out("Fig6_FAPROTAX_ASV_function_assignment.csv"))

# functional categories as laid out in the paper's Fig. 6 (microeco's FAPROTAX has no grouping)
CATS <- list(
  "Energy source" = c("photoheterotrophy", "anaerobic_chemoheterotrophy", "aerobic_chemoheterotrophy",
                      "chemoheterotrophy", "photoautotrophy", "chemoautotrophy"),
  "C-cycle" = c("xylanolysis", "oil_bioremediation", "methylotrophy", "methanotrophy",
                "hydrocarbon_degradation", "fermentation", "chitinolysis", "cellulolysis",
                "methanogenesis", "ligninolysis"),
  "N-cycle" = c("ureolysis", "nitrous_oxide_denitrification", "nitrogen_respiration", "nitrogen_fixation",
                "nitrite_respiration", "nitrite_denitrification", "nitrate_respiration", "nitrate_reduction",
                "nitrate_denitrification", "denitrification", "nitrification", "aerobic_ammonia_oxidation",
                "aerobic_nitrite_oxidation", "nitrite_ammonification", "nitrate_ammonification"))
cat_of <- function(f) {
  hit <- names(CATS)[sapply(CATS, function(v) f %in% v)]
  if (length(hit)) hit[1] else "Other"
}
md <- mt$sample_table
d <- data.frame(sample = rownames(fp), fp, check.names = FALSE) %>%
  pivot_longer(-sample, names_to = "func", values_to = "perc") %>%
  left_join(md %>% mutate(sample = rownames(md)) %>% select(sample, cultivar, health), by = "sample") %>%
  mutate(category = sapply(func, cat_of))

# canker vs healthy test per function (Wilcoxon, BH-adjusted)
st <- d %>% group_by(category, func) %>%
  summarise(mean_canker = mean(perc[health == "Canker"]), mean_healthy = mean(perc[health == "Healthy"]),
            p = tryCatch(wilcox.test(perc ~ health)$p.value, error = function(e) NA_real_), .groups = "drop") %>%
  mutate(p_adj = p.adjust(p, "BH")) %>% arrange(p)
write.csv(st, tab_out("Fig6_FAPROTAX_canker_vs_healthy_wilcoxon.csv"), row.names = FALSE)
print(head(st, 15))

# heatmap: mean % per cultivar x health
h <- d %>% group_by(category, func, cultivar, health) %>% summarise(perc = mean(perc), .groups = "drop") %>%
  mutate(sig = st$p_adj[match(func, st$func)] < 0.05,
         func = paste0(func, ifelse(!is.na(sig) & sig, " *", "")),
         func = factor(func, rev(sort(unique(func)))))
p6 <- ggplot(h, aes(cultivar, func, fill = perc)) + geom_tile(colour = "white") +
  facet_grid(category ~ health, scales = "free_y", space = "free_y") +
  scale_fill_gradientn(colours = c("#f7fbff", "#6baed6", "#08306b", "#cb181d"), name = "% of ASVs") +
  labs(x = NULL, y = NULL,
       title = "Fig. 6: FAPROTAX functional profiles of canker and healthy wood",
       subtitle = "Mean % of ASVs assigned to each function; * = canker vs healthy BH-adjusted P < 0.05") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.text.y = element_text(size = 6),
        strip.text.y = element_text(angle = 0))
pdf(pdf_out("Fig6_FAPROTAX_functions.pdf"), width = 10, height = 12)
print(p6)
dev.off()
