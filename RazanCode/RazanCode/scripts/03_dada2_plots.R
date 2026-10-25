# Step 3 plots: read tracking through DADA2 + ASVs per sample -> results/pdf/03_*.pdf
# Run from the project root directory (called by scripts/03_dada2_denoise.sh)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

d    <- "results/03_dada2/export"
meta <- read.delim("metadata/metadata.tsv", comment.char = "#")
meta <- meta[meta$sample.id != "#q2:types", ]
names(meta)[1] <- "sample"
names(meta)[names(meta) == "label"] <- "sample_name"
status_col <- c(Canker = "#B2182B", Healthy = "#2166AC")

st <- read.delim(file.path(d, "stats.tsv"), comment.char = "#")
st <- st[st[[1]] != "#q2:types", ]
names(st)[1] <- "sample"
num <- setdiff(names(st), "sample")
st[num] <- lapply(st[num], as.numeric)

tab <- read.delim(file.path(d, "feature-table.tsv"), skip = 1, check.names = FALSE, row.names = 1)
asv <- data.frame(sample = colnames(tab), n_asv = colSums(tab > 0), n_reads = colSums(tab))

st <- merge(merge(st, asv, by = "sample"), meta, by = "sample")
write.table(st, "results/03_dada2/read_tracking.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

steps <- intersect(c("input", "primer.removed", "filtered", "denoised", "non.chimeric"), names(st))
long <- st |> select(sample, sample_name, health_status, all_of(steps)) |>
  pivot_longer(all_of(steps), names_to = "step", values_to = "reads") |>
  mutate(step = factor(step, levels = steps))

pdf("results/pdf/03_dada2_read_tracking.pdf", width = 9, height = 5)
print(ggplot(long, aes(step, reads / 1000, group = sample, colour = health_status)) +
  geom_line(alpha = 0.6) + geom_point(size = 0.8) + scale_colour_manual(values = status_col) +
  labs(x = "DADA2 step", y = "Reads (x1000)", colour = "Health status",
       title = "Reads retained through DADA2 (CCS mode)") + theme_bw())
st$sample_name <- factor(st$sample_name, levels = st$sample_name[order(st$health_status, st$non.chimeric)])
print(ggplot(st, aes(sample_name, 100 * non.chimeric / input, fill = health_status)) +
  geom_col() + scale_fill_manual(values = status_col) +
  labs(x = NULL, y = "% of input reads kept (non-chimeric)", fill = "Health status") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)))
dev.off()

pdf("results/pdf/03_asvs_per_sample.pdf", width = 7, height = 5)
print(ggplot(st, aes(cultivar, n_asv, fill = health_status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1) +
  scale_fill_manual(values = status_col) +
  labs(x = NULL, y = "ASVs per sample (unrarefied)", fill = "Health status",
       title = "ASV richness per sample") + theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)))
dev.off()

cat("\nTotal non-chimeric reads:", sum(st$non.chimeric), " (paper: 1,817,968)\n")
cat("Total ASVs:", nrow(tab), " (paper: 2,322)\n")
for (g in c("Canker", "Healthy")) {
  x <- st$n_asv[st$health_status == g]
  cat(g, "ASVs/sample: range", min(x), "-", max(x), " mean", round(mean(x)), "\n")
}
cat("% input kept: min", round(min(100 * st$non.chimeric / st$input), 1),
    " median", round(median(100 * st$non.chimeric / st$input), 1), "\n")
