# Step 6a: read-quality summary (FastQC) and DADA2 read tracking -> PDF
source("scripts/R_common.R")

meta <- read.delim("metadata/sample_metadata.tsv", check.names = FALSE)
meta <- meta[!grepl("^#", meta[[1]]), ]

# ---- parse FastQC modules from each report zip ----
read_module <- function(zipf, module) {
  txt <- readLines(unz(zipf, file.path(sub("\\.zip$", "", basename(zipf)), "fastqc_data.txt")))
  s <- grep(paste0("^>>", module), txt); e <- grep("^>>END_MODULE", txt); e <- e[e > s][1]
  block <- txt[(s + 1):(e - 1)]
  block <- block[!grepl("^#", block) | grepl("^#Base|^#Length|^#Quality", block)]
  df <- read.delim(text = sub("^#", "", block), check.names = FALSE)
  df$sample <- sub("_fastqc\\.zip$", "", basename(zipf)); df
}
zips <- list.files("results/01_fastqc", "_fastqc\\.zip$", full.names = TRUE)
qual <- bind_rows(lapply(zips, read_module, "Per base sequence quality"))
lens <- bind_rows(lapply(zips, read_module, "Sequence Length Distribution"))
sq   <- bind_rows(lapply(zips, read_module, "Per sequence quality scores"))

qual <- qual %>% mutate(pos = as.numeric(sub("-.*", "", Base))) %>%
  left_join(meta, by = c(sample = "sample-id"))
lens <- lens %>% mutate(len = as.numeric(sub("-.*", "", Length))) %>%
  left_join(meta, by = c(sample = "sample-id"))
sq   <- sq %>% left_join(meta, by = c(sample = "sample-id"))

p1 <- ggplot(qual, aes(pos, Mean, group = sample, colour = health)) +
  geom_line(alpha = .6) + scale_colour_manual(values = HEALTH_COL) +
  labs(x = "Position in read (bp)", y = "Mean Phred quality", title = "FastQC: per-base quality")
p2 <- ggplot(lens, aes(len, Count, group = sample, colour = health)) +
  geom_line(alpha = .6) + scale_colour_manual(values = HEALTH_COL) +
  coord_cartesian(xlim = c(1000, 1800)) +
  labs(x = "Read length (bp)", y = "Reads", title = "FastQC: read-length distribution")
p3 <- ggplot(sq, aes(Quality, Count, group = sample, colour = health)) +
  geom_line(alpha = .6) + scale_colour_manual(values = HEALTH_COL) +
  labs(x = "Mean read quality", y = "Reads", title = "FastQC: per-read mean quality")

# ---- DADA2 read tracking ----
st <- read.delim(file.path(EXP, "dada2_stats/stats.tsv"), check.names = FALSE)
st <- st[!grepl("^#", st[[1]]), ]
names(st)[1] <- "sample"
num <- setdiff(names(st), "sample"); st[num] <- lapply(st[num], as.numeric)
write.csv(st, tab_out("dada2_read_tracking.csv"), row.names = FALSE)
cnt <- intersect(c("input", "primer-removed", "filtered", "denoised", "non-chimeric"), names(st))
stl <- st %>% select(sample, all_of(cnt)) %>% pivot_longer(-sample) %>%
  mutate(name = factor(name, cnt)) %>% left_join(meta, by = c(sample = "sample-id"))
p4 <- ggplot(stl, aes(name, value, group = sample, colour = health)) +
  geom_line(alpha = .6) + geom_point(size = .8) + scale_colour_manual(values = HEALTH_COL) +
  labs(x = NULL, y = "Reads", title = "DADA2 denoise-ccs read tracking") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf(pdf_out("S0_read_QC.pdf"), width = 11, height = 8)
print((p1 | p2) / (p3 | p4) + plot_layout(guides = "collect"))
dev.off()
cat(sprintf("Total non-chimeric reads: %s\n", format(sum(st[["non-chimeric"]]), big.mark = ",")))
