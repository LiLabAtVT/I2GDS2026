# Step 2d: plots from seqkit output -> results/pdf/02_*.pdf
# Run from the project root directory (called by scripts/02_read_qc.sh)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

qc   <- "results/02_read_qc"
meta <- read.delim("metadata/metadata.tsv", comment.char = "#") |>
  dplyr::filter(sample.id != "#q2:types")
names(meta)[1] <- "sample"
names(meta)[names(meta) == "label"] <- "sample_name"
status_col <- c(Canker = "#B2182B", Healthy = "#2166AC")

# ---- per-sample summary ----
st <- read.delim(file.path(qc, "seqkit_stats.tsv"))
st$sample <- sub("\\.fastq\\.gz$", "", basename(st$file))
st <- merge(st, meta, by = "sample")
st$sample_name <- factor(st$sample_name, levels = st$sample_name[order(st$health_status, st$cultivar, st$num_seqs)])
write.table(st[, c("sample", "sample_name", "cultivar", "health_status", "instrument",
                   "num_seqs", "min_len", "avg_len", "max_len", "Q1", "Q2", "Q3",
                   "AvgQual", "Q20...", "Q30...")],
            file.path(qc, "per_sample_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf("results/pdf/02_reads_per_sample.pdf", width = 11, height = 5)
print(ggplot(st, aes(sample_name, num_seqs / 1000, fill = health_status)) +
  geom_col() + scale_fill_manual(values = status_col) +
  labs(x = NULL, y = "Raw HiFi reads (x1000)", fill = "Health status",
       title = "Raw reads per sample (PRJNA1348523)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)))
dev.off()

# ---- per-read length / quality ----
pr <- read.delim(gzfile(file.path(qc, "per_read_len_qual.tsv.gz")), header = FALSE,
                 col.names = c("sample", "length", "avgQ"))
pr <- merge(pr, meta[, c("sample", "sample_name", "health_status", "instrument")], by = "sample")

pdf("results/pdf/02_read_length_distribution.pdf", width = 8, height = 5)
print(ggplot(pr, aes(length, fill = health_status)) +
  geom_histogram(binwidth = 10, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = c(1000, 1600), linetype = 2) +
  coord_cartesian(xlim = c(0, 2500)) + scale_fill_manual(values = status_col) +
  labs(x = "Read length (bp)", y = "Reads", fill = "Health status",
       title = "Raw read length (dashed: 1,000-1,600 bp)") + theme_bw())
print(ggplot(pr, aes(sample_name, length, fill = health_status)) +
  geom_boxplot(outlier.size = 0.2) + scale_fill_manual(values = status_col) +
  coord_cartesian(ylim = c(1000, 1800)) +
  labs(x = NULL, y = "Read length (bp)", fill = "Health status") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)))
dev.off()

pdf("results/pdf/02_read_quality_distribution.pdf", width = 8, height = 5)
print(ggplot(pr, aes(avgQ, fill = instrument)) +
  geom_histogram(binwidth = 1, alpha = 0.7, position = "identity") +
  labs(x = "Mean read Phred quality", y = "Reads", fill = "Instrument",
       title = "Mean per-read quality") + theme_bw())
print(ggplot(pr, aes(sample_name, avgQ, fill = health_status)) +
  geom_boxplot(outlier.size = 0.2) + scale_fill_manual(values = status_col) +
  labs(x = NULL, y = "Mean read Phred quality", fill = "Health status") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)))
dev.off()

# ---- primer check ----
pc <- read.delim(file.path(qc, "primer_check.tsv"))
pc$pct_fwd <- 100 * pc$fwd_27F_at_start / pc$total
pc$pct_rev <- 100 * pc$rev_1492R_at_start / pc$total
pc$pct_with_primer <- pc$pct_fwd + pc$pct_rev
write.table(pc, file.path(qc, "primer_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# summary for the log
cat("\nReads total:", sum(st$num_seqs), "\n")
cat("Reads/sample: min", min(st$num_seqs), "median", median(st$num_seqs), "max", max(st$num_seqs), "\n")
cat("Length quantiles:\n"); print(quantile(pr$length, c(0, .001, .01, .05, .5, .95, .99, .999, 1)))
cat("Fraction of reads 1000-1600 bp:", mean(pr$length >= 1000 & pr$length <= 1600), "\n")
cat("Mean-Q quantiles:\n"); print(quantile(pr$avgQ, c(0, .01, .05, .5, .95, 1)))
cat("Primer check (% reads with a primer in first 40 bp): fwd",
    round(mean(pc$pct_fwd), 1), " rev", round(mean(pc$pct_rev), 1),
    " total", round(mean(pc$pct_with_primer), 1), "\n")
