# Step 3d: plots from seqkit output -> results/pdf/03_*.pdf
# Called by scripts/03_fastqc.sh
source("scripts/R_common.R")

qc   <- file.path(BASE, "results/01_read_qc")
meta <- read.delim("metadata/sample_metadata.tsv", check.names = FALSE)
meta <- meta[!grepl("^#", meta[[1]]), ]
names(meta)[1] <- "sample"

# ---- per-sample summary ----
st <- read.delim(file.path(qc, "seqkit_stats.tsv"), check.names = FALSE)
st$sample <- sub("\\.fastq\\.gz$", "", basename(st$file))
st <- merge(st, meta[, c("sample", "cultivar", "health")], by = "sample")
st$sample <- factor(st$sample, levels = st$sample[order(st$health, st$cultivar, st$num_seqs)])
write.table(st[, c("sample", "cultivar", "health", "num_seqs", "min_len", "avg_len", "max_len",
                   "Q1", "Q2", "Q3", "AvgQual", "Q20(%)", "Q30(%)")],
            file.path(qc, "per_sample_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(pdf_out("03_reads_per_sample.pdf"), width = 11, height = 5)
print(ggplot(st, aes(sample, num_seqs / 1000, fill = health)) +
  geom_col() + scale_fill_manual(values = HEALTH_COL) +
  labs(x = NULL, y = "Raw HiFi reads (x1000)", fill = "Health status",
       title = "Raw reads per sample (PRJNA1348523)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)))
dev.off()

# ---- per-read length / quality ----
pr <- read.delim(gzfile(file.path(qc, "per_read_len_qual.tsv.gz")), header = FALSE,
                 col.names = c("sample", "length", "avgQ"))
pr <- merge(pr, meta[, c("sample", "health")], by = "sample")
pr$sample <- factor(pr$sample, levels = levels(st$sample))

pdf(pdf_out("03_read_length_distribution.pdf"), width = 8, height = 5)
print(ggplot(pr, aes(length, fill = health)) +
  geom_histogram(binwidth = 10, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = c(1000, 1600), linetype = 2) +
  coord_cartesian(xlim = c(0, 2500)) + scale_fill_manual(values = HEALTH_COL) +
  labs(x = "Read length (bp)", y = "Reads", fill = "Health status",
       title = "Raw read length (dashed: 1,000-1,600 bp)"))
print(ggplot(pr, aes(sample, length, fill = health)) +
  geom_boxplot(outlier.size = 0.2) + scale_fill_manual(values = HEALTH_COL) +
  coord_cartesian(ylim = c(1000, 1800)) +
  labs(x = NULL, y = "Read length (bp)", fill = "Health status") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)))
dev.off()

pdf(pdf_out("03_read_quality_distribution.pdf"), width = 8, height = 5)
print(ggplot(pr, aes(avgQ, fill = health)) +
  geom_histogram(binwidth = 1, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = HEALTH_COL) +
  labs(x = "Mean read Phred quality", y = "Reads", fill = "Health status",
       title = "Mean per-read quality"))
print(ggplot(pr, aes(sample, avgQ, fill = health)) +
  geom_boxplot(outlier.size = 0.2) + scale_fill_manual(values = HEALTH_COL) +
  labs(x = NULL, y = "Mean read Phred quality", fill = "Health status") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)))
dev.off()

# ---- primer / orientation check ----
pc <- read.delim(file.path(qc, "primer_check.tsv"))
pc$pct_fwd <- 100 * pc$fwd_27F_at_start / pc$total
pc$pct_rev <- 100 * pc$rev_1492R_at_start / pc$total
pc$pct_with_primer <- pc$pct_fwd + pc$pct_rev
write.table(pc, file.path(qc, "primer_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pcl <- tidyr::pivot_longer(pc, c(pct_fwd, pct_rev), names_to = "orientation", values_to = "pct")
pcl$orientation <- c(pct_fwd = "27F at start (forward)", pct_rev = "1492R at start (reverse)")[pcl$orientation]
pcl$sample <- factor(pcl$sample, levels = levels(st$sample))
pdf(pdf_out("03_primer_orientation.pdf"), width = 11, height = 5)
print(ggplot(pcl, aes(sample, pct, fill = orientation)) +
  geom_col() + geom_hline(yintercept = 100, linetype = 2) +
  labs(x = NULL, y = "% of reads", fill = NULL,
       title = "Primer found in first 40 bp (IUPAC-aware, + strand)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        legend.position = "top"))
dev.off()

# summary for the log
cat("\nReads total:", sum(st$num_seqs), "\n")
cat("Reads/sample: min", min(st$num_seqs), "median", median(st$num_seqs), "max", max(st$num_seqs), "\n")
cat("Length quantiles:\n"); print(quantile(pr$length, c(0, .001, .01, .05, .5, .95, .99, .999, 1)))
cat("Fraction of reads 1000-1600 bp:", mean(pr$length >= 1000 & pr$length <= 1600), "\n")
cat("Mean-Q quantiles:\n"); print(quantile(pr$avgQ, c(0, .01, .05, .5, .95, 1)))
cat("Primer check (% reads with a primer in first 40 bp): fwd",
    round(mean(pc$pct_fwd), 1), " rev", round(mean(pc$pct_rev), 1),
    " total", round(mean(pc$pct_with_primer), 1), "\n")
