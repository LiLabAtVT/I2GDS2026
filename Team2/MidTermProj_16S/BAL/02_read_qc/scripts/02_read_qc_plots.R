#!/usr/bin/env Rscript
# Step 2c: read QC figures and summary numbers (length, quality, depth,
# the ~1,820 bp read class, primer orientation).
# Usage: Rscript 02_read_qc_plots.R <project_dir>
#
# Inputs (written earlier by 02_read_qc.sh):
#   results/02_qc/read_tables/<run>.tsv   length + mean quality of every read
#   results/02_qc/primer_orientation.tsv  primer counts per sample (02_primer_check.py)
#   results/02_qc/samples.tsv             which sample is which cultivar / health status (02a_prepare_samples.py)
# Outputs:
#   results/02_qc/per_sample_summary.tsv  one row of QC numbers per sample
#   results/02_qc/qc_summary.txt          the headline numbers, in plain text
#   results/02_qc/length_top_bins.tsv     the tallest bars of the length histogram
#   results/pdf/02_*.pdf (+ results/png/) the figures

# Load packages quietly: data.table (fast tables), ggplot2 (plots),
# patchwork (combining plots), scales (nicer axis labels, e.g. 1,000)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})
# The project folder is passed in from the command line
proj <- commandArgs(trailingOnly = TRUE)[1]
# Shared colours, plot theme and save_fig() helper
source(file.path(proj, "scripts", "vt_colors.R"))
out <- file.path(proj, "results", "02_qc")
pdf_dir <- file.path(proj, "results", "pdf")

# Clean sample sheet written by 02a_prepare_samples.py (sample-id, run,
# cultivar, health_status). Fix the order of cultivars so every plot lists them
# the same way: known cultivars first (vt_colors.R), then any others A-Z.
meta <- fread(file.path(out, "samples.tsv"))
meta[, cultivar := factor(cultivar, levels = c(intersect(names(cultivar_cols), cultivar),
                                               sort(setdiff(cultivar, names(cultivar_cols)))))]

# per-read table: length + mean quality
# Read every sample's table and stack them into one big table
# (one row per read, ~3.5 million rows), labelled with the run it came from
reads <- rbindlist(lapply(meta$run, function(r) {
  d <- fread(file.path(out, "read_tables", paste0(r, ".tsv")), header = FALSE,
             col.names = c("length", "qual"))
  d[, `:=`(length = as.numeric(length), run = r)]
}))
# Attach sample name, cultivar and health status to every read
reads <- merge(reads, meta[, .(run, `sample-id`, cultivar, health_status)], by = "run")

# The ~1,820 bp class: locate the peak above 1,700 bp, then take +/- 75 bp
# These reads are too long to be bacterial 16S (likely apple mitochondria).
# Find the most common length between 1,700 and 2,100 bp and call everything
# within 75 bp of it the "long class".
long_peak <- reads[length > 1700 & length < 2100, as.numeric(names(which.max(table(length))))]
long_lo <- long_peak - 75; long_hi <- long_peak + 75
# Small yes/no helpers: is a read inside the 1,000-1,600 bp window DADA2 keeps
# in step 3? Is it in the long class?
in_amp  <- function(x) x >= 1000 & x <= 1600
in_long <- function(x) x >= long_lo & x <= long_hi

# --- overall numbers ---------------------------------------------------------------
# The same five numbers, calculated three ways: all reads together, by health
# status (canker vs healthy), and for each sample.
# mean() of a yes/no column = fraction of "yes"; x 100 gives a percentage.
overall <- reads[, .(
  reads = .N, median_length = median(length), median_qual = median(qual),
  pct_1000_1600 = 100 * mean(in_amp(length)), pct_long = 100 * mean(in_long(length)))]
by_status <- reads[, .(
  reads = .N, median_length = median(length), median_qual = median(qual),
  pct_1000_1600 = 100 * mean(in_amp(length)), pct_long = 100 * mean(in_long(length))),
  by = health_status]

per_sample <- reads[, .(
  reads = .N, median_length = median(length), median_qual = median(qual),
  pct_1000_1600 = 100 * mean(in_amp(length)), pct_long = 100 * mean(in_long(length))),
  by = .(`sample-id`, run, cultivar, health_status)]
# Add the primer counts from 02_primer_check.py and turn them into percentages
primers <- fread(file.path(out, "primer_orientation.tsv"))
per_sample <- merge(per_sample, primers[, .(`sample-id`, forward, reverse, none)], by = "sample-id")
per_sample[, `:=`(pct_forward = 100 * forward / reads, pct_reverse = 100 * reverse / reads,
                  pct_none = 100 * none / reads)]
setorder(per_sample, cultivar, health_status, `sample-id`)
fwrite(per_sample, file.path(out, "per_sample_summary.tsv"), sep = "\t")

# Do healthy samples have a higher share of long reads than cankers?
# Wilcoxon rank-sum test: compares two groups without assuming a bell curve,
# with each sample counted once (so big samples don't dominate).
wt <- wilcox.test(pct_long ~ health_status, data = per_sample)
# sink() sends everything printed into qc_summary.txt until sink() is called again
sink(file.path(out, "qc_summary.txt"))
cat("Overall\n"); print(overall)
cat("\nBy health status\n"); print(by_status)
cat(sprintf("\nReads per sample: %d - %d\n", min(per_sample$reads), max(per_sample$reads)))
cat(sprintf("Long-read peak: %d bp (class %d-%d bp)\n", long_peak, long_lo, long_hi))
cat(sprintf("Long-read class, %% per sample: canker median %.1f, healthy median %.1f, Wilcoxon P = %.2g\n",
            per_sample[health_status == "Canker", median(pct_long)],
            per_sample[health_status == "Healthy", median(pct_long)], wt$p.value))
cat(sprintf("Primer orientation (all reads): forward %.1f%%, reverse %.1f%%, none %.1f%%\n",
            100 * sum(primers$forward) / sum(primers$total), 100 * sum(primers$reverse) / sum(primers$total),
            100 * sum(primers$none) / sum(primers$total)))
sink()
# Also print the summary into the job log
cat(readLines(file.path(out, "qc_summary.txt")), sep = "\n")

# --- figures ---------------------------------------------------------------------
# Shortcut: save a plot into results/pdf/ (save_fig also writes a PNG copy)
save_pdf <- function(p, name, w = 8, h = 5) save_fig(p, file.path(pdf_dir, name), w, h)

# 1. read length distribution
# Histogram of read lengths in 10 bp bars, canker and healthy overlaid.
# "Density" scales each group to the same total area, so the shapes can be
# compared even though cankers have more reads. The grey band marks the
# 1,000-1,600 bp window that DADA2 keeps in step 3.
p_len <- ggplot(reads[length <= 2500], aes(length, fill = health_status, colour = health_status)) +
  annotate("rect", xmin = 1000, xmax = 1600, ymin = -Inf, ymax = Inf, fill = vt[["smoke"]], alpha = 0.7) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 10, alpha = 0.55, position = "identity", linewidth = 0.1) +
  annotate("text", x = 1300, y = Inf, vjust = 1.5, label = "1,000-1,600 bp\n(DADA2 window)", size = 3) +
  annotate("text", x = long_peak, y = Inf, vjust = 1.5, label = paste0("~", comma(long_peak), " bp"), size = 3) +
  scale_fill_manual(values = status_cols) + scale_colour_manual(values = status_cols) +
  scale_x_continuous(labels = comma) +
  labs(title = "Read length distribution", x = "Read length (bp)", y = "Density",
       fill = NULL, colour = NULL) + theme_vt()
save_pdf(p_len, "02_read_length_distribution.pdf")

# 2. per-read mean quality
# Histogram of each read's average Phred quality score (higher = fewer errors;
# Q30 = 1 error in 1,000 bases). The dashed line is the overall median.
p_q <- ggplot(reads, aes(qual, fill = health_status)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.5, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = overall$median_qual, linetype = 2, colour = vt[["stone"]]) +
  scale_fill_manual(values = status_cols) +
  labs(title = "Per-read mean quality", subtitle = sprintf("Median Q%.1f", overall$median_qual),
       x = "Mean Phred quality", y = "Density", fill = NULL) + theme_vt()
save_pdf(p_q, "02_read_quality.pdf")

# 3. reads per sample
# One bar per sample showing how many reads it has (sequencing depth),
# grouped into panels by cultivar. Very low bars would flag weak samples.
# Lock the sample order to the sorted table so bars appear in that order
per_sample[, `sample-id` := factor(`sample-id`, levels = `sample-id`)]
p_depth <- ggplot(per_sample, aes(`sample-id`, reads, fill = health_status)) +
  geom_col() + facet_grid(~ cultivar, scales = "free_x", space = "free_x", labeller = label_wrap_gen(10)) +
  scale_fill_manual(values = status_cols) + scale_y_continuous(labels = comma) +
  labs(title = "Reads per sample", x = NULL, y = "HiFi reads", fill = NULL) +
  theme_vt() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
save_pdf(p_depth, "02_reads_per_sample.pdf", w = 11, h = 5)

# 4. ~1,820 bp read class by health status
# Box plot of the % of long reads in each sample, canker vs healthy.
# Each dot is one sample, coloured by cultivar; the subtitle gives the
# Wilcoxon test result from above.
p_long <- ggplot(per_sample, aes(health_status, pct_long, fill = health_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(colour = cultivar), width = 0.15, size = 2) +
  scale_fill_manual(values = status_cols, guide = "none") +
  scale_colour_manual(values = cultivar_cols) +
  labs(title = sprintf("Reads of %d-%d bp (likely host mitochondria)", long_lo, long_hi),
       subtitle = sprintf("Wilcoxon P = %.2g", wt$p.value),
       x = NULL, y = "% of reads in sample", colour = NULL) + theme_vt()
save_pdf(p_long, "02_long_read_class.pdf", w = 6, h = 5)

# 5. primer orientation
# One stacked bar per sample: % of reads starting with 27F, with 1492R, or
# with neither. melt() reshapes the three % columns into one column so
# ggplot can stack them.
po <- melt(per_sample[, .(`sample-id`, cultivar, pct_forward, pct_reverse, pct_none)],
           id.vars = c("sample-id", "cultivar"), variable.name = "orientation", value.name = "pct")
po[, orientation := factor(orientation, c("pct_forward", "pct_reverse", "pct_none"),
                           c("Starts with 27F", "Starts with 1492R", "No primer found"))]
p_po <- ggplot(po, aes(`sample-id`, pct, fill = orientation)) +
  geom_col() + facet_grid(~ cultivar, scales = "free_x", space = "free_x", labeller = label_wrap_gen(10)) +
  scale_fill_manual(values = c(vt[["maroon"]], vt[["orange"]], vt[["grey"]])) +
  labs(title = "Read orientation (primer at read start)", x = NULL, y = "% of reads", fill = NULL) +
  theme_vt() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
save_pdf(p_po, "02_primer_orientation.pdf", w = 11, h = 5)
cat("figures written to", pdf_dir, "\n")

# --- read-count versions (to compare directly with the course demo step 2 figure) ---
# Same layout as the demo figure: stacked counts, 10 bp bins, 0-2,600 bp, dashed 1,000/1,600 bp.
# Here the y-axis is the actual number of reads (not density), and the canker
# and healthy bars are stacked, so each bar's height is the total for that length.
p_len_n <- ggplot(reads[length <= 2600], aes(length, fill = health_status)) +
  geom_histogram(binwidth = 10, boundary = 0) +
  geom_vline(xintercept = c(1000, 1600), linetype = 2) +
  scale_fill_manual(values = status_cols) +
  scale_x_continuous(labels = comma, limits = c(0, 2600)) + scale_y_continuous(labels = comma) +
  labs(title = "Raw read length (dashed: 1,000-1,600 bp)", x = "Read length (bp)", y = "Reads",
       fill = "Health status") + theme_vt()
save_pdf(p_len_n, "02_read_length_counts.pdf")

# Quality histogram, also as stacked read counts
p_q_n <- ggplot(reads, aes(qual, fill = health_status)) +
  geom_histogram(binwidth = 1, boundary = 0) +
  geom_vline(xintercept = overall$median_qual, linetype = 2, colour = vt[["stone"]]) +
  scale_fill_manual(values = status_cols) + scale_y_continuous(labels = comma) +
  labs(title = "Per-read mean quality", subtitle = sprintf("Median Q%.1f", overall$median_qual),
       x = "Mean Phred quality", y = "Reads", fill = "Health status") + theme_vt()
save_pdf(p_q_n, "02_read_quality_counts.pdf")

# tallest 10 bp bins, for comparing bar heights with the course demo figure
# ggplot bins are right-closed: (1440, 1450] -> reads 1,441-1,450 bp
# Label every read with its 10 bp bin, count reads per bin and group, and
# save the 6 tallest bins (with canker and healthy counts side by side).
reads[, bin := paste0(ceiling(length / 10) * 10 - 9, "-", ceiling(length / 10) * 10)]
top_bins <- reads[, .N, by = .(bin, health_status)][
  , total := sum(N), by = bin][order(-total)]
top_bins <- dcast(top_bins, bin + total ~ health_status, value.var = "N", fill = 0)[order(-total)][1:6]
fwrite(top_bins, file.path(out, "length_top_bins.tsv"), sep = "\t")
cat("\nTallest 10 bp length bins (reads):\n"); print(top_bins)
