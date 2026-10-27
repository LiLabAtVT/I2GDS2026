# Step 2b: summarise FastQC (per-base quality, length distribution, read counts) into one PDF.
suppressPackageStartupMessages(library(ggplot2))
setwd("/projects/intro2gds/I2GDS2026/Team1/Kamal/Claude_Code")
dir.create("results/pdf", TRUE, FALSE)
meta <- read.delim("metadata/sample-metadata.tsv", row.names = 1)
for (z in Sys.glob("results/fastqc/*_fastqc.zip")) unzip(z, exdir = "results/fastqc", overwrite = TRUE)
module <- function(f, name) {  # read one ">>module" block of fastqc_data.txt into a data.frame
  l <- readLines(f); s <- grep(paste0("^>>", name), l); e <- grep("^>>END_MODULE", l); e <- e[e > s][1]
  b <- l[(s + 1):(e - 1)]; b <- b[!grepl("^#Total", b)]
  read.delim(text = sub("^#", "", b), check.names = FALSE)
}
files <- Sys.glob("results/fastqc/*_fastqc/fastqc_data.txt")
id <- sub("_fastqc$", "", basename(dirname(files)))
qual <- do.call(rbind, Map(function(f, s) cbind(sample = s, module(f, "Per base sequence quality")), files, id))
qual$pos <- as.numeric(sub("-.*", "", qual$Base)); qual$status <- meta[qual$sample, "status"]
len <- do.call(rbind, Map(function(f, s) cbind(sample = s, module(f, "Sequence Length Distribution")), files, id))
len$len <- as.numeric(sub("-.*", "", len$Length)); len$status <- meta[len$sample, "status"]
basic <- do.call(rbind, Map(function(f, s) { b <- module(f, "Basic Statistics"); data.frame(sample = s,
  reads = as.numeric(b$Value[b$Measure == "Total Sequences"]), length = b$Value[b$Measure == "Sequence length"]) }, files, id))
basic <- cbind(basic, meta[basic$sample, c("sample_name", "status")])
write.csv(basic, "results/tables/fastqc_read_counts.csv", row.names = FALSE)
pdf("results/pdf/Fig0_FastQC_summary.pdf", 11, 7)
print(ggplot(qual, aes(pos, Mean, group = sample, colour = status)) + geom_line(alpha = .5) + theme_bw() +
  geom_hline(yintercept = c(20, 30), linetype = 2) + labs(x = "Position (bp)", y = "Mean Phred quality", title = "FastQC per-base quality (48 HiFi samples)"))
print(ggplot(len, aes(len, Count, group = sample, colour = status)) + geom_line(alpha = .5) + theme_bw() + xlim(1000, 1800) +
  labs(x = "Read length (bp)", y = "Reads", title = "Read length distribution"))
print(ggplot(basic, aes(reorder(sample_name, reads), reads, fill = status)) + geom_col() + coord_flip() + theme_bw() +
  labs(x = NULL, y = "Raw HiFi reads", title = sprintf("Raw reads per sample (total %s)", format(sum(basic$reads), big.mark = ","))))
dev.off()
