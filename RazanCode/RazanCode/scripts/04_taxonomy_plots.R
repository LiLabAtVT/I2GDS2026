# Step 4 plots: what was filtered + class-level composition (paper Fig. 2A)
# Run from the project root directory (called by scripts/04_taxonomy_tree.sh)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

d    <- "results/04_taxonomy/export"
meta <- read.delim("metadata/metadata.tsv", comment.char = "#")
names(meta)[1] <- "sample"
names(meta)[names(meta) == "label"] <- "sample_name"

read_tab <- function(f) as.matrix(read.delim(f, skip = 1, check.names = FALSE, row.names = 1))
tab_all <- read_tab(file.path(d, "unfiltered/feature-table.tsv"))
tab     <- read_tab(file.path(d, "feature-table.tsv"))

tax <- read.delim(file.path(d, "taxonomy.tsv"), check.names = FALSE)
names(tax)[1:2] <- c("asv", "taxon")
ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
sp <- strsplit(tax$taxon, ";\\s*")
for (i in seq_along(ranks)) tax[[ranks[i]]] <- sapply(sp, function(x) if (length(x) >= i) x[i] else NA)
write.table(tax, file.path(d, "taxonomy_ranks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- what was removed ----
cat_of <- function(t) ifelse(grepl("mitochondria", t, ignore.case = TRUE), "Mitochondria",
                      ifelse(grepl("chloroplast", t, ignore.case = TRUE), "Chloroplast",
                      ifelse(!grepl("p__", t), "Unassigned / no phylum", "Kept (bacteria/archaea)")))
tax$category <- cat_of(tax$taxon)
rm_s <- tax |> filter(asv %in% rownames(tab_all)) |> group_by(category) |>
  summarise(ASVs = n(), reads = sum(rowSums(tab_all[asv, , drop = FALSE])), .groups = "drop")
print(rm_s)
write.table(rm_s, "results/04_taxonomy/filtering_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

per_s <- lapply(unique(tax$category), function(k) {
  a <- intersect(tax$asv[tax$category == k], rownames(tab_all))
  data.frame(sample = colnames(tab_all), category = k, reads = colSums(tab_all[a, , drop = FALSE]))
}) |> bind_rows() |> left_join(meta, by = "sample")
per_s$sample_name <- factor(per_s$sample_name, levels = meta$sample_name)

pdf("results/pdf/04_filtered_reads_per_sample.pdf", width = 11, height = 5)
print(ggplot(per_s, aes(sample_name, reads, fill = category)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("Kept (bacteria/archaea)" = "grey60", Mitochondria = "#E08214",
                               Chloroplast = "#1B7837", "Unassigned / no phylum" = "#762A83")) +
  facet_grid(~ health_status, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Fraction of DADA2 reads", fill = NULL,
       title = "Reads removed by taxonomy filter (mitochondria, chloroplast, unassigned)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)))
dev.off()

# ---- class-level relative abundance (Fig. 2A) ----
cls <- tax$Class[match(rownames(tab), tax$asv)]
cls[is.na(cls) | cls %in% c("c__", "")] <- "Unclassified"
ctab <- rowsum(tab, cls)
rel  <- sweep(ctab, 2, colSums(ctab), "/")
top  <- names(sort(rowMeans(rel), decreasing = TRUE))[1:min(21, nrow(rel))]  # paper shows 21 classes
long <- as.data.frame(rel) |> mutate(Class = rownames(rel)) |>
  pivot_longer(-Class, names_to = "sample", values_to = "rel") |>
  mutate(Class = ifelse(Class %in% top, Class, "Other")) |>
  group_by(sample, Class) |> summarise(rel = sum(rel), .groups = "drop") |>
  left_join(meta, by = "sample")
long$Class <- factor(long$Class, levels = c("Other", rev(top)))
ord <- meta |> arrange(desc(health_status == "Healthy"), cultivar, as.numeric(replicate))
long$sample_name <- factor(long$sample_name, levels = ord$sample_name)
long$health_status <- factor(long$health_status, levels = c("Healthy", "Canker"))

pal <- c("#1F4E79", "#E07B39", "#1B7837", "#8E44AD", "#2E86C1", "#C0392B", "#F1C40F", "#16A085",
         "#7F8C8D", "#D35400", "#27AE60", "#5B2C6F", "#A93226", "#117864", "#B7950B", "#2874A6",
         "#935116", "#1D8348", "#6C3483", "#B03A2E", "#0E6655", "grey80")
names(pal) <- c(top, "Other")[seq_along(pal)]

pdf("results/pdf/04_class_relative_abundance.pdf", width = 13, height = 6.5)
print(ggplot(long, aes(sample_name, 100 * rel, fill = Class)) +
  geom_col(width = 0.9) + scale_fill_manual(values = pal) +
  facet_grid(~ health_status + cultivar, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Relative abundance (%) (Class level)", fill = "Class") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                     strip.text = element_text(size = 7), legend.text = element_text(size = 7)))
dev.off()
write.table(round(100 * rel, 4), "results/04_taxonomy/class_relative_abundance.tsv", sep = "\t", quote = FALSE)

# ---- console summary ----
cat("\nASVs after filter:", nrow(tab), " reads:", sum(tab), "\n")
cat("Reads/sample after filter: min", min(colSums(tab)), "median", median(colSums(tab)), "\n")
cat("Mean class % by health status:\n")
print(long |> group_by(health_status, Class) |> summarise(mean_pct = round(100 * mean(rel), 1), .groups = "drop") |>
        arrange(health_status, desc(mean_pct)) |> group_by(health_status) |> slice_head(n = 6))
g <- tax$Genus[match(rownames(tab), tax$asv)]
ea <- grepl("Erwinia", g)
cat("Erwinia ASVs:", sum(ea), " species labels:\n")
print(table(tax$Species[match(rownames(tab)[ea], tax$asv)], useNA = "ifany"))
