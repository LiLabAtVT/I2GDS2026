# Shared settings for the R re-analysis of BakhshiGanje et al. 2026 (PRJNA1348523)
suppressPackageStartupMessages({
  library(phyloseq); library(ggplot2); library(dplyr); library(tidyr)
  library(reshape2); library(magrittr); library(patchwork)
})

BASE <- "/projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code"
setwd(BASE)
# Rscript opens a default PDF device named Rplots.pdf. If no pages are drawn
# on that device (plots go to explicit pdf() files instead), R writes a
# 0-page PDF that viewers report as damaged. pdf(NULL) disables it.
if (!interactive()) pdf(NULL)
EXP  <- file.path(BASE, "results/02_qiime2/export")   # QIIME2 exports
RDS  <- file.path(BASE, "results/03_R/rds")           # intermediate R objects
TAB  <- file.path(BASE, "results/03_R/tables")        # CSV tables
PDF  <- file.path(BASE, "results/pdf")                # all PDF figures
for (d in c(RDS, TAB, PDF)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

SEED <- 2026
CULTIVARS <- c("Aztec Fuji", "Cripps Pink", "Gala", "Granniwinkle", "Hewe's Crab", "Unknown cultivar")
CULT_COL  <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#a6761d"), CULTIVARS)
HEALTH_COL <- c(Canker = "#c0392b", Healthy = "#2e86c1")

theme_set(theme_bw(base_size = 10) +
            theme(panel.grid.minor = element_blank(),
                  strip.background = element_rect(fill = "grey92")))

pdf_out <- function(name, w = 8, h = 6) file.path(PDF, name)
tab_out <- function(name) file.path(TAB, name)

# relative abundance (%) at a taxonomic rank; returns long data frame
rank_relab <- function(ps, rank) {
  p <- tax_glom(ps, rank, NArm = FALSE)
  p <- transform_sample_counts(p, function(x) 100 * x / sum(x))
  psmelt(p)
}

# phyloseq -> microeco microtable, with k__/p__/... rank prefixes microeco expects
to_meco <- function(ps) {
  tt <- as(tax_table(ps), "matrix")
  pre <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  for (j in seq_len(ncol(tt))) tt[, j] <- paste0(pre[j], ifelse(is.na(tt[, j]), "", tt[, j]))
  tax_table(ps) <- tax_table(tt)
  mt <- file2meco::phyloseq2meco(ps)
  mt$tidy_dataset()
  mt
}

# nicer taxon label: fill NA with the lowest known rank
tax_label <- function(tt, rank) {
  tt <- as.data.frame(tt, stringsAsFactors = FALSE)
  ranks <- colnames(tt)
  out <- tt[[rank]]
  for (i in seq_along(out)) {
    if (is.na(out[i])) {
      known <- ranks[seq_len(match(rank, ranks))]
      v <- unlist(tt[i, known]); v <- v[!is.na(v)]
      out[i] <- if (length(v)) paste0("Unclassified ", tail(v, 1)) else "Unassigned"
    }
  }
  out
}

# pheatmap() switches the active device back to the pdf(NULL) above, so the
# heatmap (and every plot after it) never reaches the open PDF file.
# Build it silently, return to the PDF device, then draw it there.
draw_pheatmap <- function(...) {
  d <- dev.cur()
  ph <- pheatmap::pheatmap(..., silent = TRUE)
  dev.set(d); grid::grid.newpage(); grid::grid.draw(ph$gtable)
  invisible(ph)
}
