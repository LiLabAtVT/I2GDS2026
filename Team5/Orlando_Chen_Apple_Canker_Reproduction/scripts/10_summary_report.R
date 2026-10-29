#!/usr/bin/env Rscript
# Step 10: one summary report of the reproduction (PDF + Markdown) - pipeline, paper-vs-ours comparison, deviations, paths.
suppressPackageStartupMessages({ library(grid); library(gridExtra) })
B <- "/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code"
PDF <- file.path(B, "pdf/00_reproduction_summary.pdf"); MD <- file.path(B, "REPRODUCTION_SUMMARY.md")

title <- "Reproduction of BakhshiGanje et al. 2026, Phytobiomes J. 10:277-292"
subtitle <- "Fire blight canker microbiome - PacBio full-length 16S, BioProject PRJNA1348523 (48 samples: 32 canker, 16 healthy bark)"
S <- list()

S$overview <- list(h = "1. Summary", p = c(
  "The public reads were re-analysed following the paper's Methods (QIIME2 2024.10, DADA2, SILVA 138, rarefaction, Bray-Curtis/PERMANOVA, Random Forest, LEfSe, Dirichlet-multinomial mixtures, FAPROTAX via microeco, SparCC networks). All computation ran on TinkerCliffs compute nodes (Slurm).",
  "Reproduced well: read retention and ASV counts (1,839,176 reads / 2,333 ASVs vs 1,817,968 / 2,322), the canker community (top-5 genera Erwinia, Sphingomonas, Massilia, Pantoea, Novosphingobium in the paper's order), the canker-vs-healthy separation (PERMANOVA P = 0.001), LEfSe canker markers, DMM drivers (4 E. amylovora ASVs, Curtobacterium, Novosphingobium, Methylobacterium), the cultivar ranking of Erwinia ASVs (Cripps Pink > Gala > Aztec Fuji), FAPROTAX fermentation/plant-pathogen enrichment in cankers and methanol oxidation in healthy Gala > Aztec Fuji > unknown, and Erwinia network partners (+Pantoea, +Novosphingobium, -Sphingomonas roseiflava).",
  "Main discrepancy: healthy bark reads are 99.8% apple chloroplast 16S (canker 53%); 8 of 16 healthy samples contain no bacterial reads after chloroplast removal. The paper does not mention chloroplast filtering and reports a Bacillus-dominated healthy microbiome; no Bacillus is found in healthy samples here. Consequently the healthy-bark composition, LEfSe healthy markers (Bacillus, Methylobacterium) and the healthy network (51 nodes in the paper, 2 here) cannot be reproduced from the public data.",
  "Other differences: with chloroplast kept, cultivar is significant and the health x cultivar interaction is not (paper: the reverse); DMM model fit prefers k = 1 (paper presents k = 4); SparCC networks are much sparser than the paper's MENAP networks (canker T5: 251 vs 1,586 edges); FAPROTAX ureolysis does not differ (paper: the only N-cycle function that differs)."))

S$pipeline <- list(h = "2. Pipeline (scripts in Claude_Code/scripts/, PDFs in Claude_Code/pdf/)", tab = data.frame(check.names = FALSE,
  Step = c("1", "2", "3", "4", "5", "6", "7", "8", "8b", "9", "10"),
  What = c("SRA/ENA metadata, sample sheet", "Download 48 HiFi FASTQs (md5 checked)", "FastQC + read length/quality summary",
           "QIIME2 import; DADA2 denoise-ccs (27F/1492R, 1,000-1,600 bp, maxEE 2, consensus chimeras)",
           "SILVA 138 99% naive Bayes classifier (sklearn 1.4.2); chloroplast/mitochondria check",
           "Phylogeny (MAFFT/FastTree); rarefaction; alpha; Bray-Curtis PCoA; PERMANOVA (A: all, 14,881; B: canker bacteria, 505)",
           "Composition, top-20 heatmap, Random Forest, LEfSe, DMM k = 1-7, Erwinia ASVs (bacteria-only, not rarefied)",
           "FAPROTAX v1.2.6 via microeco (genus-level taxonomy)", "BLAST vs NCBI 16S RefSeq -> species; FAPROTAX again",
           "SparCC (FastSpar, 1,000 bootstraps, |r| >= 0.2, P < 0.01); topology, Zi-Pi, Erwinia edges", "This summary"),
  Scripts = c("01_get_metadata.sh, 01b_make_sample_sheet.py", "02_download_reads.sh", "03_fastqc.sh, 03_read_summary.py",
              "04_dada2.sh, 04b_dada2_summaries.sh, 04_dada2_report.py", "05_taxonomy.sh, 05_taxonomy_report.py, 05b_bacterial_read_check.py",
              "06_diversity.sh, 06_diversity_report.py", "07_composition_ml.sh / .R", "08_faprotax.sh / .R",
              "08b_blast_species.sh, 08b_assign_species.py", "09_networks.sh / .R", "10_summary_report.R"),
  PDF = c("-", "-", "03_read_QC.pdf", "04_dada2_stats.pdf", "05_taxonomy.pdf", "06_diversity.pdf", "07_composition_ml.pdf",
          "08_faprotax.pdf", "08b_faprotax_species.pdf", "09_networks.pdf", "00_reproduction_summary.pdf")))

S$compare1 <- list(h = "3. Paper vs this re-analysis (1/2)", tab = data.frame(check.names = FALSE,
  Result = c("Reads after chimera removal", "ASVs", "ASVs per canker / healthy sample", "Canker top-5 genera (Fig. 4A)",
             "Canker classes (Fig. 2A)", "Healthy composition (Fig. 2A)", "PERMANOVA health status", "PERMANOVA cultivar",
             "PERMANOVA health x cultivar", "PCoA axes (Fig. 3)", "Canker Simpson / Shannon (Fig. 2B)"),
  Paper = c("1,817,968", "2,322", "16-470 (mean 130) / 1-22 (mean 5)", "Erwinia, Sphingomonas, Massilia, Pantoea, Novosphingobium",
            "Gammaproteobacteria (up to 96%), Alphaproteobacteria", "Bacilli, Acidobacteria", "R2 17.2%, F 9.35, P 0.001",
            "n.s. (P 0.223)", "R2 12.4%, P 0.016", "23.7% / 6.9%", "~0.85-0.95 / ~1.5-5"),
  Ours = c("1,839,176", "2,333", "19-471 (mean 138) / 1-12 (mean 4)", "Same five, same order (46.6, 6.4, 4.1, 4.0, 4.0%)",
           "Gammaproteobacteria 69%, Alphaproteobacteria 17%", "99.8% chloroplast; bacteria: Alphaproteobacteria, Gammaproteobacteria, Bacteroidia, Acidobacteriae; no Bacilli",
           "A: R2 21.5%, F 16.3, P 0.001", "A: P 0.002 (chloroplast share varies by cultivar); B: P 0.47", "A: R2 8.3%, P 0.263",
           "A: 62.7% / 10.3%", "B: 0.71-0.99 (mean 0.89) / 1.5-4.8 (ln)"),
  Verdict = c("match", "match", "match", "match", "match", "NOT reproduced", "match", "differs (A) / match (B)", "differs", "differs", "match (B)")))

S$compare2 <- list(h = "3. Paper vs this re-analysis (2/2)", tab = data.frame(check.names = FALSE,
  Result = c("Random Forest by cultivar (Fig. 4B)", "LEfSe canker (Fig. 5A)", "LEfSe healthy (Fig. 5A)", "DMM clusters (Fig. 5B)",
             "DMM drivers (Fig. 5D)", "Erwinia ASV ranking (Fig. 5C)", "Ea species identity", "FAPROTAX (Fig. 6)",
             "Canker network T5 (Fig. 7A, Table 1)", "Canker network T0 (Fig. 7C)", "Erwinia edges (Fig. 7B/D)", "Healthy network (Fig. 8)"),
  Paper = c("Friedmanniella, Rhizobium grp, Cutibacterium, Nocardioides...", "Erwinia LDA > 5; Pantoea, Massilia, Kosakonia, Rhizobium grp...",
            "Bacillus, Methylobacterium, Sphingomonas, Hymenobacter, Pseudomonas", "k = 4; sizes 4/5/11/12", "4 x E. amylovora, Curtobacterium, N. fluoreni, Methylobacterium",
            "Cripps Pink > Gala > Aztec Fuji > Granniwinkle > Hewe's Crab > unknown", "BLAST: E. amylovora",
            "Canker: fermentation, plant pathogen; healthy: methanol oxidation (Gala > Fuji > unknown); ureolysis differs",
            "160 nodes, 1,586 links; hub Mucilaginibacter; avgCC 0.267", "17 nodes; hub E. amylovora ASV1032", "+Pantoea, +Novosphingobium, -S. roseiflava ...",
            "51 nodes, 698 links, Bacillus hubs"),
  Ours = c("Rhizobium grp, Friedmanniella, Cutibacterium, Duganella... (OOB error 69%)", "Erwinia LDA 5.7; 36 canker markers, largely the same list (stable over 10 seeds)",
           "Sphingomonas, Pseudomonas, Hymenobacter, Comamonadaceae; no Bacillus/Methylobacterium", "k = 4 sizes 3/4/9/16; best fit k = 1",
           "Same (ASV0001-0004 E. amylovora, Curtobacterium, Novosphingobium, Methylobacterium)", "Cripps Pink (10) > Gala (8) > Aztec Fuji (6) > others (4)",
           "10/13 Erwinia ASVs E. amylovora (99.7-99.9%, margin >= 0.9)",
           "Fermentation 56% vs 0; plant pathogen 47% vs 0 (with species); methanol oxidation Gala 62 > Fuji 25 > unknown 18; ureolysis n.s.",
           "138 nodes, 251 links; max degree Mucilaginibacter; avgCC 0.267", "6 nodes; hub = 4 E. amylovora ASVs", "+Pantoea (0.50), +Novosphingobium (0.45), -S. roseiflava (-0.40)",
           "2 nodes, 1 link (8 samples, 7-177 bacterial reads)"),
  Verdict = c("partly (weak model)", "match", "partly", "partly", "match", "match (top 3)", "match", "mostly", "partly (sparser)", "partly", "match", "NOT reproduced")))

S$notes <- list(h = "4. Decisions, deviations and caveats", p = c(
  "DADA2 parameters are not given in the paper; standard PacBio CCS settings were used, with max length 1,600 bp (user's choice). Reads of ~1,800 bp (~19% of raw reads, 23% in healthy samples) are removed by this length filter; they are probably host mitochondrial rRNA but this was not tested.",
  "Chloroplast: the paper's ASV counts and healthy alpha values match only with chloroplast kept, while its canker alpha values, class plots (no Cyanobacteria) and cultivar result match only with chloroplast removed. Step 6 therefore reports both (A: chloroplast kept, all samples; B: canker bacteria only). Steps 7-9 use the bacteria-only table.",
  "Hypothesis (untested): the paper's healthy-bark 'Bacillus' may be mis-labelled chloroplast; our 14 chloroplast ASVs are assigned with confidence > 0.9997.",
  "Species: SILVA 138 assigns no species; species were added by BLAST vs NCBI 16S RefSeq only when unambiguous (>= 99.5% identity, same genus as SILVA, >= 0.2 points ahead of the next species on full-length hits). E. pyrifoliae reaches 99.5% only on partial reference sequences.",
  "Networks: MENAP could not be run; FastSpar (SparCC) with 1,000 bootstraps was used. 'Majority selection' was implemented as prevalence: T5 = ASV in >= 5 of 32 cankers, T0 = in more than half of cankers (as the paper defines 0), healthy T1 = in >= 1 healthy sample.",
  "ASV IDs (ASV0001...) are ours, ranked by abundance; the paper's numbering (e.g. ASV1032) cannot be recovered. Lookup: results/07_composition_ml/asv_lookup.tsv.",
  "Randomness: rarefaction (QIIME2, no seed option) may shift alpha values slightly between runs; R steps use fixed seeds (2026).",
  "Failed/superseded runs are kept in logs/ with _FAILED_ or _superseded_ in the name (read-only home in container, metadata header, LEfSe empty result, marker_table class, broken pipe, species margin rule, bind_rows)."))

S$paths <- list(h = "5. Where everything is", p = c(
  paste0("Project folder: ", B),
  "  raw_data/     48 FASTQs (raw_data/fastq/), SraRunInfo.csv, ena_filereport.tsv, md5_check.txt, sample_metadata.tsv",
  "  scripts/      every script, numbered by step (01-10)",
  "  results/      QIIME2 artifacts (.qza/.qzv) and tables, one folder per step (03_fastqc ... 09_networks)",
  "  pdf/          00_reproduction_summary.pdf (this report) and one PDF per step (03-09)",
  "  reference/    SILVA 138 99% classifier (sklearn 1.4.2), NCBI 16S rRNA BLAST database",
  "  logs/         Slurm logs of every job",
  "  REPRODUCTION_SUMMARY.md   Markdown copy of this report",
  "",
  "Software: QIIME2 2024.10 (Apptainer, /apps/common/software/QIIME2/2024.10/); R 4.3.3 conda env ~/miniconda3/envs/canker_r (phyloseq 1.46, randomForest 4.7-1.2, DirichletMultinomial 1.44, microbiomeMarker 1.8, microeco 1.1, igraph 2.1.4); FastSpar env ~/miniconda3/envs/fastspar; BLAST+ 2.16.0.",
  "QIIME2 .qzv files can be viewed at https://view.qiime2.org."))

# ---------------- Markdown ----------------
md <- c(paste("#", title), "", paste0("*", subtitle, "*"), "")
for (s in S) {
  md <- c(md, paste("##", s$h), "")
  if (!is.null(s$p)) md <- c(md, if (s$h == S$paths$h) c("```", s$p, "```") else paste0(s$p, "\n"), "")
  if (!is.null(s$tab)) { t <- s$tab
    md <- c(md, paste0("| ", paste(names(t), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(t)), collapse = "|"), "|"),
            apply(t, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "") }
}
writeLines(md, MD)

# ---------------- PDF ----------------
wrap_tab <- function(t, widths) { for (j in seq_along(t)) t[[j]] <- sapply(t[[j]], function(x) paste(strwrap(x, widths[j]), collapse = "\n")); t }
page <- function(h, body_grob = NULL, paras = NULL, mono = FALSE) {
  grid.newpage()
  grid.text(h, x = 0.04, y = 0.96, just = c("left", "top"), gp = gpar(fontsize = 15, fontface = "bold", col = "#1f3b5a"))
  grid.lines(x = c(0.04, 0.96), y = c(0.935, 0.935), gp = gpar(col = "#1f3b5a"))
  if (!is.null(paras)) {
    y <- 0.91
    for (p in paras) {
      txt <- if (mono && nchar(p) <= 165) p else paste(strwrap(p, if (mono) 165 else 150, exdent = if (mono) 2 else 0), collapse = "\n")
      grid.text(txt, x = 0.05, y = y, just = c("left", "top"),
                gp = gpar(fontsize = if (mono) 9 else 10, fontfamily = if (mono) "mono" else "sans", lineheight = 1.15))
      y <- y - (length(strsplit(txt, "\n")[[1]]) * (if (mono) 0.021 else 0.023) + (if (mono) 0.004 else 0.016))
    }
  }
  if (!is.null(body_grob)) { pushViewport(viewport(y = 0.46, height = 0.86, width = 0.94)); grid.draw(body_grob); popViewport() }
}
th <- ttheme_default(base_size = 10, core = list(fg_params = list(hjust = 0, x = 0.02)), colhead = list(fg_params = list(hjust = 0, x = 0.02)))
col_verdict <- function(g, v) {   # colour the Verdict column
  cols <- ifelse(grepl("^match", v), "#d9f0d3", ifelse(grepl("NOT", v), "#f4cccc", "#fff2cc"))
  nc <- ncol(g); for (i in seq_along(v)) {
    idx <- which(g$layout$t == i + 1 & g$layout$l == nc & g$layout$name == "core-bg"); g$grobs[[idx]]$gp$fill <- cols[i] }
  g }

pdf(PDF, width = 14, height = 8.5)
grid.newpage()
grid.text(title, y = 0.62, gp = gpar(fontsize = 20, fontface = "bold", col = "#1f3b5a"))
grid.text(subtitle, y = 0.55, gp = gpar(fontsize = 12))
grid.text(paste0("Re-analysis of public data on VT ARC TinkerCliffs - ", format(Sys.Date(), "%Y-%m-%d")), y = 0.49, gp = gpar(fontsize = 11, col = "grey30"))
grid.text(paste("Files:", B), y = 0.40, gp = gpar(fontsize = 10, fontfamily = "mono"))
page(S$overview$h, paras = S$overview$p)
page(S$pipeline$h, tableGrob(wrap_tab(S$pipeline$tab, c(4, 70, 45, 26)), rows = NULL, theme = th))
for (k in c("compare1", "compare2")) {
  t <- S[[k]]$tab; g <- tableGrob(wrap_tab(t, c(32, 50, 60, 18)), rows = NULL, theme = th)
  page(S[[k]]$h, col_verdict(g, t$Verdict))
}
page(S$notes$h, paras = S$notes$p)
page(S$paths$h, paras = S$paths$p, mono = TRUE)
dev.off()
cat("Wrote", PDF, "and", MD, "\n")
