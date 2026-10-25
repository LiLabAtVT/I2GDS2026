#!/usr/bin/env python3
"""Step 13: PDF report comparing the re-analysis with BakhshiGanje et al. 2026 (Phytobiomes J.)

Output: Report_Reanalysis_vs_Paper.pdf (in the project root directory)
Figures are re-rendered from results/pdf/ into results/report_figs/ with pdftoppm.
Run from the project root directory (env fb_r has reportlab):
    python scripts/13_make_report.py
"""
import os
import subprocess

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table,
                                TableStyle)

# ---------- fonts (DejaVu covers >=, ~, x, minus, Greek) ----------
FD = "/usr/share/fonts/dejavu-sans-fonts"
pdfmetrics.registerFont(TTFont("DV", f"{FD}/DejaVuSans.ttf"))
pdfmetrics.registerFont(TTFont("DV-B", f"{FD}/DejaVuSans-Bold.ttf"))
pdfmetrics.registerFont(TTFont("DV-I", f"{FD}/DejaVuSans-Oblique.ttf"))
pdfmetrics.registerFont(TTFont("DV-BI", f"{FD}/DejaVuSans-BoldOblique.ttf"))
pdfmetrics.registerFontFamily("DV", normal="DV", bold="DV-B", italic="DV-I", boldItalic="DV-BI")

ss = getSampleStyleSheet()
INK, MUTED, ACCENT = colors.HexColor("#1f2933"), colors.HexColor("#52606d"), colors.HexColor("#1f4e79")
body = ParagraphStyle("body", parent=ss["Normal"], fontName="DV", fontSize=9.5, leading=13.5, textColor=INK, spaceAfter=6)
small = ParagraphStyle("small", parent=body, fontSize=8, leading=10.5, spaceAfter=0)
cap = ParagraphStyle("cap", parent=body, fontSize=8, leading=10.5, textColor=MUTED, spaceAfter=10)
h1 = ParagraphStyle("h1", parent=body, fontName="DV-B", fontSize=14, leading=18, textColor=ACCENT, spaceBefore=10, spaceAfter=6)
h2 = ParagraphStyle("h2", parent=body, fontName="DV-B", fontSize=11, leading=14, textColor=INK, spaceBefore=8, spaceAfter=4)
title = ParagraphStyle("title", parent=body, fontName="DV-B", fontSize=18, leading=23, textColor=INK, spaceAfter=4)
sub = ParagraphStyle("sub", parent=body, fontSize=10, leading=14, textColor=MUTED, spaceAfter=12)
bullet = ParagraphStyle("bullet", parent=body, leftIndent=12, bulletIndent=2, spaceAfter=3)
cell = ParagraphStyle("cell", parent=body, fontSize=7.8, leading=10, spaceAfter=0, alignment=TA_LEFT)
cellb = ParagraphStyle("cellb", parent=cell, fontName="DV-B")

MATCH = {"Yes": colors.HexColor("#d8efdc"), "Partly": colors.HexColor("#fdf0c9"), "No": colors.HexColor("#f8d7d4")}

FIG_DIR = "results/report_figs"
os.makedirs(FIG_DIR, exist_ok=True)


def fig(pdf, page, name):
    out = os.path.join(FIG_DIR, name)
    subprocess.run(["pdftoppm", "-png", "-r", "130", "-f", str(page), "-l", str(page), "-singlefile",
                    os.path.join("results/pdf", pdf), out], check=True)
    return out + ".png"


def image(path, width):
    from reportlab.lib.utils import ImageReader
    iw, ih = ImageReader(path).getSize()
    return Image(path, width=width, height=width * ih / iw)


def P(t, s=body):
    return Paragraph(t, s)


def bullets(items):
    return [Paragraph(i, bullet, bulletText="•") for i in items]


def table(rows, widths, header=True, match_col=None):
    data = [[P(c, cellb if (header and r == 0) else cell) for c in row] for r, row in enumerate(rows)]
    t = Table(data, colWidths=widths, repeatRows=1 if header else 0)
    st = [("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#c7ccd1")),
          ("VALIGN", (0, 0), (-1, -1), "TOP"),
          ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4),
          ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]
    if header:
        st.append(("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e4e8ec")))
    if match_col is not None:
        for r, row in enumerate(rows[1:], start=1):
            key = row[match_col].split()[0] if row[match_col] else ""
            if key in MATCH:
                st.append(("BACKGROUND", (match_col, r), (match_col, r), MATCH[key]))
    t.setStyle(TableStyle(st))
    return t


def footer(canv, doc):
    canv.saveState()
    canv.setFont("DV", 7.5)
    canv.setFillColor(MUTED)
    canv.drawString(0.75 * inch, 0.5 * inch, "Re-analysis of PRJNA1348523 vs BakhshiGanje et al. 2026, Phytobiomes J.")
    canv.drawRightString(letter[0] - 0.75 * inch, 0.5 * inch, f"Page {doc.page}")
    canv.restoreState()


W = letter[0] - 1.5 * inch
S = []

# ---------------------------------------------------------------- title + summary
S += [P("Fire blight canker microbiome: re-analysis versus the published results", title),
      P("Raw data: NCBI SRA BioProject PRJNA1348523 (48 PacBio HiFi full-length 16S runs). Paper: BakhshiGanje M. et al. (2026) "
        "<i>Fire Blight Canker Microbiome: Insights into the Key Ecological Niche for Erwinia amylovora</i>. Phytobiomes J., "
        "doi:10.1094/PBIOMES-10-25-0086-R. Analyst: Razan Alajoleen. Re-analysis run 24–25 September 2026 on Virginia Tech ARC (TinkerCliffs).", sub)]

S += [P("1. Summary", h1),
      P("Most of the paper's main findings reproduce from the raw reads. Both the read and ASV totals match within 1%. "
        "Cankers are much richer than healthy bark and are dominated by <i>Erwinia amylovora</i> (Ea). "
        "The top 20 canker genera are the same set as the paper's Fig. 4A. The per-sample Ea profile (Fig. 5C) "
        "matches almost sample for sample, and the bacteria correlated with Ea (Fig. 7B) are the same genera."),
      P("The disagreements have three sources:"),
      *bullets([
          "<b>The healthy-bark samples contain almost no bacterial DNA.</b> 99.5% of their reads are one apple "
          "chloroplast sequence. This one fact explains the largest mismatch: the paper's Bacillus-dominated healthy "
          "microbiome cannot be reproduced. There are 0 Bacillus reads in any healthy sample.",
          "<b>Processing choices that the Methods do not state.</b> The paper's numbers are reproduced only when "
          "chloroplast/mitochondria are removed and the data are <i>not</i> rarefied. The paper describes rarefying and "
          "never mentions organelle removal. Tool settings are also not given (DADA2 settings, which SparCC "
          "implementation, the DMM feature filter).",
          "<b>Some published conclusions are not statistically supported by the same data.</b> The Random Forest "
          "cannot predict cultivar (72% error). DMM model fit prefers 1 cluster, not 4. Cultivar has no significant effect "
          "on canker communities once chloroplast is removed."]),
      P("<b>Conclusions that hold:</b> Ea dominance of cankers; the healthy-vs-canker difference; the main canker "
        "genera; Ea's positive and negative partner taxa; plant-pathogen/fermentation functions in cankers. "
        "<b>Conclusions that do not hold or are fragile:</b> the healthy-bark (Bacillus) microbiome and its network; "
        "cultivar-specific canker taxa; the choice of four DMM community types; network hubs and density.")]

# ---------------------------------------------------------------- methods
S += [P("2. How the re-analysis was done", h1),
      P("Every step followed the paper's Methods where they were specific. Where they were not, standard settings were "
        "used and are listed below. All scripts are in <font name='DV-B'>scripts/</font> and every step is logged in "
        "<font name='DV-B'>README.md</font>."),
      table([
          ["Step", "Paper", "This re-analysis"],
          ["Data", "PRJNA1348523, 32 canker + 16 healthy, 6 cultivars", "Same 48 runs, 3,529,407 reads (matches NCBI). Hewe's Crab runs were sequenced on a Sequel IIe; all others on a Sequel"],
          ["Read QC", "FastQC, quality cut-off", "FastQC + MultiQC + seqkit; primer and length checks"],
          ["Denoising", "QIIME2 2024.10 + DADA2; settings not given", "QIIME2 2024.10 dada2 denoise-ccs; 27F/1492R primers, length 1,000–1,600 bp, maxEE 2, consensus chimeras, minFold 3.5"],
          ["Taxonomy", "SILVA 138", "SILVA 138 99% Naive Bayes (sklearn 1.4.2); BLAST vs NCBI 16S type strains for species names"],
          ["Diversity", "QIIME2 alpha metrics, rarefied; Bray–Curtis PCoA; PERMANOVA", "Same, on 3 tables: paper-style (all ASVs), filtered (organelles removed, rarefied) and relaxed (organelles removed, not rarefied)"],
          ["Random Forest / LEfSe / DMM / FAPROTAX", "randomForest, LEfSe, DirichletMultinomial (k = 1–7), microeco FAPROTAX", "Same packages (microeco 2.4.0, FAPROTAX 1.2.12)"],
          ["Networks", "SparCC in the MENAP web pipeline, |r| ≥ 0.2, P < 0.01", "FastSpar 1.0 (C++ SparCC), 1,000 bootstraps, |r| ≥ 0.2, P < 0.01 and P < 0.05; igraph"],
      ], [0.9 * inch, 2.2 * inch, W - 3.1 * inch])]

# ---------------------------------------------------------------- comparison table
S += [P("3. Result-by-result comparison", h1),
      P("Match: <b>Yes</b> = same result; <b>Partly</b> = same direction or partly the same values; <b>No</b> = not reproduced. "
        "The \"Why\" column is expanded in Section 4.")]
rows = [["Paper result (figure)", "Paper", "Re-analysis", "Match", "Why (section)"],
        ["Reads / ASVs after DADA2", "1,817,968 reads, 2,322 ASVs", "1,839,144 reads, 2,334 ASVs", "Yes", "Small differences from unreported DADA2 settings (4.5)"],
        ["ASVs per sample", "Canker 16–470 (mean 130); healthy 1–22 (mean 5)", "Canker 16–469 (mean 134); healthy 1–12 (mean 4)", "Yes", "—"],
        ["Healthy bark composition (Fig. 2A)", "Bacilli and Acidobacteria dominate", "99.5% chloroplast; 0 Bacilli reads; bacteria left are Methylobacterium, Sphingomonas, Pseudomonas, Hymenobacter, Terriglobus", "No", "Chloroplast (4.1)"],
        ["Canker composition (Fig. 2A)", "Gammaproteobacteria up to 96%, then Alphaproteobacteria", "Gammaproteobacteria 69% mean, Alphaproteobacteria 17%", "Yes", "—"],
        ["Alpha diversity (Fig. 2B)", "Canker > healthy; Cripps Pink canker highest (\"a\"); other cultivars not different", "Canker > healthy (P ≤ 1e-5); Cripps Pink canker \"a\", other cankers \"ab\"; 4 of 6 healthy groups \"b\" (2 are \"ab\")", "Partly", "Richness matches; healthy evenness depends on organelle handling (4.2)"],
        ["Canker Simpson / Shannon", "~0.9 / ~3", "Paper-style 0.57 / 2.0 nats; organelles removed 0.89 / 3.1 nats", "Partly", "Matches only after chloroplast removal (4.2)"],
        ["PCoA axes (Fig. 3)", "23.7% / 6.9%", "Relaxed table 24.6 / 7.6 (paper-style 62.8 / 10.4)", "Yes", "4.2"],
        ["PERMANOVA", "Status R<super>2</super> 17.2% P 0.001; cultivar alone P 0.223; status × cultivar R<super>2</super> 12.4% P 0.016", "Relaxed table: status 15.3% P 0.001; cultivar alone P 0.51; interaction 12.7% P 0.024", "Yes", "4.2; dispersion also differs (4.4)"],
        ["Top-20 canker genera (Fig. 4A)", "Erwinia, Sphingomonas, Massilia, Pantoea, Novosphingobium, …", "Identical 20 genera, top 5 in the same order (organelles removed); paper-style ranks Chloroplast first (53%)", "Yes", "Fig. 4A has no chloroplast, so the paper removed it (4.2)"],
        ["Random Forest by cultivar (Fig. 4B)", "8 cultivar-specific genera", "Top 3 match (Rhizobium, Cutibacterium, Friedmanniella); others rank 14–97; out-of-bag error 72%", "Partly", "Model cannot classify cultivar (4.4)"],
        ["LEfSe (Fig. 5A)", "Canker: Erwinia + ~20 genera. Healthy: Bacillus, Methylobacterium, Sphingomonas, Hymenobacter, Pseudomonas", "Canker: Erwinia (LDA 5.37) + 19 of the paper's canker genera. Healthy: Sphingomonas, Pseudomonas, Hymenobacter. Missing: Bacillus, Methylobacterium", "Partly", "Chloroplast (4.1); LEfSe implementation (4.3)"],
        ["Ea ASVs (Fig. 5C)", "4 Ea ASVs; most in Cripps Pink, then Gala", "10 BLAST-confirmed Ea ASVs, top 4 dominate; Cripps Pink 60%, Gala 53%; per-sample order ≈ identical (Nic1 highest, ~95%)", "Yes", "Ea cultivar differences not significant (P 0.64) (4.4)"],
        ["DMM clusters (Fig. 5B)", "k = 4 by lowest Laplace; clusters 12/11/5/4", "Lowest Laplace k = 1 (k = 1–2 for all filters); at k = 4: 11/8/7/6", "No", "4.4"],
        ["DMM drivers (Fig. 5D)", "4 Ea ASVs, C. flaccumfaciens, N. fluoreni, Methylobacterium, Microbacterium", "4 Ea ASVs, C. flaccumfaciens (100%), Variovorax, N. fluoreni, Methylobacterium", "Yes", "—"],
        ["FAPROTAX (Fig. 6)", "Plant pathogen, fermentation, chemoheterotrophy in cankers; methanol oxidation in healthy Gala / Aztec Fuji / unknown; ureolysis differs", "Plant pathogen 49% vs 0 and fermentation 56% vs 0 (both P = 0.0003); methanol oxidation only in the same healthy samples (not significant); ureolysis n.s.", "Partly", "Species labels needed (4.3); n = 8 healthy (4.4)"],
        ["Ea network partners (Fig. 7B)", "Negative: Pedobacter, Plantibacter, S. roseiflava, Pigmentiphaga, Kineococcus… (r −0.24 to −0.45). Positive: Friedmanniella, Curtobacterium, Novosphingobium, Nocardioides, Methylobacterium (0.26–0.92), Pantoea", "Same negative genera (r −0.32 to −0.44) and same positive genera (0.30–0.93) incl. Pantoea (P < 0.05)", "Yes", "Only at P < 0.05, as quoted for Fig. 7B (4.3)"],
        ["Network size (Table 1)", "Canker T5: 160 nodes / 1,586 edges; T0: 17 / 27; healthy: 51 / 698", "Input ASVs 159 / 15 / 24 (match). Edges at P < 0.01: 230 / 7 / 4; at P < 0.05: 1,011 / 10 / 12", "Partly", "SparCC implementation (4.3); healthy data too sparse (4.1)"],
        ["Network hubs", "Module hubs Novosphingobium, Variovorax, Nocardioides; Ea ASV1032 central in T0", "No module hubs; the 4 Ea ASVs are the most connected nodes in T0", "Partly", "4.3"]]
S += [table(rows, [1.15 * inch, 1.6 * inch, 1.95 * inch, 0.6 * inch, W - 5.3 * inch], match_col=3)]

# ---------------------------------------------------------------- why
S += [P("4. Why the results differ", h1)]

S += [P("4.1 The healthy-bark samples are almost entirely chloroplast DNA", h2),
      P("27F/1492R primers also amplify plant chloroplast 16S, which is bacterial in origin. In this dataset, chloroplast "
        "makes up 56% of all reads (1.02 M of 1.84 M). A single chloroplast ASV makes up about 99.5% of every healthy "
        "sample. Once chloroplast and mitochondria are removed:"),
      *bullets(["8 of the 16 healthy samples have <b>zero</b> bacterial reads;",
                "the other 8 have only 7–177 bacterial reads each, and several cankers are also mostly chloroplast "
                "(e.g. HW2 canker keeps 1%).",
                "Not a single Bacillus read occurs in any healthy sample: all 39 Bacillus ASVs come from cankers."]),
      P("So healthy bark carries very little bacterial DNA compared with plant DNA, while necrotic cankers are full of "
        "bacteria. This is biologically expected. However, it means the data cannot support the paper's healthy-bark "
        "results: a Bacillus-dominated community (Fig. 2A, LEfSe) and a 51-node healthy network (Fig. 8). We cannot tell "
        "from the public data why the paper reports Bacillus. Possible explanations are a different taxonomy step, a "
        "labelling error in the figures, or data that differ from what was deposited in the SRA. The paper does not "
        "report chloroplast handling, negative controls or PNA/blocking primers."),
      image(fig("04_filtered_reads_per_sample.pdf", 1, "f04_chloroplast"), W),
      P("Figure 1. Share of each sample's reads that is chloroplast (green), mitochondria (orange) or bacterial (grey). "
        "Healthy samples are almost entirely chloroplast.", cap)]

S += [P("4.2 Processing steps the Methods do not state", h2),
      P("The Methods say that the ASV table was rarefied, and they do not mention removing organelles. Rarefying "
        "while keeping chloroplast (\"paper-style\", 14,881 reads) reproduces the richness numbers but not the others. "
        "Four independent results instead point to chloroplast removal <i>without</i> rarefaction (\"relaxed\" table, "
        "32 cankers + 8 healthy):"),
      *bullets(["Fig. 4A contains no chloroplast, yet chloroplast would be the #1 canker genus (53%) if it were kept.",
                "PCoA axes: paper 23.7/6.9%; relaxed 24.6/7.6%; paper-style 62.8/10.4%.",
                "PERMANOVA: the paper's significant status × cultivar interaction (R<super>2</super> 12.4%, P 0.016) appears only in the "
                "relaxed table (12.7%, P 0.024). With rarefaction it disappears (P 0.23).",
                "Healthy Simpson 0–0.75 and Shannon 0–2 in the paper; relaxed 0–0.82 and 0–1.9. The paper-style values are "
                "near 0, because each healthy sample is one chloroplast ASV."]),
      P("Rarefying to any depth that keeps healthy samples is impossible once chloroplast is removed (maximum 177 reads). "
        "The paper therefore most likely analysed relative abundances without rarefying. As a result its healthy-bark "
        "statistics rest on a few dozen reads per sample, and half of the healthy samples cannot be included at all. "
        "The PCoA supports this: roughly 9–10 healthy points appear in Fig. 3, not 16."),
      KeepTogether([image(fig("07b_sensitivity_relaxed.pdf", 1, "f07b_pcoa"), W * 0.8),
                    P("Figure 2. Bray–Curtis PCoA of the relaxed table (8 healthy triangles, 32 cankers). The axis "
                      "percentages (24.6/7.6%) match the paper's Fig. 3 (23.7/6.9%).", cap)])]

S += [P("4.3 Different tools or implementations", h2),
      *bullets([
          "<b>LEfSe.</b> microeco labels each genus by the group with the higher <i>median</i>. The original LEfSe cited "
          "by the paper uses the higher <i>mean</i>. In sparse healthy samples the median is 0, so microeco's default "
          "put Sphingomonas, Pseudomonas and Hymenobacter on the canker side. With the original mean rule they are "
          "healthy-enriched, as in the paper. Methylobacterium has a high healthy mean (26% vs 2%) but occurs in "
          "only 4 of 8 healthy samples, so its Kruskal–Wallis test is not significant.",
          "<b>FAPROTAX.</b> <i>plant_pathogen</i> is defined at species level (E. amylovora, Curtobacterium "
          "flaccumfaciens, Pantoea agglomerans…). SILVA 138 gives no species names, so no plant pathogen is detected at "
          "all until species names are added. Here they were added by BLAST (≥ 99% identity and the same genus: 673 ASVs, "
          "70% of reads). The paper must have had species names from some source it does not report.",
          "<b>SparCC networks.</b> MENAP (a web pipeline) cannot be re-run from a script. FastSpar implements the same "
          "algorithm, but its bootstrap P-values are more conservative. The ASV sets going into the networks match the "
          "paper (159 vs 160; 15 vs 17), but at P < 0.01 the canker network has 230 edges versus 1,586. At P < 0.05, "
          "the threshold quoted for Fig. 7B, it has 1,011 edges, and Ea's partner genera and correlation ranges match the paper.",
          "<b>ASV numbers.</b> The paper's ASV IDs (e.g. ASV1032) are not published with sequences, so they cannot be "
          "mapped to ours. Here ASVs are numbered by abundance (ASV0001 = most abundant = Ea).",
          "<b>Shannon.</b> QIIME2 reports log<sub>2</sub>; the paper's values look like natural logs (divide by 1.443)."]),
      KeepTogether([image(fig("12_network_canker_T5_p0.05.pdf", 3, "f12_ea_partners"), W * 0.82),
                    P("Figure 3. SparCC partners of the four main E. amylovora ASVs in the canker network (P < 0.05). "
                      "Negative (red): S. roseiflava, Pedobacter, Plantibacter, Pigmentiphaga, Kineococcus. Positive "
                      "(teal): Curtobacterium, Friedmanniella, Novosphingobium, Methylobacterium, Nocardioides and Pantoea, "
                      "as in the paper's Fig. 7B.", cap)])]

S += [P("4.4 Conclusions that the data do not support statistically", h2),
      *bullets([
          "<b>Cultivar-specific canker taxa.</b> Random Forest has an out-of-bag error of 72% (chance ≈ 80%), so "
          "cultivar cannot be predicted from canker composition. Among the 32 cankers, cultivar is not significant by "
          "PERMANOVA (organelles removed P = 0.35; relaxed P = 0.51) or for any alpha metric. The genera the paper "
          "ranks as cultivar markers (Deinococcus, Aurantisolimonas, Gibbsiella, Shimwellia) are each ≤ 0.09% abundant.",
          "<b>Four DMM community types.</b> The lowest Laplace value is at k = 1 (k = 1–2 with 50, 100 or 200 ASVs or at "
          "genus level). k = 4 is not the best-fitting model, and its clusters are not cultivar-specific.",
          "<b>Ea load by cultivar.</b> Cripps Pink > Gala > others matches the paper, but the differences are not "
          "significant (Kruskal–Wallis P = 0.64; 6 samples per cultivar).",
          "<b>Health-status tests.</b> Canker vs healthy is highly significant, but group dispersions differ strongly "
          "(betadisper P ≤ 0.015). The healthy group rests on 8 samples with very few reads, so the PERMANOVA R<super>2</super> "
          "mixes differences in location and in spread."]),
      KeepTogether([image(fig("10_dmm_clustering.pdf", 1, "f10_laplace"), W * 0.55),
                    P("Figure 4. DMM model fit. The Laplace approximation is lowest at k = 1 (red), not k = 4.", cap)])]

S += [P("4.5 Small, expected differences", h2),
      P("Read and ASV totals differ by about 1% (1,839,144 vs 1,817,968 reads; 2,334 vs 2,322 ASVs). The expected "
        "sources are different DADA2 filter and chimera settings (not reported), software versions and classifier builds. "
        "About 18% of the raw reads are ~1,820 bp host SSU amplicons; the 1,000–1,600 bp length filter removes them. "
        "A different length cut-off in the paper would shift the totals slightly. None of these affect the conclusions.")]

# ---------------------------------------------------------------- what reproduced well (figures)
S += [P("5. Results that reproduce well", h1),
      KeepTogether([image(fig("07_genus_heatmap_cankers_filtered.pdf", 1, "f07_heatmap"), W),
                    P("Figure 5. Top 20 genera in canker samples. They are exactly the 20 genera of the paper's Fig. 4A, with "
                      "Erwinia, Sphingomonas, Massilia, Pantoea and Novosphingobium in the same order.", cap)]),
      KeepTogether([image(fig("10_erwinia_amylovora_asvs_per_canker.pdf", 1, "f10_ea"), W),
                    P("Figure 6. E. amylovora ASVs as % of bacterial reads per canker (BLAST 99.3–99.9% to E. amylovora type "
                      "strains). The sample order and the maximum (Nic1, ~95%) closely match the paper's Fig. 5C. "
                      "E. amylovora has 7 rRNA operon copies, so one genome yields several ASVs.", cap)]),
      KeepTogether([image(fig("05_alpha_diversity_paper.pdf", 3, "f05_observed"), W * 0.9),
                    P("Figure 7. Observed ASVs (rarefied to 14,881 reads, all ASVs). Canker ≫ healthy, and Cripps Pink "
                      "canker is highest (\"a\"), as in Fig. 2B.", cap)]),
      KeepTogether([image(fig("09_lefse_canker_vs_healthy.pdf", 1, "f09_lefse"), W * 0.62),
                    P("Figure 8. LEfSe (LDA > 2, Kruskal–Wallis P < 0.05, original mean-based class assignment). "
                      "Erwinia is the top canker marker. Sphingomonas, Pseudomonas and Hymenobacter are healthy-enriched; "
                      "Bacillus is absent.", cap)])]

# ---------------------------------------------------------------- recommendations + files
S += [P("6. Implications and recommendations", h1),
      *bullets([
          "<b>Robust:</b> Ea dominates canker tissue (median 53% of bacterial reads, up to 95%). Canker and healthy "
          "tissue differ strongly. The canker community (Sphingomonas, Massilia, Pantoea, Novosphingobium, Curtobacterium…) "
          "and Ea's positive and negative partners reproduce.",
          "<b>Treat with caution:</b> all healthy-bark results (only 8 of 16 samples usable, 7–177 reads each); "
          "cultivar-specific differences; the four DMM community types; network density, hubs and the healthy network.",
          "<b>For future sampling of healthy bark:</b> use PNA clamps or chloroplast-blocking primers, sequence more "
          "deeply, include extraction and PCR negative controls, and report organelle filtering and normalisation explicitly.",
          "<b>For reporting:</b> state DADA2 settings, the organelle filter, whether and how data were rarefied, the "
          "LEfSe implementation, the DMM feature filter and model-fit values, the Random Forest error, and deposit ASV "
          "sequences with IDs."]),
      P("7. Files", h1),
      table([["What", "Where (inside the project root)"],
             ["Step-by-step log and re-run commands", "README.md"],
             ["Scripts (steps 00–13)", "scripts/"],
             ["Raw reads, metadata", "raw_data/fastq/ (48 × .fastq.gz), metadata/"],
             ["All figures (42 PDFs, numbered by step); all merged", "results/pdf/; results/ALL_FIGURES_combined.pdf"],
             ["Tables (ASV table, taxonomy, statistics, BLAST, networks)", "results/03_dada2 … results/12_networks"],
             ["This report and its figure images", "Report_Reanalysis_vs_Paper.pdf; results/report_figs/"]],
            [2.6 * inch, W - 2.6 * inch])]

doc = SimpleDocTemplate("Report_Reanalysis_vs_Paper.pdf", pagesize=letter, leftMargin=0.75 * inch, rightMargin=0.75 * inch,
                        topMargin=0.7 * inch, bottomMargin=0.75 * inch,
                        title="Fire blight canker microbiome: re-analysis vs published results", author="Razan Alajoleen")
doc.build(S, onFirstPage=footer, onLaterPages=footer)
print("wrote Report_Reanalysis_vs_Paper.pdf")
