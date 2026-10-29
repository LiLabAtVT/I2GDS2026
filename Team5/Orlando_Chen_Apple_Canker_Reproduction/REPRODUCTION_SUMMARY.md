# Reproduction of BakhshiGanje et al. 2026, Phytobiomes J. 10:277-292

*Fire blight canker microbiome - PacBio full-length 16S, BioProject PRJNA1348523 (48 samples: 32 canker, 16 healthy bark)*

## 1. Summary

The public reads were re-analysed following the paper's Methods (QIIME2 2024.10, DADA2, SILVA 138, rarefaction, Bray-Curtis/PERMANOVA, Random Forest, LEfSe, Dirichlet-multinomial mixtures, FAPROTAX via microeco, SparCC networks). All computation ran on TinkerCliffs compute nodes (Slurm).

Reproduced well: read retention and ASV counts (1,839,176 reads / 2,333 ASVs vs 1,817,968 / 2,322), the canker community (top-5 genera Erwinia, Sphingomonas, Massilia, Pantoea, Novosphingobium in the paper's order), the canker-vs-healthy separation (PERMANOVA P = 0.001), LEfSe canker markers, DMM drivers (4 E. amylovora ASVs, Curtobacterium, Novosphingobium, Methylobacterium), the cultivar ranking of Erwinia ASVs (Cripps Pink > Gala > Aztec Fuji), FAPROTAX fermentation/plant-pathogen enrichment in cankers and methanol oxidation in healthy Gala > Aztec Fuji > unknown, and Erwinia network partners (+Pantoea, +Novosphingobium, -Sphingomonas roseiflava).

Main discrepancy: healthy bark reads are 99.8% apple chloroplast 16S (canker 53%); 8 of 16 healthy samples contain no bacterial reads after chloroplast removal. The paper does not mention chloroplast filtering and reports a Bacillus-dominated healthy microbiome; no Bacillus is found in healthy samples here. Consequently the healthy-bark composition, LEfSe healthy markers (Bacillus, Methylobacterium) and the healthy network (51 nodes in the paper, 2 here) cannot be reproduced from the public data.

Other differences: with chloroplast kept, cultivar is significant and the health x cultivar interaction is not (paper: the reverse); DMM model fit prefers k = 1 (paper presents k = 4); SparCC networks are much sparser than the paper's MENAP networks (canker T5: 251 vs 1,586 edges); FAPROTAX ureolysis does not differ (paper: the only N-cycle function that differs).


## 2. Pipeline (scripts in Claude_Code/scripts/, PDFs in Claude_Code/pdf/)

| Step | What | Scripts | PDF |
|---|---|---|---|
| 1 | SRA/ENA metadata, sample sheet | 01_get_metadata.sh, 01b_make_sample_sheet.py | - |
| 2 | Download 48 HiFi FASTQs (md5 checked) | 02_download_reads.sh | - |
| 3 | FastQC + read length/quality summary | 03_fastqc.sh, 03_read_summary.py | 03_read_QC.pdf |
| 4 | QIIME2 import; DADA2 denoise-ccs (27F/1492R, 1,000-1,600 bp, maxEE 2, consensus chimeras) | 04_dada2.sh, 04b_dada2_summaries.sh, 04_dada2_report.py | 04_dada2_stats.pdf |
| 5 | SILVA 138 99% naive Bayes classifier (sklearn 1.4.2); chloroplast/mitochondria check | 05_taxonomy.sh, 05_taxonomy_report.py, 05b_bacterial_read_check.py | 05_taxonomy.pdf |
| 6 | Phylogeny (MAFFT/FastTree); rarefaction; alpha; Bray-Curtis PCoA; PERMANOVA (A: all, 14,881; B: canker bacteria, 505) | 06_diversity.sh, 06_diversity_report.py | 06_diversity.pdf |
| 7 | Composition, top-20 heatmap, Random Forest, LEfSe, DMM k = 1-7, Erwinia ASVs (bacteria-only, not rarefied) | 07_composition_ml.sh / .R | 07_composition_ml.pdf |
| 8 | FAPROTAX v1.2.6 via microeco (genus-level taxonomy) | 08_faprotax.sh / .R | 08_faprotax.pdf |
| 8b | BLAST vs NCBI 16S RefSeq -> species; FAPROTAX again | 08b_blast_species.sh, 08b_assign_species.py | 08b_faprotax_species.pdf |
| 9 | SparCC (FastSpar, 1,000 bootstraps, |r| >= 0.2, P < 0.01); topology, Zi-Pi, Erwinia edges | 09_networks.sh / .R | 09_networks.pdf |
| 10 | This summary | 10_summary_report.R | 00_reproduction_summary.pdf |

## 3. Paper vs this re-analysis (1/2)

| Result | Paper | Ours | Verdict |
|---|---|---|---|
| Reads after chimera removal | 1,817,968 | 1,839,176 | match |
| ASVs | 2,322 | 2,333 | match |
| ASVs per canker / healthy sample | 16-470 (mean 130) / 1-22 (mean 5) | 19-471 (mean 138) / 1-12 (mean 4) | match |
| Canker top-5 genera (Fig. 4A) | Erwinia, Sphingomonas, Massilia, Pantoea, Novosphingobium | Same five, same order (46.6, 6.4, 4.1, 4.0, 4.0%) | match |
| Canker classes (Fig. 2A) | Gammaproteobacteria (up to 96%), Alphaproteobacteria | Gammaproteobacteria 69%, Alphaproteobacteria 17% | match |
| Healthy composition (Fig. 2A) | Bacilli, Acidobacteria | 99.8% chloroplast; bacteria: Alphaproteobacteria, Gammaproteobacteria, Bacteroidia, Acidobacteriae; no Bacilli | NOT reproduced |
| PERMANOVA health status | R2 17.2%, F 9.35, P 0.001 | A: R2 21.5%, F 16.3, P 0.001 | match |
| PERMANOVA cultivar | n.s. (P 0.223) | A: P 0.002 (chloroplast share varies by cultivar); B: P 0.47 | differs (A) / match (B) |
| PERMANOVA health x cultivar | R2 12.4%, P 0.016 | A: R2 8.3%, P 0.263 | differs |
| PCoA axes (Fig. 3) | 23.7% / 6.9% | A: 62.7% / 10.3% | differs |
| Canker Simpson / Shannon (Fig. 2B) | ~0.85-0.95 / ~1.5-5 | B: 0.71-0.99 (mean 0.89) / 1.5-4.8 (ln) | match (B) |

## 3. Paper vs this re-analysis (2/2)

| Result | Paper | Ours | Verdict |
|---|---|---|---|
| Random Forest by cultivar (Fig. 4B) | Friedmanniella, Rhizobium grp, Cutibacterium, Nocardioides... | Rhizobium grp, Friedmanniella, Cutibacterium, Duganella... (OOB error 69%) | partly (weak model) |
| LEfSe canker (Fig. 5A) | Erwinia LDA > 5; Pantoea, Massilia, Kosakonia, Rhizobium grp... | Erwinia LDA 5.7; 36 canker markers, largely the same list (stable over 10 seeds) | match |
| LEfSe healthy (Fig. 5A) | Bacillus, Methylobacterium, Sphingomonas, Hymenobacter, Pseudomonas | Sphingomonas, Pseudomonas, Hymenobacter, Comamonadaceae; no Bacillus/Methylobacterium | partly |
| DMM clusters (Fig. 5B) | k = 4; sizes 4/5/11/12 | k = 4 sizes 3/4/9/16; best fit k = 1 | partly |
| DMM drivers (Fig. 5D) | 4 x E. amylovora, Curtobacterium, N. fluoreni, Methylobacterium | Same (ASV0001-0004 E. amylovora, Curtobacterium, Novosphingobium, Methylobacterium) | match |
| Erwinia ASV ranking (Fig. 5C) | Cripps Pink > Gala > Aztec Fuji > Granniwinkle > Hewe's Crab > unknown | Cripps Pink (10) > Gala (8) > Aztec Fuji (6) > others (4) | match (top 3) |
| Ea species identity | BLAST: E. amylovora | 10/13 Erwinia ASVs E. amylovora (99.7-99.9%, margin >= 0.9) | match |
| FAPROTAX (Fig. 6) | Canker: fermentation, plant pathogen; healthy: methanol oxidation (Gala > Fuji > unknown); ureolysis differs | Fermentation 56% vs 0; plant pathogen 47% vs 0 (with species); methanol oxidation Gala 62 > Fuji 25 > unknown 18; ureolysis n.s. | mostly |
| Canker network T5 (Fig. 7A, Table 1) | 160 nodes, 1,586 links; hub Mucilaginibacter; avgCC 0.267 | 138 nodes, 251 links; max degree Mucilaginibacter; avgCC 0.267 | partly (sparser) |
| Canker network T0 (Fig. 7C) | 17 nodes; hub E. amylovora ASV1032 | 6 nodes; hub = 4 E. amylovora ASVs | partly |
| Erwinia edges (Fig. 7B/D) | +Pantoea, +Novosphingobium, -S. roseiflava ... | +Pantoea (0.50), +Novosphingobium (0.45), -S. roseiflava (-0.40) | match |
| Healthy network (Fig. 8) | 51 nodes, 698 links, Bacillus hubs | 2 nodes, 1 link (8 samples, 7-177 bacterial reads) | NOT reproduced |

## 4. Decisions, deviations and caveats

DADA2 parameters are not given in the paper; standard PacBio CCS settings were used, with max length 1,600 bp (user's choice). Reads of ~1,800 bp (~19% of raw reads, 23% in healthy samples) are removed by this length filter; they are probably host mitochondrial rRNA but this was not tested.

Chloroplast: the paper's ASV counts and healthy alpha values match only with chloroplast kept, while its canker alpha values, class plots (no Cyanobacteria) and cultivar result match only with chloroplast removed. Step 6 therefore reports both (A: chloroplast kept, all samples; B: canker bacteria only). Steps 7-9 use the bacteria-only table.

Hypothesis (untested): the paper's healthy-bark 'Bacillus' may be mis-labelled chloroplast; our 14 chloroplast ASVs are assigned with confidence > 0.9997.

Species: SILVA 138 assigns no species; species were added by BLAST vs NCBI 16S RefSeq only when unambiguous (>= 99.5% identity, same genus as SILVA, >= 0.2 points ahead of the next species on full-length hits). E. pyrifoliae reaches 99.5% only on partial reference sequences.

Networks: MENAP could not be run; FastSpar (SparCC) with 1,000 bootstraps was used. 'Majority selection' was implemented as prevalence: T5 = ASV in >= 5 of 32 cankers, T0 = in more than half of cankers (as the paper defines 0), healthy T1 = in >= 1 healthy sample.

ASV IDs (ASV0001...) are ours, ranked by abundance; the paper's numbering (e.g. ASV1032) cannot be recovered. Lookup: results/07_composition_ml/asv_lookup.tsv.

Randomness: rarefaction (QIIME2, no seed option) may shift alpha values slightly between runs; R steps use fixed seeds (2026).

Failed/superseded runs are kept in logs/ with _FAILED_ or _superseded_ in the name (read-only home in container, metadata header, LEfSe empty result, marker_table class, broken pipe, species margin rule, bind_rows).


## 5. Where everything is

```
Project folder: /projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code
  raw_data/     48 FASTQs (raw_data/fastq/), SraRunInfo.csv, ena_filereport.tsv, md5_check.txt, sample_metadata.tsv
  scripts/      every script, numbered by step (01-10)
  results/      QIIME2 artifacts (.qza/.qzv) and tables, one folder per step (03_fastqc ... 09_networks)
  pdf/          00_reproduction_summary.pdf (this report) and one PDF per step (03-09)
  reference/    SILVA 138 99% classifier (sklearn 1.4.2), NCBI 16S rRNA BLAST database
  logs/         Slurm logs of every job
  REPRODUCTION_SUMMARY.md   Markdown copy of this report

Software: QIIME2 2024.10 (Apptainer, /apps/common/software/QIIME2/2024.10/); R 4.3.3 conda env ~/miniconda3/envs/canker_r (phyloseq 1.46, randomForest 4.7-1.2, DirichletMultinomial 1.44, microbiomeMarker 1.8, microeco 1.1, igraph 2.1.4); FastSpar env ~/miniconda3/envs/fastspar; BLAST+ 2.16.0.
QIIME2 .qzv files can be viewed at https://view.qiime2.org.
```

