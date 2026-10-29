# Apple fire blight canker microbiome: reproduction project (ALS 5224, Fall 2026)

Re-analysis of **BakhshiGanje et al. 2026, "Fire Blight Canker Microbiome: Insights into the Key Ecological Niche for *Erwinia amylovora*", Phytobiomes Journal 10:277–292** ([doi:10.1094/PBIOMES-10-25-0086-R](https://doi.org/10.1094/PBIOMES-10-25-0086-R)), using the public PacBio HiFi full-length 16S reads (NCBI BioProject **PRJNA1348523**: 48 samples, 32 canker and 16 healthy bark, six apple cultivars).

Orlando Chen, Team 5, I2GDS / ALS 5224, Virginia Tech, Fall 2026. All computation ran as Slurm jobs on VT ARC TinkerCliffs.

**Start here:** [`REPRODUCTION_SUMMARY.md`](REPRODUCTION_SUMMARY.md), or [`pdf/00_reproduction_summary.pdf`](pdf/00_reproduction_summary.pdf), which compares every result with the paper.

## Headline results

| | Paper | This re-analysis |
|---|---|---|
| Non-chimeric reads / ASVs | 1,817,968 / 2,322 | 1,839,176 / 2,333 |
| Canker top-5 genera | Erwinia, Sphingomonas, Massilia, Pantoea, Novosphingobium | same five, same order |
| PERMANOVA, health status | R² 17.2%, P 0.001 | R² 21.5%, P 0.001 |
| Healthy bark microbiome | Bacillus-dominated | 99.8% apple chloroplast reads; not reproducible |

## Pipeline

| Step | What | Script(s) | Output PDF |
|---|---|---|---|
| 1–2 | SRA metadata, sample sheet, download 48 FASTQs | `01_get_metadata.sh`, `01b_make_sample_sheet.py`, `02_download_reads.sh` | – |
| 3 | FastQC and read length/quality summary | `03_fastqc.sh`, `03_read_summary.py` | `03_read_QC.pdf` |
| 4 | QIIME 2 2024.10 import + DADA2 `denoise-ccs` (27F/1492R, 1,000–1,600 bp, maxEE 2) | `04_dada2.sh`, `04b_dada2_summaries.sh`, `04_dada2_report.py` | `04_dada2_stats.pdf` |
| 5 | SILVA 138 99% naive Bayes taxonomy; chloroplast/mitochondria check | `05_taxonomy.sh`, `05_taxonomy_report.py`, `05b_bacterial_read_check.py` | `05_taxonomy.pdf` |
| 6 | Rarefaction, alpha/beta diversity, PERMANOVA | `06_diversity.sh`, `06_diversity_report.py` | `06_diversity.pdf` |
| 7 | Composition, Random Forest, LEfSe, Dirichlet-multinomial, Erwinia ASVs | `07_composition_ml.sh` / `.R` | `07_composition_ml.pdf` |
| 8 | FAPROTAX (microeco) | `08_faprotax.sh` / `.R` | `08_faprotax.pdf` |
| 8b | BLAST vs NCBI 16S → species; FAPROTAX with species | `08b_blast_species.sh`, `08b_assign_species.py` | `08b_faprotax_species.pdf` |
| 9 | SparCC networks (FastSpar) | `09_networks.sh` / `.R` | `09_networks.pdf` |
| 10 | Summary report | `10_summary_report.R` | `00_reproduction_summary.pdf` |

## Repository layout

```
scripts/        every analysis script, numbered by step (*.sh are Slurm job scripts)
pdf/            one PDF report per step + the summary report
results/        result tables (.tsv/.csv), ASV sequences (.fasta), QIIME 2 visualizations (.qzv)
raw_data/       sample metadata only (sample_metadata.tsv, SraRunInfo.csv, SRR_Acc_List.txt, ...)
presentation/   DADA2 Q&A presentation: transcripts (short + full) and slide sources
REPRODUCTION_SUMMARY.md
```

`.qzv` files open at <https://view.qiime2.org>. ASV IDs (ASV0001…) map to QIIME 2 feature hashes in `results/07_composition_ml/asv_lookup.tsv`.

## Not in this repository

Kept on TinkerCliffs only, because of size or because they can be regenerated:
- **Raw reads:** download them from SRA with `scripts/02_download_reads.sh`.
- **Reference databases:** the SILVA 138 classifier and the NCBI 16S BLAST database are fetched by the step 5 and 8b scripts.
- **Large intermediate files:** QIIME 2 data artifacts (`.qza`, including the 2 GB `demux.qza`), BIOM tables, SparCC bootstrap archives and the raw BLAST hit table.
- **Slurm logs.**

## Software

QIIME 2 2024.10 (Apptainer image on TinkerCliffs); R 4.3.3 conda env (phyloseq 1.46, randomForest 4.7-1.2, DirichletMultinomial 1.44, microbiomeMarker 1.8, microeco 1.1 with FAPROTAX 1.2.6, igraph 2.1.4); FastSpar; BLAST+ 2.16.0. Scripts use absolute TinkerCliffs paths (`/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code`); change `B=` at the top of each script to run elsewhere.
