#!/bin/bash
# Step 10: write README.md (pipeline summary + comparison with the paper; network table appended from step 9)
P=$HOME/fireblight_canker_reanalysis
cat > $P/README.md <<'MD'
# Reanalysis: Fire Blight Canker Microbiome (BakhshiGanje et al. 2026, Phytobiomes J. 10:277-292)

Data: NCBI BioProject PRJNA1348523 (48 PacBio HiFi full-length 16S runs: 32 canker, 16 healthy wood).
Run everything on a compute node with `bash scripts/run_all.sh`. All figures: `results/pdf/`, tables: `results/tables/`, logs: `logs/`.

## Steps
| Step | Script | What it does | Output |
|---|---|---|---|
| 0 | `00_install_R_packages.R` | installs microeco, file2meco, SpiecEasi (Oct-2023 commit; newer needs R>=4.5) | - |
| 1 | `01_get_metadata.sh` | ENA run list -> `metadata/sample_metadata.tsv` | 48 samples |
| 2 | `02_download_raw.sh` | downloads FASTQs, MD5 check; SRR36091178 is corrupt on ENA -> pulled from NCBI SRA | `data/raw/` (1.9 GB) |
| 3 | `03_fastqc.sh` | FastQC | `results/fastqc/` |
| 4 | `04_qiime2_dada2_taxonomy.sh` | QIIME2 2024.10: DADA2 denoise-ccs (27F/1492R), SILVA 138 99% NB classifier, MAFFT/FastTree | `qiime/` |
| 5 | `05_build_phyloseq_qc.R` | phyloseq object, rarefaction, organelle removal, QC | `Step05_QC_reads_and_ASVs.pdf` |
| 6 | `06_composition_alpha_beta.R` | Fig 2A/2B, Fig S2, Fig 3 + PERMANOVA | `Fig2_*.pdf`, `FigS2_*.pdf`, `Fig3_*.pdf` |
| 7 | `07_taxa_RF_LEfSe_DMM.R` | Fig 4A heatmap, 4B Random Forest, 5A LEfSe, 5B-D DMM + E. amylovora ASVs | `Fig4A_*`, `Fig4B_*`, `Fig5A_*`, `Fig5BCD_*` |
| 8 | `08_FAPROTAX.R` | Fig 6 FAPROTAX (microeco) | `Fig6_FAPROTAX_heatmap.pdf` |
| 9 | `09_SparCC_networks.R` | Figs 7-8, Table 1: SparCC (|r|>=0.2, bootstrap P<0.01), igraph topology, Zi/Pi | `Fig7_*`, `Fig8_*`, `Table1_*` |

## Choices not specified in the paper
- DADA2 cutoffs: PacBio HiFi defaults (1000-1600 bp, maxEE 2, minQ 0).
- Rarefaction: full table (incl. chloroplast) rarefied to the minimum depth (14,872), *then* chloroplast/mitochondria removed. Healthy bark is ~99.8% chloroplast; this is the only order consistent with the paper's ASV counts. 8 healthy samples have 0 bacterial reads afterwards and drop out of Bray-Curtis/LEfSe/FAPROTAX.
- E. amylovora ASVs: SILVA gives no species; Erwinia ASVs with >=99.5% identity to the Ea type-strain 16S (NR_041970) are called Ea (also used so FAPROTAX can assign `plant_pathogen`).
- Networks: MENAP (web) replaced by SparCC in R (SpiecEasi, 100 bootstraps). "Majority selection": canker T5 = ASV in >=5 samples; canker T0 = ASV in >half of samples; healthy T1 = ASV in >=1 sample (>=10 reads).
- ASV numbers (ASV0001...) are ranked by abundance here and do not correspond to the paper's ASV IDs.

## Comparison with the paper
| Result | Paper | This reanalysis |
|---|---|---|
| Non-chimeric reads | 1,817,968 | 1,838,653 |
| ASVs | 2,322 | 2,333 (2,220 bacterial after organelle removal) |
| ASVs per canker sample | 16-470, mean 130 | 15-454, mean 131 |
| ASVs per healthy sample | 1-22, mean 5 | 0-8, mean 2 |
| PERMANOVA status | R2 17.2%, P=0.001 | R2 15.3%, P=0.001 |
| PERMANOVA cultivar | P=0.223 | P=0.353 |
| PERMANOVA status x cultivar | R2 12.4%, P=0.016 | R2 12.7%, P=0.023 |
| PCoA axes | 23.7% / 6.9% | 24.7% / 7.6% |
| Top canker genera | Erwinia, Sphingomonas, Massilia, Pantoea, Novosphingobium | Erwinia, Sphingomonas, Pantoea, Novosphingobium, Massilia |
| LEfSe top canker taxon | Erwinia (LDA >5) | Erwinia (LDA 5.37) |
| Random Forest top genera | Rhizobia group, Friedmanniella, Cutibacterium ... | Rhizobia group, Friedmanniella, Duganella, Cutibacterium ... |
| Ea ASVs | 4 ASVs, highest in Nic1; lowest HW5, Nic3, Nic4, CP6, UN2 | 4 top ASVs (99.6-99.7% to Ea), highest in Nic1; lowest Nic4, CP6, UN2, HW5, Nic3 |
| DMM top drivers | 4 Ea ASVs, then Curtobacterium, Novosphingobium ... | 4 Ea ASVs, then Novosphingobium, Variovorax, Pseudomonas, Curtobacterium |
| DMM k | k=4 (lowest Laplace) | Laplace favours k=1 (k=3 with prevalence>=8); k=4 shown as in the paper |
| FAPROTAX canker | fermentation, chemoheterotrophy, plant pathogen | fermentation 56%, anaerobic chemoheterotrophy 58%, plant_pathogen 47% (0% healthy) |
| FAPROTAX healthy | methanol oxidation, methylotrophy, hydrocarbon degradation, ureolysis | same functions enriched in healthy bark |

| Canker network T5 | 160 nodes, 1,586 links | 152 nodes, 425 links (198 +, 227 -) |
| Ea positive partners (T5) | Friedmanniella, Curtobacterium, Novosphingobium, Nocardioides, Methylobacterium-Methylorubrum; Pantoea (T0) | Friedmanniella, Methylobacterium-Methylorubrum, Novosphingobium, Pantoea |
| Ea negative partners (T5) | Pedobacter, Plantibacter, Sphingomonas, Pigmentiphaga, Pseudomonas, ... | Pigmentiphaga, Plantibacter, Sphingomonas |
| Canker network T0 | 17 nodes, 27 links; Ea ASV1032 has max degree | 10 nodes, 10 links; Ea ASVs have max degree |
| Module hubs (T5) | Novosphingobium, Variovorax, Nocardioides | one module hub (Oxalobacteraceae ASV0023) |
| Healthy network (Fig 8) | 51 nodes, 698 links (Bacillus) | **not reproducible**: no significant SparCC edges (healthy bark has <450 bacterial reads in total) |

Known differences: healthy-bark class composition (paper: Bacilli/Acidobacteria; here Alphaproteobacteria dominate) is based on very few bacterial reads (350 in total across all healthy samples after rarefaction) and is not robust. For the same reason the healthy-bark network (Fig 8) could not be rebuilt, and no Fig8 PDF is produced. The paper's dense Bacillus-dominated healthy network cannot be derived from the public reads as processed here.
SparCC is very slow under multi-threaded OpenBLAS on ARC nodes; scripts set OPENBLAS_NUM_THREADS=1.
MD
T=$P/results/tables/Table1_network_topology.csv
if [ -f $T ]; then { echo; echo "## Network topology (Table 1 equivalent, from step 9)"; echo; echo '```'; column -s, -t < $T | cut -c1-200; echo '```'; } >> $P/README.md; fi
echo "README written $(date)"
