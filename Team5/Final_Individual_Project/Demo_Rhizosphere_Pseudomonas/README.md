# Demo: what does the soybean root recruit from the soil?

**From "which bacteria are enriched?" to "which strain wins in each plot?"**

Orlando Chen · I2GDS 2026 final individual project · Team 5

This demo reproduces one sub-story of a field study of soybean rhizosphere bacteria at two Virginia sites. You start from processed 16S tables (the DADA2 part is already done) and run two connected analyses:

1. **Part 1, differential abundance.** Which bacterial genera does the root enrich or deplete relative to the surrounding soil, at each site? Three methods are compared on the same data: ANCOM-BC2, DESeq2 and a CLR Wilcoxon test. The strongest signal shared by both sites is the genus *Pseudomonas*.
2. **Part 2, *Pseudomonas* variants.** Part 1 says *Pseudomonas* is enriched, but which *Pseudomonas*? You name the sequence variants by BLAST against type strains, find the winning variant in every plot, check whether it was already in the bulk soil, and test whether soybean genotype decides the winner.

Methods that go beyond the class 16S module: compositional differential abundance with bias correction (ANCOM-BC2), a comparison of three differential-abundance methods, a paired (blocked) design, strain-level naming with BLAST and tie handling, and PERMANOVA with restricted permutations.

**Run time:** building the software environment once takes about 15–25 minutes. The demo itself takes about 6 minutes on 8 cores (most of it in ANCOM-BC2).

---

## 1. Study background

- **Sites:**
  - **Kentland Farm** (Blacksburg, VA): Hayter loam.
  - **Eastern Virginia AREC** (Warsaw, VA): Kempsville loam.
- **Plants:** 6 soybean Plant Introduction lines, labelled A–F. A–C have a high root:shoot ratio and D–F a low one. They were chosen from a hydroponic screen of 200 genotypes and grown in 4 randomized blocks at each site in 2025.
- **Samples:** from each plot, at pod stage (R3–R4):
  - **Rhizosphere:** soil attached to the roots.
  - **Bulk soil:** soil between the rows.
- **Sequencing:** 16S rRNA gene V3–V4 region (primers 341F / 805R), Illumina 2 × 250 bp.
- **After quality control:** 94 samples. Kentland has 24 plots × 2 compartments; Warsaw has 23 rhizosphere and 23 bulk samples, of which 22 plots have both.

**Data use:** these are unpublished data from the Gan lab (Virginia Tech, School of Plant and Environmental Sciences), shared as derived tables only, for I2GDS 2026 peer review. Please do not reuse them for other purposes without permission. Raw reads are not included.

---

## 2. Folder contents

```
.
├── README.md
├── environment.yml                  software (conda)
├── scripts/
│   ├── 00_setup_env.sh              builds ./env from environment.yml (run once)
│   ├── run_demo.sh                  runs Part 1 then Part 2
│   ├── 01_differential_abundance.R  Part 1
│   └── 02_pseudomonas_variants.R    Part 2
├── input_dataset/                   everything the demo reads (see section 3)
└── expected_output/                 our results, for checking yours
```

---

## 3. Input dataset (`input_dataset/`)

| File | Format | Contents | Used in |
|---|---|---|---|
| `sample_metadata.tsv` | TSV, 94 rows | `sample_id`, `site` (Kentland / Warsaw), `compartment` (Bulk / Rhizosphere), `plot`, `block` (1–4), `site_plot` (pairs a plot's two samples), `genotype` (PI line), `genotype_label` (A–F), `root_shoot_class` (High / Low) | Parts 1 and 2 |
| `genus_counts.tsv` | TSV, 1,398 genera × 94 samples | **Unrarefied** read counts summed to genus. Reads without a named genus are grouped as `<lowest named rank> (unclassified)`; placeholder names such as "uncultured" are treated as unclassified | Part 1 |
| `genus_taxonomy.tsv` | TSV, 1,398 rows | `genus`, `Phylum`, `Class`, `Order`, `Family` | Part 1 |
| `pseudomonas_asv_counts_rarefied80k.tsv` | TSV, 406 ASVs × 94 samples | Read counts of every *Pseudomonas* ASV in the table **rarefied to 80,000 reads per sample** (one rarefaction, seed 1). Each sample column of the full table sums to 80,000, so `count / 80000` is the share of all reads | Part 2 |
| `pseudomonas_asv_info.tsv` | TSV, 406 rows | `asv` (ASV1, ASV2 … numbered by mean rhizosphere abundance), `asv_md5` (DADA2 sequence ID), `blasted_in_study` (TRUE for the 69 ASVs named in the study: all rhizosphere-enriched ASVs plus the most abundant) | Part 2 |
| `pseudomonas_asv_sequences.fasta` | FASTA, 406 sequences | The *Pseudomonas* ASV sequences (merged V3–V4 reads, about 427 bp), named ASV1, ASV2 … | Part 2 |
| `reference/ncbi16S_typestrains_subset.fasta` | FASTA, 426 sequences | Subset of the NCBI 16S RefSeq type-strain database (Bacteria and Archaea type strains, BLAST database of 3 Oct 2026): all 356 *Pseudomonas* type strains plus every other type strain hit by the study's full-database search. The subset keeps the database order, so tied hits come out in the same order as in the study | Part 2 |
| `reference/source_database_info.txt` | text | Version information of the full database the subset came from | — |

**How these tables were made upstream** (already done; not part of the demo):
1. **Primer removal:** cutadapt removed the primers (QIIME 2 2024.10).
2. **Denoising:** DADA2 with truncation at 280 / 200 bases, one error model for the run, read merging (overlap ≥ 12 bp, no mismatches), consensus chimera removal, and a length filter of 380–460 bp.
3. **Taxonomy:** naive-Bayes classifier on SILVA 138 (99%); only Bacteria and Archaea kept.
4. **Prevalence filter:** ASVs found in at least 2 samples kept.
5. **Rarefaction:** one rarefaction to 80,000 reads (seed 1) for the *Pseudomonas* table.

---

## 4. How to run

On Virginia Tech ARC (replace `<your_allocation>` with your Slurm account):

```bash
git clone https://github.com/LiLabAtVT/I2GDS2026.git          # or copy this folder
cd I2GDS2026/Team5/Final_Individual_Project/Demo_Rhizosphere_Pseudomonas

sbatch --account=<your_allocation> --partition=normal_q scripts/00_setup_env.sh   # once; creates ./env
sbatch --account=<your_allocation> --partition=normal_q scripts/run_demo.sh       # Part 1 + Part 2
```

On another Linux machine with conda: `bash scripts/00_setup_env.sh`, then `bash scripts/run_demo.sh`.

Results are written to `output/`. The log is `run_demo_<jobid>.out`.

**Check your results:**

```bash
diff output/1_summary.txt expected_output/1_summary.txt
diff output/2_summary.txt expected_output/2_summary.txt
```

Both should print nothing, apart from the PERMANOVA p-values in step 2e (see the note there). Each run of this demo gives the same numbers; we tested this by running it twice from a clean copy.

**Software** (`environment.yml`): R 4.4, ANCOMBC 2.8, DESeq2 1.46, CVXR 1.0-15 (ANCOMBC 2.8 needs this older CVXR), vegan, dplyr, tidyr, ggplot2, BLAST+ 2.16.

---

## 5. Part 1: Which genera does the root enrich or deplete? (`01_differential_abundance.R`)

**Input:** `sample_metadata.tsv`, `genus_counts.tsv`, `genus_taxonomy.tsv`

**Design:**
- Each site is analysed separately.
- Each rhizosphere sample is compared with the bulk soil **from the same plot**, so only plots with both samples are used: 24 at Kentland and 22 at Warsaw.
- Pairing removes plot-to-plot differences in soil, so the comparison isolates the effect of the root.

### Step 1A: ANCOM-BC2 (the study's method)

Sequencing gives proportions, not true amounts: when one genus blooms, all others look smaller even if they did not change. ANCOM-BC2 estimates each sample's bias and removes it before testing each genus, so a single bloom does not distort the others.

**Parameters:**

| Parameter | Value | Meaning |
|---|---|---|
| `fix_formula` | `compartment + site_plot` | Rhizosphere vs bulk soil, with plot as a fixed effect (the pairing) |
| reference level | `Bulk` | Fold changes are rhizosphere relative to bulk soil |
| `prv_cut` | 0.10 | A genus must be present in ≥ 10% of the site's samples |
| `lib_cut` | 1000 | Samples with fewer than 1,000 reads are dropped (none are) |
| `s0_perc` | 0.05 | Small variance offset that stabilises the test statistic |
| `struc_zero`, `neg_lb` | TRUE | Detects genera present in only one compartment |
| `p_adj_method` | `"BH"` | Benjamini–Hochberg correction for testing many genera |
| `alpha` | 0.05 | Significance threshold on the adjusted p-value (q) |
| `pseudo_sens` | TRUE | Sensitivity test for the pseudo-count. A genus counts as significant only if q < 0.05 **and** it passes this test |
| seed | 123 | — |

Effect size: `ancombc2_lfc`, the natural-log fold change.

### Step 1B: DESeq2 (comparison)

A negative-binomial count model, widely used for RNA-seq and microbiome data. It does not model compositional bias.
- **Model:** `~ site_plot + compartment`, contrast Rhizosphere vs Bulk.
- **Size factors:** `type = "poscounts"`, needed for data with many zeros.
- **Significance:** BH-adjusted p < 0.05.
- **Genera tested:** the same genera ANCOM-BC2 tested.
- **Effect size:** `deseq2_log2fc`, the log2 fold change.

### Step 1C: Paired Wilcoxon test on CLR values (comparison)

A simple compositional check:
- **Transform:** counts become centred log-ratio (CLR) values, `log(count + 0.5)` minus the sample's mean log value, computed over all genera and all 94 samples.
- **Test:** for each genus, the rhizosphere-minus-bulk differences across plots are tested with a Wilcoxon signed-rank test, then BH-adjusted (q < 0.05).
- **Effect size:** `clr_diff`.

### Step 1D: Combining the sites

The same genus is compared between the two sites and given one of five ANCOM-BC2 categories:
- enriched at both sites;
- depleted at both sites;
- opposite directions;
- significant at one site only;
- not significant.

**Outputs:**

| File | Contents |
|---|---|
| `1_DA_Kentland.tsv`, `1_DA_Warsaw.tsv` | One row per genus: results of all three methods, phylum, direction |
| `1_method_comparison.tsv` | Significant genera per method and site; genera significant by all three; correlation of effect sizes between methods |
| `1_rhizosphere_effect_both_sites.tsv` | ANCOM-BC2 fold change at both sites and the category |
| `1_fig_quadrant_kentland_vs_warsaw.png` | **Rhizosphere effect on each genus at Kentland (x) vs Warsaw (y).** Top right = enriched at both sites |
| `1_fig_ancombc2_vs_deseq2.png` | **Effect sizes from ANCOM-BC2 vs DESeq2 for the same genera** |
| `1_summary.txt` | All of the above as text |

**What you should find:**

| | Kentland | Warsaw |
|---|---|---|
| Genera tested | 618 | 727 |
| ANCOM-BC2 enriched / depleted | 32 / 324 | 220 / 149 |
| DESeq2 enriched / depleted | 87 / 251 | 257 / 221 |
| CLR Wilcoxon enriched / depleted | 46 / 398 | 232 / 249 |
| Correlation of effect sizes, ANCOM-BC2 vs DESeq2 | 0.80 | 0.88 |

- **Agreement across sites:** only **4 genera are enriched at both sites**: *Pseudomonas*, unclassified *Bacillales*, *Psychrobacillus* and *Clostridium sensu stricto 1*. 68 genera are depleted at both sites, and 116 change in opposite directions, so most of the root's effect depends on the site.
- ***Pseudomonas*** is the most strongly enriched named genus at both sites, by every method. ANCOM-BC2 gives a natural-log fold change of 5.79 at Kentland and 4.29 at Warsaw. **This leads to Part 2.**
- **Agreement among methods:** the three methods agree on the direction of most effects but disagree on how many genera are significant. DESeq2 calls more genera enriched; this is the point made by Nearing et al. (2022).

---

## 6. Part 2: Which *Pseudomonas* wins? (`02_pseudomonas_variants.R`)

**Input:** `sample_metadata.tsv`, `pseudomonas_asv_counts_rarefied80k.tsv`, `pseudomonas_asv_info.tsv`, `pseudomonas_asv_sequences.fasta`, `reference/ncbi16S_typestrains_subset.fasta`

### Step 2a: Naming the variants with BLAST

The 69 variants named in the study are searched against the type-strain reference.

**Parameters:**
- **Database:** `makeblastdb -dbtype nucl` built from the reference subset.
- **Search:** `blastn -max_target_seqs 100 -outfmt "6 qseqid sseqid pident length qlen mismatch gapopen evalue bitscore stitle"`.
- **Filter:** a hit counts only if it covers ≥ 95% of the query (`length / qlen ≥ 0.95`).
- **Ties:** every hit within **0.25% identity** of the best hit is a tie, and `n_tied` counts the tied species.
- **Naming:** a variant is named after the first-listed tied species, e.g. "*P. corrugata* (ASV3)".

Names are **putative**: the V3–V4 region often cannot separate closely related species (ASV10 ties with 12 species).

### Step 2b: How much *Pseudomonas*?

*Pseudomonas* reads ÷ 80,000 × 100, per compartment and site.

### Step 2c: The winner in every rhizosphere plot

- **Winner:** the variant with the most reads in each rhizosphere sample.
- **Also recorded:** its share of the sample's *Pseudomonas* reads, and whether it was detected in the bulk soil of the same plot.

### Step 2d: The dominant variant in bulk soil

The most abundant *Pseudomonas* variant in each bulk-soil sample, and how often ASV3 is detected there.

### Step 2e: Does genotype decide the variants?

- **Test:** PERMANOVA (`vegan::adonis2`) on Bray–Curtis dissimilarity of each rhizosphere sample's *Pseudomonas* composition (variant counts ÷ total *Pseudomonas* reads).
- **Model:** `~ genotype`, with 999 permutations kept **within blocks** (`how(blocks = block)`), because plots in the same block share field conditions.
- **Samples:** only those with ≥ 200 *Pseudomonas* reads (all 47 qualify).
- **Seed:** 1.

### Step 2f: Figure

The 8 most abundant variants in every rhizosphere plot, arranged by genotype (A–F) and block (1–4).

**Outputs:**

| File | Contents |
|---|---|
| `2a_blast_hits.tsv`, `2a_variant_names.tsv` | Raw BLAST hits; best identity, number of tied species and closest type strains per variant |
| `2b_pseudomonas_share.tsv` | *Pseudomonas* % of reads per site and compartment |
| `2c_winner_per_plot.tsv` | Winner, its share, and its reads in the plot's own bulk soil |
| `2d_bulk_dominant_per_plot.tsv` | Dominant variant and ASV3 reads in every bulk-soil sample |
| `2e_genotype_permanova.tsv` | PERMANOVA R² and p per site |
| `2_fig_pseudomonas_variants_per_plot.png` | ***Pseudomonas* variants in each rhizosphere plot, by genotype and block** |
| `2_summary.txt` | All of the above as text |

**What you should find:**

| | Kentland | Warsaw |
|---|---|---|
| *Pseudomonas*, % of reads: bulk soil → rhizosphere (mean) | 0.41 → 33.4 | 0.48 → 19.8 |
| Different winning variants among rhizosphere plots | **9 in 24 plots** | **7 in 23 plots** |
| Most frequent winner | ASV1 (*P. silesiensis*), 7 plots | **ASV3 (*P. corrugata*), 12 plots** |
| Median share of *Pseudomonas* reads held by the winner | 69% | 24% |
| Winner detected in its own plot's bulk soil | 1 of 24 plots | 4 of 22 plots |
| Dominant variant in bulk soil | **ASV3 in all 24 plots** | No single one (ASV36: 7, ASV41: 7, ASV9: 3, ASV3: 2, …) |
| ASV3 detected in bulk soil | 24 of 24 | 5 of 23 |
| Genotype effect on variant composition (PERMANOVA) | R² = 0.214, p ≈ 0.5 | R² = 0.224, p ≈ 0.5 |

**Note on p-values:** these vary by about ±0.03 between computers, because permutation tests use random draws. R² does not change. In our run, p = 0.521 (Kentland) and 0.547 (Warsaw).

**The sub-story:**
- **Kentland:** a different *Pseudomonas* variant takes over each plot, usually one that was undetected in that plot's bulk soil. ASV3 dominates Kentland bulk soil, yet it wins no Kentland plot.
- **Warsaw:** ASV3 is rare in bulk soil, yet it wins half of the rhizosphere plots.
- **Genotype:** it does not decide the winner at either site.

In the full study, a null-model analysis (iCAMP) labelled the Kentland pattern "dispersal limitation", which fits priority effects (whichever strain arrives first takes over), and the Warsaw pattern "homogeneous selection".

---

## 7. Things to try

1. **Ignore the pairing.** In Part 1, change the formulas to `compartment` only (ANCOM-BC2) and `~ compartment` (DESeq2). How many genera become significant, and why does pairing matter here?
2. **Prevalence filter.** Change `PRV` from 0.10 to 0.25. Which genera drop out, and does *Pseudomonas* stay on top?
3. **Tie window.** In Part 2, change `TIE` from 0.25 to 0 or 1. How do the names and `n_tied` change for ASV1 and ASV10?
4. **Unrestricted permutations.** Remove `blocks = factor(d$block)` from the PERMANOVA. Does the genotype p-value change?
5. **Look at the figures.** Is the plot-to-plot variation larger between genotypes or within a genotype?

---

## 8. References

- Anderson MJ (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecol* 26:32–46. <https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x>
- Callahan BJ et al. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. *Nat Methods* 13:581–583. <https://doi.org/10.1038/nmeth.3869>
- Camacho C et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics* 10:421. <https://doi.org/10.1186/1471-2105-10-421>
- Debray R et al. (2022). Priority effects in microbiome assembly. *Nat Rev Microbiol* 20:109–121. <https://doi.org/10.1038/s41579-021-00604-w>
- Lin H, Peddada SD (2024). Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures. *Nat Methods* 21:83–91. <https://doi.org/10.1038/s41592-023-02092-7>
- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15:550. <https://doi.org/10.1186/s13059-014-0550-8>
- Nearing JT et al. (2022). Microbiome differential abundance methods produce different results across 38 datasets. *Nat Commun* 13:342. <https://doi.org/10.1038/s41467-022-28034-z>
- Ning D et al. (2020). A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. *Nat Commun* 11:4717. <https://doi.org/10.1038/s41467-020-18560-z>
- Quast C et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Res* 41:D590–D596. <https://doi.org/10.1093/nar/gks1219>
