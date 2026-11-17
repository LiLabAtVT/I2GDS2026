# Pipeline plan: assembly of soybean rhizosphere bacterial communities at two field sites

**Orlando Chen** · School of Plant and Environmental Sciences, Virginia Tech · Final individual project, I2GDS 2026

---

## 1. Research question and reference paper

**Research question.** How does the soybean rhizosphere bacterial community differ from bulk soil across two Virginia field sites and six genotypes, and what is the relative importance of selection, dispersal and ecological drift in assembling it?

**Reference paper.** Mendes LW, Kuramae EE, Navarrete AA, van Veen JA, Tsai SM (2014). Taxonomical and functional microbial community selection in soybean rhizosphere. *ISME J* 8:1577–1587. <https://doi.org/10.1038/ismej.2014.17>

Mendes et al. showed, at one site, that the soybean rhizosphere is a niche-selected subset of the bulk-soil community. This project asks whether that holds across two soils and six genotypes, and measures how much selection, dispersal and drift each contribute to rhizosphere assembly.

---

## 2. Input data

| Data | Format | Description |
|---|---|---|
| 16S rRNA gene amplicons | Paired-end Illumina FASTQ (`.fastq.gz`), 121 samples × 2 reads | V3–V4 region, primers 341F `CCTACGGGNGGCWGCAG` / 805R `GACTACHVGGGTATCTAATCC`, 2 × 250 bp, one sequencing run. Rhizosphere (soil attached to roots) and bulk soil (between rows) from every plot |
| Sample metadata | TSV | Site, compartment, plot, block, genotype (Plant Introduction line) and root:shoot class for each sample |
| Plot-level field data | CSV | N₂O and CO₂ flux, soil moisture and temperature, shoot and root dry weight, root:shoot ratio, nodule count |
| Taxonomy reference | QIIME 2 classifier (`.qza`) | SILVA 138 (99%), full-length 16S |
| Naming reference | BLAST database | NCBI 16S RefSeq type material (used to name *Pseudomonas* variants) |

- **Study design:** 2 sites (Kentland Farm, Blacksburg, VA, on Hayter loam; Eastern Virginia AREC, Warsaw, VA, on Kempsville loam) × 6 Plant Introduction (PI) lines (3 with high and 3 with low root:shoot ratio, chosen from a hydroponic screen of 200 genotypes) × 4 blocks.
- **Sampling:** pod stage (R3–R4), 2025. After quality control, 94 samples are compared (24 plots × 2 compartments at Kentland; 23 × 2 at Warsaw).
- **Source:** field trial of the Gan lab, Virginia Tech. Earlier per-site processed files are kept only to compare with the new pipeline. Raw data are unpublished and are not stored in this repository.

---

## 3. Two analytical approaches considered

| | **A. ASV + compositional + bin-level null model (selected)** | **B. OTU + count-based + community-level null model** |
|---|---|---|
| **Methods** | DADA2 exact sequence variants (ASVs); ANCOM-BC2 paired by plot; split-plot PERMANOVA; iCAMP (phylogenetic bins) | 97% OTU clustering (UPARSE); DESeq2; PERMANOVA; βNTI / RC-Bray on the whole community |
| **Scientific relevance** | Single-nucleotide resolution lets us follow individual *Pseudomonas* variants plot by plot; assembly processes are estimated for each lineage | One assembly answer per community; closely related *Pseudomonas* variants merge into one OTU |
| **Feasibility** | Computationally heavy (about 1,000 null communities per bin per plot pair); needs HPC, about one day per site | Faster and widely used; runs on a workstation |
| **Expected outcome** | Process shares per site, compartment and lineage; enrichment corrected for compositional bias | Community-level process shares; enrichment lists more sensitive to library size |
| **Limitations** | Processes are inferred, not measured; the V3–V4 region limits species-level names | Lower accuracy on simulated communities; OTU labels are not reproducible across studies; differential-abundance results depend on the tool |

## 4. Selected approach and justification

We selected **Approach A**, for four reasons:

1. **Accuracy and resolution of the assembly analysis.** On simulated communities, iCAMP reached 10–160% higher accuracy than the whole-community null model (Ning et al. 2020). It also reports processes for each phylogenetic bin, which is how we could see that one *Pseudomonas* lineage drove the Kentland pattern.
2. **Exact, reproducible units.** ASVs are exact sequences that mean the same thing in every study (Callahan et al. 2016, 2017). We need them to track which variant won each plot.
3. **Reliable differential abundance.** Differential-abundance tools disagree widely (Nearing et al. 2022), and ANCOM-type methods keep false discoveries low (Weiss et al. 2017). ANCOM-BC2 also corrects for compositional bias and handles paired samples within a plot (Lin & Peddada 2024).
4. **A test that matches the design.** Site varies between plots, while compartment varies within plots, so the PERMANOVA permutes at those two levels (Anderson 2001).

The main cost, computing time, is handled on Virginia Tech's Advanced Research Computing (ARC) cluster.

---

## 5. Methods explained

### 5.1 From raw reads to an ASV table

- **Primer removal (cutadapt; Martin 2011).** cutadapt finds the 341F and 805R primer sequences at the start of each read and cuts them off. Read pairs missing a primer are discarded. Primers are synthetic DNA, not part of the bacteria, and leftover primer variants create false sequence variants. In the earlier Kentland processing, untrimmed primers inflated the variant count from 33,690 to 47,480.
- **Denoising (DADA2; Callahan et al. 2016).** DADA2 learns the error pattern of the sequencing run, then decides which reads are true biological sequences and which are copies with sequencing errors. The result is exact amplicon sequence variants (ASVs) instead of clusters.
  - **Truncation (280 forward / 200 reverse bases):** removes the low-quality read ends while keeping enough overlap to merge the two reads. It was chosen with a pilot run.
  - **Error model:** one model for the whole run, since all samples were sequenced together, with a monotone fit suited to the instrument's binned quality scores.
  - **Merging and cleanup:** paired reads are merged (overlap ≥ 12 bp, no mismatches), chimeras are removed by consensus, and merged lengths outside 380–460 bp are dropped.
- **Taxonomy (naive-Bayes classifier on SILVA 138; Bokulich et al. 2018; Quast et al. 2013).** The classifier breaks each ASV into short words (k-mers) and finds the most probable taxon at each rank (phylum down to genus), with a confidence score. Only Bacteria and Archaea are kept; chloroplast and mitochondrial sequences (plant DNA) are removed.
- **Prevalence filter.** ASVs found in only one sample are mostly noise. Keeping ASVs present in at least 2 samples kept 96.4% of reads.
- **Phylogenetic tree (MAFFT and FastTree; Katoh & Standley 2013; Price et al. 2010).** MAFFT aligns all ASV sequences, and FastTree builds an approximate maximum-likelihood tree. The tree is needed for phylogenetic diversity, NRI/NTI and the iCAMP bins.
- **One analysis object (phyloseq; McMurdie & Holmes 2013).** The ASV table, taxonomy, tree and sample metadata are stored together in one R object, so every later step uses the same data. All steps run inside QIIME 2 (Bolyen et al. 2019) and R.

### 5.2 Diversity within samples (alpha diversity)

- **Rarefaction.** Every sample is randomly subsampled to 80,000 reads so that sequencing depth does not drive the comparison (Weiss et al. 2017). For diversity, this is repeated 100 times and averaged.
- **Hill numbers (Chao et al. 2014).** One family of diversity measures, each in "effective number of ASVs". q = 0 counts all ASVs (richness); q = 1 weights ASVs by abundance (the exponential of Shannon diversity, "common" ASVs); q = 2 counts mainly the dominant ones (inverse Simpson). Looking at all three shows whether a community has lost rare members or become dominated by a few.
- **Linear mixed models (lme4; Bates et al. 2015).** Diversity is modelled with site, compartment and root:shoot class as fixed effects, and block and plot as random effects. Plot must be random because the rhizosphere and bulk samples of a plot are not independent.
- **Phylogenetic structure (NRI and NTI; Webb et al. 2002; picante, Kembel et al. 2010).** These measure whether a sample's abundant ASVs are more closely related than random draws from the tree. Values above 2 mean clustering, a sign that the environment filters for related lineages.

### 5.3 Differences in community composition (beta diversity)

- **Bray–Curtis dissimilarity and PCoA.** Bray–Curtis compares two samples by the abundance of each genus: 0 means identical and 1 means nothing in common. Principal coordinates analysis (PCoA) draws all samples in 2D so that similar samples sit close together. Five other distance measures (ASV-level Bray–Curtis, Jaccard, robust Aitchison, weighted and unweighted UniFrac) serve as checks.
- **Split-plot PERMANOVA (Anderson 2001).** PERMANOVA tests how much of the variation in distances is explained by each factor, comparing it with random shuffles of the labels. Here the shuffles respect the design: site is tested by shuffling whole plots, while compartment and the site × compartment interaction are tested by swapping the two samples within a plot. Genotype is tested within each site and compartment, with shuffles kept within blocks.

### 5.4 Which bacteria the rhizosphere enriches or depletes

- **ANCOM-BC2 (Lin & Peddada 2024).** Sequencing data only give proportions, so an increase in one taxon makes all others look smaller. ANCOM-BC2 estimates and removes this sample-specific bias, then tests each genus.
  - **Model:** each rhizosphere sample is compared with the bulk soil of its own plot (plot as a fixed effect), separately at each site.
  - **Significance:** P values are adjusted for multiple testing (Benjamini–Hochberg, q < 0.05), and a genus counts only if it also passes a pseudo-count sensitivity test.
  - **Safety check:** results are compared with a second, simpler method (paired Wilcoxon test on centred log-ratio values), as Nearing et al. (2022) recommend.
- **Core and site-specific taxa.** A "core" ASV is one found in at least 75% of the samples of a group. Comparing the cores of the two rhizospheres shows how much they share. ASVs found only in one genotype are compared with a null distribution from 999 label shuffles to see whether they exceed chance.

### 5.5 How the community is assembled (iCAMP)

Community ecology describes assembly with four processes (Vellend 2010). **Selection** is ecological filtering, where the environment favours certain traits; it is not natural selection in the evolutionary sense. **Dispersal** is the movement of organisms, **drift** is random change in abundance through chance births and deaths, and **diversification** is the origin of new variants. Null models estimate these from patterns: if two communities are more alike or more different than random communities would be, a deterministic process is at work (Stegen et al. 2013).

iCAMP (Ning et al. 2020) applies this to small groups of related bacteria rather than to the whole community:

1. **Bins:** ASVs are grouped into bins of close relatives on the phylogenetic tree (at least 24 ASVs each).
2. **Plot pairs:** within each site and compartment, every plot is compared with every other plot.
3. **Two tests per bin per pair:**
   - **Phylogenetic test:** the bin's mean phylogenetic distance between the two plots (βMPD) is compared with 1,000 null values. Lower than 97.5% of them means homogeneous selection; higher means heterogeneous selection.
   - **Taxonomic test:** if the phylogenetic test is in the middle 95%, the Bray–Curtis distance is compared with 1,000 random communities drawn from the site's pool. Lower means homogenizing dispersal; higher means dispersal limitation; in the middle means drift.
4. **Adding up:** each bin's labels are weighted by its abundance and averaged over all plot pairs. This gives the share of each process for each group and for each lineage.
5. **Uncertainty:** 95% confidence intervals come from 1,000 bootstrap resamples of plots.

### 5.6 Naming the dominant *Pseudomonas* variants

- **BLAST (Camacho et al. 2009).** Each *Pseudomonas* ASV is searched against NCBI 16S type-strain sequences, requiring at least 95% query coverage; hits within 0.25% of the best identity count as ties. Names are putative, because the V3–V4 region cannot separate closely related species.
- **Winner per plot.** The most abundant variant in each rhizosphere plot is recorded and compared between sites and genotypes (PERMANOVA on *Pseudomonas*-only composition, with permutations kept within blocks).

---

## 6. Pipeline outline and expected outputs

| # | Step | Tools | Intermediate result | Output |
|---|---|---|---|---|
| 1 | Primer removal and QC | cutadapt (QIIME 2) | Trimmed reads (97% kept) | QC summaries |
| 2 | Denoising | DADA2 | Per-sample ASVs | ASV table: 80,049 ASVs, median 334,768 reads per sample |
| 3 | Taxonomy and tree | SILVA 138 classifier; MAFFT, FastTree | Taxonomy, rooted tree | phyloseq object (`.rds`) |
| 4 | Alpha diversity | vegan, iNEXT, lme4, picante | Rarefied tables | Hill numbers, NRI/NTI, model tables |
| 5 | Beta diversity | vegan (Bray–Curtis, PCoA, PERMANOVA) | Distance matrices | Variance explained by site, compartment, genotype |
| 6 | Enrichment | ANCOM-BC2; occupancy analysis | Genus- and ASV-level tests | Enriched and depleted genera; shared core |
| 7 | Assembly | iCAMP; bootstrap | Process share per bin and plot pair | Selection / dispersal / drift shares with 95% CIs |
| 8 | *Pseudomonas* variants | BLAST+; vegan | Typed ASVs | Putative names, winner per plot |

**Final outputs:** result tables and figures for each step, a detailed methods document, a conference poster and a 5-minute talk.

**Computing environment:**
- **Cluster:** Virginia Tech ARC TinkerCliffs, with all jobs submitted through Slurm.
- **Software:** the QIIME 2 2024.10 container; R 4.4.3 with phyloseq, vegan, ANCOMBC 2.8.1, iCAMP 1.9.1, iNEXT, picante and ggplot2; BLAST+.

---

## 7. Validation strategy

- **Against established processing:** compare with the earlier per-site pipeline. Check that taxonomy from forward reads alone agrees with taxonomy from merged reads (r ≥ 0.996).
- **Agreement between methods:** ANCOM-BC2 vs a centred log-ratio Wilcoxon test (consensus, as Nearing et al. 2022 recommend); PERMANOVA on six distance measures.
- **Sensitivity and uncertainty:** unfiltered vs prevalence-filtered tables; iCAMP rerun with a stricter prevalence filter; bootstrap 95% confidence intervals.
- **Null models and power:** genotype-specific ASV counts compared with 999 label permutations; a power simulation showing which genotype effect sizes the design could detect.
- **Cross-check with Approach B and the reference paper:** βNTI / RC-Bray on the same data (planned); and whether the rhizosphere is again a selected subset of bulk soil, as in Mendes et al. (2014).

---

## 8. References

1. Anderson MJ (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecol* 26:32–46. <https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x>
2. Bates D, Mächler M, Bolker B, Walker S (2015). Fitting linear mixed-effects models using lme4. *J Stat Softw* 67(1). <https://doi.org/10.18637/jss.v067.i01>
3. Bokulich NA et al. (2018). Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2's q2-feature-classifier plugin. *Microbiome* 6:90. <https://doi.org/10.1186/s40168-018-0470-z>
4. Bolyen E et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. *Nat Biotechnol* 37:852–857. <https://doi.org/10.1038/s41587-019-0209-9>
5. Callahan BJ et al. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. *Nat Methods* 13:581–583. <https://doi.org/10.1038/nmeth.3869>
6. Callahan BJ, McMurdie PJ, Holmes SP (2017). Exact sequence variants should replace operational taxonomic units in marker-gene data analysis. *ISME J* 11:2639–2643. <https://doi.org/10.1038/ismej.2017.119>
7. Camacho C et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics* 10:421. <https://doi.org/10.1186/1471-2105-10-421>
8. Chao A et al. (2014). Rarefaction and extrapolation with Hill numbers. *Ecol Monogr* 84:45–67. <https://doi.org/10.1890/13-0133.1>
9. Edgar RC (2013). UPARSE: highly accurate OTU sequences from microbial amplicon reads. *Nat Methods* 10:996–998. <https://doi.org/10.1038/nmeth.2604>
10. Katoh K, Standley DM (2013). MAFFT multiple sequence alignment software version 7. *Mol Biol Evol* 30:772–780. <https://doi.org/10.1093/molbev/mst010>
11. Kembel SW et al. (2010). Picante: R tools for integrating phylogenies and ecology. *Bioinformatics* 26:1463–1464. <https://doi.org/10.1093/bioinformatics/btq166>
12. Lin H, Peddada SD (2024). Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures. *Nat Methods* 21:83–91. <https://doi.org/10.1038/s41592-023-02092-7>
13. Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15:550. <https://doi.org/10.1186/s13059-014-0550-8>
14. Martin M (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal* 17:10–12. <https://doi.org/10.14806/ej.17.1.200>
15. McMurdie PJ, Holmes S (2013). phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS One* 8:e61217. <https://doi.org/10.1371/journal.pone.0061217>
16. Mendes LW et al. (2014). Taxonomical and functional microbial community selection in soybean rhizosphere. *ISME J* 8:1577–1587. <https://doi.org/10.1038/ismej.2014.17>
17. Nearing JT et al. (2022). Microbiome differential abundance methods produce different results across 38 datasets. *Nat Commun* 13:342. <https://doi.org/10.1038/s41467-022-28034-z>
18. Ning D et al. (2020). A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. *Nat Commun* 11:4717. <https://doi.org/10.1038/s41467-020-18560-z>
19. Price MN, Dehal PS, Arkin AP (2010). FastTree 2: approximately maximum-likelihood trees for large alignments. *PLoS One* 5:e9490. <https://doi.org/10.1371/journal.pone.0009490>
20. Quast C et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Res* 41:D590–D596. <https://doi.org/10.1093/nar/gks1219>
21. Stegen JC et al. (2013). Quantifying community assembly processes and identifying features that impose them. *ISME J* 7:2069–2079. <https://doi.org/10.1038/ismej.2013.93>
22. Vellend M (2010). Conceptual synthesis in community ecology. *Q Rev Biol* 85:183–206. <https://doi.org/10.1086/652373>
23. Webb CO, Ackerly DD, McPeek MA, Donoghue MJ (2002). Phylogenies and community ecology. *Annu Rev Ecol Syst* 33:475–505. <https://doi.org/10.1146/annurev.ecolsys.33.010802.150448>
24. Weiss S et al. (2017). Normalization and microbial differential abundance strategies depend upon data characteristics. *Microbiome* 5:27. <https://doi.org/10.1186/s40168-017-0237-y>
