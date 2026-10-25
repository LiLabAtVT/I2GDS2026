# Re-analysis of the fire blight canker microbiome (PRJNA1348523)

Reproduction of the 16S analysis in BakhshiGanje et al. (2026), *Fire Blight Canker Microbiome: Insights into the
Key Ecological Niche for Erwinia amylovora*, Phytobiomes J. 10:277-292, https://doi.org/10.1094/PBIOMES-10-25-0086-R

- **Data:** NCBI SRA BioProject [PRJNA1348523](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1348523): 48 PacBio HiFi full-length
  16S runs (32 canker, 16 healthy bark; 6 apple cultivars). Public, so nothing needs to be shared separately.
- **Main finding:** canker results reproduce well; healthy-bark results do not (about 99.5% of healthy reads are one
  chloroplast ASV). See `docs/Report_Reanalysis_vs_Paper.pdf`.

## Quick start

```bash
# 1. export SBATCH_ACCOUNT=<your_account> and add the reference files (see "Configuration")
# 2. Create the conda environments once
#    sra, qiime2-2024.10, fb_r (R 4.4), fastspar  (see envs/*.yml)
# 3. Run the steps in order, from this directory
python3 scripts/00_make_metadata.py
sbatch scripts/01_download_sra.sh
sbatch scripts/02_read_qc.sh
sbatch scripts/03_dada2_denoise.sh        # about 4 h; calls 03b_dada2_summary.sh
sbatch scripts/04_taxonomy_tree.sh        # needs the SILVA 138 classifier in ref/
sbatch scripts/05_alpha_diversity.sh
# R steps: conda activate fb_r first (see "R environment note")
Rscript scripts/06_beta_diversity.R
Rscript scripts/07_genus_heatmap.R
Rscript scripts/07b_sensitivity_relaxed.R
Rscript scripts/08_random_forest.R
Rscript scripts/09_lefse.R
sbatch scripts/10_erwinia_dmm.sh
bash scripts/11a_blast_all_asvs.sh && Rscript scripts/11_faprotax.R
sbatch scripts/12_sparcc_networks.sh
```

## Pipeline overview

| Step | Script(s) | What it does | Tool |
|---|---|---|---|
| 0-1 | `00_`, `01_` | Build metadata/manifest; download 48 runs | sra-tools |
| 2 | `02_read_qc*` | FastQC, MultiQC, seqkit, primer and length checks | FastQC, seqkit |
| 3 | `03_`, `03b_` | Denoise to ASVs (PacBio CCS mode) | QIIME 2 2024.10, DADA2 |
| 4 | `04_` | SILVA 138 taxonomy, remove chloroplast/mitochondria, tree | classify-sklearn, MAFFT, FastTree |
| 5 | `05_` | Alpha diversity, ANOVA + Tukey HSD | QIIME 2, R |
| 6-7b | `06_`, `07_`, `07b_` | Bray-Curtis PCoA, PERMANOVA, genus heatmap, sensitivity table | vegan |
| 8-9 | `08_`, `09_` | Random Forest by cultivar; LEfSe | randomForest, microeco |
| 10 | `10_` | *E. amylovora* ASVs (BLAST) and Dirichlet multinomial mixtures | DirichletMultinomial |
| 11 | `11_`, `11a_` | FAPROTAX functions (needs species labels from BLAST) | microeco |
| 12 | `12_`, `12b_` | SparCC networks | FastSpar, igraph |
| 13 | `13_make_report.py` | Comparison report PDF | reportlab |

Three tables are compared throughout: **paper-style** (all ASVs, rarefied), **filtered** (organelles removed, rarefied)
and **relaxed** (organelles removed, not rarefied). The relaxed table best matches the paper's numbers.

## Configuration

The scripts no longer contain a personal path or Slurm account. Before running on a cluster:

```bash
export SBATCH_ACCOUNT=<your_slurm_account>   # once per session; replaces the old hard-coded account
```

Every script is run from the project root (the folder that contains `scripts/`, `metadata/`, `results/`).
The `#SBATCH --partition=normal_q` lines are VT ARC defaults; change them for another cluster.

Reference files that are not in the repo (put them in `ref/`):

| File | Used by | How to get it |
|---|---|---|
| `ref/silva-138-99-nb-classifier.qza` | step 4 | QIIME 2 data resources page (must match the QIIME 2 / scikit-learn version) |
| `ref/ncbi_16S/16S_ribosomal_RNA` (BLAST db) | steps 10a, 11a | NCBI 16S rRNA BLAST database, e.g. `update_blastdb.pl --decompress 16S_ribosomal_RNA` in `ref/ncbi_16S/` |

## R environment note

On ARC, a personal `~/.Renviron` can crash the conda R with "illegal operation". Run R steps with:

```bash
conda activate fb_r
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=$CONDA_PREFIX/lib/R/library
```

## Key decisions and known caveats

- The paper gives no DADA2 settings. Used the standard HiFi settings: front 27F, adapter 1492R, length 1,000-1,600 bp, maxEE 2.
- The paper does not state organelle removal or normalisation; results are shown for all three tables above.
- FastSpar replaces the paper's MENAP web pipeline, so networks are sparser at the same P cutoff.
- 8 of 16 healthy samples have zero bacterial reads after the organelle filter; the other 8 have 7-177. Treat all
  healthy-bark results as weak.
- Bug fixed: distance-based network measures (path length, efficiency, closeness...) now use unweighted edges.



## Repository layout

```
scripts/           one numbered script per step
metadata/          runinfo.csv, SRR_Acc_List.txt, metadata.tsv, manifest.tsv
```

