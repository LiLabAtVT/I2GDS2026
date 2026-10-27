# Step 2: Read quality control (PacBio HiFi full-length 16S)

Part of the Team 2 re-analysis of BakhshiGanje et al. 2026 (*Phytobiomes Journal*), fire blight canker microbiome, BioProject **PRJNA1348523** (48 samples: 32 canker, 16 healthy, 6 apple cultivars).

This step checks the raw reads **before** denoising: read length, quality, depth per sample, and primer orientation. It is a diagnostic step. No reads are removed here; filtering happens in step 3 (DADA2 `denoise-ccs`).

## Quick start

Run from the project folder you used for step 1 (the one containing `raw_data/fastq/` and `metadata/metadata.tsv`):

```bash
git clone https://github.com/lilabatvt/i2gds2026.git ~/i2gds2026        # once, to get the scripts
mkdir -p scripts logs
cp ~/i2gds2026/Team2/MidTermProj_16S/BAL/02_read_qc/scripts/* scripts/
sbatch --account=introtogds --partition=normal_q --qos=tc_normal_short \
       --export=ALL,CONDA_MODULE=Miniforge3/24.11.3-0 scripts/02_read_qc.sh
```

- The settings above are for the I2GDS 2026 class on VT ARC. Elsewhere, use your own `--account` and `--partition`, and drop `CONDA_MODULE` if `conda` already works in your terminal.
- Takes ~20 min the first time (it installs the software), a few minutes after that. Check on it with `squeue -u $USER`; the log is in `logs/`.
- Output: tables in `results/02_qc/`, figures in `results/pdf/` and `results/png/`. See [FIGURES.md](FIGURES.md) for how to read them.
- If step 1 put things somewhere else, see [Using the output of step 1](#using-the-output-of-step-1).

Everything below is reference.

## Where the data comes from

The reads are public, so nothing needs to be shared separately:

- **Paper:** BakhshiGanje et al. (2026). Fire Blight Canker Microbiome: Insights into the Key Ecological Niche for *Erwinia amylovora*. *Phytobiomes Journal* 10:277–292. https://doi.org/10.1094/PBIOMES-10-25-0086-R
- **Reads:** NCBI SRA BioProject [PRJNA1348523](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1348523): 48 PacBio HiFi runs (one SRR run per sample), amplified with the 27F/1492R full-length 16S primers. About 2 GB as `.fastq.gz`.
- **Download:** not part of this folder. Step 1 (downloading the reads and making the sample sheet) is uploaded separately by another team in this repository. Step 2 only needs its output; see [Using the output of step 1](#using-the-output-of-step-1). To download a run by hand instead: `prefetch SRR36091169 && fasterq-dump SRR36091169 && gzip SRR36091169.fastq` (sra-tools).

## Using the output of step 1

Run step 2 from the **same project folder** as step 1. It needs two things step 1 makes:

| Step 1 output | Default location | What step 2 needs from it |
|---|---|---|
| FASTQ files | `raw_data/fastq/` | one `<run>.fastq.gz` per sample, named by SRR number |
| Sample sheet | `metadata/metadata.tsv` | tab-separated; one row per sample with its run number, a sample name, the cultivar and the health status |

You do not need to rename columns. The first script, `02a_prepare_samples.py`, finds them automatically and writes a clean copy to `results/02_qc/samples.tsv`, which the rest of step 2 uses. It recognises:

| Needed | Column names it looks for |
|---|---|
| run number | `run`, or any column whose values are all SRR/ERR/DRR numbers |
| sample name | `sample-id` (or `label`, when `sample-id` holds the run number) |
| cultivar | `cultivar` |
| health status | `health_status` or `status` |

QIIME 2 directive rows such as `#q2:types` are skipped. Tested with three different sample-sheet layouts, including the course demo's and Team 3's.

If step 1 saved files somewhere else, or used other column names, say so when submitting:

```bash
sbatch ... --export=ALL,FQ=/path/to/fastq,METADATA=/path/to/sheet.tsv,STATUS_COL=disease scripts/02_read_qc.sh
```

The log starts by printing which columns were used. If anything is missing (the sheet, a column, or a sample's FASTQ file), the job stops within seconds with a message saying what to fix.

## Files

```
02_read_qc/
├── README.md                   this file
├── FIGURES.md                  how to read each figure, with results
├── scripts/
│   ├── 02_read_qc.sh           Slurm job: builds envs, runs everything below
│   ├── 02a_prepare_samples.py  reads the step 1 sample sheet, checks FASTQs, writes samples.tsv
│   ├── 02_primer_check.py      counts reads starting with 27F / 1492R / neither
│   ├── 02_read_qc_plots.R      summary tables, Wilcoxon test, all figures
│   ├── 02d_insilico_pcr.py     predicts 27F-1492R product sizes in reference genomes
│   └── vt_colors.R             shared colour palette (sourced by the R script)
└── figures/                    PNG previews (view on GitHub); pdf/ has the originals
```

## How the scripts connect

```
sbatch scripts/02_read_qc.sh          (bash: Slurm job, runs each tool in order)
  ├─ python 02a_prepare_samples.py    -> results/02_qc/samples.tsv   (from metadata/metadata.tsv)
  ├─ fastqc + multiqc                 -> results/02_qc/fastqc/, multiqc/
  ├─ seqkit fx2tab (per read)         -> results/02_qc/read_tables/<run>.tsv
  ├─ seqkit stats                     -> results/02_qc/seqkit_stats.tsv
  ├─ python 02_primer_check.py        -> results/02_qc/primer_orientation.tsv
  └─ Rscript 02_read_qc_plots.R       reads the TSVs above
                                      -> per_sample_summary.tsv, qc_summary.txt, results/pdf/02_*.pdf
```

The languages pass data to each other only through TSV files, so the R plots can be re-run on their own without repeating the slow steps.

## Settings and other commands

The basic run is in [Quick start](#quick-start). Optional settings, passed with `--export=ALL,NAME=value`:

| Variable | Default | Meaning |
|---|---|---|
| `PROJ` | directory `sbatch` was run from | project directory |
| `FQ` | `$PROJ/raw_data/fastq` | folder of `.fastq.gz` files |
| `METADATA` | `$PROJ/metadata/metadata.tsv` | sample sheet from step 1 |
| `RUN_COL`, `SAMPLE_COL`, `CULTIVAR_COL`, `STATUS_COL` | detected | column names in the sample sheet, only if not detected |
| `ENVS` | `$PROJ/envs` | where conda environments are built |
| `CONDA_MODULE` | empty | module that provides conda/mamba, if it is not already on `PATH` |

Re-run only the plots (after a full run):

```bash
sbatch --account=<your_account> --partition=<partition> --output=logs/%x_%j.out \
       --wrap "envs/fb_r/bin/Rscript scripts/02_read_qc_plots.R $PWD"
```

In-silico PCR (run separately, takes seconds; genomes from NCBI: NC_018554.1, NC_061549.1, NZ_CP055227.2):

```bash
envs/qc/bin/python scripts/02d_insilico_pcr.py reference/*.fa
```

## Software

| Environment | Contents |
|---|---|
| `envs/qc` | FastQC, MultiQC, seqkit, Python (the Python scripts use the standard library only) |
| `envs/fb_r` | R 4.4: data.table, ggplot2, patchwork, scales (plus packages used in later steps) |

Versions used for the results below: FastQC 0.13.0, MultiQC 1.35, seqkit 2.14.0, Python 3.12, R 4.4.3.

## Key settings

- Primers (PacBio degenerate): 27F `AGRGTTYGATYMTGGCTCAG`, 1492R `RGYTACCTTGTTACGACTT`
- Primer check: exact IUPAC match within the **first 60 bp** of each read
- Long-read class: tallest peak between 1,700 and 2,100 bp, +/- 75 bp
- Amplicon window shown on plots: 1,000-1,600 bp (the DADA2 length filter used in step 3)

## Results

| Metric | Result |
|---|---|
| Reads per sample | 35,672 to 117,612 |
| Median per-read mean quality | Q37.4 (canker 38.6, healthy 34.2) |
| Median read length | 1,456 bp |
| Reads 1,000-1,600 bp | 82.0% |
| Long-read class (~1,819 bp) | 17.9% of reads; likely apple mitochondrial SSU rRNA (in-silico PCR predicts 1,819 bp) |
| Long class, canker vs healthy (per-sample median) | 19.1% vs 21.7%, Wilcoxon P = 0.14 |
| Orientation | 27F first 34.2%, 1492R first 47.5%, no primer found 18.3% |

**What this means for step 3:**
- Reads come in both orientations, so DADA2 `denoise-ccs` is needed (it finds primers on either strand and flips reads).
- The ~1,819 bp reads fall outside the 1,000-1,600 bp window and are removed by the length filter.
- Reads with no detectable primer are removed at the primer-trimming stage (step 3 dropped 17.1% of reads there).
- Chloroplast 16S (~1,451 bp) sits inside the bacterial size range, so it must be removed by taxonomy in step 4, not by length.

See [FIGURES.md](FIGURES.md) for each figure and its interpretation.

## Points for reviewers

- The primer check is deliberately strict (exact match, first 60 bp). It is a diagnostic, not a filter; DADA2 uses its own primer search with mismatch tolerance.
- The length plots exist in two versions: density (`02_read_length_distribution`, compares group shapes) and stacked counts (`02_read_length_counts`, comparable to the course demo figure).
- Colours and the cultivar order are set for this dataset in `vt_colors.R`. Other cultivars still work (they are added after the known ones), but the Wilcoxon test and summary expect the health status groups to be called `Canker` and `Healthy`.
- Testing: the original version of `02_read_qc.sh` produced the results above on ARC. This generic version was tested locally with three sample-sheet layouts and stand-in tools (real Python and R), not yet with a fresh conda install on a cluster; newer tool versions may differ slightly from those listed under Software.
- Large intermediate files (`read_tables/`, FastQC output, MultiQC HTML) and result tables are not uploaded; running the script recreates them.
