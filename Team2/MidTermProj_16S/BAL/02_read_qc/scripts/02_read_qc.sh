#!/bin/bash
# The #SBATCH lines below are instructions to Slurm, the job scheduler on a
# computing cluster. They reserve resources for this job on a compute node.
# (To bash they are just comments; only sbatch reads them.)
#SBATCH --job-name=02_read_qc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=logs/%x_%j.out
# Above: 1 machine, 8 CPU cores, 32 GB memory, up to 1.5 hours.
# Everything the job prints goes to logs/02_read_qc_<job number>.out.

# =============================================================================
# Step 2: read quality control (QC) for PacBio HiFi full-length 16S reads
# =============================================================================
# Before turning reads into ASVs (step 3, DADA2), this checks the raw reads are
# what we expect: ~1,450-1,500 bp long, good quality, enough reads per sample,
# and carrying the PCR primers. No reads are removed; it only measures and plots.
#
# Steps (numbered the same in the code below):
#   1. clean sample sheet (02a_prepare_samples.py)   4. per-read length + quality (seqkit)
#   2. software setup, first run only (conda)         5. primer orientation (02_primer_check.py)
#   3. FastQC + MultiQC reports                       6. tables + figures (02_read_qc_plots.R)
#
# Inputs from step 1: raw_data/fastq/<run>.fastq.gz and metadata/metadata.tsv
#
# Usage, from the project folder (I2GDS 2026 class settings on VT ARC):
#   mkdir -p logs
#   sbatch --account=introtogds --partition=normal_q --qos=tc_normal_short \
#          --export=ALL,CONDA_MODULE=Miniforge3/24.11.3-0 scripts/02_read_qc.sh
# Elsewhere, use your own --account/--partition. Other settings (FQ, METADATA,
# ENVS, column names) are listed in the README.

# Safety switches: stop at the first error (-e), treat a misspelled/unset
# variable as an error (-u), and fail if any command in a pipe "|" fails.
# Without this, the script would carry on and make plots from broken data.
set -euo pipefail

# --- settings ----------------------------------------------------------------
# ${NAME:-default} means "use NAME if the user set it, otherwise the default".
# SLURM_SUBMIT_DIR is the folder you were in when you typed sbatch.
PROJ=${PROJ:-${SLURM_SUBMIT_DIR:-$PWD}}
FQ=${FQ:-$PROJ/raw_data/fastq}
METADATA=${METADATA:-$PROJ/metadata/metadata.tsv}
ENVS=${ENVS:-$PROJ/envs}
CONDA_MODULE=${CONDA_MODULE:-}
OUT=$PROJ/results/02_qc
# Use as many CPU cores as Slurm gave us (8 by default)
THREADS=${SLURM_CPUS_PER_TASK:-8}
QC_ENV=$ENVS/qc
R_ENV=$ENVS/fb_r
# The clean sample sheet written by 02a_prepare_samples.py; everything below reads this
SAMPLES=$OUT/samples.tsv

# --- check the inputs exist before doing any work ------------------------------
# If the sheet is missing, print a clear message and stop (exit 1 = failure).
[[ -f $METADATA ]] || { echo "missing sample sheet $METADATA (set METADATA=...)" >&2; exit 1; }
# Make the output folders (-p: no error if they already exist)
mkdir -p "$OUT"/{fastqc,multiqc,read_tables} "$PROJ/results/pdf"

# Record where and when the job ran, to help with troubleshooting later
echo "host: $(hostname)  job: ${SLURM_JOB_ID:-none}  start: $(date)"
echo "project: $PROJ  fastq: $FQ  metadata: $METADATA  envs: $ENVS"

# --- 1. clean sample sheet ------------------------------------------------------
# Done before building software so that a wrong path or column name fails in
# seconds, not after the environments are built. Uses the system python3 if
# there is one, otherwise waits for the qc environment below.
rm -f "$SAMPLES"
prepare_samples() { "$1" "$PROJ/scripts/02a_prepare_samples.py" "$METADATA" "$FQ" "$SAMPLES"; }
if command -v python3 > /dev/null; then prepare_samples python3; fi

# --- 2. environments --------------------------------------------------------------
# A conda "environment" is a self-contained folder of software at fixed
# versions. We build two: "qc" (QC tools + Python) and "fb_r" (R + packages).
# They are only built if missing, so later runs skip straight to the analysis.

# On clusters, conda is often provided as a "module" that has to be loaded first
[[ -n $CONDA_MODULE ]] && module load "$CONDA_MODULE"
# Use mamba (a faster conda) if available, otherwise plain conda
MAMBA=$(command -v mamba || command -v conda) || { echo "mamba/conda not found; set CONDA_MODULE" >&2; exit 1; }

# -p puts the environment inside the project folder (not your home folder);
# -c lists the software channels to download from (bioconda has bio tools).
if [[ ! -x $QC_ENV/bin/multiqc ]]; then
    # No qc environment yet: create it with everything this step needs
    "$MAMBA" create -y -p "$QC_ENV" -c conda-forge -c bioconda fastqc multiqc seqkit python
elif [[ ! -x $QC_ENV/bin/seqkit ]]; then
    # An older qc environment without seqkit: add the missing tools to it
    "$MAMBA" install -y -p "$QC_ENV" -c conda-forge -c bioconda seqkit python
fi
# fb_r holds R for all later plotting/statistics steps as well;
# step 2 itself only needs data.table, ggplot2, patchwork and scales.
if [[ ! -x $R_ENV/bin/Rscript ]]; then
    "$MAMBA" create -y -p "$R_ENV" -c conda-forge -c bioconda \
        r-base=4.4 r-data.table r-ggplot2 r-patchwork r-scales r-dplyr r-tidyr \
        r-ggrepel r-vegan r-agricolae r-randomforest r-igraph \
        bioconductor-phyloseq bioconductor-dirichletmultinomial
fi
# Put the qc tools first on PATH, so "fastqc", "seqkit", "python" below
# mean the versions inside the qc environment
export PATH=$QC_ENV/bin:$PATH
[[ -f $SAMPLES ]] || prepare_samples python
# The FASTQ files to analyse: one per sample in the clean sheet (column 2 = run)
mapfile -t FILES < <(awk -F'\t' -v d="$FQ" 'NR > 1 { print d "/" $2 ".fastq.gz" }' "$SAMPLES")
echo "FASTQ files: ${#FILES[@]}"
# Print the version of every tool into the log, so results can be reproduced
echo "fastqc: $(fastqc --version)"
echo "multiqc: $(multiqc --version)"
echo "seqkit: $(seqkit version)"
echo "python: $(python --version)"
echo "R: $("$R_ENV/bin/Rscript" --version 2>&1)"

# --- 3. FastQC + MultiQC ----------------------------------------------------------
# FastQC writes a quality report (an HTML page) for each FASTQ file.
# MultiQC then combines all the reports into one page for easy comparison.
fastqc -t "$THREADS" --memory 2048 -q --outdir "$OUT/fastqc" "${FILES[@]}"
multiqc -f -q -o "$OUT/multiqc" "$OUT/fastqc"

# --- 4. per-read length and mean quality --------------------------------------------
# For every sample, write a 2-column table: length and average quality of each read.
# seqkit fx2tab -n -i -l -q prints read ID, length and mean quality; cut keeps
# columns 2-3. xargs -P runs several samples at the same time (one per core);
# inside sh -c, $1 is the FASTQ file and $2 the output folder.
printf '%s\n' "${FILES[@]}" | xargs -P "$THREADS" -I{} sh -c \
    'r=$(basename "$1" .fastq.gz); seqkit fx2tab -n -i -l -q "$1" | cut -f2,3 > "$2/read_tables/$r.tsv"' _ {} "$OUT"
# One-line summary per sample: number of reads, min/mean/max length, etc.
seqkit stats -j "$THREADS" -T -a "${FILES[@]}" > "$OUT/seqkit_stats.tsv"

# --- 5. primer orientation ----------------------------------------------------------
# Count, per sample, how many reads start with the forward primer (27F), the
# reverse primer (1492R), or neither. PacBio reads either strand, so both are expected.
python "$PROJ/scripts/02_primer_check.py" "$SAMPLES" "$FQ" "$OUT/primer_orientation.tsv"

# --- 6. plots + summary -------------------------------------------------------------
# R reads the tables written above and makes the summary numbers and figures.
"$R_ENV/bin/Rscript" "$PROJ/scripts/02_read_qc_plots.R" "$PROJ"

echo "done: $(date)"
