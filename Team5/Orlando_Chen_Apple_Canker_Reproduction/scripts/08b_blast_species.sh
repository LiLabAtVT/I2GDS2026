#!/bin/bash
#SBATCH --job-name=blast16S
#SBATCH --account=introtogds
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code/logs/08b_blast_species_%j.log
# Step 8b: BLAST bacterial ASVs against NCBI 16S rRNA RefSeq (the paper identified its Erwinia / DMM-driver ASVs
# with NCBI BLAST, Fig. 5D), assign species only where unambiguous, then rerun FAPROTAX with species (Step 8).
set -euo pipefail
B=/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code
OUT=$B/results/08b_blast; DB=$B/reference/ncbi_16S_ribosomal_RNA
mkdir -p $OUT $DB
module load apptainer/1.4.0
SIF=/apps/common/software/QIIME2/2024.10/qiime2-amplicon-2024.10.sif
export TMPDIR=${TMPDIR:-/tmp}
Q() { apptainer exec --bind /projects --bind $TMPDIR $SIF "$@"; }
Q blastn -version > $OUT/blast_version.txt; cat $OUT/blast_version.txt

# 1. NCBI 16S rRNA RefSeq database (Bacteria + Archaea type material)
if [ ! -s $DB/16S_ribosomal_RNA.nsq ] && [ ! -s $DB/16S_ribosomal_RNA.00.nsq ]; then
  wget -q -O $DB/16S_ribosomal_RNA.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
  wget -q -O $DB/16S_ribosomal_RNA.tar.gz.md5 https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz.md5
  (cd $DB && md5sum -c 16S_ribosomal_RNA.tar.gz.md5 && tar -xzf 16S_ribosomal_RNA.tar.gz && rm 16S_ribosomal_RNA.tar.gz)
fi
Q blastdbcmd -db $DB/16S_ribosomal_RNA -info > $OUT/blastdb_info.txt; head -6 $OUT/blastdb_info.txt

# 2. BLAST all ASV sequences (bacterial ones are selected in the Python step)
Q blastn -query $B/results/04_dada2/export/dna-sequences.fasta -db $DB/16S_ribosomal_RNA -task megablast \
  -outfmt "6 qseqid sacc pident length qlen slen mismatch gapopen evalue bitscore stitle" \
  -max_target_seqs 50 -perc_identity 90 -num_threads $SLURM_CPUS_PER_TASK -out $OUT/blast_all_asvs.tsv
wc -l $OUT/blast_all_asvs.tsv

# 3. Species calls + Erwinia / DMM-driver report
Q python $B/scripts/08b_assign_species.py

# 4. FAPROTAX again with species added
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate canker_r
Rscript $B/scripts/08_faprotax.R $OUT/species_assignments.tsv 08b_faprotax_species
echo "Step 8b done"
