#!/bin/bash
# Step 1: fetch run list for BioProject PRJNA1348523 from ENA and build a QIIME2 metadata file
set -euo pipefail
P=$HOME/fireblight_canker_reanalysis
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1348523&result=read_run&fields=run_accession,sample_accession,sample_alias,library_name,instrument_model,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp&format=tsv" \
  > $P/metadata/ena_runinfo.tsv
# sample_alias looks like "Cripps Pink_1_Healthy_HiFi.fastq" -> cultivar / replicate / status
awk -F'\t' 'BEGIN{OFS="\t"; print "sample-id","run","cultivar","replicate","status","cultivar_code","group"}
NR>1{ split($3,a,"_"); cv=a[1]; if(cv ~ /Unknown/) cv="Unknown cultivar";
  code=(cv=="Cripps Pink")?"CP":(cv=="Aztec Fuji")?"Fuji":(cv=="Gala")?"Gala":(cv=="Hewe'"'"'s Crab")?"HW":(cv=="Granniwinkle")?"Nic":"UN";
  sid=code a[2] "_" a[3];
  print sid,$1,cv,a[2],a[3],code,a[3] "_" code }' $P/metadata/ena_runinfo.tsv \
  | (read h; echo "$h"; sort) > $P/metadata/sample_metadata.tsv
echo "samples: $(($(wc -l < $P/metadata/sample_metadata.tsv)-1))"
