#!/bin/bash
# Step 2: download raw PacBio HiFi FASTQs (PRJNA1348523) from ENA and verify MD5 checksums
set -euo pipefail
P=$HOME/fireblight_canker_reanalysis
cd $P/data/raw
tail -n +2 $P/metadata/ena_runinfo.tsv | awk -F'\t' '{print $1"\t"$9"\t"$10}' > urls.tsv
dl() { run=$1; md5=$2; url=$3
  f=${run}.fastq.gz
  if [[ -s $f ]] && [[ $(md5sum $f | cut -d' ' -f1) == $md5 ]]; then echo "$run ok (cached)"; return; fi
  curl -s --retry 5 -o $f "https://$url"
  if [[ $(md5sum $f | cut -d' ' -f1) == $md5 ]]; then echo "$run ok"; return; fi
  # ENA copy of SRR36091178 is corrupt at source -> fall back to NCBI SRA (sratoolkit in refs/)
  echo "$run ENA MD5 FAIL -> NCBI fasterq-dump"
  SRA=$(ls -d $P/refs/sratoolkit.*/bin 2>/dev/null || true)
  if [[ -z $SRA ]]; then (cd $P/refs && curl -sL https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz | tar xz); SRA=$(ls -d $P/refs/sratoolkit.*/bin); fi
  T=$(mktemp -d /scratch/$USER/sra_XXXX); (cd $T && $SRA/prefetch -q $run && $SRA/fasterq-dump -q $run -O .)
  gzip -c $T/$run.fastq > $f && rm -rf $T && gzip -t $f && echo "$run ok (NCBI)"; }
export -f dl; export P
xargs -P 4 -L 1 bash -c 'dl "$@"' _ < urls.tsv
rm urls.tsv
