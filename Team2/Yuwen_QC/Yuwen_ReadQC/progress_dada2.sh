#!/bin/bash
# Helper: count files in the DADA2 temp directories of a running job (run on its compute node)
# usage: bash progress_dada2.sh JOBID
J=$1
for d in $(find /localscratch/"$J" /tmp -maxdepth 4 -type d -user "$USER" 2>/dev/null); do
  printf "%5d files  %s\n" "$(ls "$d" 2>/dev/null | wc -l)" "$d"
done | sort -k3 | head -30
