#!/bin/bash
# Helper: show my Slurm jobs and the tail of their logs
cd /projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code
squeue -u "$USER" -o '%.10i %.14j %.9T %.10M %R'
sacct -u "$USER" -S today -X -o JobID,JobName%16,State,Elapsed,ExitCode | tail -8
for j in "$@"; do echo "== log $j"; tail -n 15 logs/*_"$j".out 2>/dev/null; done
