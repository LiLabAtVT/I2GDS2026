#!/bin/bash
# Syntax-check all R scripts (parse only, nothing is executed)
cd /projects/intro2gds/I2GDS2026/Team2/Yuwen/Claude_Code
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2 >/dev/null 2>&1
for f in scripts/*.R; do
  Rscript -e "invisible(parse('$f')); cat('OK  $f\n')" 2>&1 | tail -3
done
Rscript -e 'library(microeco); tf <- trans_func$new; cat("func_group_list in trans_func:", "func_group_list" %in% names(trans_func$public_fields) || grepl("func_group_list", paste(deparse(trans_func$public_methods$initialize), collapse="")), "\n")'
