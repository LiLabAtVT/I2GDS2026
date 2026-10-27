library(dada2)
fs <- Sys.glob("/projects/intro2gds/I2GDS2026/Team1/Kamal/Claude_Code/raw_data/*.fastq.gz")
out <- file.path(Sys.getenv("TMPDIR"), basename(fs))
for (i in seq_along(fs)) {
  r <- tryCatch(removePrimers(fs[i], out[i], "AGRGTTYGATYMTGGCTCAG", "AAGTCGTAACAAGGTARCY", max.mismatch = 0, orient = TRUE, verbose = FALSE), error = function(e) conditionMessage(e))
  cat(basename(fs[i]), if (is.character(r)) paste("ERROR", r) else paste(r[1], r[2]), "\n")
}
