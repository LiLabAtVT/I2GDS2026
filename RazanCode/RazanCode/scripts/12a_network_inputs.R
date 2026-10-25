# Step 12a: FastSpar input tables for the three SparCC networks (paper Figs. 7-8, Table 1)
#   paper (MENAP/iNAP): "majority selection" 5 and 0 for canker, 1 for healthy; SparCC |r| threshold 0.2; P < 0.01.
#   Interpretation used here: T = k > 0 -> keep ASVs present in >= k samples; T = 0 -> keep ASVs present in more than
#   half of the samples (the paper's definition of 0).
#   Table (user decision): relaxed = bacteria only; cankers = 32 samples, healthy = 8 samples with >= 1 bacterial read.
# Run from the project root directory (env fb_r):
#   R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_LIBS_USER=$CONDA_PREFIX/lib/R/library Rscript scripts/12a_network_inputs.R
meta <- read.delim("metadata/metadata.tsv", comment.char = "#"); names(meta)[1] <- "sample"
map  <- read.delim("results/10_erwinia_dmm/asv_id_map.tsv")
tab  <- as.matrix(read.delim("results/04_taxonomy/export/feature-table.tsv", skip = 1, check.names = FALSE, row.names = 1))
rownames(tab) <- map$asv_num[match(rownames(tab), map$feature_id)]

nets <- list(canker_T5 = list(status = "Canker", T = 5),
             canker_T0 = list(status = "Canker", T = 0),
             healthy_T1 = list(status = "Healthy", T = 1))
for (nm in names(nets)) {
  s <- intersect(meta$sample[meta$health_status == nets[[nm]]$status], colnames(tab))
  x <- tab[, s, drop = FALSE]; x <- x[, colSums(x) > 0, drop = FALSE]
  prev <- rowSums(x > 0)
  keep <- if (nets[[nm]]$T == 0) prev > ncol(x) / 2 else prev >= nets[[nm]]$T
  x <- x[keep, , drop = FALSE]
  d <- file.path("results/12_networks", nm); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  out <- data.frame(`#OTU ID` = rownames(x), x, check.names = FALSE)
  write.table(out, file.path(d, "input.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%-11s samples %2d  ASVs kept %4d\n", nm, ncol(x), nrow(x)))
}
