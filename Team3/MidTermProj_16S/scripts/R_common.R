# Shared settings for all R steps
suppressPackageStartupMessages({library(phyloseq); library(ggplot2); library(dplyr); library(tidyr); library(reshape2); library(magrittr)})
P   <- file.path(Sys.getenv("HOME"), "fireblight_canker_reanalysis")
PDF <- file.path(P, "results/pdf"); TAB <- file.path(P, "results/tables")
dir.create(PDF, recursive=TRUE, showWarnings=FALSE); dir.create(TAB, recursive=TRUE, showWarnings=FALSE)
cv_levels <- c("Cripps Pink","Aztec Fuji","Gala","Hewe's Crab","Granniwinkle","Unknown cultivar")
cv_cols   <- c("Cripps Pink"="#B5A000","Aztec Fuji"="#F8766D","Gala"="#00BA38","Hewe's Crab"="#619CFF","Granniwinkle"="#00BFC4","Unknown cultivar"="#F564E3")
theme_set(theme_bw(base_size=10))
load_ps <- function() readRDS(file.path(P, "results/phyloseq.rds"))
# strip SILVA rank prefixes, fill unassigned with higher-rank label
clean_tax <- function(tt) { tt <- as.matrix(tt); tt <- gsub("^[a-z]__", "", tt); tt[is.na(tt) | tt %in% c("", "uncultured", "metagenome")] <- NA; tt }
