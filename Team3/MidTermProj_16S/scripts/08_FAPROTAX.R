# Step 8: Fig. 6 -- FAPROTAX functional prediction (microeco), canker vs healthy per cultivar
source(file.path(Sys.getenv("HOME"), "fireblight_canker_reanalysis/scripts/R_common.R"))
suppressPackageStartupMessages(library(microeco))
ps <- load_ps()
tt <- data.frame(tax_table(ps))
# SILVA gives no species for Erwinia; label ASVs matching the E. amylovora type-strain 16S (>=99.5%, step 7)
# as E. amylovora so FAPROTAX can assign "plant_pathogen" (paper identified Ea ASVs by BLAST)
ea <- read.csv(file.path(TAB, "Erwinia_ASVs_vs_Ea_16S.csv")); tt[ea$asv[ea$pct_id_Ea >= 99.5], "Species"] <- "Erwinia amylovora"
 for (r in colnames(tt)) tt[[r]] <- paste0(tolower(substr(r, 1, 1)), "__", ifelse(is.na(tt[[r]]), "", tt[[r]]))
sd <- data.frame(sample_data(ps)); sd$Sample <- paste0(sd$cultivar_code, sd$replicate, "_", sd$status)
me <- microtable$new(otu_table=data.frame(otu_table(ps), check.names=FALSE), tax_table=tt, sample_table=sd); me$tidy_dataset()
fn <- trans_func$new(me); fn$cal_spe_func(prok_database="FAPROTAX"); fn$cal_func_FR(abundance_weighted=TRUE)
res <- fn$res_func_FR; sd <- sd[rownames(res), ]   # microeco drops samples with 0 bacterial reads
res <- 100 * res[, colSums(res) > 0, drop=FALSE];   # fraction -> % write.csv(res, file.path(TAB, "Fig6_FAPROTAX_function_percent.csv"))
# order samples: Canker then Healthy, by cultivar
ord <- rownames(sd)[order(sd$status, sd$cultivar, sd$replicate)]
fn$res_func_FR <- res[ord, ] / 100
# Kruskal/Wilcoxon canker vs healthy per function
wt <- sapply(colnames(res), function(f) suppressWarnings(wilcox.test(res[sd$status == "Canker", f], res[sd$status == "Healthy", f])$p.value))
write.csv(data.frame(func=names(wt), p=wt, padj=p.adjust(wt, "BH"), mean_canker=colMeans(res[sd$status == "Canker", ]), mean_healthy=colMeans(res[sd$status == "Healthy", ])),
          file.path(TAB, "Fig6_FAPROTAX_wilcoxon.csv"), row.names=FALSE)
pdf(file.path(PDF, "Fig6_FAPROTAX_heatmap.pdf"), width=14, height=9)
print(fn$plot_func_FR(order_x=ord) + ggtitle("Fig. 6: FAPROTAX functional groups (% of abundance), samples ordered Canker | Healthy by cultivar"))
dev.off()
