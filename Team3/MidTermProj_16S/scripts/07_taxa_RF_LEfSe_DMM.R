# Step 7: Fig. 4A (top-20 genera heatmap, cankers), Fig. 4B (Random Forest, cultivar in cankers),
#         Fig. 5A (LEfSe canker vs healthy), Fig. 5B/D (Dirichlet-multinomial mixtures k=1..7, show k=4),
#         Fig. 5C (E. amylovora ASVs per canker sample)
source(file.path(Sys.getenv("HOME"), "fireblight_canker_reanalysis/scripts/R_common.R"))
suppressPackageStartupMessages({library(randomForest); library(DirichletMultinomial); library(microeco); library(pheatmap)})
set.seed(2026)
ps <- load_ps(); can <- prune_samples(sample_data(ps)$status == "Canker", ps); can <- prune_taxa(taxa_sums(can) > 0, can)
md <- data.frame(sample_data(can)); md$cultivar <- droplevels(md$cultivar)
rel_genus <- function(x) { g <- tax_glom(transform_sample_counts(x, function(v) v / sum(v)), "Genus", NArm=TRUE)
  m <- as(otu_table(g), "matrix"); rownames(m) <- tax_table(g)[, "Genus"]; rowsum(m, rownames(m)) }

## Fig. 4A: heatmap of top 20 genera in canker samples, log10 % relative abundance
gm <- rel_genus(can) * 100; top20 <- names(sort(rowMeans(gm), decreasing=TRUE))[1:20]
ord <- rownames(md)[order(md$cultivar, md$replicate)]
ann <- data.frame(Cultivar=md[ord, "cultivar"], row.names=ord)
write.csv(gm[top20, ord], file.path(TAB, "Fig4A_top20_genera_canker_pct.csv"))
pdf(file.path(PDF, "Fig4A_top20_genera_heatmap.pdf"), width=10, height=6)
pheatmap(log10(gm[top20, ord] + 0.01), cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=ann, annotation_colors=list(Cultivar=cv_cols[levels(md$cultivar)]),
  gaps_col=cumsum(table(md[ord, "cultivar"])), color=colorRampPalette(c("navy","#4575B4","#FFFFBF","#F46D43","#A50026"))(100),
  breaks=seq(-2, 2, length.out=101), legend_breaks=c(-2, -1, 0, 1, 2), legend_labels=c("0.01","0.1","1","10","100"),
  main="Fig. 4A: top 20 genera in cankers (%, log scale)", fontsize_col=6)
dev.off()

## Fig. 4B: Random Forest predicting cultivar from genus relative abundance (canker samples)
X <- t(rel_genus(can)); X <- X[, colSums(X > 0) >= 2, drop=FALSE]; colnames(X) <- make.names(colnames(X))
rf <- randomForest(x=X, y=md[rownames(X), "cultivar"], ntree=1000, importance=TRUE); print(rf)
capture.output(print(rf), file=file.path(TAB, "Fig4B_randomForest_summary.txt"))
imp <- data.frame(genus=colnames(X), Gini=importance(rf, type=2)[, 1]) %>% arrange(desc(Gini))
write.csv(imp, file.path(TAB, "Fig4B_RF_importance.csv"), row.names=FALSE)
top8 <- imp$genus[1:8]
ra <- data.frame(X[, top8], cultivar=md[rownames(X), "cultivar"]) %>% pivot_longer(-cultivar, names_to="genus") %>%
  group_by(genus, cultivar) %>% summarise(mean=mean(value), se=sd(value) / sqrt(n()), .groups="drop") %>% mutate(genus=factor(genus, rev(top8)))
p1 <- ggplot(imp[1:8, ], aes(reorder(genus, Gini), Gini)) + geom_col(fill="grey40") + coord_flip() + labs(x=NULL, y="Mean Decrease Gini")
p2 <- ggplot(ra, aes(genus, mean, fill=cultivar)) + geom_col(position=position_dodge(0.9)) + geom_errorbar(aes(ymin=mean, ymax=mean + se), position=position_dodge(0.9), width=0.3) +
  coord_flip() + scale_fill_manual(values=cv_cols) + labs(x=NULL, y="Relative abundance")
pdf(file.path(PDF, "Fig4B_RandomForest_cultivar.pdf"), width=11, height=5)
gridExtra::grid.arrange(p1, p2, ncol=2, widths=c(1, 1.4), top=sprintf("Fig. 4B: Random Forest (OOB error %.1f%%)", 100 * tail(rf$err.rate[, 1], 1)))
dev.off()

## Fig. 5A: LEfSe (canker vs healthy, genus level) via microeco
tt <- data.frame(tax_table(ps)); otu <- data.frame(otu_table(ps), check.names=FALSE); sd <- data.frame(sample_data(ps))
for (r in colnames(tt)) tt[[r]] <- paste0(tolower(substr(r, 1, 1)), "__", ifelse(is.na(tt[[r]]), "", tt[[r]]))
me <- microtable$new(otu_table=otu, tax_table=tt, sample_table=sd); me$tidy_dataset()
lf <- trans_diff$new(dataset=me, method="lefse", group="status", taxa_level="Genus", alpha=0.05, lefse_subgroup=NULL, p_adjust_method="none")
write.csv(lf$res_diff, file.path(TAB, "Fig5A_LEfSe_genus.csv"), row.names=FALSE)
pdf(file.path(PDF, "Fig5A_LEfSe_canker_vs_healthy.pdf"), width=9, height=9)
lr <- lf$res_diff %>% mutate(genus=sub(".*\\|", "", Taxa)) %>% mutate(genus=sub("^g__", "", genus)) %>% arrange(desc(LDA)) %>% slice_head(n=30)
print(ggplot(lr, aes(LDA, reorder(genus, LDA), colour=Group, size=LDA)) + geom_point() + scale_colour_manual(values=c(Canker="#E07B7B", Healthy="#4C9F70"), name="Enriched in") +
  scale_size(range=c(2, 6), name="LDA magnitude") + geom_vline(xintercept=c(-2, 2), linetype=2, colour="grey60") + xlim(-2, 6) +
  labs(x="LDA score (log10)", y="Bacterial genera", title="Fig. 5A: LEfSe, canker vs healthy (genus level)"))
dev.off()

## E. amylovora ASVs: Erwinia ASVs checked against the E. amylovora type-strain 16S (NR_041970)
ref <- Biostrings::readDNAStringSet(file.path(P, "refs/Ea_16S_NR_041970.fasta"))
erw <- taxa_names(can)[which(tax_table(can)[, "Genus"] == "Erwinia")]
idn <- sapply(erw, function(a) { s <- refseq(can)[[a]]
  pa <- pwalign::pairwiseAlignment(s, ref[[1]], type="global-local"); pwalign::pid(pa, "PID1") })
ea_tab <- data.frame(asv=erw, pct_id_Ea=round(idn, 2), reads=taxa_sums(can)[erw], silva_species=tax_table(can)[erw, "Species"]) %>% arrange(desc(reads))
write.csv(ea_tab, file.path(TAB, "Erwinia_ASVs_vs_Ea_16S.csv"), row.names=FALSE); print(head(ea_tab, 10))
ea <- ea_tab$asv[ea_tab$pct_id_Ea >= 99.5]; ea4 <- head(ea, 4)
rc <- transform_sample_counts(can, function(x) 100 * x / sum(x))
d5c <- psmelt(prune_taxa(ea4, rc)) %>% mutate(Sample=paste0(cultivar_code, replicate))
lvl <- d5c %>% group_by(Sample) %>% summarise(t=sum(Abundance)) %>% arrange(t) %>% pull(Sample)
p5c <- ggplot(d5c, aes(factor(Sample, lvl), Abundance, fill=OTU)) + geom_col() + labs(x=NULL, y="Relative abundance (%)", fill=NULL, title="Fig. 5C: Erwinia amylovora ASVs per canker sample") +
  theme(axis.text.x=element_text(angle=90, vjust=0.5))

## Fig. 5B/D: Dirichlet-multinomial mixtures on canker ASV counts (k = 1..7)
cnt <- t(as(otu_table(can), "matrix")); cnt <- cnt[, colSums(cnt) > 0]
fits <- lapply(1:7, function(k) dmn(cnt, k, verbose=FALSE, seed=2026))
lap <- sapply(fits, laplace); best <- which.min(lap); cat("Laplace:", round(lap, 1), " best k =", best, "\n")
k4 <- fits[[4]]; cl <- apply(mixture(k4), 1, which.max)
write.csv(data.frame(sample=names(cl), cluster=cl, md[names(cl), c("cultivar")]), file.path(TAB, "Fig5B_DMM_k4_clusters.csv"), row.names=FALSE)
# driver ASVs: difference of cluster-1 fitted from single-component fit (DirichletMultinomial vignette)
p0 <- fitted(fits[[1]], scale=TRUE); p4 <- fitted(k4, scale=TRUE); colnames(p4) <- paste0("m", 1:4)
diff <- rowSums(abs(p4 - as.vector(p0))); drv <- head(sort(diff, decreasing=TRUE), 10)
drv_df <- data.frame(ASV=names(drv), value=drv, taxon=apply(tax_table(can)[names(drv), c("Genus","Species")], 1, function(x) paste(na.omit(x), collapse=" ")))
write.csv(drv_df, file.path(TAB, "Fig5D_DMM_top10_drivers.csv"), row.names=FALSE)
pdf(file.path(PDF, "Fig5BCD_DMM_and_Ea_ASVs.pdf"), width=11, height=7)
plot(1:7, lap, type="b", xlab="Number of Dirichlet components (k)", ylab="Model fit (Laplace)", main="DMM model selection (paper presents k = 4)")
topA <- names(sort(colSums(cnt), decreasing=TRUE))[1:40]
so <- names(sort(cl)); a2 <- data.frame(Cluster=factor(cl[so]), Cultivar=md[so, "cultivar"], row.names=so)
pheatmap(t(sqrt(cnt[so, topA])), cluster_cols=FALSE, gaps_col=cumsum(table(cl[so])), annotation_col=a2, fontsize_row=6, fontsize_col=6,
  main=sprintf("Fig. 5B: DMM k=4 clusters (sizes %s), top 40 ASVs (sqrt counts)", paste(table(cl), collapse="/")))
print(p5c)
print(ggplot(drv_df, aes(reorder(ASV, value), value)) + geom_col(fill="grey30") + geom_text(aes(y=max(value) * 1.02, label=taxon), hjust=0, size=3) +
  coord_flip(ylim=c(0, max(drv_df$value) * 1.9)) + labs(x="ASV", y="Value", title="Fig. 5D: top 10 drivers, DMM community types"))
dev.off()
