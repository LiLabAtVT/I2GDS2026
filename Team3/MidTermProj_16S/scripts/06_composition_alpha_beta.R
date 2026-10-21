# Step 6: Fig. 2A (class-level composition), Fig. 2B (alpha diversity + ANOVA/Tukey letters),
#         Supplementary Fig. S2 (ACE, Jackknife, Faith PD), Fig. 3 (Bray-Curtis PCoA + PERMANOVA)
source(file.path(Sys.getenv("HOME"), "fireblight_canker_reanalysis/scripts/R_common.R"))
suppressPackageStartupMessages({library(vegan); library(agricolae)})
set.seed(2026)
ps <- load_ps(); md <- data.frame(sample_data(ps))

## Fig. 2A: relative abundance at class level, samples grouped by status/cultivar
rel <- transform_sample_counts(ps, function(x) x / sum(x))
cl <- psmelt(rel) %>% mutate(Class=ifelse(is.na(Class), "Unassigned", Class)) %>% group_by(Sample, status, cultivar, replicate, Class) %>% summarise(Abundance=sum(Abundance), .groups="drop") %>% filter(!is.na(Abundance))
top <- cl %>% group_by(Class) %>% summarise(a=sum(Abundance)) %>% arrange(desc(a)) %>% slice_head(n=21) %>% pull(Class)
cl$Class <- factor(ifelse(cl$Class %in% top, cl$Class, "Other"), c(rev(top), "Other"))
cl$Sample <- factor(cl$Sample, md %>% arrange(desc(status == "Healthy"), cultivar, replicate) %>% rownames())
pal <- setNames(c(colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(top)), "grey80"), c(top, "Other"))
p2a <- ggplot(cl, aes(Sample, Abundance * 100, fill=Class)) + geom_col(width=0.9) + facet_grid(~status + cultivar, scales="free_x", space="free_x") +
  scale_fill_manual(values=pal) + labs(y="Relative abundance (%) (Class level)", x=NULL, title="Fig. 2A") +
  theme(axis.text.x=element_text(angle=90, size=5), strip.text=element_text(size=6), legend.text=element_text(size=6), legend.key.size=unit(3, "mm"))
write.csv(cl %>% group_by(status, Class) %>% summarise(mean_rel=mean(Abundance) * 100, .groups="drop"), file.path(TAB, "Fig2A_class_relabund_by_status.csv"), row.names=FALSE)

## Alpha diversity on rarefied table (paper: rarefied to equal depth)
psr <- ps   # already rarefied (before organelle removal) in step 5
al <- estimate_richness(psr, measures=c("Observed","Chao1","ACE","Shannon","Simpson"))
al$Jackknife <- estimateR_jack <- vegan::specpool(as(t(otu_table(psr)), "matrix"), pool=rownames(al))$jack1
al$FaithPD <- picante_free_pd <- sapply(sample_names(psr), function(s) { keep <- taxa_names(psr)[otu_table(psr)[, s] > 0]
  if (length(keep) < 2) return(NA); sum(ape::keep.tip(phy_tree(psr), keep)$edge.length) })
al <- cbind(md[rownames(al), c("cultivar","status")], al); al$sample <- rownames(al)
al$grp <- interaction(al$status, al$cultivar, sep=": ")
write.csv(al, file.path(TAB, "alpha_diversity.csv"), row.names=FALSE)

# ANOVA (status vs metric) and Tukey HSD letters across status x cultivar groups
stats <- lapply(c("Simpson","Shannon","Observed","Chao1","ACE","Jackknife","FaithPD"), function(m) {
  a <- summary(aov(as.formula(paste(m, "~ status")), al))[[1]]
  data.frame(metric=m, F=a$`F value`[1], P=a$`Pr(>F)`[1]) })
stats <- do.call(rbind, stats); print(stats); write.csv(stats, file.path(TAB, "alpha_anova_status.csv"), row.names=FALSE)
letters_for <- function(m) { h <- HSD.test(aov(as.formula(paste(m, "~ grp")), al), "grp", group=TRUE)$groups
  data.frame(grp=rownames(h), letter=h$groups) }
alpha_plot <- function(m) {
  L <- letters_for(m); y <- al %>% group_by(grp) %>% summarise(y=max(.data[[m]], na.rm=TRUE)) %>% left_join(L, by="grp")
  y <- y %>% separate(grp, c("status","cultivar"), sep=": ", remove=FALSE)
  ggplot(al, aes(cultivar, .data[[m]], colour=cultivar)) + geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.15, size=1) +
    geom_text(data=y, aes(cultivar, y * 1.05, label=letter), colour="black", size=3) + facet_grid(~status) +
    scale_colour_manual(values=cv_cols) + labs(x=NULL, y=m) + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") }

pdf(file.path(PDF, "Fig2_composition_alpha_diversity.pdf"), width=12, height=7)
print(p2a)
print(patchwork_free <- gridExtra::grid.arrange(grobs=lapply(c("Simpson","Shannon","Observed","Chao1"), function(m) ggplotGrob(alpha_plot(m))), ncol=2, top="Fig. 2B: alpha diversity (Tukey HSD letters, P<0.05)"))
dev.off()
pdf(file.path(PDF, "FigS2_alpha_ACE_Jackknife_FaithPD.pdf"), width=12, height=4)
gridExtra::grid.arrange(grobs=lapply(c("ACE","Jackknife","FaithPD"), function(m) ggplotGrob(alpha_plot(m))), ncol=3)
dev.off()

## Fig. 3: Bray-Curtis PCoA + PERMANOVA (status * cultivar) and pairwise status test
pb <- prune_samples(sample_sums(psr) > 0, psr)   # Bray-Curtis undefined for samples with 0 bacterial reads
cat("Samples in PCoA:", nsamples(pb), " dropped (0 bacterial reads):", paste(setdiff(sample_names(psr), sample_names(pb)), collapse=", "), "\n")
relr <- transform_sample_counts(pb, function(x) x / sum(x))
bc <- phyloseq::distance(relr, "bray"); ord <- ordinate(relr, "PCoA", distance=bc)
mdr <- data.frame(sample_data(relr))
perm <- adonis2(bc ~ status * cultivar, data=mdr, permutations=999, by="terms"); print(perm)
capture.output(perm, file=file.path(TAB, "Fig3_PERMANOVA.txt"))
ev <- round(100 * ord$values$Relative_eig[1:2], 1)
p3 <- plot_ordination(relr, ord, color="cultivar", shape="status") + geom_point(size=2.5) + stat_ellipse(aes(group=status), linetype=2, colour="black") +
  scale_colour_manual(values=cv_cols) + labs(x=paste0("PCoA1 [", ev[1], "%]"), y=paste0("PCoA2 [", ev[2], "%]"),
  title="Fig. 3: PCoA (Bray-Curtis)", subtitle=sprintf("PERMANOVA status R2=%.3f P=%.3f; cultivar P=%.3f; status:cultivar R2=%.3f P=%.3f",
  perm$R2[1], perm$`Pr(>F)`[1], perm$`Pr(>F)`[2], perm$R2[3], perm$`Pr(>F)`[3]))
pdf(file.path(PDF, "Fig3_PCoA_BrayCurtis.pdf"), width=7.5, height=6); print(p3); dev.off()
