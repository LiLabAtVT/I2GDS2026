# Step 5: assemble QIIME2 outputs into a phyloseq object; QC plots
#   - DADA2 read tracking, ASVs per sample (Supplementary Fig. S1), rarefaction curves
source(file.path(Sys.getenv("HOME"), "fireblight_canker_reanalysis/scripts/R_common.R"))
E <- file.path(P, "qiime/export")
otu <- read.delim(file.path(E, "feature-table.tsv"), skip=1, row.names=1, check.names=FALSE)
tax <- read.delim(file.path(E, "taxonomy/taxonomy.tsv"), row.names=1)
ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tt <- t(sapply(strsplit(tax$Taxon, ";\\s*"), function(x) { length(x) <- 7; x })); dimnames(tt) <- list(rownames(tax), ranks)
tt <- clean_tax(tt)
meta <- read.delim(file.path(P, "metadata/sample_metadata.tsv"), row.names=1, check.names=FALSE)
meta$cultivar <- factor(meta$cultivar, cv_levels); meta$status <- factor(meta$status, c("Canker","Healthy"))
tree <- ape::read.tree(file.path(E, "rooted-tree/tree.nwk"))
ps <- phyloseq(otu_table(as.matrix(otu), taxa_are_rows=TRUE), tax_table(tt), sample_data(meta), phy_tree(tree),
               refseq(Biostrings::readDNAStringSet(file.path(E, "rep-seqs/dna-sequences.fasta"))))
# Rarefy the FULL table (incl. chloroplast) to the minimum depth, THEN remove host organelles.
# Healthy bark is ~99.8% chloroplast; rarefying after organelle removal would leave <10 reads/sample.
# This order reproduces the paper's numbers (canker 16-470 ASVs, mean ~130; healthy 1-22, mean ~5).
depth <- min(sample_sums(ps)); cat("Rarefying full table to", depth, "reads/sample\n")
psr <- rarefy_even_depth(ps, sample.size=depth, rngseed=2026, replace=FALSE, trimOTUs=FALSE, verbose=FALSE)
filt <- function(x) { x <- subset_taxa(x, Kingdom == "Bacteria" & (is.na(Order) | Order != "Chloroplast") & (is.na(Family) | Family != "Mitochondria"))
  prune_taxa(taxa_sums(x) > 0, x) }
ps <- filt(ps); psr <- filt(psr)
# readable ASV labels (ASV0001.. by total abundance)
o <- order(taxa_sums(ps), decreasing=TRUE); ids <- sprintf("ASV%04d", seq_len(ntaxa(ps)))
map <- data.frame(feature_id=taxa_names(ps)[o], asv=ids); write.csv(map, file.path(TAB, "asv_id_map.csv"), row.names=FALSE)
taxa_names(ps) <- ids[match(taxa_names(ps), map$feature_id)]; taxa_names(psr) <- ids[match(taxa_names(psr), map$feature_id)]
saveRDS(ps, file.path(P, "results/phyloseq_unrarefied.rds")); saveRDS(psr, file.path(P, "results/phyloseq.rds"))
cat("Rarefied+filtered: reads", sum(sample_sums(psr)), " ASVs", ntaxa(psr), "\n")
ps_unrar <- ps; ps <- psr   # downstream QC summaries use the rarefied table (as in the paper)
cat("Reads:", sum(sample_sums(ps)), " ASVs:", ntaxa(ps), "\n")

st <- read.delim(file.path(E, "dada2-stats/stats.tsv"), skip=0, comment.char="#", check.names=FALSE)
st <- st[st[[1]] != "#q2:types", ]; write.csv(st, file.path(TAB, "dada2_read_tracking.csv"), row.names=FALSE)
nasv <- data.frame(sample=sample_names(ps), ASVs=colSums(otu_table(ps) > 0), reads=sample_sums(ps), meta[sample_names(ps), ])
write.csv(nasv, file.path(TAB, "asvs_per_sample.csv"), row.names=FALSE)
print(nasv %>% group_by(status) %>% summarise(min=min(ASVs), mean=round(mean(ASVs)), max=max(ASVs), reads=sum(reads)))

pdf(file.path(PDF, "Step05_QC_reads_and_ASVs.pdf"), width=10, height=6)
org <- data.frame(sample=sample_names(ps_unrar), bacterial=sample_sums(ps_unrar))
sl <- st %>% select(1, input, filtered, denoised, `non-chimeric`) %>% pivot_longer(-1, names_to="step", values_to="reads") %>%
  mutate(step=factor(step, c("input","filtered","denoised","non-chimeric")), reads=as.numeric(reads))
print(ggplot(sl, aes(.data[[names(sl)[1]]], reads, fill=step)) + geom_col(position="dodge") + coord_flip() + labs(x=NULL, title="DADA2 read tracking per sample"))
print(ggplot(nasv, aes(reorder(sample, ASVs), ASVs, fill=cultivar)) + geom_col() + facet_grid(~status, scales="free_x", space="free") +
  scale_fill_manual(values=cv_cols) + theme(axis.text.x=element_text(angle=90, hjust=1, size=6)) + labs(x=NULL, title="Number of ASVs per sample (cf. Supplementary Fig. S1)"))
m <- as(t(otu_table(ps)), "matrix"); vegan::rarecurve(m, step=500, label=FALSE, col=cv_cols[as.character(meta[rownames(m),"cultivar"])], main="Rarefaction curves")
dev.off()
