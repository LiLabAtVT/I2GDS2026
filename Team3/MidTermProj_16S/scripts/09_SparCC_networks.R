# Step 9: Figs. 7-8 + Table 1 -- SparCC co-occurrence networks (Friedman & Alm 2012), MENAP-like settings:
#   |r| >= 0.2 and bootstrap pseudo-P < 0.01 (100 bootstraps). "Majority selection" from the paper:
#   canker T=5 (ASV present in >=5 canker samples), canker T=0 (ASV present in >half of samples),
#   healthy T=1 (ASV present in >=1 healthy sample, ASVs with >=10 reads).
#   Network topology (Table 1) with igraph; module roles (Zi/Pi) with greedy modularity.
source(file.path(Sys.getenv("HOME"), "fireblight_canker_reanalysis/scripts/R_common.R"))
suppressPackageStartupMessages({library(SpiecEasi); library(igraph)})
set.seed(2026)
ps <- load_ps(); NB <- as.integer(Sys.getenv("SPARCC_BOOT", "100"))
ea <- read.csv(file.path(TAB, "Erwinia_ASVs_vs_Ea_16S.csv")); ea_ids <- ea$asv[ea$pct_id_Ea >= 99.5]
TT <- as(tax_table(ps), "matrix")
lab <- function(x, ids) { g <- unname(TT[ids, "Genus"]); f <- unname(TT[ids, "Family"]); g[is.na(g)] <- f[is.na(g)]; g[is.na(g)] <- "unclassified"
  ifelse(ids %in% ea_ids, "Erwinia amylovora", g) }

run_net <- function(x, min_prev, name, min_reads=10) {
  ck <- file.path(P, "results", paste0("sparcc_", name, ".rds")); if (file.exists(ck)) { cat(name, ": loaded checkpoint\n"); return(readRDS(ck)) }
  m <- t(as(otu_table(x), "matrix")); keep <- colSums(m > 0) >= min_prev & colSums(m) >= min_reads; m <- m[, keep]
  cat(name, ": samples", nrow(m), " ASVs kept", ncol(m), "\n")
  sp <- sparcc(m); bt <- sparccboot(m, R=NB, ncpus=8); pv <- pval.sparccboot(bt)
  P <- diag(ncol(m)); P[upper.tri(P)] <- pv$pvals; P[lower.tri(P)] <- t(P)[lower.tri(P)]
  R <- sp$Cor; dimnames(R) <- dimnames(P) <- list(colnames(m), colnames(m))
  R <- (R + t(R)) / 2   # remove floating-point asymmetry
  A <- (abs(R) >= 0.2) & (P < 0.01); diag(A) <- FALSE
  g <- graph_from_adjacency_matrix(A * R, mode="undirected", weighted=TRUE, diag=FALSE); g <- delete_vertices(g, degree(g) == 0)
  V(g)$taxon <- lab(x, V(g)$name); E(g)$sign <- ifelse(E(g)$weight > 0, "positive", "negative")
  el <- as_data_frame(g); el$from_taxon <- lab(x, el$from); el$to_taxon <- lab(x, el$to)
  write.csv(el, file.path(TAB, paste0("network_edges_", name, ".csv")), row.names=FALSE)
  out <- list(g=g, R=R, P=P, name=name); saveRDS(out, ck); cat(name, ": done", format(Sys.time()), "\n"); out }

topology <- function(g) { ga <- g; E(ga)$weight <- abs(E(ga)$weight); d <- degree(ga)
  dd <- table(d); fit <- tryCatch(summary(lm(log(as.numeric(dd)) ~ log(as.numeric(names(dd)))))$r.squared, error=function(e) NA)
  bt <- betweenness(ga, weights=NA); ev <- eigen_centrality(ga, weights=NA)$vector; cl <- closeness(ga, weights=NA)
  sp <- distances(ga, weights=NA); fin <- sp[upper.tri(sp) & is.finite(sp)]
  nm <- function(v) paste(paste0(names(v)[v == max(v)], "-", V(g)$taxon[v == max(v)]), collapse="; ")
  data.frame(Nodes=vcount(g), Links=ecount(g), Positive=sum(E(g)$weight > 0), Negative=sum(E(g)$weight < 0),
    R2_power_law=fit, avgK=mean(d), avgCC=transitivity(ga, "average"), GD=mean(fin), E=mean(1 / fin), HD=1 / mean(1 / fin),
    MaxDegree=max(d), NodesMaxDegree=nm(d), CD=centr_degree(ga)$centralization, MaxBetweenness=max(bt), NodesMaxBetweenness=nm(bt),
    CB=centr_betw(ga)$centralization, NodesMaxEigen=nm(ev), CE=centr_eigen(ga)$centralization, NodesMaxCloseness=nm(cl),
    CCL=centr_clo(ga, mode="all")$centralization, Density=edge_density(ga), Transitivity=transitivity(ga, "global"),
    Connectedness=max(components(ga)$csize) / vcount(ga), Modules=length(unique(membership(cluster_fast_greedy(ga)))),
    Modularity=modularity(cluster_fast_greedy(ga)), check.names=FALSE) }

zipi <- function(g) { ga <- g; E(ga)$weight <- abs(E(ga)$weight); mod <- membership(cluster_fast_greedy(ga)); A <- as_adjacency_matrix(ga, sparse=FALSE) > 0
  k_in <- sapply(seq_along(mod), function(i) sum(A[i, mod == mod[i]])); z <- ave(k_in, mod, FUN=function(v) if (sd(v) > 0) (v - mean(v)) / sd(v) else 0)
  k <- rowSums(A); p <- 1 - sapply(seq_along(mod), function(i) sum((tapply(A[i, ], mod, sum) / k[i])^2))
  data.frame(ASV=V(ga)$name, taxon=V(ga)$taxon, module=as.integer(mod), Zi=z, Pi=p,
    role=ifelse(z > 2.5 & p > 0.62, "Network hub", ifelse(z > 2.5, "Module hub", ifelse(p > 0.62, "Connector", "Peripheral")))) }

plot_net <- function(nt, title) { g <- nt$g; tx <- V(g)$taxon; top <- names(sort(table(tx), decreasing=TRUE))[1:min(11, length(unique(tx)))]
  grp <- ifelse(tx %in% top, tx, "Other"); cols <- setNames(c(RColorBrewer::brewer.pal(max(3, length(top)), "Set3")[seq_along(top)], "grey85"), c(top, "Other"))
  V(g)$color <- cols[grp]; V(g)$color[tx == "Erwinia amylovora"] <- "red"
  E(g)$color <- ifelse(E(g)$weight > 0, adjustcolor("#2CA25F", 0.5), adjustcolor("#E34A33", 0.5))
  plot(g, layout=layout_with_fr(g, weights=abs(E(g)$weight)), vertex.size=3 + 1.2 * sqrt(degree(g)), vertex.label=if (vcount(g) <= 40) tx else NA, vertex.label.cex=0.6,
       vertex.label.color="black", edge.width=0.5 + 2 * abs(E(g)$weight), main=title)
  legend("bottomleft", legend=c(names(cols), "E. amylovora", "positive", "negative"), pch=c(rep(19, length(cols) + 1), NA, NA), lty=c(rep(NA, length(cols) + 1), 1, 1),
         col=c(cols, "red", "#2CA25F", "#E34A33"), cex=0.55, bty="n") }

ea_dot <- function(nt, title) { e <- read.csv(file.path(TAB, paste0("network_edges_", nt$name, ".csv")))
  e <- rbind(e, transform(e, from=to, to=from, from_taxon=to_taxon, to_taxon=from_taxon)); e <- e[e$from %in% ea_ids, ]
  if (!nrow(e)) return(NULL)
  e$ea <- paste("E. amylovora", e$from); e$partner <- paste(e$to_taxon, e$to)
  ggplot(e, aes(partner, ea, size=abs(weight), colour=ifelse(weight > 0, "Positive", "Negative"))) + geom_point(alpha=0.8) +
    scale_colour_manual(values=c(Positive="#2CA25F", Negative="#E34A33"), name="Interaction") + scale_size(name="|Correlation|", range=c(2, 8)) +
    labs(x=NULL, y=NULL, title=title) + theme(axis.text.x=element_text(angle=45, hjust=1, size=7)) }

can <- prune_samples(sample_data(ps)$status == "Canker", ps); hea <- prune_samples(sample_data(ps)$status == "Healthy", ps)
# Healthy bark: after rarefaction only ~350 bacterial reads remain in all 16 samples (bark is ~99.8% chloroplast),
# which yields no significant SparCC edges. Use the unrarefied bacterial counts for healthy samples (still very sparse).
hea_u <- readRDS(file.path(P, "results/phyloseq_unrarefied.rds")); hea_u <- prune_samples(sample_data(hea_u)$status == "Healthy", hea_u)
hea_u <- prune_samples(sample_sums(hea_u) > 0, hea_u); hea_u <- prune_taxa(taxa_sums(hea_u) > 0, hea_u)
nets <- list(run_net(can, 5, "canker_T5"), run_net(can, floor(nsamples(can) / 2) + 1, "canker_T0"), run_net(hea_u, 1, "healthy_T1_unrarefied", min_reads=2))
nets <- Filter(function(n) { if (vcount(n$g) == 0) { cat("WARNING:", n$name, "has no significant edges - skipped\n"); FALSE } else TRUE }, nets)
tab1 <- t(do.call(rbind, lapply(nets, function(n) topology(n$g)))); colnames(tab1) <- sapply(nets, `[[`, "name")
write.csv(tab1, file.path(TAB, "Table1_network_topology.csv")); print(tab1)
zp <- do.call(rbind, lapply(nets, function(n) cbind(network=n$name, zipi(n$g)))); write.csv(zp, file.path(TAB, "network_ZiPi_roles.csv"), row.names=FALSE)
print(subset(zp, role != "Peripheral"))

pdf(file.path(PDF, "Fig7_canker_SparCC_networks.pdf"), width=11, height=11)
plot_net(nets[[1]], "Fig. 7A: canker network (majority selection T = 5)"); p <- ea_dot(nets[[1]], "Fig. 7B: E. amylovora ASV correlations (T = 5)"); if (!is.null(p)) print(p)
plot_net(nets[[2]], "Fig. 7C: canker network (majority selection T = 0, >half of samples)"); p <- ea_dot(nets[[2]], "Fig. 7D: E. amylovora ASV correlations (T = 0)"); if (!is.null(p)) print(p)
dev.off()
if (length(nets) >= 3) {
pdf(file.path(PDF, "Fig8_healthy_SparCC_network.pdf"), width=11, height=11)
plot_net(nets[[3]], "Fig. 8A: healthy bark network (majority selection T = 1)")
g <- nets[[3]]$g; top40 <- names(sort(degree(g), decreasing=TRUE))[1:min(40, vcount(g))]; nt <- nets[[3]]; nt$g <- induced_subgraph(g, top40)
plot(nt$g, layout=layout_in_circle(nt$g), vertex.label=paste(V(nt$g)$name, V(nt$g)$taxon), vertex.label.cex=0.55, vertex.size=2 + degree(nt$g) / 3,
     vertex.color="grey20", edge.color=ifelse(E(nt$g)$weight > 0, adjustcolor("#2CA25F", 0.6), adjustcolor("#E34A33", 0.6)), main="Fig. 8B: ASV interaction network (top 40 ASVs by degree)")
dev.off() }
pdf(file.path(PDF, "Table1_network_topology.pdf"), width=11, height=10)
gridExtra::grid.table(cbind(Parameter=rownames(tab1), apply(tab1, 2, function(v) { n <- suppressWarnings(as.numeric(v)); ifelse(is.na(n), strtrim(v, 60), signif(n, 4)) })),
                      theme=gridExtra::ttheme_minimal(base_size=7))
dev.off()
