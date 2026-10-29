# Step 6g: Figs. 7-8 / Table 1 / Supplementary Fig. S3 - SparCC co-occurrence networks
# The paper ran SparCC in the MENAP web pipeline and drew networks in Cytoscape/Gephi.
# Here the same logic is scripted in R: SpiecEasi::sparcc + bootstrap pseudo P-values,
# |r| >= 0.2 and P < 0.01 for all networks (as in the paper's Methods), network properties and Zi-Pi roles with igraph.
#   "majority selection" T = keep ASVs present in >= T samples; T = 0 -> present in > half of samples.
source("scripts/R_common.R")
library(SpiecEasi); library(igraph); library(pheatmap)
ps <- readRDS(file.path(RDS, "ps_filtered.rds"))
NCPU  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 4))
NBOOT <- as.integer(Sys.getenv("SPARCC_BOOT", 200))
R_MIN <- 0.2; P_MAX <- 0.01
asv_tax <- read.csv(tab_out("ASV_taxonomy.csv"), row.names = 1)
genus_of <- function(ids) { g <- asv_tax[ids, "Genus"]; ifelse(is.na(g), "Unclassified", g) }

run_sparcc <- function(ps_sub, T, name, p_max = P_MAX) {
  X <- t(as(otu_table(ps_sub), "matrix"))
  min_n <- if (T == 0) floor(nrow(X) / 2) + 1 else T             # T = 1 for healthy bark, as in the paper
  X <- X[, colSums(X > 0) >= min_n, drop = FALSE]
  cat(sprintf("[%s] %d samples, %d ASVs present in >= %d samples\n", name, nrow(X), ncol(X), min_n))
  if (ncol(X) < 4 || nrow(X) < 4) { cat(sprintf("[%s] too few ASVs/samples for SparCC - skipped\n", name)); return(NULL) }
  set.seed(SEED)
  sp <- sparcc(X)
  bt <- sparccboot(X, R = NBOOT, ncpus = NCPU)
  pv <- pval.sparccboot(bt)
  n <- ncol(X); P <- matrix(1, n, n); P[upper.tri(P)] <- pv$pvals; P[lower.tri(P)] <- t(P)[lower.tri(P)]
  C <- sp$Cor; dimnames(C) <- dimnames(P) <- list(colnames(X), colnames(X))
  saveRDS(list(X = X, cor = C, p = P), file.path(RDS, paste0("sparcc_", name, ".rds")))

  A <- C * (abs(C) >= R_MIN & P < p_max); diag(A) <- 0
  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- delete_vertices(g, degree(g) == 0)
  if (vcount(g) < 3) { cat(sprintf("[%s] < 3 connected nodes after filtering - skipped\n", name)); return(NULL) }
  E(g)$sign <- ifelse(E(g)$weight > 0, "positive", "negative")
  V(g)$genus <- genus_of(V(g)$name)

  # modules (fast greedy on |r|) and Zi-Pi roles (Guimera & Amaral; MENAP thresholds 2.5 / 0.62)
  gw <- g; E(gw)$weight <- abs(E(gw)$weight)
  mod <- cluster_fast_greedy(gw); V(g)$module <- membership(mod)
  k <- degree(g); zi <- pi_ <- numeric(vcount(g))
  adj <- as_adjacency_matrix(g, sparse = FALSE) != 0
  for (m in unique(V(g)$module)) {
    idx <- which(V(g)$module == m)
    kin <- rowSums(adj[idx, idx, drop = FALSE])
    zi[idx] <- if (length(idx) > 1 && sd(kin) > 0) (kin - mean(kin)) / sd(kin) else 0
  }
  for (i in seq_len(vcount(g))) {
    kt <- tapply(adj[i, ], V(g)$module, sum); pi_[i] <- 1 - sum((kt / k[i])^2)
  }
  role <- ifelse(zi > 2.5 & pi_ > 0.62, "Network hub", ifelse(zi > 2.5, "Module hub",
                 ifelse(pi_ > 0.62, "Connector", "Peripheral")))
  nodes <- data.frame(ASV = V(g)$name, genus = V(g)$genus, species = asv_tax[V(g)$name, "Species"],
                      module = V(g)$module, degree = k,
                      betweenness = betweenness(g, normalized = TRUE, weights = NA),
                      closeness = closeness(g, normalized = TRUE, weights = NA),
                      eigenvector = eigen_centrality(g, weights = NA)$vector,
                      Zi = zi, Pi = pi_, role = role) %>% arrange(-degree)
  write.csv(nodes, tab_out(paste0("Network_", name, "_nodes.csv")), row.names = FALSE)
  el <- as_data_frame(g, "edges"); el$genus_from <- genus_of(el$from); el$genus_to <- genus_of(el$to)
  el$p <- P[cbind(el$from, el$to)]
  write.csv(el, tab_out(paste0("Network_", name, "_edges.csv")), row.names = FALSE)
  erw <- el %>% filter(genus_from == "Erwinia" | genus_to == "Erwinia")
  write.csv(erw, tab_out(paste0("Network_", name, "_Erwinia_edges.csv")), row.names = FALSE)

  props <- data.frame(network = name, samples = nrow(X), ASVs_input = ncol(X), threshold_T = T, r_min = R_MIN, p_max = p_max,
                      nodes = vcount(g), edges = ecount(g),
                      positive = sum(E(g)$weight > 0), negative = sum(E(g)$weight < 0),
                      density = edge_density(g), avg_degree = mean(k),
                      avg_clustering = transitivity(g, "average"),
                      avg_path_length = mean_distance(g, weights = NA),
                      modularity = modularity(mod), modules = length(unique(V(g)$module)),
                      mean_closeness = mean(nodes$closeness, na.rm = TRUE),
                      mean_betweenness = mean(nodes$betweenness),
                      mean_eigenvector = mean(nodes$eigenvector),
                      max_degree_node = paste(nodes$ASV[1], nodes$genus[1]),
                      max_betweenness_node = with(nodes[which.max(nodes$betweenness), ], paste(ASV, genus)),
                      module_hubs = paste(with(nodes[nodes$role %in% c("Module hub", "Network hub"), ],
                                               paste0(ASV, " (", genus, ")")), collapse = "; "),
                      Erwinia_edges_pos = sum(erw$weight > 0), Erwinia_edges_neg = sum(erw$weight < 0))
  list(g = g, nodes = nodes, props = props, C = C, P = P, name = name, p_max = p_max)
}

plot_net <- function(r, title) {
  g <- r$g
  topg <- names(sort(table(V(g)$genus), decreasing = TRUE))
  topg <- head(setdiff(topg, "Unclassified"), 11)
  pal <- setNames(RColorBrewer::brewer.pal(12, "Paired")[-6][seq_along(topg)], topg)   # [-6]: drop red, reserved for Erwinia
  vcol <- ifelse(V(g)$genus == "Erwinia", "red", ifelse(V(g)$genus %in% topg, pal[V(g)$genus], "grey75"))
  set.seed(SEED)
  lay <- layout_with_fr(g, weights = abs(E(g)$weight))
  par(mar = c(1, 1, 3, 12), xpd = TRUE)
  plot(g, layout = lay, vertex.size = ifelse(V(g)$genus == "Erwinia", 7, 3 + 4 * degree(g) / max(degree(g))),
       vertex.color = vcol, vertex.frame.color = NA,
       vertex.label = ifelse(V(g)$genus == "Erwinia" | r$nodes$role[match(V(g)$name, r$nodes$ASV)] != "Peripheral",
                             paste0(V(g)$name, "\n", V(g)$genus), NA),
       vertex.label.cex = .5, vertex.label.color = "black",
       edge.color = ifelse(E(g)$weight > 0, adjustcolor("#2ca02c", .5), adjustcolor("#d62728", .5)),
       edge.width = .4 + 1.5 * abs(E(g)$weight), main = title)
  legend("topright", inset = c(-0.28, 0), bty = "n", cex = .65, pt.cex = 1.2, pch = 16,
         col = c("red", pal, "grey75"), legend = c("Erwinia", names(pal), "other"))
  legend("bottomright", inset = c(-0.28, 0), bty = "n", cex = .65, lwd = 2,
         col = c("#2ca02c", "#d62728"), legend = c("positive", "negative"))
  mtext(sprintf("%d nodes, %d edges (%d +, %d -), density %.3f, modularity %.2f",
                r$props$nodes, r$props$edges, r$props$positive, r$props$negative,
                r$props$density, r$props$modularity), side = 3, line = 0, cex = .7)
}

can <- prune_samples(sample_data(ps)$health == "Canker", ps)
hea <- prune_samples(sample_data(ps)$health == "Healthy", ps)
nets <- list(run_sparcc(can, 5, "canker_T5"),
             run_sparcc(can, 0, "canker_T0"),
             # healthy bark: 8 samples with 7-177 bacterial reads (99.8% of reads are chloroplast);
             # no pair reaches P < 0.01, so this network is empty and Fig. 8 is not reproducible
             run_sparcc(hea, 1, "healthy_T1"))
titles <- c("Fig. 7A: canker network, majority selection T = 5",
            "Fig. 7C: canker network, majority selection T = 0 (> half of samples)",
            "Fig. 8A: healthy bark network, majority selection T = 1")
ok <- !sapply(nets, is.null); nets <- nets[ok]; titles <- titles[ok]
props <- bind_rows(lapply(nets, `[[`, "props"))
write.csv(props, tab_out("Table1_network_properties.csv"), row.names = FALSE)
print(t(props))

pdf(pdf_out("Fig7_8_SparCC_networks.pdf"), width = 11, height = 9)
for (i in seq_along(nets)) plot_net(nets[[i]], titles[i])
if (!ok[3]) {   # say so in the PDF instead of silently leaving Fig. 8 out
  plot.new(); title("Fig. 8A/B: healthy bark network - not reproducible")
  text(.5, .6, paste("8 healthy samples with bacterial reads (7-177 reads each; 99.8% of healthy reads are chloroplast).",
                     "No ASV pair passes |r| >= 0.2 and SparCC P < 0.01 (paper's Methods), so the network has no edges.",
                     sep = "\n"), cex = .9)
}
for (r in nets) {
  # Fig. 7B/D, 8B: correlation heatmap of the 30 highest-degree nodes (edges failing the filter set to 0)
  top <- head(r$nodes$ASV, 30)
  M <- r$C[top, top] * (r$P[top, top] < r$p_max & abs(r$C[top, top]) >= R_MIN); diag(M) <- 1
  draw_pheatmap(M, color = colorRampPalette(c("#d62728", "white", "#2ca02c"))(101), breaks = seq(-1, 1, length.out = 102),
           labels_row = paste(top, genus_of(top)), labels_col = paste(top, genus_of(top)),
           fontsize = 6, main = sprintf("%s: SparCC correlations |r| >= %.1f%s, top 30 nodes by degree", r$name, R_MIN,
                         if (is.finite(r$p_max)) sprintf(", P < %g", r$p_max) else ", no P-value filter"))
}
dev.off()

zp <- bind_rows(lapply(nets, function(r) mutate(r$nodes, network = r$name)))
p_zp <- ggplot(zp, aes(Pi, Zi, colour = role)) + geom_point(size = 1.5) +
  geom_hline(yintercept = 2.5, linetype = 2) + geom_vline(xintercept = 0.62, linetype = 2) +
  ggrepel::geom_text_repel(data = subset(zp, role != "Peripheral"), aes(label = paste(ASV, genus)), size = 2.5) +
  facet_wrap(~ network) + labs(title = "Supplementary Fig. S3: Zi-Pi plot of node roles",
                               x = "Among-module connectivity (Pi)", y = "Within-module connectivity (Zi)")
pdf(pdf_out("S3_ZiPi_node_roles.pdf"), width = 12, height = 4.5)
print(p_zp)
dev.off()
