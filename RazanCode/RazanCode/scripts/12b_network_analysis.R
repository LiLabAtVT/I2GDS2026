# Step 12b: network topology (paper Table 1), module roles (Zi/Pi, Suppl. Fig. S3), E. amylovora partners (Fig. 7B/D),
#           network plots (Figs. 7A/C, 8A) + GraphML for Cytoscape/Gephi
#   Edges: FastSpar SparCC |r| >= 0.2 and bootstrap P < 0.01 (paper settings). Modules: igraph fast-greedy.
#   Zi/Pi roles (Olesen et al. 2007, as in MENAP): module hub Zi > 2.5 & Pi <= 0.62; connector Zi <= 2.5 & Pi > 0.62;
#   network hub Zi > 2.5 & Pi > 0.62; else peripheral.
# Run via scripts/12_sparcc_networks.sh (env fb_r)
suppressPackageStartupMessages({ library(igraph); library(ggplot2); library(dplyr) })
set.seed(2026)

R_MIN <- 0.2
args  <- commandArgs(trailingOnly = TRUE)
P_MAX <- if (length(args)) as.numeric(args[1]) else 0.01        # paper: P < 0.01 (networks); Fig. 7B quotes P < 0.05
TAG   <- if (P_MAX == 0.01) "" else paste0("_p", P_MAX)          # outputs for other thresholds get a suffix
base <- "results/12_networks"
map  <- read.delim("results/10_erwinia_dmm/asv_id_map.tsv")
sp   <- read.delim("results/11_faprotax/species_labels_used.tsv")
map$species <- sub("^s__", "", sp$Species[match(map$feature_id, sp$feature_id)])
map$node_label <- ifelse(!is.na(map$species) & map$species != "", map$species, map$genus)
map$node_label[is.na(map$node_label) | map$node_label == ""] <- "unclassified"
ea <- map$asv_num[grepl("Erwinia", map$genus) & grepl("amylovora", map$species)]

read_sq <- function(f) { m <- as.matrix(read.delim(f, check.names = FALSE, row.names = 1)); m[rownames(m), rownames(m)] }

zipi <- function(g, memb) {
  A <- as_adjacency_matrix(g, sparse = FALSE); k <- rowSums(A)
  Z <- P <- setNames(numeric(vcount(g)), V(g)$name)
  for (m in unique(memb)) {
    idx <- which(memb == m); kin <- rowSums(A[idx, idx, drop = FALSE])
    Z[idx] <- if (length(idx) > 1 && sd(kin) > 0) (kin - mean(kin)) / sd(kin) else 0
  }
  for (i in seq_along(k)) {
    if (k[i] == 0) { P[i] <- 0; next }
    kis <- tapply(A[i, ], memb, sum); P[i] <- 1 - sum((kis / k[i])^2)
  }
  role <- ifelse(Z > 2.5 & P > 0.62, "Network hub", ifelse(Z > 2.5, "Module hub", ifelse(P > 0.62, "Connector", "Peripheral")))
  data.frame(ASV = names(Z), Zi = round(Z, 3), Pi = round(P, 3), role = role)
}

# distance-based measures count edges (hops), as Cytoscape NetworkAnalyzer does; |r| is only used for modules and plots
topology <- function(g, nm) {
  d <- degree(g); n <- vcount(g)
  dd <- table(d[d > 0]); pl <- if (length(dd) > 2) summary(lm(log(as.numeric(dd)) ~ log(as.numeric(names(dd)))))$r.squared else NA
  D <- distances(g, weights = NA); D[is.infinite(D)] <- NA; diag(D) <- NA
  cs <- components(g)$csize
  bt <- betweenness(g, weights = NA)
  # stress centrality: number of shortest paths (all pairs) passing through a node as an intermediate
  st <- setNames(numeric(n), V(g)$name)
  for (s in seq_len(n)) for (pth in all_shortest_paths(g, s, weights = NA)$vpaths)
    if (length(pth) > 2) { mid <- as.integer(pth[-c(1, length(pth))]); st[mid] <- st[mid] + 1 }
  st <- st / 2                                                    # undirected: each path counted from both ends
  ev <- eigen_centrality(g)$vector; cl <- closeness(g, weights = NA)
  top <- function(x) paste(paste0(names(x)[x == max(x)], "-", map$node_label[match(names(x)[x == max(x)], map$asv_num)]), collapse = "; ")
  mod <- cluster_fast_greedy(g, weights = E(g)$weight)
  data.frame(network = nm, nodes = n, edges = ecount(g), positive = sum(E(g)$sign > 0), negative = sum(E(g)$sign < 0),
             R2_power_law = round(pl, 3), avgK = round(mean(d), 3), avgCC = round(transitivity(g, "average", isolates = "zero"), 3),
             avg_path_GD = round(mean(D, na.rm = TRUE), 3), geodesic_efficiency_E = round(global_efficiency(g, weights = NA), 3),
             harmonic_distance_HD = round(1 / global_efficiency(g, weights = NA), 3),
             max_degree = max(d), nodes_max_degree = top(d), centralization_degree = round(centr_degree(g)$centralization, 3),
             max_betweenness = round(max(bt), 3), nodes_max_betweenness = top(bt),
             centralization_betweenness = round(centr_betw(g)$centralization, 4),
             max_stress = max(st), nodes_max_stress = top(st),
             max_eigenvector = round(max(ev), 3), nodes_max_eigenvector = top(ev),
             centralization_eigenvector = round(centr_eigen(g)$centralization, 3),
             max_closeness = round(max(cl, na.rm = TRUE), 4), nodes_max_closeness = top(replace(cl, is.na(cl), 0)),
             centralization_closeness = round(suppressWarnings(centr_clo(g)$centralization), 4),
             density = round(edge_density(g), 3), transitivity = round(transitivity(g, "global"), 3),
             connectedness = round(sum(cs * (cs - 1)) / (n * (n - 1)), 3),
             modules = length(mod), modularity = round(modularity(mod), 3))
}

tabs <- list(); roles <- list(); eas <- list()
for (nm in c("canker_T5", "canker_T0", "healthy_T1")) {
  D <- file.path(base, nm)
  if (!file.exists(file.path(D, "pvalues.tsv"))) { cat("no SparCC result for", nm, "\n"); next }
  r <- read_sq(file.path(D, "correlation.tsv")); p <- read_sq(file.path(D, "pvalues.tsv"))[rownames(r), colnames(r)]
  ut <- which(upper.tri(r) & abs(r) >= R_MIN & p < P_MAX, arr.ind = TRUE)
  ed <- data.frame(from = rownames(r)[ut[, 1]], to = colnames(r)[ut[, 2]], r = r[ut], p = p[ut])
  write.table(ed, file.path(D, paste0("edges", TAG, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  g <- graph_from_data_frame(ed, directed = FALSE, vertices = data.frame(name = unique(c(ed$from, ed$to))))
  E(g)$weight <- abs(ed$r); E(g)$sign <- sign(ed$r); E(g)$r <- ed$r
  V(g)$label <- map$node_label[match(V(g)$name, map$asv_num)]
  V(g)$genus <- map$genus[match(V(g)$name, map$asv_num)]
  V(g)$is_Ea <- V(g)$name %in% ea
  if (ecount(g) == 0) { cat(nm, ": no significant edges\n"); next }

  tabs[[nm]] <- topology(g, nm)
  mod <- cluster_fast_greedy(g, weights = E(g)$weight); V(g)$module <- membership(mod)
  zp <- zipi(g, membership(mod)); zp$label <- V(g)$label; zp$network <- nm; zp$degree <- degree(g)
  roles[[nm]] <- zp; V(g)$role <- zp$role
  write_graph(g, file.path(D, paste0(nm, TAG, ".graphml")), format = "graphml")

  # E. amylovora partners (Fig. 7B / 7D)
  e1 <- ed |> filter(from %in% ea | to %in% ea) |>
    mutate(Ea_ASV = ifelse(from %in% ea, from, to), partner = ifelse(from %in% ea, to, from),
           partner_label = paste0(partner, " ", map$node_label[match(partner, map$asv_num)]),
           interaction = ifelse(r > 0, "Positive", "Negative"), network = nm)
  if (nrow(e1)) eas[[nm]] <- e1

  # network plot
  top_gen <- names(sort(table(V(g)$genus), decreasing = TRUE))[1:min(10, length(unique(V(g)$genus)))]
  gcol <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#66C2A5", "#FFD92F", "#8DA0CB"),
                   top_gen)[seq_along(top_gen)]
  vcol <- ifelse(V(g)$genus %in% names(gcol), gcol[V(g)$genus], "grey75"); vcol[V(g)$is_Ea] <- "black"
  lab <- ifelse(degree(g) >= quantile(degree(g), 0.9) | V(g)$is_Ea | V(g)$role != "Peripheral", V(g)$label, NA)
  pdf(sprintf("results/pdf/12_network_%s%s.pdf", nm, TAG), width = 11, height = 10)
  par(mar = c(1, 1, 3, 1))
  lay <- layout_with_fr(g, weights = E(g)$weight)
  plot(g, layout = lay, vertex.size = 3 + 9 * degree(g) / max(degree(g)), vertex.color = vcol, vertex.frame.color = "grey30",
       vertex.label = lab, vertex.label.cex = 0.6, vertex.label.color = "black", vertex.label.dist = 0.6,
       edge.color = ifelse(E(g)$sign > 0, adjustcolor("#1B9E77", 0.45), adjustcolor("#D95F02", 0.45)),
       edge.width = 0.5 + 2 * E(g)$weight, edge.curved = 0.15)
  legend("topleft", legend = c("E. amylovora", names(gcol), "other"), pt.bg = c("black", gcol, "grey75"), pch = 21, cex = 0.7, bty = "n")
  legend("topright", legend = c("positive", "negative"), col = c("#1B9E77", "#D95F02"), lwd = 2, cex = 0.8, bty = "n")
  title(sprintf("%s: %d nodes, %d edges (%d+ / %d-)  |r| >= %.1f, P < %.2f", nm, vcount(g), ecount(g),
                sum(E(g)$sign > 0), sum(E(g)$sign < 0), R_MIN, P_MAX))
  if (nm == "healthy_T1") mtext("CAUTION: 8 healthy samples with 7-177 bacterial reads each - correlations are not reliable", col = "red")
  print(ggplot(zp, aes(Pi, Zi, colour = role)) + geom_point(size = 2) +
          geom_hline(yintercept = 2.5, linetype = 2) + geom_vline(xintercept = 0.62, linetype = 2) +
          geom_text(data = subset(zp, role != "Peripheral"), aes(label = paste(ASV, label)), size = 2.5, vjust = -0.8, show.legend = FALSE) +
          labs(title = paste("Zi-Pi module roles -", nm), x = "Among-module connectivity (Pi)", y = "Within-module connectivity (Zi)") + theme_bw())
  if (nrow(e1)) print(ggplot(e1, aes(partner_label, Ea_ASV, colour = interaction, size = abs(r))) + geom_point() +
          scale_colour_manual(values = c(Negative = "#F08080", Positive = "#20B2AA")) +
          scale_size_continuous(name = "|Correlation|", range = c(2, 7)) +
          labs(x = NULL, y = NULL, title = paste("E. amylovora ASV partners -", nm)) + theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)))
  dev.off()
}

T1 <- bind_rows(tabs); write.table(T1, file.path(base, paste0("table1_topology", TAG, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
RL <- bind_rows(roles); write.table(RL, file.path(base, paste0("zipi_roles", TAG, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
EA <- bind_rows(eas); write.table(EA, file.path(base, paste0("erwinia_amylovora_edges", TAG, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

paper <- data.frame(network = c("canker_T5", "canker_T0", "healthy_T1"), nodes_p = c(160, 17, 51), edges_p = c(1586, 27, 698),
                    avgK_p = c(19.825, 3.176, 27.372), avgCC_p = c(0.267, 0.439, 0.701), density_p = c(0.124, 0.198, 0.547),
                    maxdeg_node_p = c("ASV451-Mucilaginibacter", "ASV1032-E. amylovora", "ASV1308-Bacillus"))
cat("\n==== Table 1 (this analysis vs paper) ====\n")
print(merge(T1[, c("network", "nodes", "edges", "positive", "negative", "avgK", "avgCC", "density", "modularity",
                   "nodes_max_degree", "nodes_max_betweenness")], paper, by = "network"), row.names = FALSE)
cat("\n==== Hubs / connectors (non-peripheral) ====\n"); print(subset(RL, role != "Peripheral"), row.names = FALSE)
cat("\n==== E. amylovora edges ====\n")
print(EA |> group_by(network, interaction) |> summarise(n = n(), min_r = round(min(r), 2), max_r = round(max(r), 2), .groups = "drop"))
print(EA |> select(network, Ea_ASV, partner_label, r) |> arrange(network, r) |> mutate(r = round(r, 2)), row.names = FALSE)
