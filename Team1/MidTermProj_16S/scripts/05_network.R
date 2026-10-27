# Step 5: SparCC co-occurrence networks for canker and healthy samples (paper used the MENAP web pipeline;
# this reproduces its settings locally: majority selection, |SparCC r| >= 0.2, pseudo P < 0.01).
suppressPackageStartupMessages({ library(SpiecEasi); library(igraph) })
setwd("/projects/intro2gds/I2GDS2026/Team1/Kamal/Claude_Code")
set.seed(123)
ex <- "results/qiime2/export"; out <- "results/pdf"; tab <- "results/tables"
dir.create(out, TRUE, FALSE); dir.create(tab, TRUE, FALSE)
meta <- read.delim("metadata/sample-metadata.tsv", row.names = 1)
otu <- as.matrix(read.delim(file.path(ex, "rarefied_table.tsv"), skip = 1, row.names = 1, check.names = FALSE))
tax <- read.delim(file.path(ex, "taxonomy/taxonomy.tsv"), row.names = 1)
genus <- sub(".*g__", "", ifelse(grepl("g__", tax$Taxon), tax$Taxon, "Unassigned")); genus <- sub(";.*", "", genus)
names(genus) <- rownames(tax)

# ponytail: MENAP "majority selection" read as "ASV present in >= N samples"; nboot=200 so P<0.01 is resolvable.
build <- function(status, majority, nboot = 200) {
  x <- t(otu[, rownames(meta)[meta$status == status]])
  x <- x[, colSums(x > 0) >= max(majority, 1), drop = FALSE]
  cat(status, ": ", ncol(x), " ASVs in ", nrow(x), " samples\n", sep = "")
  bt <- sparccboot(x, R = nboot, ncpus = 16)
  pv <- pval.sparccboot(bt)
  n <- ncol(x); r <- p <- matrix(0, n, n)
  r[upper.tri(r)] <- pv$cors; p[upper.tri(p)] <- pv$pvals
  r <- r + t(r); p <- p + t(p); diag(p) <- 1
  adj <- ifelse(abs(r) >= 0.2 & p < 0.01, r, 0); dimnames(adj) <- list(colnames(x), colnames(x))
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- delete_vertices(g, degree(g) == 0)
  V(g)$genus <- genus[V(g)$name]; E(g)$sign <- ifelse(E(g)$weight > 0, "positive", "negative")
  e <- as_data_frame(g); e$genus_from <- genus[e$from]; e$genus_to <- genus[e$to]
  write.csv(e, file.path(tab, sprintf("network_%s_majority%d_edges.csv", status, majority)), row.names = FALSE)
  g
}
props <- function(g, label) {
  wc <- cluster_fast_greedy(as_undirected(g, mode = "each") |> (\(h) { E(h)$weight <- abs(E(h)$weight); h })())
  data.frame(network = label, nodes = vcount(g), edges = ecount(g), positive = sum(E(g)$weight > 0), negative = sum(E(g)$weight < 0),
             avg_degree = mean(degree(g)), avg_clustering = transitivity(g, "average"), avg_path = mean_distance(g),
             density = edge_density(g), modularity = modularity(wc), n_modules = length(wc),
             erwinia_nodes = sum(V(g)$genus == "Erwinia"))
}
plot_net <- function(g, label) {
  cols <- ifelse(V(g)$genus == "Erwinia", "red", "grey70")
  plot(g, layout = layout_with_fr(g), vertex.size = ifelse(V(g)$genus == "Erwinia", 6, 3), vertex.label = NA, vertex.color = cols,
       edge.color = ifelse(E(g)$weight > 0, adjustcolor("forestgreen", .5), adjustcolor("firebrick", .5)), edge.width = .6,
       main = sprintf("%s: %d nodes, %d edges (red = Erwinia ASVs)", label, vcount(g), ecount(g)))
  legend("bottomleft", c("positive", "negative"), col = c("forestgreen", "firebrick"), lty = 1, bty = "n")
  # Erwinia ego network: direct partners of Ea ASVs, labelled by genus
  ea <- V(g)[genus == "Erwinia"]
  if (length(ea)) {
    h <- induced_subgraph(g, unique(unlist(ego(g, 1, ea))))
    plot(h, layout = layout_with_fr(h), vertex.label = V(h)$genus, vertex.label.cex = .6, vertex.size = 5,
         vertex.color = ifelse(V(h)$genus == "Erwinia", "red", "lightblue"),
         edge.color = ifelse(E(h)$weight > 0, "forestgreen", "firebrick"), main = paste(label, "- Erwinia neighbourhood"))
  }
}
nets <- list(Canker_maj5 = build("Canker", 5), Canker_maj0 = build("Canker", 0), Healthy_maj1 = build("Healthy", 1))
pr <- do.call(rbind, Map(props, nets, names(nets))); print(pr)
write.csv(pr, file.path(tab, "network_global_properties.csv"), row.names = FALSE)
pdf(file.path(out, "Fig8_SparCC_networks.pdf"), 10, 10)
for (n in names(nets)) plot_net(nets[[n]], n)
dev.off()
cat("05_network.R done\n")
