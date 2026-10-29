#!/usr/bin/env Rscript
# Step 9: SparCC co-occurrence networks (paper: SparCC via MENAP; |r| >= 0.2, P < 0.01; Figs. 7-8, Table 1).
# Usage: Rscript 09_networks.R prep      -> write FastSpar input tables
#        Rscript 09_networks.R analyze   -> build networks from FastSpar correlations / p-values
# "Majority selection" thresholds as described in the paper:
#   canker T5  = ASVs present in >= 5 of the 32 canker samples
#   canker T0  = ASVs present in more than half of the canker samples (paper: "0: keep only ASVs detected in more than half")
#   healthy T1 = ASVs present in >= 1 healthy sample (only healthy samples with bacterial reads)
suppressPackageStartupMessages({ library(dplyr); library(igraph); library(ggplot2); library(ggraph) })
mode <- commandArgs(trailingOnly = TRUE)[1]
B <- "/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code"; O <- file.path(B, "results/09_networks")
dir.create(O, showWarnings = FALSE, recursive = TRUE)
NETS <- c("canker_T5", "canker_T0", "healthy_T1")
R_MIN <- 0.2; P_MAX <- 0.01

lk <- read.delim(file.path(B, "results/07_composition_ml/asv_lookup.tsv"), check.names = FALSE)   # bacterial ASVs
meta <- read.delim(file.path(B, "raw_data/sample_metadata.tsv"), check.names = FALSE)
meta <- meta[meta[["sample-id"]] != "#q2:types", ]; rownames(meta) <- meta[["sample-id"]]

if (mode == "prep") {
  ft <- as.matrix(read.delim(file.path(B, "results/04_dada2/export/feature-table.tsv"), skip = 1,
                             row.names = 1, check.names = FALSE, comment.char = ""))
  ft <- ft[lk$feature_id, ]; rownames(ft) <- lk$ASV
  can <- ft[, meta[colnames(ft), "health_status"] == "Canker"]
  hea <- ft[, meta[colnames(ft), "health_status"] == "Healthy"]; hea <- hea[, colSums(hea) > 0]
  sets <- list(canker_T5 = can[rowSums(can > 0) >= 5, ],
               canker_T0 = can[rowSums(can > 0) > ncol(can) / 2, ],
               healthy_T1 = hea[rowSums(hea > 0) >= 1, ])
  for (n in names(sets)) {
    d <- file.path(O, n); dir.create(d, showWarnings = FALSE)
    m <- sets[[n]]
    write.table(data.frame(`#OTU ID` = rownames(m), m, check.names = FALSE), file.path(d, "counts.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(n, ":", nrow(m), "ASVs x", ncol(m), "samples\n")
  }
  quit(save = "no")
}

# ---------------- analyze ----------------
sp <- tryCatch(read.delim(file.path(B, "results/08b_blast/species_assignments.tsv")), error = function(e) NULL)
lab <- setNames(paste0(lk$ASV, "-", lk$Genus), lk$ASV)
if (!is.null(sp)) { s <- sp[sp$assigned == "yes", ]; lab[s$ASV] <- paste0(s$ASV, "-", s$species) }
genus <- setNames(lk$Genus, lk$ASV)

read_mat <- function(f) { m <- as.matrix(read.delim(f, row.names = 1, check.names = FALSE, comment.char = "")); m }
net_stats <- function(g) {
  k <- degree(g); n <- vcount(g); d <- distances(g); d <- d[upper.tri(d)]
  fin <- is.finite(d)
  dd <- table(k[k > 0]); pl <- if (length(dd) > 2) summary(lm(log(as.numeric(dd)) ~ log(as.numeric(names(dd)))))$r.squared else NA
  mods <- cluster_fast_greedy(g, weights = abs(E(g)$rho))
  data.frame(
    Nodes = n, Links = ecount(g), Positive = sum(E(g)$rho > 0), Negative = sum(E(g)$rho < 0),
    `R square of power-law` = pl, `Average degree (avgK)` = mean(k),
    `Average clustering coefficient (avgCC)` = mean(transitivity(g, type = "local"), na.rm = TRUE),
    `Average path distance (GD)` = mean(d[fin]), `Geodesic efficiency (E)` = mean(ifelse(fin, 1 / d, 0)),
    `Harmonic geodesic distance (HD)` = 1 / mean(ifelse(fin, 1 / d, 0)),
    `Maximal degree` = max(k), `Nodes with maximal degree` = paste(lab[names(k)[k == max(k)]], collapse = "; "),
    `Centralization of degree (CD)` = centr_degree(g)$centralization,
    `Maximal betweenness` = max(betweenness(g)), `Nodes with maximal betweenness` = paste(lab[names(which.max(betweenness(g)))], collapse = "; "),
    `Centralization of betweenness (CB)` = centr_betw(g)$centralization,
    `Centralization of closeness (CCL)` = tryCatch(centr_clo(g)$centralization, error = function(e) NA),
    `Nodes with maximal eigenvector centrality` = lab[names(which.max(eigen_centrality(g)$vector))],
    `Centralization of eigenvector centrality (CE)` = centr_eigen(g)$centralization,
    `Density (D)` = edge_density(g), `Transitivity (Trans)` = transitivity(g, type = "global"),
    `Connectedness (Con)` = mean(fin), Modules = length(mods), Modularity = modularity(mods),
    check.names = FALSE)
}
zipi <- function(g) {  # within-module degree z (Zi) and among-module connectivity (Pi); Olesen et al. / MENAP roles
  mods <- membership(cluster_fast_greedy(g, weights = abs(E(g)$rho)))
  A <- as_adjacency_matrix(g, sparse = FALSE) != 0
  k <- rowSums(A)
  kin <- sapply(seq_along(k), function(i) sum(A[i, mods == mods[i]]))
  z <- ave(kin, mods, FUN = function(x) if (length(x) > 1 && sd(x) > 0) (x - mean(x)) / sd(x) else 0)
  P <- 1 - sapply(seq_along(k), function(i) if (k[i] == 0) 1 else sum((tapply(A[i, ], mods, sum) / k[i])^2))
  data.frame(ASV = names(k), label = lab[names(k)], module = as.integer(mods), degree = k, Zi = z, Pi = P,
             role = ifelse(z > 2.5 & P > 0.62, "network hub", ifelse(z > 2.5, "module hub", ifelse(P > 0.62, "connector", "peripheral"))))
}
paper <- data.frame(check.names = FALSE,
  canker_T5 = c("160", "1586", "0.033", "19.825", "0.267", "2.032", "0.5361", "1.865", "39", "ASV451-Mucilaginibacter", "0.120", "ASV345-Nocardioides", "0.0160", "0.124", "0.276", "1"),
  canker_T0 = c("17", "27", "0.913", "3.176", "0.439", "1.990", "0.588", "1.700", "11", "ASV1032-E. amylovora", "0.489", "ASV1032-E. amylovora", "0.411", "0.198", "0.381", "0.780"),
  healthy_T1 = c("51", "698", "0.0007", "27.372", "0.701", "1.456", "0.773", "1.293", "45", "ASV1308-Bacillus", "0.352", "ASV1308-Bacillus", "0.0241", "0.547", "0.677", "1"),
  row.names = c("Nodes", "Links", "R square of power-law", "Average degree (avgK)", "Average clustering coefficient (avgCC)",
                "Average path distance (GD)", "Geodesic efficiency (E)", "Harmonic geodesic distance (HD)", "Maximal degree",
                "Nodes with maximal degree", "Centralization of degree (CD)", "Nodes with maximal betweenness",
                "Centralization of betweenness (CB)", "Density (D)", "Transitivity (Trans)", "Connectedness (Con)"))

pdf(file.path(B, "pdf/09_networks.pdf"), width = 14, height = 11)
allstats <- list(); erw_edges <- list()
for (n in NETS) {
  d <- file.path(O, n)
  r <- read_mat(file.path(d, "median_correlation.tsv")); p <- read_mat(file.path(d, "pvalues.tsv"))
  p <- p[rownames(r), colnames(r)]
  ut <- which(upper.tri(r) & abs(r) >= R_MIN & p < P_MAX, arr.ind = TRUE)
  ed <- data.frame(from = rownames(r)[ut[, 1]], to = colnames(r)[ut[, 2]], rho = r[ut], p = p[ut])
  write.table(ed %>% mutate(from_label = lab[from], to_label = lab[to]), file.path(O, paste0(n, "_edges.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  g <- graph_from_data_frame(ed, directed = FALSE, vertices = data.frame(name = unique(c(ed$from, ed$to))))
  cat(n, ": ASVs tested", nrow(r), "-> nodes", vcount(g), "edges", ecount(g), "\n")
  if (ecount(g) == 0) { allstats[[n]] <- data.frame(Nodes = 0, Links = 0, check.names = FALSE); next }
  st <- net_stats(g); allstats[[n]] <- st
  zp <- zipi(g); write.table(zp, file.path(O, paste0(n, "_zi_pi.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  hubs/connectors:\n"); print(zp[zp$role != "peripheral", c("label", "module", "degree", "Zi", "Pi", "role")], row.names = FALSE)
  # Erwinia edges (Fig. 7B/D)
  ee <- ed %>% filter(genus[from] == "Erwinia" | genus[to] == "Erwinia") %>%
    mutate(Erwinia_ASV = ifelse(genus[from] == "Erwinia", lab[from], lab[to]), partner = ifelse(genus[from] == "Erwinia", lab[to], lab[from]),
           sign = ifelse(rho > 0, "positive", "negative"), network = n) %>% arrange(Erwinia_ASV, rho)
  erw_edges[[n]] <- ee
  # (correlation stored as "rho", not "weight", so igraph path/centrality measures stay unweighted)
  # plot network, nodes coloured by genus (top 12 genera), edges by sign
  V(g)$genus <- genus[V(g)$name]; topg <- names(sort(table(V(g)$genus), decreasing = TRUE))[1:min(12, length(unique(V(g)$genus)))]
  V(g)$gcol <- ifelse(V(g)$genus %in% topg, V(g)$genus, "other"); V(g)$deg <- degree(g)
  V(g)$lab <- ifelse(V(g)$genus == "Erwinia" | V(g)$deg >= quantile(V(g)$deg, 0.9), lab[V(g)$name], "")
  E(g)$sign <- ifelse(E(g)$rho > 0, "positive", "negative")
  set.seed(1)
  print(ggraph(g, layout = "fr") + geom_edge_link(aes(colour = sign), alpha = 0.35, width = 0.4) +
    geom_node_point(aes(fill = gcol, size = deg), shape = 21, colour = "grey20") +
    geom_node_text(aes(label = lab), size = 2.3, repel = TRUE, max.overlaps = 50) +
    scale_edge_colour_manual(values = c(positive = "#1b9e77", negative = "#d95f02")) + scale_size(range = c(1.5, 7)) +
    labs(title = sprintf("%s SparCC network: %d nodes, %d edges (|r| >= %.1f, P < %.2f)", n, vcount(g), ecount(g), R_MIN, P_MAX),
         fill = "Genus", size = "Degree", edge_colour = "Correlation") + theme_void(11))
  if (nrow(ee) > 0)
    print(ggplot(ee, aes(partner, Erwinia_ASV, colour = sign, size = abs(rho))) + geom_point() +
      scale_colour_manual(values = c(positive = "#1b9e77", negative = "#e7298a")) +
      labs(x = NULL, y = NULL, size = "|r|", title = paste0(n, ": correlations of Erwinia ASVs with other taxa (Fig. 7B/D)")) +
      theme_bw(9) + theme(axis.text.x = element_text(angle = 60, hjust = 1)))
  zp$lbl <- ifelse(zp$role != "peripheral", zp$label, "")
  print(ggplot(zp, aes(Pi, Zi, colour = role)) + geom_point() + geom_vline(xintercept = 0.62, lty = 2) + geom_hline(yintercept = 2.5, lty = 2) +
    geom_text(aes(label = lbl), size = 2.5, vjust = -0.6, show.legend = FALSE) + labs(title = paste0(n, ": Zi-Pi node roles (Supplementary Fig. S3)")) + theme_bw(11))
}
fmt <- function(x) if (is.numeric(x)) formatC(x, digits = 4, format = "fg") else as.character(x)
st <- bind_rows(allstats, .id = "network")
tab <- t(sapply(split(st, st$network)[NETS], function(r) sapply(r[-1], fmt)))
tab <- as.data.frame(t(tab), check.names = FALSE)
write.table(tab, file.path(O, "network_topology_ours.tsv"), sep = "\t", quote = FALSE)
cmp <- data.frame(Parameter = rownames(paper), check.names = FALSE)
for (n in NETS) { cmp[[paste(n, "ours")]] <- tab[rownames(paper), n]; cmp[[paste(n, "paper")]] <- paper[[n]] }
write.table(cmp, file.path(O, "table1_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 250); print(cmp, row.names = FALSE); print(tab[c("Positive", "Negative", "Modules", "Modularity"), , drop = FALSE])
ee <- bind_rows(Filter(function(x) nrow(x) > 0, erw_edges)); write.table(ee, file.path(O, "erwinia_edges.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(ee)) print(ee %>% count(network, sign))
grid::grid.newpage()
gridExtra_ok <- requireNamespace("gridExtra", quietly = TRUE)
if (gridExtra_ok) gridExtra::grid.table(cmp, rows = NULL, theme = gridExtra::ttheme_minimal(base_size = 8)) else
  grid::grid.text(paste(capture.output(print(cmp, row.names = FALSE)), collapse = "\n"), gp = grid::gpar(fontfamily = "mono", fontsize = 7))
dev.off()
cat("Step 9 done\n")
