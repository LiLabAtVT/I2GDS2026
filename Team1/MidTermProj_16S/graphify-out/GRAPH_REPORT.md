# Graph Report - Claude_Code  (2026-10-27)

## Corpus Check
- 17 files · ~24,778 words
- Verdict: corpus is large enough that graph structure adds value.
- Unclassified: 61 file(s) not represented in the graph (top: .gz 48, .log 4, .tsv 3)

## Summary
- 156 nodes · 154 edges · 19 communities (8 shown, 11 thin omitted)
- Extraction: 94% EXTRACTED · 6% INFERRED · 0% AMBIGUOUS · INFERRED: 9 edges (avg confidence: 0.82)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- Stats & ML Analysis (04_stats.R)
- Study Design & Sequencing
- Diversity & Pathobiome Findings
- Co-occurrence Network (05_network.R)
- FastQC Summary Plotting
- Graphify Project Rules
- DADA2 Test Script
- QIIME2 Processing Script
- Data Download Script
- FastQC Run Script
- R Stage Launcher
- Taxa Aggregation Helpers
- Pipeline Orchestrator
- Claude Code Token Usage: Tools and Methods

## God Nodes (most connected - your core abstractions)
1. `Fire Blight Canker Microbiome: Insights into the Key Ecological Niche for Erwinia amylovora (Phytobiomes J 2026)` - 9 edges
2. `ASV table (1,817,968 reads, 2,322 ASVs)` - 9 edges
3. `Claude Code Token Usage: Tools and Methods` - 6 edges
4. `SparCC co-occurrence networks via MENAP (Cytoscape/Gephi)` - 6 edges
5. `2. Shrink` - 5 edges
6. `Fig0 FastQC summary (48 HiFi samples)` - 5 edges
7. `1. Measure` - 4 edges
8. `SMRT PacBio HiFi full-length 16S rRNA metabarcoding (27F/1492R, Sequel IIe)` - 4 edges
9. `Microbiome-to-pathobiome shift (dysbiosis)` - 4 edges
10. `Graphify Project Rules` - 4 edges

## Surprising Connections (you probably didn't know these)
- `Raw HiFi reads per sample (total 3,529,407; canker > healthy)` --semantically_similar_to--> `ASV table (1,817,968 reads, 2,322 ASVs)`  [INFERRED] [semantically similar]
  results/pdf/Fig0_FastQC_summary.pdf → paper.txt
- `Fig0 FastQC summary (48 HiFi samples)` --shares_data_with--> `Sampling: 32 canker + 16 healthy wood samples, Augusta County VA, Nov 2022`  [INFERRED]
  results/pdf/Fig0_FastQC_summary.pdf → paper.txt
- `Fig0 FastQC summary (48 HiFi samples)` --implements--> `FastQC read quality assessment (Andrews 2010)`  [INFERRED]
  results/pdf/Fig0_FastQC_summary.pdf → paper.txt
- `Read length distribution (peaks ~1400-1500 bp)` --conceptually_related_to--> `SMRT PacBio HiFi full-length 16S rRNA metabarcoding (27F/1492R, Sequel IIe)`  [INFERRED]
  results/pdf/Fig0_FastQC_summary.pdf → paper.txt

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **16S amplicon processing pipeline (PacBio HiFi -> FastQC -> QIIME2/DADA2 -> Silva -> ASV table)** — paper_pacbio_hifi_16s, paper_fastqc, paper_qiime2, paper_dada2, paper_silva_138, paper_asv_table [EXTRACTED 1.00]
- **Machine-learning discrimination of canker taxa (Random Forest, LEfSe, DMM)** — paper_random_forest, paper_lefse, paper_dmm_clustering, paper_canker_dominant_taxa [EXTRACTED 1.00]
- **Evidence for microbiome-to-pathobiome transition** — paper_alpha_diversity, paper_beta_diversity_bray_curtis, paper_sparcc_menap_network, paper_faprotax, paper_pathobiome_shift [INFERRED 0.85]

## Communities (19 total, 11 thin omitted)

### Community 0 - "Stats & ML Analysis (04_stats.R)"
Cohesion: 0.04
Nodes (42): aic, al, alpha, aov_res, bc, bic, can, cl (+34 more)

### Community 1 - "Study Design & Sequencing"
Cohesion: 0.13
Nodes (17): Apple cultivars (Gala, Aztec Fuji, Cripps Pink, Granniwinkle, Hewe's Crab, unknown), Authors: BakhshiGanje, Jo, Li, Kong, Park, Acimovic (Virginia Tech / Kangwon / Chungbuk), Cripps Pink cultivar, Erwinia amylovora (Ea), Fire blight disease, Fire blight canker (ecological niche), Fire Blight Canker Microbiome: Insights into the Key Ecological Niche for Erwinia amylovora (Phytobiomes J 2026), NCBI SRA BioProject PRJNA1348523 (+9 more)

### Community 2 - "Diversity & Pathobiome Findings"
Cohesion: 0.13
Nodes (10): Acute oak decline polymicrobial canker (Brenneria goodwinii), ASV table (1,817,968 reads, 2,322 ASVs), Bacillus / Amnibacterium keystone taxa in healthy bark, Canker-dominant genera: Erwinia, Sphingomonas, Massilia, Pantoea, Novosphingobium, Canker module hubs: Novosphingobium, Variovorax, Nocardioides, Pantoea (dual antagonist/opportunist role), Microbiome-to-pathobiome shift (dysbiosis), PICRUSt2 (+2 more)

### Community 3 - "Co-occurrence Network (05_network.R)"
Cohesion: 0.13
Nodes (9): ex, genus, meta, nets, otu, out, pr, tab (+1 more)

### Community 4 - "FastQC Summary Plotting"
Cohesion: 0.22
Nodes (6): basic, files, id, len, meta, qual

### Community 5 - "Graphify Project Rules"
Cohesion: 0.60
Nodes (5): GRAPH_REPORT.md, graphify-out Knowledge Graph, graphify query/path/explain commands, Graphify Project Rules, graphify update (AST-only rebuild)

### Community 6 - "DADA2 Test Script"
Cohesion: 0.40
Nodes (3): fs, out, r

### Community 13 - "Claude Code Token Usage: Tools and Methods"
Cohesion: 0.12
Nodes (15): 1.1 `/context`: what fills the context window, 1.2 `/usage`: plan limits and what drives them, 1.3 `codeburn`: your session logs as a dashboard, 1. Measure, 2.1 `/clear` and `/compact`: start fresh or summarize, 2.2 graphify: a knowledge graph of the project, 2.3 GitDiagram: an architecture diagram from a GitHub repo, 2.4 Ponytail: shorter code and shorter replies (+7 more)

## Knowledge Gaps
- **9 isolated node(s):** `TMPDIR`, `Silva 138 reference taxonomy`, `Authors: BakhshiGanje, Jo, Li, Kong, Park, Acimovic (Virginia Tech / Kangwon / Chungbuk)`, `Per-base mean Phred quality by health status`, `Acute oak decline polymicrobial canker (Brenneria goodwinii)` (+4 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 113 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **11 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Fire Blight Canker Microbiome: Insights into the Key Ecological Niche for Erwinia amylovora (Phytobiomes J 2026)` connect `Study Design & Sequencing` to `Diversity & Pathobiome Findings`?**
  _High betweenness centrality (0.024) - this node is a cross-community bridge._
- **What connects `TMPDIR`, `Silva 138 reference taxonomy`, `Authors: BakhshiGanje, Jo, Li, Kong, Park, Acimovic (Virginia Tech / Kangwon / Chungbuk)` to the rest of the system?**
  _9 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Stats & ML Analysis (04_stats.R)` be split into smaller, more focused modules?**
  _Cohesion score 0.0392156862745098 - nodes in this community are weakly interconnected._
- **Why does `ASV table (1,817,968 reads, 2,322 ASVs)` connect `Diversity & Pathobiome Findings` to `Study Design & Sequencing`?**
  _High betweenness centrality (0.022) - this node is a cross-community bridge._
- **Should `Study Design & Sequencing` be split into smaller, more focused modules?**
  _Cohesion score 0.13157894736842105 - nodes in this community are weakly interconnected._
- **Should `Diversity & Pathobiome Findings` be split into smaller, more focused modules?**
  _Cohesion score 0.1323529411764706 - nodes in this community are weakly interconnected._
- **Should `Co-occurrence Network (05_network.R)` be split into smaller, more focused modules?**
  _Cohesion score 0.13333333333333333 - nodes in this community are weakly interconnected._