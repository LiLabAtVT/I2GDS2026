#!/usr/bin/env python3
"""Step 5: taxonomy report (PDF) - assignment depth, host organelle reads, phylum/genus composition."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
B = "/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code"
O = f"{B}/results/05_taxonomy"
RANKS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
meta = pd.read_csv(f"{B}/raw_data/sample_metadata.tsv", sep="\t", skiprows=[1]).set_index("sample-id")
ft = pd.read_csv(f"{B}/results/04_dada2/export/feature-table.tsv", sep="\t", skiprows=1, index_col=0)
tax = pd.read_csv(f"{O}/export/taxonomy.tsv", sep="\t", index_col=0)

# Split "d__Bacteria; p__...; ..." into ranks; empty or bare "g__" counts as unassigned
lv = tax["Taxon"].str.split(";", expand=True).reindex(columns=range(7))
lv.columns = RANKS
lv = lv.apply(lambda c: c.str.strip().str.replace(r"^[a-z]__", "", regex=True)).replace("", np.nan)
lv["Confidence"] = tax["Confidence"]
lv.to_csv(f"{O}/taxonomy_by_rank.tsv", sep="\t")
ft = ft.loc[lv.index.intersection(ft.index)]

# Assignment depth: % ASVs and % reads classified at each rank
reads = ft.sum(axis=1)
depth = pd.DataFrame({"% ASVs": [lv[r].notna().mean() * 100 for r in RANKS],
                      "% reads": [reads[lv[r].notna()].sum() / reads.sum() * 100 for r in RANKS]}, index=RANKS).round(1)

# Host organelles / non-bacterial
cat = pd.Series("Bacteria", index=lv.index)
cat[lv["Domain"] != "Bacteria"] = "Unassigned / non-bacterial"
cat[lv["Order"].str.lower() == "chloroplast"] = "Chloroplast"
cat[lv["Family"].str.lower() == "mitochondria"] = "Mitochondria"
host = ft.groupby(cat).sum().T
host_pct = host.div(host.sum(axis=1), axis=0) * 100
host_pct = host_pct.join(meta[["cultivar", "health_status"]]).sort_values(["health_status", "cultivar"])
host_pct.to_csv(f"{O}/host_organelle_percent.tsv", sep="\t")
cats = [c for c in ["Bacteria", "Mitochondria", "Chloroplast", "Unassigned / non-bacterial"] if c in host_pct]
asv_cat = cat.value_counts()

# Relative abundance (bacterial reads only) at phylum and genus
bact = ft.loc[cat == "Bacteria"]
rel = bact.div(bact.sum(axis=0), axis=1) * 100
def collapse(rank, top):
    lab = lv.loc[rel.index, rank].fillna("Unassigned")
    g = rel.groupby(lab).sum()
    keep = g.mean(axis=1).sort_values(ascending=False).index[:top]
    g = pd.concat([g.loc[keep], g.drop(keep).sum().to_frame("Other").T])
    return g
order = host_pct.index
phy, gen = collapse("Phylum", 10), collapse("Genus", 15)
phy.to_csv(f"{O}/phylum_relative_abundance.tsv", sep="\t"); gen.to_csv(f"{O}/genus_relative_abundance.tsv", sep="\t")
grp_gen = rel.groupby(lv.loc[rel.index, "Genus"].fillna("Unassigned")).sum().T.join(meta["health_status"]).groupby("health_status").mean().T
grp_gen = grp_gen.sort_values("Canker", ascending=False)
grp_gen.round(3).to_csv(f"{O}/genus_mean_by_health_status.tsv", sep="\t")

print(depth.to_string()); print(asv_cat.to_string())
print(host_pct.groupby("health_status")[cats].mean().round(1).to_string()); print(grp_gen.head(15).round(2).to_string())

def stacked(ax, df, title):
    df = df[order]; bottom = np.zeros(df.shape[1])
    colors = list(plt.cm.tab20(np.linspace(0, 1, 20)))
    for i, (name, row) in enumerate(df.iterrows()):
        ax.bar(range(df.shape[1]), row.values, bottom=bottom, label=name,
               color="lightgrey" if name == "Other" else colors[i % 20], width=0.9)
        bottom += row.values
    ax.set_xticks(range(df.shape[1])); ax.set_xticklabels(df.columns, rotation=90, fontsize=7)
    ax.axvline((meta.loc[order, "health_status"] == "Canker").sum() - 0.5, c="k", lw=1.2)
    ax.set_ylabel("Relative abundance (%)"); ax.set_ylim(0, 100); ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7)

with PdfPages(f"{B}/pdf/05_taxonomy.pdf") as pdf:
    fig, axs = plt.subplots(1, 2, figsize=(12, 3.5))
    for ax, tbl, ttl in [(axs[0], depth, "Assignment depth (SILVA 138, naive Bayes)"),
                         (axs[1], asv_cat.to_frame("ASVs"), "ASVs by category")]:
        ax.axis("off"); t = ax.table(cellText=tbl.values, colLabels=tbl.columns, rowLabels=tbl.index, loc="center")
        t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.5); ax.set_title(ttl)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(14, 5)); bottom = np.zeros(len(host_pct))
    for c, col in zip(cats, ["#4c72b0", "#dd8452", "#55a868", "#8c8c8c"]):
        ax.bar(range(len(host_pct)), host_pct[c], bottom=bottom, color=col, label=c, width=0.9); bottom += host_pct[c].values
    ax.set_xticks(range(len(host_pct))); ax.set_xticklabels(host_pct.index, rotation=90, fontsize=7)
    ax.axvline((host_pct.health_status == "Canker").sum() - 0.5, c="k", lw=1.2)
    ax.set_ylabel("% of reads"); ax.set_title("Bacterial vs host-organelle reads per sample (left: canker, right: healthy)")
    ax.legend(fontsize=8); pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    for df, ttl in [(phy, "Top 10 phyla (bacterial reads only)"), (gen, "Top 15 genera (bacterial reads only)")]:
        fig, ax = plt.subplots(figsize=(15, 5.5)); stacked(ax, df, ttl + " - left: canker, right: healthy")
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    top = grp_gen.head(20)[::-1]
    fig, ax = plt.subplots(figsize=(8, 7)); y = np.arange(len(top))
    ax.barh(y - 0.2, top["Canker"], 0.4, label="Canker", color="#c44e52")
    ax.barh(y + 0.2, top["Healthy"], 0.4, label="Healthy", color="#55a868")
    ax.set_yticks(y); ax.set_yticklabels(top.index, fontsize=8); ax.set_xlabel("Mean relative abundance (%)")
    ax.set_title("Top 20 genera in canker tissue vs healthy bark"); ax.legend()
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
