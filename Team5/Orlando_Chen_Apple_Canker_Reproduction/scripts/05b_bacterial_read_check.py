#!/usr/bin/env python3
"""Step 5b: per-sample bacterial reads/ASVs after removing chloroplast/mitochondria/unassigned, and class-level composition."""
import pandas as pd
B = "/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code"; O = f"{B}/results/05_taxonomy"
meta = pd.read_csv(f"{B}/raw_data/sample_metadata.tsv", sep="\t", skiprows=[1]).set_index("sample-id")
ft = pd.read_csv(f"{B}/results/04_dada2/export/feature-table.tsv", sep="\t", skiprows=1, index_col=0)
lv = pd.read_csv(f"{O}/taxonomy_by_rank.tsv", sep="\t", index_col=0).reindex(ft.index)
bact = (lv.Domain == "Bacteria") & (lv.Order != "Chloroplast") & (lv.Family != "Mitochondria")
b = ft[bact]
s = pd.DataFrame({"total_reads": ft.sum(), "bacterial_reads": b.sum(), "bacterial_ASVs": (b > 0).sum(),
                  "all_ASVs": (ft > 0).sum()}).join(meta[["cultivar", "health_status"]]).sort_values(["health_status", "cultivar"])
s.to_csv(f"{O}/bacterial_reads_per_sample.tsv", sep="\t")
pd.set_option("display.width", 200)
print(s.to_string())
print(s.groupby("health_status")[["bacterial_reads", "bacterial_ASVs"]].describe().round(0).T)
cls = b.groupby(lv.loc[b.index, "Class"]).sum(); rel = cls / cls.sum() * 100
c = rel.T.join(meta.health_status).groupby("health_status").mean().T
c.round(2).to_csv(f"{O}/class_mean_by_health_status.tsv", sep="\t")
print(c.sort_values("Canker", ascending=False).head(10).round(1))
print(c.sort_values("Healthy", ascending=False).head(10).round(1))
g = pd.read_csv(f"{O}/genus_mean_by_health_status.tsv", sep="\t", index_col=0)
print(g.sort_values("Healthy", ascending=False).head(10).round(2))
