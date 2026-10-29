#!/usr/bin/env python3
"""Step 4b: DADA2 read-tracking and ASV-count report (PDF), with comparison to BakhshiGanje et al. 2026."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
B = "/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code"
E = f"{B}/results/04_dada2/export"
meta = pd.read_csv(f"{B}/raw_data/sample_metadata.tsv", sep="\t", skiprows=[1]).set_index("sample-id")
st = pd.read_csv(f"{E}/stats.tsv", sep="\t", skiprows=[1], index_col=0)
ft = pd.read_csv(f"{E}/feature-table.tsv", sep="\t", skiprows=1, index_col=0)
st = st.join(meta[["cultivar", "health_status"]])
st["n_ASVs"] = (ft > 0).sum(axis=0).reindex(st.index).fillna(0).astype(int)
st = st.sort_values(["health_status", "cultivar"])
st.to_csv(f"{B}/results/04_dada2/denoising_summary.tsv", sep="\t")
stages = [c for c in ["input", "primer-removed", "filtered", "denoised", "non-chimeric"] if c in st.columns]

tot = st["non-chimeric"].sum(); nasv = ft.shape[0]
grp = st.groupby("health_status")["n_ASVs"].agg(["min", "mean", "max"]).round(0)
cmp = pd.DataFrame({"This re-analysis": [f"{tot:,}", f"{nasv:,}",
        f"{int(grp.loc['Canker','min'])}-{int(grp.loc['Canker','max'])} (mean {int(grp.loc['Canker','mean'])})",
        f"{int(grp.loc['Healthy','min'])}-{int(grp.loc['Healthy','max'])} (mean {int(grp.loc['Healthy','mean'])})"],
    "Paper": ["1,817,968", "2,322", "16-470 (mean 130)", "1-22 (mean 5)"]},
    index=["Non-chimeric reads", "Total ASVs", "ASVs per canker sample", "ASVs per healthy sample"])
cmp.to_csv(f"{B}/results/04_dada2/comparison_with_paper.tsv", sep="\t")
print(cmp.to_string()); print(st[stages].sum())

with PdfPages(f"{B}/pdf/04_dada2_stats.pdf") as pdf:
    fig, ax = plt.subplots(figsize=(10, 3)); ax.axis("off")
    t = ax.table(cellText=cmp.values, colLabels=cmp.columns, rowLabels=cmp.index, loc="center")
    t.auto_set_font_size(False); t.set_fontsize(10); t.scale(1, 1.8)
    ax.set_title("DADA2 (denoise-ccs) output vs. paper", fontsize=12); pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(14, 5)); x = np.arange(len(st))
    for s, c in zip(stages, plt.cm.viridis(np.linspace(0, 0.9, len(stages)))):
        ax.bar(x, st[s], color=c, label=s, width=0.8)
    ax.set_xticks(x); ax.set_xticklabels(st.index, rotation=90, fontsize=7)
    ax.axvline((st.health_status == "Canker").sum() - 0.5, c="k", lw=1)
    ax.set_ylabel("Reads"); ax.set_title("Reads retained at each DADA2 stage (overlaid)"); ax.legend(fontsize=8)
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.bar(x, 100 * st["non-chimeric"] / st["input"], color=["#c0392b" if h == "Canker" else "#2e86c1" for h in st.health_status])
    ax.set_xticks(x); ax.set_xticklabels(st.index, rotation=90, fontsize=7)
    ax.set_ylabel("% of input reads retained"); ax.set_title("Final retention (red = Canker, blue = Healthy)")
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.bar(x, st["n_ASVs"], color=["#c0392b" if h == "Canker" else "#2e86c1" for h in st.health_status])
    ax.set_xticks(x); ax.set_xticklabels(st.index, rotation=90, fontsize=7)
    ax.set_ylabel("Observed ASVs"); ax.set_title("ASVs per sample (paper: canker 16-470, healthy 1-22)")
    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
