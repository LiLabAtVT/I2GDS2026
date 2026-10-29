#!/usr/bin/env python3
"""Step 6: diversity report (PDF) - alpha boxplots (as paper Fig. 2B), Bray-Curtis PCoA (Fig. 3), PERMANOVA vs paper."""
import glob
import pandas as pd, numpy as np
from scipy import stats
from skbio import OrdinationResults
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
B = "/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code"; O = f"{B}/results/06_diversity"
meta = pd.read_csv(f"{B}/raw_data/sample_metadata.tsv", sep="\t", skiprows=[1]).set_index("sample-id")
CULT = ["Cripps Pink", "Aztec Fuji", "Gala", "Hewe's Crab", "Granniwinkle", "Unknown cultivar"]
COL = dict(zip(CULT, ["#c9a227", "#e8736c", "#2ca02c", "#5b8def", "#17becf", "#e377c2"]))
METRICS = [("simpson", "Simpson"), ("shannon", "Shannon"), ("observed_features", "Observed"),
           ("chao1", "Chao1"), ("faith_pd", "Faith PD")]

def alpha(label):
    C = f"{O}/core-{label}/export"
    df = pd.concat([pd.read_csv(f"{C}/{m}/alpha-diversity.tsv", sep="\t", index_col=0).iloc[:, 0].rename(m)
                    for m, _ in METRICS], axis=1)
    return df.join(meta[["cultivar", "health_status"]])

def adonis(label):
    f = glob.glob(f"{O}/core-{label}/export/adonis/*.tsv")[0]
    return pd.read_csv(f, sep="\t", index_col=0)

def alpha_tests(df, by):
    rows = []
    for m, name in METRICS:
        groups = [g[m].values for _, g in df.groupby(by)]
        F, pa = stats.f_oneway(*groups); H, pk = stats.kruskal(*groups)
        rows.append(dict(metric=name, ANOVA_F=round(F, 2), ANOVA_P=f"{pa:.2g}", KruskalWallis_H=round(H, 2), KW_P=f"{pk:.2g}"))
    return pd.DataFrame(rows).set_index("metric")

def table_page(pdf, frames, title, h=6):
    fig, axs = plt.subplots(len(frames), 1, figsize=(11, h)); axs = np.atleast_1d(axs)
    for ax, (sub, tbl) in zip(axs, frames):
        ax.axis("off")
        t = ax.table(cellText=tbl.values, colLabels=tbl.columns, rowLabels=tbl.index, loc="center")
        t.auto_set_font_size(False); t.set_fontsize(8.5); t.scale(1, 1.4); ax.set_title(sub, fontsize=10)
    fig.suptitle(title, fontsize=12); pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

def boxplots(pdf, df, title, statuses):
    fig, axs = plt.subplots(1, len(METRICS), figsize=(4 * len(METRICS), 4.8))
    for ax, (m, name) in zip(axs, METRICS):
        pos, ticks, labs = 0, [], []
        for st in statuses:
            for c in CULT:
                v = df[(df.health_status == st) & (df.cultivar == c)][m]
                if len(v) == 0: continue
                ax.boxplot(v, positions=[pos], widths=0.6, showfliers=False,
                           medianprops=dict(color=COL[c]), boxprops=dict(color=COL[c]), whiskerprops=dict(color=COL[c]))
                ax.scatter(np.full(len(v), pos) + np.random.uniform(-0.15, 0.15, len(v)), v, s=10, color=COL[c], zorder=3)
                ticks.append(pos); labs.append(f"{c} ({st[0]})"); pos += 1
            pos += 0.8
        ax.set_xticks(ticks); ax.set_xticklabels(labs, rotation=90, fontsize=7); ax.set_title(name)
    fig.suptitle(title + "   (C = canker, H = healthy)"); fig.tight_layout(); pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

def pcoa(pdf, label, title):
    o = OrdinationResults.read(f"{O}/core-{label}/export/bray_curtis_pcoa/ordination.txt")
    s = o.samples.iloc[:, :2].join(meta[["cultivar", "health_status"]]); ev = o.proportion_explained.values[:2] * 100
    fig, ax = plt.subplots(figsize=(8, 6.5))
    for (c, st), g in s.groupby(["cultivar", "health_status"]):
        ax.scatter(g.iloc[:, 0], g.iloc[:, 1], color=COL[c], marker="o" if st == "Canker" else "^", s=45,
                   edgecolor="k", lw=0.4, label=f"{c} - {st}")
    ax.set_xlabel(f"PC1 [{ev[0]:.1f}%]"); ax.set_ylabel(f"PC2 [{ev[1]:.1f}%]"); ax.set_title(title)
    ax.legend(fontsize=7, bbox_to_anchor=(1.01, 1), loc="upper left"); pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
    return ev

aA, aB = alpha("A"), alpha("B")
aA.to_csv(f"{O}/alpha_A_all_samples.tsv", sep="\t"); aB.to_csv(f"{O}/alpha_B_canker_bacteria.tsv", sep="\t")
adA, adB = adonis("A"), adonis("B")
summ = aA.groupby("health_status")[[m for m, _ in METRICS]].agg(["min", "mean", "max"]).round(2)

paper = pd.DataFrame({"Paper R2": ["17.2%", "n.s.", "12.4%"], "Paper F": ["9.35", "-", "1.43"], "Paper P": ["0.001", "0.223", "0.016"]},
                     index=["health_status", "cultivar", "health_status:cultivar"])
rs = [c for c in adA.columns if c.lower() in ("r2", "r.squared")][0]
fs = [c for c in adA.columns if c.lower() in ("f", "f.model")][0]
ps = [c for c in adA.columns if c.lower().startswith("pr(") or c.lower() in ("p", "pr..f.")][0]
cmpA = pd.DataFrame({"Our R2": (adA[rs] * 100).round(1).astype(str) + "%", "Our F": adA[fs].round(2), "Our P": adA[ps]})
cmpA = cmpA.join(paper, how="left").fillna("")
cmpA.to_csv(f"{O}/permanova_A_vs_paper.tsv", sep="\t"); adB.to_csv(f"{O}/permanova_B_canker_cultivar.tsv", sep="\t")
tA = alpha_tests(aA, "health_status"); tAc = alpha_tests(aA[aA.health_status == "Canker"], "cultivar")
tAh = alpha_tests(aA[aA.health_status == "Healthy"], "cultivar"); tB = alpha_tests(aB, "cultivar")
for n, t in [("A_health", tA), ("A_canker_cultivar", tAc), ("A_healthy_cultivar", tAh), ("B_canker_cultivar", tB)]:
    t.to_csv(f"{O}/alpha_tests_{n}.tsv", sep="\t")
print(cmpA.to_string()); print(adB.to_string()); print(summ.to_string()); print(tA.to_string()); print(tAc.to_string()); print(tB.to_string())

with PdfPages(f"{B}/pdf/06_diversity.pdf") as pdf:
    table_page(pdf, [("A: PERMANOVA (adonis, Bray-Curtis, 999 perm.) - our result vs paper", cmpA),
                     ("B: PERMANOVA, canker bacteria only (cultivar)", adB.round(4).astype(str))],
               "Step 6 - beta diversity statistics", h=6)
    table_page(pdf, [("A: canker vs healthy", tA), ("A: among cultivars - canker", tAc),
                     ("A: among cultivars - healthy", tAh), ("B: among cultivars - canker, bacteria only", tB)],
               "Step 6 - alpha diversity tests (one-way ANOVA as in paper; Kruskal-Wallis)", h=10)
    boxplots(pdf, aA, "A: alpha diversity, all 48 samples, rarefied to 14,881 (chloroplast kept, as paper)", ["Canker", "Healthy"])
    evA = pcoa(pdf, "A", "A: Bray-Curtis PCoA, all samples (paper Fig. 3: PC1 23.7%, PC2 6.9%)")
    boxplots(pdf, aB, "B: alpha diversity, canker bacteria only, rarefied to 505", ["Canker"])
    pcoa(pdf, "B", "B: Bray-Curtis PCoA, canker samples, bacteria only")
    table_page(pdf, [("A: alpha diversity range by health status (paper: canker Observed 16-470; healthy Shannon 0-2, Simpson 0-0.75)",
                      summ.T.round(2))], "Step 6 - alpha diversity summary", h=6)
print("PCoA A explained:", evA)
