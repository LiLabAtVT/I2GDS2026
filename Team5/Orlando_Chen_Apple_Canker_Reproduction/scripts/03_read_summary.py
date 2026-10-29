#!/usr/bin/env python3
"""Step 3b: Per-sample read summary + QC plots (PDF) from raw PacBio HiFi FASTQs and FastQC output.
Run inside the QIIME2 2024.10 container (has numpy/pandas/matplotlib)."""
import gzip, glob, os, zipfile
from multiprocessing import Pool
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

B = "/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code"
FQ, OUT, PDF = f"{B}/raw_data/fastq", f"{B}/results/03_fastqc", f"{B}/pdf/03_read_QC.pdf"
meta = pd.read_csv(f"{B}/raw_data/sample_metadata.tsv", sep="\t", comment=None, skiprows=[1]).set_index("sample-id")
CV_COL = {"Aztec Fuji":"#1b9e77","Cripps Pink":"#d95f02","Gala":"#7570b3","Granniwinkle":"#e7298a",
          "Hewe's Crab":"#66a61e","Unknown cultivar":"#a6761d"}

def scan(path):
    sid = os.path.basename(path).replace(".fastq.gz", "")
    L, Q = [], []
    with gzip.open(path, "rb") as f:
        for i, line in enumerate(f):
            k = i % 4
            if k == 1: L.append(len(line) - 1)
            elif k == 3: Q.append(np.frombuffer(line.rstrip(), dtype=np.uint8).mean() - 33)
    return sid, np.array(L), np.array(Q)

if __name__ == "__main__":
    with Pool(int(os.environ.get("SLURM_CPUS_PER_TASK", 4))) as p:
        res = p.map(scan, sorted(glob.glob(f"{FQ}/*.fastq.gz")))
    rows, lens = [], {}
    for sid, L, Q in res:
        lens[sid] = L
        rows.append(dict(sample_id=sid, reads=len(L), mean_len=L.mean(), median_len=np.median(L),
                         min_len=L.min(), max_len=L.max(), pct_1000_1600=100*np.mean((L>=1000)&(L<=1600)),
                         pct_gt1600=100*np.mean(L>1600), mean_readQ=Q.mean(), pct_readQ_ge30=100*np.mean(Q>=30)))
    df = pd.DataFrame(rows).set_index("sample_id").join(meta[["cultivar","health_status","platform"]])
    df = df.sort_values(["health_status","cultivar","sample_id"])
    df.round(2).to_csv(f"{OUT}/read_summary.tsv", sep="\t")

    # FastQC module PASS/WARN/FAIL
    fq = {}
    for z in sorted(glob.glob(f"{OUT}/*_fastqc.zip")):
        sid = os.path.basename(z).replace("_fastqc.zip", "")
        with zipfile.ZipFile(z) as zz:
            txt = zz.read(f"{sid}_fastqc/summary.txt").decode()
        fq[sid] = {l.split("\t")[1]: l.split("\t")[0] for l in txt.strip().split("\n")}
    fqdf = pd.DataFrame(fq).T.reindex(df.index)
    fqdf.to_csv(f"{OUT}/fastqc_module_summary.tsv", sep="\t")

    with PdfPages(PDF) as pdf:
        # 1. read-length histogram, canker vs healthy
        fig, ax = plt.subplots(figsize=(9, 5))
        bins = np.arange(0, 3001, 20)
        for st, c in [("Canker", "#c0392b"), ("Healthy", "#2e86c1")]:
            allL = np.concatenate([lens[s] for s in df.index[df.health_status == st]])
            ax.hist(allL, bins=bins, histtype="step", lw=1.5, color=c, label=f"{st} (n={len(allL):,} reads)")
        for x in (1000, 1600): ax.axvline(x, ls="--", c="grey", lw=0.8)
        ax.set_yscale("log"); ax.set_xlabel("Read length (bp)"); ax.set_ylabel("Reads (log)")
        ax.set_title("Raw HiFi read-length distribution (dashed: 1,000-1,600 bp)"); ax.legend()
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # 2. reads per sample
        fig, ax = plt.subplots(figsize=(14, 5))
        cols = [CV_COL[c] for c in df.cultivar]
        hatch = ["//" if p == "Sequel IIe" else "" for p in df.platform]
        bars = ax.bar(range(len(df)), df.reads, color=cols, edgecolor="k", lw=0.4)
        for b, h in zip(bars, hatch): b.set_hatch(h)
        ax.set_xticks(range(len(df))); ax.set_xticklabels(df.index, rotation=90, fontsize=7)
        ax.axvline((df.health_status == "Canker").sum() - 0.5, c="k", lw=1)
        ax.set_ylabel("Raw HiFi reads"); ax.set_title("Reads per sample (left: Canker | right: Healthy; hatched = Sequel IIe)")
        ax.legend(handles=[plt.Rectangle((0,0),1,1,color=v) for v in CV_COL.values()], labels=list(CV_COL), fontsize=8, ncol=3)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # 3. per-sample length and quality boxplots
        fig, axs = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
        order = list(df.index)
        bp = axs[0].boxplot([lens[s] for s in order], showfliers=False, patch_artist=True)
        for p, c in zip(bp["boxes"], cols): p.set_facecolor(c)
        axs[0].set_ylabel("Read length (bp)"); axs[0].set_title("Per-sample read length (no outliers)")
        qd = {sid: Q for sid, L, Q in res}
        bp = axs[1].boxplot([qd[s] for s in order], showfliers=False, patch_artist=True)
        for p, c in zip(bp["boxes"], cols): p.set_facecolor(c)
        axs[1].axhline(30, ls="--", c="grey", lw=0.8)
        axs[1].set_ylabel("Mean read Phred Q"); axs[1].set_title("Per-sample mean read quality (dashed: Q30)")
        axs[1].set_xticks(range(1, len(order)+1)); axs[1].set_xticklabels(order, rotation=90, fontsize=7)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # 4. FastQC module status heatmap
        if len(fqdf.dropna(how="all")):
            code = fqdf.replace({"PASS": 0, "WARN": 1, "FAIL": 2}).astype(float)
            fig, ax = plt.subplots(figsize=(8, 14))
            from matplotlib.colors import ListedColormap
            ax.imshow(code.values, cmap=ListedColormap(["#27ae60", "#f1c40f", "#e74c3c"]), vmin=0, vmax=2, aspect="auto")
            ax.set_xticks(range(code.shape[1])); ax.set_xticklabels(code.columns, rotation=90, fontsize=8)
            ax.set_yticks(range(code.shape[0])); ax.set_yticklabels(code.index, fontsize=7)
            ax.set_title("FastQC modules (green PASS / yellow WARN / red FAIL)\nNote: FastQC thresholds are designed for short reads")
            pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # 5. summary table
        fig, ax = plt.subplots(figsize=(14, 13)); ax.axis("off")
        t = df[["cultivar","health_status","platform","reads","mean_len","pct_1000_1600","pct_gt1600","mean_readQ","pct_readQ_ge30"]].round(1)
        tb = ax.table(cellText=t.values, colLabels=t.columns, rowLabels=t.index, loc="center", fontsize=6)
        tb.auto_set_font_size(False); tb.set_fontsize(6); tb.scale(1, 1.1)
        ax.set_title("Per-sample raw read summary")
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
    print(df[["reads","mean_len","pct_1000_1600","pct_gt1600","mean_readQ"]].describe().round(1))
    print(df.groupby("health_status")[["reads","mean_len","mean_readQ"]].mean().round(1))
