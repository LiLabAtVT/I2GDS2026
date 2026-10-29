#!/usr/bin/env python3
"""Step 8b: species calls from BLAST vs NCBI 16S RefSeq.
A species is assigned only if: best hit >= 99.5% identity over >= 95% of the ASV, its genus equals the SILVA genus,
and no hit of a different species is within 0.2 identity points of the best (otherwise 'ambiguous')."""
import pandas as pd, re
B = "/projects/intro2gds/I2GDS2026/Team5/orlando/Claude_Code"; O = f"{B}/results/08b_blast"
MIN_ID, MIN_COV, MARGIN = 99.5, 0.95, 0.2
cols = "qseqid sacc pident length qlen slen mismatch gapopen evalue bitscore stitle".split()
bl = pd.read_csv(f"{O}/blast_all_asvs.tsv", sep="\t", names=cols)
lk = pd.read_csv(f"{B}/results/07_composition_ml/asv_lookup.tsv", sep="\t")   # bacterial ASVs only, ASV0001...
lk = lk.set_index("feature_id")

def species_of(title):
    w = title.replace("[", "").replace("]", "").split()
    return " ".join(w[:3]) if w[0] == "Candidatus" else " ".join(w[:2])

bl = bl[bl.qseqid.isin(lk.index)].copy()
bl["cov"] = bl.length / bl.qlen
bl["species"] = bl.stitle.map(species_of)
bl = bl.sort_values(["qseqid", "bitscore", "pident"], ascending=[True, False, False])

# Per ASV and species: best identity among near-full-length hits (cov >= MIN_COV); compare best species vs the next species
full = bl[bl["cov"] >= MIN_COV]
bysp = full.groupby(["qseqid", "species"]).agg(pident=("pident", "max"), sacc=("sacc", "first")).reset_index()
rows = []
for q in bl.qseqid.unique():
    g = bysp[bysp.qseqid == q].sort_values("pident", ascending=False)
    silva_genus = str(lk.loc[q, "Genus"])
    if g.empty:
        top = bl[bl.qseqid == q].iloc[0]
        rows.append(dict(ASV=lk.loc[q, "ASV"], feature_id=q, silva_genus=silva_genus, total_reads=lk.loc[q, "total_reads"],
                         species=top.species, best_acc=top.sacc, best_pident=top.pident, second_species=None, second_pident=None,
                         margin=None, assigned="no", reason="low_identity_or_coverage"))
        continue
    best = g.iloc[0]; sp = best.species
    second = g.iloc[1] if len(g) > 1 else None
    margin = best.pident - (second.pident if second is not None else 0)
    gen = sp.split()[1] if sp.startswith("Candidatus") else sp.split()[0]
    ok_id = best.pident >= MIN_ID
    ok_gen = gen in silva_genus.split("-") or gen == silva_genus
    ok_margin = margin > MARGIN
    rows.append(dict(ASV=lk.loc[q, "ASV"], feature_id=q, silva_genus=silva_genus, total_reads=lk.loc[q, "total_reads"],
                     species=sp, best_acc=best.sacc, best_pident=best.pident,
                     second_species=None if second is None else second.species,
                     second_pident=None if second is None else second.pident, margin=round(margin, 2),
                     assigned="yes" if (ok_id and ok_gen and ok_margin) else "no",
                     reason="ok" if (ok_id and ok_gen and ok_margin) else
                            ";".join(r for r, c in [("low_identity", not ok_id), ("genus_mismatch", not ok_gen),
                                                    ("ambiguous", not ok_margin)] if c)))
res = pd.DataFrame(rows).sort_values("ASV")
nohit = lk.index.difference(bl.qseqid.unique())
res.to_csv(f"{O}/species_assignments.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20)
print(f"Bacterial ASVs: {len(lk)}; with BLAST hit >=90%: {len(res)}; no hit: {len(nohit)}")
print(f"Species assigned: {(res.assigned == 'yes').sum()} ASVs "
      f"({res.loc[res.assigned == 'yes', 'total_reads'].sum() / lk.total_reads.sum() * 100:.1f}% of bacterial reads)")
print(res.reason.value_counts().to_string())

show = ["ASV", "total_reads", "species", "best_pident", "second_species", "second_pident", "margin", "assigned"]
erw = res[res.silva_genus == "Erwinia"]
print("\n== Erwinia ASVs =="); print(erw[show].to_string(index=False))
# top hits per species for Erwinia ASVs (how close are E. amylovora vs relatives?)
eb = full[full.qseqid.isin(erw.feature_id)].groupby(["qseqid", "species"]).pident.max().reset_index()
eb["ASV"] = eb.qseqid.map(lk.ASV)
piv = eb.pivot_table(index="ASV", columns="species", values="pident")
piv = piv[piv.max().sort_values(ascending=False).index[:8]]
piv.round(2).to_csv(f"{O}/erwinia_asv_identity_by_species.tsv", sep="\t")
print("\n== Erwinia ASVs: best identity to each of the closest species =="); print(piv.round(2).to_string())

drv = pd.read_csv(f"{B}/results/07_composition_ml/dmm_k4_top10_drivers.tsv", sep="\t").ASV.unique()
print("\n== DMM drivers (Step 7) =="); print(res[res.ASV.isin(drv)][show + ["silva_genus"]].to_string(index=False))
