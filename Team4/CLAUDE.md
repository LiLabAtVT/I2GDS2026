# CLAUDE.md — Botrytis cinerea mycovirome re-analysis + VHAMSTeR host prediction

## 1. Project goal

Re-run (as closely as practical) the metatranscriptomic mycovirus discovery pipeline of

> Ruiz-Padilla A, Rodríguez-Romero J, Gómez-Cid I, Pacifico D, Ayllón MA. 2021. *Novel mycoviruses discovered in the mycovirome of a necrotrophic fungus.* mBio 12:e03705-20. https://doi.org/10.1128/mBio.03705-20

then add a second layer: run **VHAMSTeR** (Grigson, Neri & Roux 2026, bioRxiv, doi 10.64898/2026.10.23.761443; code at github.com/susiegriggo/VHAMSTeR) on the recovered viral sequences and test whether its host assignments agree with the known biology (all sequences come from pure *B. cinerea* cultures, so the expected host class is **fungi**).

Outputs: reproduced contig/virus tables, presence/absence + read-count matrices, phylogenies, VHAMSTeR host calls per segment, and a concordance report with figures.

## 2. Ground rules for Claude Code

- Never fabricate results, accession numbers, tool flags or file contents. If unsure, run `--help`, read the README, or ask.
- **Do not guess VHAMSTeR CLI flags.** Run `vhamster --help` (and `vhamster-features --help`) and read the installed package README first; write the verified command into `docs/vhamster_commands.md`.
- Never print, log, echo or commit `HF_TOKEN` or any credential. Read it from the environment only. Add `.env`, `*.token`, `models/` to `.gitignore`.
- No heavy compute on login nodes. Everything beyond file inspection and tiny tests goes through Slurm (`sbatch`). Ask before submitting any job that requests > 8 h, > 128 GB or a GPU.
- Keep every deliverable. Write final tables/figures to `results/` (never only to `scratch/`), and write an md5 manifest (`results/MANIFEST.md5`) after each stage. Copy key outputs off scratch when a stage finishes.
- Log everything: each Slurm script writes `logs/<stage>_<jobid>.out/.err`; each stage appends tool versions and parameters to `docs/run_manifest.tsv` (date, stage, tool, version, command, input checksum).
- Record every deviation from the paper in `docs/deviations.md` (e.g., newer tool versions, different database snapshot).
- Make scripts idempotent and resumable (skip if output exists and checksum matches).
- Work stage by stage. At the end of each stage, summarize counts versus the paper's expected values (Section 5) and stop for review before the next stage.

## 3. Environment (fill in once)

- Cluster: Slurm. Account `<ACCOUNT>`, CPU partition `<CPU_PARTITION>`, GPU partition `<GPU_PARTITION>` (confirm with `sinfo`, `sacctmgr show assoc user=$USER`).
- Project root: `<PROJECT_ROOT>`; bulk data/scratch: `<SCRATCH_ROOT>`.
- Check whether compute nodes have internet. If not, do all downloads (SRA, GenBank, Zenodo, HuggingFace models) on a login/data-transfer node, then run the rest offline.
- Use separate conda/mamba envs (pin versions with `conda env export --no-builds > envs/<name>.yml`):
  - `mycoviro`: sra-tools, bbtools, trinity, blast, diamond, bwa, samtools, seqkit, mmseqs2 or cd-hit, mafft, clustalo, muscle, iqtree, trimal, python (pandas, biopython, matplotlib, seaborn)
  - `vhamster`: vhamster (see Section 7). Keep it separate because of PyTorch/CUDA pins.

### Directory layout

```
project/
  CLAUDE.md
  config/            pools.tsv, samples.tsv, paths.env (no secrets)
  envs/              exported conda envs
  scripts/           NN_stage_name.sh|.slurm|.py  (numbered, one stage each)
  data/raw/          FASTQ (scratch), data/ref/ databases, data/supp/ paper supplements
  work/              intermediates (scratch)
  results/           final tables, figures, MANIFEST.md5
  logs/
  docs/              run_manifest.tsv, deviations.md, vhamster_commands.md, report.md
```

## 4. Data

- Reads: BioProject **PRJNA632510**; BioSamples SAMN14911182–SAMN14911210; SRA experiments SRX8335942–SRX8335970 (29 libraries = 29 pools).
- Pools: 17 Spanish (BCS1–BCS17) and 12 Italian (BCI1–BCI12); 248 isolates total (150 Spain, 98 Italy), sampled 2018 from grapevine.
- Sequencing: Illumina NovaSeq 6000, paired-end 150 bp, > 100 M read pairs per library (~2,696 M reads total, ~93 M average).
- **Do not assume SRR/SRX order equals pool order.** Build `config/pools.tsv` (SRX/SRR → BioSample → pool name BCSn/BCIn) from the SRA run table / BioSample titles, and ask the user to confirm any pool that cannot be resolved unambiguously.
- Paper supplements (mBio online only, not in the PDF) are needed: Table S1 (reads per virus per pool), Table S2 (primers), Fig S1 (pool/region/sample counts). Ask the user to place them in `data/supp/`. Reference sequences: the viral accessions in Table 1 of the paper (MN6xxxxx, MN9xxxxx, MT1xxxxx etc.) — fetch from GenBank into `data/ref/paper_viruses.fasta` and use them as the ground-truth set.
- Wet-lab work (qPCR/RT-PCR validation, RACE, single-spore isolates) cannot be reproduced computationally. Skip it, but use Table 2/3 as annotation of which viruses were experimentally validated.

## 5. Expected values from the paper (for sanity checks, not exact targets)

| Item | Paper value |
|---|---|
| Pools / isolates | 29 pools (17 BCS, 12 BCI) / 248 isolates |
| Contigs > 1,000 bp assigned to mycoviruses | 1,269 (670 BCS, 599 BCI) |
| After 90% nt identity / 90% length dereplication | 158 representatives |
| Unique sequences with complete CDS | 109 (79 Spain, 30 Italy); 19 already in databases |
| Mycoviruses | 92 total: 57 ssRNA(+), 19 dsRNA, 15 ssRNA(−), 1 ssDNA |
| New species | 62 |
| Geographic split | 55 shared, 32 Spain only, 5 Italy only |
| ssRNA(+) breakdown | 11 mitoviruses, 5 narnaviruses (BcNV4 + BcBNV1/2/3/5), 19 ourmia-like (17 new BcOLV1–17), 3 umbra-like, 5 fusariviruses, 6 hypoviruses, 2 endornaviruses, 2 deltaflexiviruses, 1 mycoflexivirus, 1 flexivirus, 1 alpha-like, 1 tymo-like |
| dsRNA breakdown | 2 botybirnaviruses, 1 quadrivirus, 11 victoriviruses, 2 partitiviruses, 3 unclassified |
| Presence threshold | > 500 reads mapped per virus per pool (log10 = 2.7 in Figs 1–2) |

Because the original used an NCBI nr snapshot from October 2019 and today's nr is much larger, expect **taxonomic labels and first-hit identities to drift**. The core quantity to reproduce is the set of RdRp-bearing viral contigs, not exact BLAST hit names.

### Known inconsistencies inside the paper (do not "fix" silently; note in `docs/deviations.md`)

- BcssDV1 length: text says 1,426 nt, Table 1/Fig 12 say 1,694 nt (the 1,694-nt sequence with a 381-nt putative intron is the one deposited, MN625247).
- BcFlV1 protein: 4,619 aa in Table 1 vs 3,986 aa in text; BcUV1 RdRp: 519 aa (Table 1) vs 619 aa (text).
- Names vary: BcAV1 vs BcALV1; BcBV1 sometimes typed "BcBYV1"; "BCI16" for PoOLV2 variant where the pool is probably BCI6; BcNSRV7 variants listed with BCI7/BCS13 isolate labels.
- "Quadrivirus" is spelled Quadriviridae/Quadraviridae in different places.

## 6. Workflow (one numbered script per stage)

**Stage 00 — setup.** Create envs, directory layout, `config/paths.env`, `.gitignore`. Confirm Slurm account/partitions. Check disk quota: raw FASTQ for 29 libraries is several hundred GB.

**Stage 01 — download.** `prefetch` + `fasterq-dump --split-files` (or ENA FASTQ) for the 29 runs; gzip; md5 into `results/MANIFEST.md5`; `seqkit stats` into `results/01_read_stats.tsv`. Verify read counts are in the expected range (~93 M pairs average).

**Stage 02 — cleaning (BBTools).** Remove adapters, artifacts/PhiX, short reads, and ribosomal sequences (bbduk with an rRNA reference, e.g., SILVA). Record reads before/after per pool. Paper cites BBTools without parameters; use defaults documented in `docs/deviations.md`.

**Stage 03 — de novo assembly.** Trinity per pool (paper used v2.3.2; a newer version is acceptable if recorded). Trinity on ~100 M pairs is memory- and time-heavy: test on one small pool or a subsample (`seqtk sample`/`reformat.sh samplereads`) first, then submit one array job per pool with `--max_memory` and `--CPU` matched to the Slurm request. Keep only `Trinity.fasta` and the gene-trans map; delete Trinity's intermediate directories after checking outputs.

**Stage 04 — viral contig identification.**
1. BLASTn against a viral nucleotide database (paper: custom viral DB, NCBI BLAST 2.8; record whichever DB you use, e.g., RefSeq viral).
2. DIAMOND blastx of candidates (and, separately, all contigs for MEGAN-style taxonomic profiling) against NCBI nr; paper used the October 2019 snapshot, so record the date of whatever nr you download.
3. Keep contigs > 1,000 nt, manually inspect borderline ones (identity, alignment length, query length), and separate real viruses from endogenized viral elements/host sequence (mapping to the *B. cinerea* genome helps).
4. Export a table of per-contig best hits, genome type guess (dsRNA, ssRNA+, ssRNA−, ssDNA), and pool of origin.

**Stage 05 — dereplication and curation.** Within each pool, cluster at 90% nt identity over 90% length and keep the longest contig per cluster (mmseqs2 `easy-cluster` or cd-hit-est, document settings). Then pick sequences with complete coding sequences, predict ORFs (`getorf`/`orfipy`; use the right genetic code, mitoviruses use the mitochondrial code), and group variants (> 95% aa identity and identical genome organization). Pair multi-segment viruses by conserved 5′/3′ termini and co-occurrence across pools (the paper's bisegmented binarnaviruses BcBNV1/2/3/5 and the trisegmented BcBV1 are the key test cases; BcMyV4 is tetrasegmented). Output `results/05_curated_viruses.tsv` and `results/05_curated_viruses.fasta` with a `virus_id`, `segment`, `pool`, `genome_type` column set.

**Stage 06 — abundance and presence/absence.** Map cleaned reads of each pool to the curated set (bwa mem; samtools; the paper reports bwa). Count mapped reads per virus per pool; apply the > 500-read presence threshold. Reproduce the Fig 1 and Fig 2 style bubble charts (log10 reads) and the geographic summary (55 shared / 32 Spain / 5 Italy).

**Stage 07 — taxonomy and phylogeny.** For each virus family/clade, extract RdRp proteins, add reference sequences from the paper's trees (GenBank accessions are printed in Figs 4, 5, 7, 9, 11, 12), align (Clustal Omega, MUSCLE or MAFFT), trim, and run IQ-TREE (paper: VT+F+I+G4 / LG+F+I+G4 / WAG+F+I+G5 depending on the tree, 1,000 bootstrap replicates). Standard bootstrap is slow; ultrafast bootstrap (`-B 1000`) is an acceptable documented deviation. Check that the binarnavirus clade (BcBNV1/2/3/5) is separate from true narnaviruses and Matryoshka RNA virus 1.

**Stage 08 — VHAMSTeR host prediction.** See Section 7.

**Stage 09 — integration, figures, report.** Merge Stage 05–08 into `results/master_table.tsv` (one row per virus segment: pool, genome type, taxonomy, reads, validated?, VHAMSTeR host, probability/score). Write `docs/report.md` with methods, deviations, reproduction-fidelity table (observed versus Section 5), and VHAMSTeR concordance analysis.

## 7. VHAMSTeR integration

### What it is (as documented at time of writing)

- Classifies viral sequences into five host categories: **prokaryote, animal, plant, fungi, protist**.
- Combines fine-tuned genomic language-model embeddings (NTv3 base model) with explicit features for genome architecture and protein homology; uses mmseqs2 and a geNomad marker database; 5-fold ensemble (`fold_0`–`fold_4`).
- Preprint reports MCC > 0.91 for DNA viruses and > 0.69 for RNA viruses. It is a **new, not yet peer-reviewed preprint** (posted Oct 2026): treat predictions as hypotheses and pin the version.
- Python package `vhamster` (PyPI 1.4.1 at time of writing; GitHub release v1.4.0) with commands `vhamster`, `vhamster-install-models`, `vhamster-features`.

### Install (on a login node with internet)

```bash
mamba create -n vhamster -c conda-forge -c bioconda vhamster      # preferred, pulls mmseqs2
# or: pip install vhamster  (then install mmseqs2 separately)
# GPU optional: if using the repo's environment.yml, edit the pytorch-cuda version to match the cluster.
conda activate vhamster
export HF_TOKEN=...            # headless HPC: NTv3 needs a HuggingFace licence acceptance; use env var or `huggingface-cli login`
vhamster-install-models -o <SCRATCH_ROOT>/vhamster_models      # weights + geNomad DB + NTv3
```

Ask the user to accept the NTv3 licence on HuggingFace and supply the token interactively; never write it into a file. Record `pip show vhamster` and the Zenodo/HF model versions in `docs/run_manifest.tsv`. Smoke-test on the example data shipped with the repo before running the real inputs.

### Inputs

Build three nucleotide FASTA sets from the Stage 05 output, **one record per genome segment**, headers `virusID|segment|pool`:

1. `paper_curated_109.fasta` — the 109 deposited GenBank sequences from Table 1 (ground-truth set).
2. `my_curated.fasta` — your own Stage 05 contigs.
3. `controls.fasta` — reference viruses with known hosts, taken from accessions in the paper's trees: Ourmia melon virus EU770623 (plant), Sindbis virus J02363 (animal), Enterobacteria phage MS2 V00642 (prokaryote), Sclerotinia sclerotiorum hypovirulence-associated DNA virus 1 GQ365709 (fungi), Bemisia-associated genomovirus NfO KY230625 (animal host expected), Leishmania RNA virus 1-1 M92355 (protist). Verify each accession and host by fetching the GenBank record before using it as a control.

Notes: mitoviruses are best treated by their RNA sequence as deposited; the ssRNA(−) genomes are deposited as the genome/antigenome orientation given in GenBank, so keep the sequence orientation as deposited and note it in the table; for dsRNA and ssRNA genomes use cDNA sequence as deposited. Run `vhamster-features` separately only if needed to inspect features.

### Run

Use the commands verified from `vhamster --help` (Section 2). Submit as a Slurm job (GPU if available, otherwise CPU; record runtime), output to `results/08_vhamster/<set>/`. Keep raw outputs untouched; parse them with `scripts/08_parse_vhamster.py` into `results/08_vhamster_calls.tsv` (virus, segment, predicted host, per-class scores/probabilities, any geNomad/homology flags the tool reports).

### Analysis questions

1. **Accuracy against expectation.** What fraction of the 92 viruses (and of all 109 sequences) are called *fungi*? Break down by genome type (ssRNA+, dsRNA, ssRNA−, ssDNA) and by family. Report MCC/confusion tables only for the control set (the *B. cinerea* viruses are a one-class test, so report recall of "fungi" rather than MCC).
2. **Informative disagreements.** List every non-fungal call and compare with the paper's own cross-kingdom/horizontal-transfer discussion. Pay particular attention to:
   - BcBV1 (trisegmented bocivirus; RdRp and NCP cluster with plant coguvirus-like viruses; RNA2 movement-protein-like ORF) — the paper's strongest cross-kingdom hypothesis.
   - BcssDV1 (genomovirus; nearest non-fungal neighbor is Bemisia-associated genomovirus NfO).
   - Ourmia-like/botourmia, umbra-like, alpha-like, tymo-like and deltaflexi-like viruses (relatives of plant or animal viruses).
   - Binarnaviruses BcBNV1/2/3/5 and the RdRp-less "HP" second segments (protein homology features are weak here; compare segment 1 versus segment 2 calls) and BcNV4 (classic narnavirus).
   - ssRNA(−) bunya-like viruses (BcNSRV2, 6, 8–11, BcOBV1).
3. **Segment consistency.** For multi-segment viruses, do all segments get the same host? Flag discordance.
4. **Training-set leakage check.** Many of these sequences were deposited in 2020–2021 and may be in VHAMSTeR's training or benchmark data. Download the authors' supplementary code/data (Zenodo record 23670375, repo susiegriggo/VHAMSTeR) and test by accession whether the paper's MN/MT accessions are in the training/benchmark lists. Report performance separately for seen versus unseen sequences. If many are seen, say so prominently: agreement then says little about generalization.
5. **Score calibration.** Check how confidence scores relate to genome completeness (full-length versus partial contigs) and sequence length.
6. **Optional extension.** Run VHAMSTeR on all viral contigs > 1 kb from Stage 04 that were *excluded* during curation (partial, dereplicated, low-confidence) to see whether any non-fungal calls point to plant/animal/prokaryotic contaminants or carry-over viruses (e.g., from grapevine, mites, or bacteria in field-derived isolates).

### Interpretation rules

- A predicted host class is a statistical inference from sequence, not proof of infection. Cultures here are axenic *B. cinerea* isolates, so "fungi" is the biologically expected label; a plant/animal/protist call on a sequence recovered from such a culture is a hypothesis about evolutionary origin or relatedness, not evidence that the virus infects that host.
- Do not present agreement with VHAMSTeR as validation of the paper's novelty claims (new genera/families such as Binarnaviridae, Mybuviridae, Bocivirus). Those rest on phylogeny of RdRp proteins.
- Report deviations from the preprint's benchmark conditions (input completeness, orientation, segment-wise versus genome-wise input).

## 8. Definition of done

- `results/` contains: read stats, assembly stats, curated virus table/FASTA, presence/absence and read-count matrices, bubble-chart figures, phylogenies (treefiles + PDF/SVG), `08_vhamster_calls.tsv`, `master_table.tsv`, and `MANIFEST.md5`.
- `docs/report.md` states: what reproduced (counts versus Section 5), what did not and why, all deviations, VHAMSTeR concordance with leakage-stratified results, and limitations (pooled samples, database drift, new preprint, no wet-lab validation).
- All scripts run from a clean clone via `scripts/` in numbered order, with versions recorded in `docs/run_manifest.tsv`.
