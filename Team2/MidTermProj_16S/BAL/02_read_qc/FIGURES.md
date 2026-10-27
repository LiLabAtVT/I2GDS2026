# Step 2 figures: read quality control

These five figures (plus two read-count versions at the end) describe the raw PacBio HiFi reads for all 48 samples (32 canker, 16 healthy, 6 cultivars) **before** any filtering or denoising. They answer one question: is the raw data good enough for DADA2 (step 3), and is there anything odd we should know about first?

- "Course demo" means the instructor's re-analysis of the same dataset, presented in class (I2GDS 2026); it is used here as the reference to compare against.
- Colors: canker = Chicago Maroon, healthy = Burnt Orange.
- In the per-sample plots, panels are cultivars.
- Numbers come from `results/02_qc/qc_summary.txt` and `per_sample_summary.tsv`, which `02_read_qc_plots.R` writes when the step is run.
- Scripts: `scripts/02_read_qc.sh`, `02_primer_check.py`, `02_read_qc_plots.R`.

**Bottom line:**
- The data are high quality.
- The read length, quality and depth match the course demo re-analysis.
- About 18% of reads are a non-target ~1,820 bp product. DADA2's length window will remove these.

---

## 1. `02_read_length_distribution.pdf`: read length

**What it shows:** a histogram of the length of every read, drawn as a density so the canker and healthy groups are comparable even though canker samples have more reads. The grey band is the 1,000–1,600 bp window that DADA2 keeps in step 3.

**What to look for:**
- **The main peaks sit at 1,450–1,500 bp.** This is the expected size of a full-length 16S amplicon (27F to 1492R), so the libraries contain what they should.
- **There are two main peaks, and they differ by group.**
  - Healthy reads pile up almost entirely at **~1,451 bp** (median 1,451 bp in every healthy sample).
  - Canker reads have a second peak at **~1,500 bp** (canker medians 1,451–1,502 bp).
  - A likely explanation is that healthy bark is dominated by apple **chloroplast** 16S (the course demo found 99.8% chloroplast in healthy bark), while cankers contain real bacteria such as *Erwinia*. Step 4 (taxonomy) will confirm or rule this out.
- **There's a separate spike at ~1,819 bp.** It holds 17.9% of all reads, matching the course demo's "~18% at ~1,820 bp". It is **apple mitochondrial DNA**, not bacteria (see the in-silico PCR check below).
  - These reads fall outside the DADA2 window, so they are dropped in step 3 and do not affect the downstream results.
- **Almost nothing is short (<1,000 bp)**, so there is no sign of degraded DNA or truncated reads.

**In-silico PCR check** (`scripts/02d_insilico_pcr.py`): we searched three reference genomes from NCBI for the 27F and 1492R primer sites and measured the product each would give.

| Genome (NCBI accession) | Predicted 27F–1492R product | Peak in our reads |
|---|---|---|
| Apple mitochondrion (NC_018554.1) | **1,819 bp** (0 primer mismatches) | 1,819 bp spike |
| Apple chloroplast (NC_061549.1) | **1,451 bp** (0 mismatches; 2 copies) | 1,451 bp peak dominating healthy samples |
| *Erwinia amylovora* Ea1189 (NZ_CP055227.2) | **1,500–1,501 bp** (7 rRNA copies, 0 mismatches) | 1,501 bp peak, almost all in cankers |

**Why the mitochondrial product is longer:**
- Plant mitochondria descend from bacteria, so their small-subunit rRNA gene still has perfect binding sites for "universal" bacterial primers.
- The plant gene carries extra variable domains, which make it longer than bacterial 16S (Spencer et al. 1984).
- *E. coli* 16S is 1,541 nt (Brosius et al. 1978), so bacterial 27F–1492R products are ~1,450–1,510 bp. An 1,819 bp product is outside that range.

The chloroplast product (1,451 bp) sits inside the bacterial size range, so chloroplast cannot be removed by length. It has to be removed by taxonomy (step 4).

References:
- Brosius J, Palmer ML, Kennedy PJ, Noller HF (1978). Complete nucleotide sequence of a 16S ribosomal RNA gene from *Escherichia coli*. *PNAS* 75:4801–4805. doi:10.1073/pnas.75.10.4801
- Spencer DF, Schnare MN, Gray MW (1984). Pronounced structural similarities between the small subunit ribosomal RNA genes of wheat mitochondria and *Escherichia coli*. *PNAS* 81:493–497.
- Goremykin VV et al. (2012). The mitochondrial genome of *Malus domestica* and the import-driven hypothesis of mitochondrial genome expansion in seed plants. *Plant J* 71:615–626 (source of NC_018554).

**Ours vs course demo:** median length 1,456 bp (demo: 1,456), and 82.0% of reads within 1,000–1,600 bp (demo: 82%).

---

## 2. `02_read_quality.pdf`: per-read mean quality

**What it shows:** a histogram of each read's mean Phred quality score, from seqkit. The dashed line marks the overall median.

**How to read Phred scores:**
- Q20 means about 1 error per 100 bases.
- Q30 means about 1 error per 1,000 bases.
- Q40 means about 1 error per 10,000 bases.

**What to look for:**
- **Nothing falls below Q20.** That is by design: PacBio only calls a read "HiFi" if it reaches Q20 or better, so the floor is a property of the data type, not a problem.
- **The median is Q37.4**, about 99.98% accuracy per base, which matches the course demo's Q37. That is excellent for amplicon work.
- **The spike at about Q93** is the maximum score PacBio assigns. These are reads where every pass agreed perfectly, which is common for shorter, well-covered molecules.
- **Healthy reads are shifted lower** (median Q34.2 vs Q38.6 for canker) and have fewer Q93 reads. That is still high quality.
  - The probable reason is fewer sequencing passes per molecule in healthy samples, but this is not verified.
  - The course demo does not split quality by group, so this is an extra detail rather than a disagreement.
- **The comb/saw-tooth pattern is a plotting artifact.** Mean quality values cluster at certain numbers, and the 0.5-wide bins catch them unevenly. Ignore it.

**Takeaway:** quality is not a concern for any sample. DADA2's `maxEE = 2` filter will remove the few error-prone reads.

---

## 3. `02_reads_per_sample.pdf`: sequencing depth

**What it shows:** the number of HiFi reads in each sample, grouped by cultivar.

**What to look for:**
- **The range is 35,672 to 117,612 reads**, identical to the course demo's numbers. No sample failed.
- **Healthy samples consistently have fewer reads**: about 36k–52k, versus about 50k–118k for cankers.
- **Why it matters:** read count affects how many ASVs you can detect. This is why alpha diversity (step 5) uses **rarefaction** (subsampling every sample to the same depth). Keep the depth difference in mind when comparing canker and healthy diversity.
- **Within a cultivar, canker depth varies** (e.g. Cripps Pink 3 at ~53k vs Cripps Pink 4 at ~117k). That is normal library-to-library variation.

---

## 4. `02_long_read_class.pdf`: the ~1,820 bp reads, per sample

**What it shows:** for each sample, the percentage of reads in the 1,744–1,894 bp class (the spike in figure 1). Boxes summarize each group, dots are samples, and dot color is cultivar. The P-value is from a Wilcoxon rank-sum test comparing canker and healthy samples.

**What to look for:**
- **Healthy samples are very consistent**, at 16–26% (median 21.7%).
- **Canker samples are all over the place**, from 0.7% to 49% (median 19.1%).
  - **Cripps Pink cankers are low**: 5 of 6 are under 5%.
  - Some Aztec Fuji, Hewe's Crab and Gala cankers are at 25–49%.
- **Canker vs healthy: Wilcoxon P = 0.14, not significant.**
  - Pooled over all reads, healthy bark does have a higher share (22.1% vs 16.9%), which agrees with the course demo's statement that these reads are "more common in healthy bark".
  - Sample by sample, though, the canker spread is too wide for a clear difference.
- **One way to read this:** if the reads are host mitochondria, the share reflects how much host DNA was amplified relative to bacterial DNA.
  - In cankers, where *Erwinia* can dominate, bacterial template can swamp the host signal (Cripps Pink cankers had the highest *E. amylovora* load in the course demo analysis).
  - This is a hypothesis to keep in mind, not a result.

---

## 5. `02_primer_orientation.pdf`: read orientation

**What it shows:** for each sample, the share of reads that begin (within the first 60 bp) with the forward primer **27F** (maroon), the reverse primer **1492R** (orange), or neither (grey).

**What to look for:**
- **Reads come in both orientations**, about 34% 27F-first and 48% 1492R-first. This is normal for PacBio, which sequences each molecule's strand in whichever direction it happens to load.
  - It is why step 3 uses DADA2's **`denoise-ccs`**, which finds the primers on either strand and flips reads into one orientation.
  - This matches the course demo.
- **"No primer found" is 18.3% overall**: about 15–21% in cankers and 23–27% in healthy samples.
  - These reads probably have a primer with sequencing errors, a partial primer, or a primer slightly further in than 60 bp.
  - In step 3, DADA2 discards reads where it cannot find the primers, so expect some loss there.
- **Caution:** the 18.3% "no primer" rate happens to be close to the 17.9% long-read share, but **they are not the same reads.**
  - Cripps Pink cankers have only ~1–2% long reads yet still ~17% no-primer reads.
  - The two numbers are a coincidence, not one population.

---

## Read-count versions (for comparing with the course demo slide)

**Files:** `02_read_length_counts.pdf` and `02_read_quality_counts.pdf` are the same data as figures 1 and 2, but with **number of reads** on the y-axis instead of density.

**What's different:**
- The groups are **stacked**: the healthy bar sits on the bottom, canker on top, and the bar's full height is the total reads in that bin.
- The length plot copies the course demo's layout: 10 bp bins, 0–2,600 bp, and dashed lines at 1,000 and 1,600 bp. Only the colors differ (VT maroon/orange instead of its red/blue).

**Why use them:** counts show how much data sits where. For example, most of the 1,450 bp peak is canker reads simply because cankers have more reads in total. Density plots hide that; they are better for comparing the *shape* of the two groups.

**Comparing bar heights with the course demo figure:**
- Small differences in individual bar heights are expected. The peaks are very sharp (most reads are within a few bp of 1,451, 1,501 and 1,819 bp), so where the bin edges fall decides which bar a read lands in.
  - Our bins are 1,441–1,450, 1,451–1,460 and so on. The demo's bin edges are not stated.
  - Shifting the edges by a few bp can move a couple of hundred thousand reads between neighboring bars.
- The totals are identical (3,529,407 reads), so compare the overall pattern, not the exact bar heights:

| Peak | Ours (tallest bar) | Course demo (approx. from slide) |
|---|---|---|
| ~1,450 bp | 1,186,611 (1,451–1,460 bp), with 610,368 in the 1,441–1,450 bar next to it | ~1,140,000 |
| ~1,500 bp | 549,434 (1,501–1,510 bp), with 249,091 in the 1,491–1,500 bar | ~680,000 |
| ~1,820 bp | 507,769 (1,811–1,820 bp) | ~430,000 |

The pattern matches: a large mixed peak at ~1,450 bp, a canker-only peak at ~1,500 bp, and the ~1,820 bp spike. Bar heights are within the range expected from different bin edges. The tallest-bin table is written to `results/02_qc/length_top_bins.tsv` when the step is run.

**One thing the stacked count plot makes clear:** the **~1,500 bp peak is almost entirely canker** (549,197 canker vs 237 healthy reads in the tallest bar). The ~1,450 bp peak has both groups. This supports the idea that 1,500 bp reads are bacterial (e.g. *Erwinia*) and 1,450 bp reads include chloroplast; step 4 will test that.

---

## How step 2 sets up step 3

| Observation | Consequence in step 3 (DADA2) |
|---|---|
| Reads in both orientations | Use `denoise-ccs` with both primers |
| Median Q37, all reads ≥ Q20 | Standard `maxEE = 2` is appropriate |
| ~18% of reads at ~1,820 bp | Removed by the 1,000–1,600 bp length filter |
| ~18% of reads without a detectable primer | Some loss at the primer-removal step |
| Healthy samples are shallower and likely chloroplast-heavy | Expect few bacterial ASVs in healthy samples (course demo: 1–12 per sample) |

Together these explain why only about half of the raw reads are expected to survive denoising (course demo: 1,839,144 of 3,529,407 non-chimeric reads).
