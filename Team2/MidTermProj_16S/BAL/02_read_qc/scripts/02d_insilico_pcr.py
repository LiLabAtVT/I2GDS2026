# Step 2d: in-silico PCR with 27F/1492R on reference genomes, to identify the read-length peaks.
# Genomes (NCBI): NC_018554.1 Malus x domestica mitochondrion; NC_061549.1 Malus domestica chloroplast;
# NZ_CP055227.2 Erwinia amylovora Ea1189 chromosome. Result: 1,819 / 1,451 / 1,500-1,501 bp.
"""In-silico PCR: 27F / 1492R product sizes in a genome FASTA.

"In-silico PCR" means doing the PCR on the computer: we search a known
genome for places where both primers would stick, and measure the length of
DNA between them. That is the length of read we would expect to see from
that organism. Comparing these predicted lengths with the peaks in our read
length plot tells us which organism each peak comes from.

A primer site must match its 3'-most 8 bases exactly (IUPAC-aware) and have
<= MAXMM mismatches over the full primer. Product length = 27F start .. 1492R end
(primers included), which is what a HiFi read of the amplicon spans.
Usage: python 02d_insilico_pcr.py genome.fa [genome2.fa ...]
Each FASTA file must hold a single sequence (one genome or chromosome).
"""
import re
import sys

F = "AGRGTTYGATYMTGGCTCAG"     # 27F
R = "RGYTACCTTGTTACGACTT"      # 1492R
# Allow up to 3 mismatched bases per primer; real PCR tolerates a few
MAXMM = 3
# Which real bases each degenerate letter can match (e.g. R = A or G)
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "M": "AC", "K": "GT",
         "S": "CG", "W": "AT", "N": "ACGT", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}
# Pairing rules for the opposite DNA strand: A<->T, C<->G (and degenerate codes)
COMP = str.maketrans("ACGTRYMKSWNBDHV", "TGCAYRKMSWNVHDB")


def rc(s):
    # Reverse complement: the sequence of the opposite DNA strand, read 5'->3'
    return s.translate(COMP)[::-1]


def sites(seq, primer):
    """Start positions where primer (5'->3') binds the given strand, 3' anchored."""
    # The 3' end of a primer (its last 8 bases) must match perfectly, because
    # that is where the polymerase starts copying. Find those spots first.
    anchor = "".join(f"[{IUPAC[b]}]" for b in primer[-8:])
    out = []
    for m in re.finditer(f"(?={anchor})", seq):
        start = m.start() - (len(primer) - 8)
        if start < 0:
            continue
        # Then count mismatches over the whole primer at that spot
        win = seq[start:start + len(primer)]
        mm = sum(c not in IUPAC[p] for c, p in zip(win, primer))
        if mm <= MAXMM:
            out.append((start, mm))
    return out


# Run once for each genome file given on the command line
for path in sys.argv[1:]:
    # FASTA: first line is the name (">..."), the rest is the sequence
    lines = open(path).read().split("\n")
    name, seq = lines[0][1:], "".join(lines[1:]).upper()
    print(f"\n{name}  ({len(seq):,} bp)")
    # Check both strands of the genome, since the gene can be on either
    for strand, s in (("+", seq), ("-", rc(seq))):
        fwd = sites(s, F)
        rev = sites(s, rc(R))           # 1492R binds the top strand as its reverse complement
        # Pair every forward site with every reverse site downstream of it
        for fs, fmm in fwd:
            for rs, rmm in rev:
                size = rs + len(R) - fs
                # Keep only products of a plausible amplicon size
                if 1000 <= size <= 3000:
                    print(f"  strand {strand}: product {size:,} bp  (27F mismatches {fmm}, 1492R mismatches {rmm})")
