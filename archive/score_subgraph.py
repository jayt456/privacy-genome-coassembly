#!/usr/bin/env python3
"""
score_subgraph.py: Score every S-line in an existing subgraph GFA against
a gene reference, reporting both contig-coverage and reference-coverage metrics.

Usage:
    python score_subgraph.py <gene.fa> <subgraph.gfa> [--sort contig|ref]
    python score_subgraph.py --kmers <unique.txt> <subgraph.gfa> [--sort contig|ref]

  --kmers FILE   use a pre-computed k-mer file (one k-mer string per line, e.g.
                 from ffind.py --output-kmers) instead of a gene FASTA

Output TSV (stdout):
    segment_id  length  hits  contig_kmers  contig_pct  ref_pct

  contig_pct = hits / contig_kmers * 100  — appropriate for short unitig GFAs
  ref_pct    = hits / ref_n * 100         — appropriate for chromosome-scale contigs

A histogram and summary are printed to stderr.
"""

import sys
import math

K = 15
MASK = (1 << (2 * K)) - 1

_ENC = [None] * 256
for _c, _v in zip('ACGTacgt', [0, 1, 2, 3, 0, 1, 2, 3]):
    _ENC[ord(_c)] = _v


def kmerize(seq):
    fwd = rev = run = 0
    for ch in seq:
        b = _ENC[ord(ch)] if ord(ch) < 256 else None
        if b is None:
            fwd = rev = run = 0
            continue
        fwd = ((fwd << 2) | b) & MASK
        rev = ((rev >> 2) | ((3 - b) << (2 * (K - 1)))) & MASK
        run += 1
        if run >= K:
            yield fwd if fwd <= rev else rev


def load_reference_kmers(fasta_path):
    kmers = set()
    buf = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if buf:
                    for km in kmerize(''.join(buf)):
                        kmers.add(km)
                    buf = []
            else:
                buf.append(line)
    if buf:
        for km in kmerize(''.join(buf)):
            kmers.add(km)
    return kmers


def load_kmers_file(path):
    """Load pre-computed k-mer strings (one per line) as canonical 30-bit ints."""
    kmers = set()
    with open(path) as fh:
        for line in fh:
            seq = line.strip().upper()
            if not seq or seq.startswith('#'):
                continue
            for km in kmerize(seq):
                kmers.add(km)
    return kmers


def histogram(values, bins=20, width=50):
    if not values:
        return
    lo, hi = 0.0, 100.0
    step = (hi - lo) / bins
    counts = [0] * bins
    for v in values:
        idx = min(int((v - lo) / step), bins - 1)
        counts[idx] += 1
    maxc = max(counts) if counts else 1
    sys.stderr.write(f"\n  {'pct range':<14} {'count':>6}  bar\n")
    sys.stderr.write(f"  {'-'*14}  {'-'*6}  {'-'*width}\n")
    for i, c in enumerate(counts):
        lo_b = lo + i * step
        hi_b = lo_b + step
        bar = '#' * int(c / maxc * width)
        sys.stderr.write(f"  {lo_b:5.1f}–{hi_b:5.1f}%   {c:>6}  {bar}\n")


def main():
    args = sys.argv[1:]
    kmers_file = None
    sort_by = "contig"
    positional = []
    i = 0
    while i < len(args):
        if args[i] == '--kmers' and i + 1 < len(args):
            kmers_file = args[i + 1]; i += 2
        elif args[i] == '--sort' and i + 1 < len(args):
            sort_by = args[i + 1]; i += 2
        else:
            positional.append(args[i]); i += 1

    if kmers_file:
        if len(positional) < 1:
            sys.exit(f"Usage: {sys.argv[0]} --kmers <file> <subgraph.gfa> [--sort contig|ref]")
        gfa_path = positional[0]
    else:
        if len(positional) < 2:
            sys.exit(f"Usage: {sys.argv[0]} <gene.fa> <subgraph.gfa> [--sort contig|ref]")
        gene_fa  = positional[0]
        gfa_path = positional[1]

    if kmers_file:
        sys.stderr.write(f"Loading genome-unique k-mers from {kmers_file} ...\n")
        ref_kmers = load_kmers_file(kmers_file)
        ref_n = len(ref_kmers)
        sys.stderr.write(f"  {ref_n:,} genome-unique {K}-mers\n\n")
    else:
        sys.stderr.write(f"Loading reference k-mers from {gene_fa} ...\n")
        ref_kmers = load_reference_kmers(gene_fa)
        ref_n = len(ref_kmers)
        sys.stderr.write(f"  {ref_n:,} unique {K}-mers\n\n")

    sys.stderr.write(f"Scoring S-lines in {gfa_path} ...\n")

    rows = []
    total = 0
    with open(gfa_path) as fh:
        for line in fh:
            if not line.startswith('S\t'):
                continue
            parts = line.split('\t', 3)
            if len(parts) < 3:
                continue
            seg_id = parts[1]
            seq = parts[2].rstrip()
            if seq == '*' or len(seq) < K:
                continue
            total += 1
            if total % 50 == 0:
                sys.stderr.write(f"  processed {total} segments...\r")

            hits = n_contig = 0
            for km in kmerize(seq):
                n_contig += 1
                if km in ref_kmers:
                    hits += 1

            contig_pct = (hits / n_contig * 100) if n_contig else 0.0
            ref_pct    = (hits / ref_n    * 100) if ref_n    else 0.0
            rows.append((seg_id, len(seq), hits, n_contig, contig_pct, ref_pct))

    sys.stderr.write(f"\n  scored {total} segments\n")

    # Sort
    if sort_by == "ref":
        rows.sort(key=lambda r: r[5], reverse=True)
    else:
        rows.sort(key=lambda r: r[4], reverse=True)

    # Print TSV
    print("segment_id\tlength\thits\tcontig_kmers\tcontig_pct\tref_pct")
    for seg_id, length, hits, n_contig, contig_pct, ref_pct in rows:
        print(f"{seg_id}\t{length}\t{hits}\t{n_contig}\t{contig_pct:.2f}\t{ref_pct:.4f}")

    # Histograms
    contig_pcts = [r[4] for r in rows]
    ref_pcts    = [r[5] for r in rows]

    sys.stderr.write("\n=== contig_pct distribution (hits / contig_kmers) ===")
    histogram(contig_pcts)

    sys.stderr.write("\n=== ref_pct distribution (hits / ref_kmers) ===")
    histogram(ref_pcts)

    # Summary stats
    def stats(vals):
        vals = sorted(vals)
        n = len(vals)
        mean = sum(vals) / n if n else 0
        median = vals[n // 2] if n else 0
        return mean, median, vals[0] if vals else 0, vals[-1] if vals else 0

    sys.stderr.write("\n=== Summary ===\n")
    sys.stderr.write(f"  {'metric':<14} {'mean':>8} {'median':>8} {'min':>8} {'max':>8}\n")
    for label, vals in [("contig_pct", contig_pcts), ("ref_pct", ref_pcts)]:
        mn, md, lo, hi = stats(vals)
        sys.stderr.write(f"  {label:<14} {mn:>8.2f} {md:>8.2f} {lo:>8.2f} {hi:>8.2f}\n")

    # Threshold table
    sys.stderr.write("\n=== Segments surviving various contig_pct thresholds ===\n")
    for t in [0, 10, 25, 50, 75, 90, 95, 99]:
        n = sum(1 for v in contig_pcts if v > t)
        sys.stderr.write(f"  >{t:3d}%:  {n:4d} / {total}\n")

    sys.stderr.write("\n=== Segments surviving various ref_pct thresholds ===\n")
    for t in [0, 25, 50, 75, 90, 95, 99, 100]:
        n = sum(1 for v in ref_pcts if v > t)
        sys.stderr.write(f"  >{t:3d}%:  {n:4d} / {total}\n")


if __name__ == '__main__':
    main()
