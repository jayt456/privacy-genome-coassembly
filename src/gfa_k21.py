#!/usr/bin/env python3
"""
gfa_k21.py: genome-unique k-mer subgraph scoring at k=21.

Identical to gfa.py except K=21. Use with hg38<GENE>.unique_k21.txt files
produced by ffind.py -k 21.

Usage:
    python gfa_k21.py <assembly.bp.p_utg.gfa> \\
        --kmers reference/hg38CYP2D6.unique_k21.txt \\
        [--threshold 1.0]

Output: one segment ID per line → gfatools view -l @segments.txt -r 0
"""

import sys

K = 21  # only change from gfa.py
MASK = (1 << (2 * K)) - 1

_ENC = [None] * 256
for _c, _v in zip('ACGTacgt', [0, 1, 2, 3, 0, 1, 2, 3]):
    _ENC[ord(_c)] = _v


def kmerize(seq):
    """Yield canonical 21-mer integers from seq using a rolling hash."""
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


def load_kmers_file(path):
    """Load k-mer strings (one per line) as canonical ints."""
    kmers = set()
    with open(path) as fh:
        for line in fh:
            seq = line.strip().upper()
            if not seq or seq.startswith('#'):
                continue
            for km in kmerize(seq):
                kmers.add(km)
    return kmers


def score_gfa(gfa_path, ref_kmers, threshold):
    """Print segment IDs where hits/len(ref_kmers)*100 > threshold."""
    ref_n = len(ref_kmers)
    kept = total = 0
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
            hits = sum(1 for km in kmerize(seq) if km in ref_kmers)
            if ref_n and (hits / ref_n * 100) > threshold:
                print(seg_id)
                kept += 1
    sys.stderr.write(f"  {kept}/{total} segments retained (>{threshold:.1f}% ref_pct)\n")


if __name__ == '__main__':
    _args = sys.argv[1:]
    kmers_file = None
    threshold_override = None
    positional = []
    _i = 0
    while _i < len(_args):
        if _args[_i] == '--kmers' and _i + 1 < len(_args):
            kmers_file = _args[_i + 1]; _i += 2
        elif _args[_i] == '--threshold' and _i + 1 < len(_args):
            threshold_override = float(_args[_i + 1]); _i += 2
        else:
            positional.append(_args[_i]); _i += 1

    if len(positional) < 1 or not kmers_file:
        sys.exit(f"Usage: {sys.argv[0]} <assembly.gfa> --kmers <unique_k21.txt> [--threshold pct]")

    asm_path = positional[0]
    threshold = threshold_override if threshold_override is not None else 1.0

    sys.stderr.write(f"Loading genome-unique k-mers from {kmers_file} ...\n")
    ref_kmers = load_kmers_file(kmers_file)
    n = len(ref_kmers)
    sys.stderr.write(f"  {n:,} genome-unique {K}-mers\n")
    sys.stderr.write(f"  threshold: >{threshold:.1f}% ref_pct\n")

    score_gfa(asm_path, ref_kmers, threshold)
