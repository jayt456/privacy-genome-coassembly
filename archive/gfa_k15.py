#!/usr/bin/env python3
"""
gfa.py: identify GFA segments overlapping a gene locus by genome-unique k-mer hits.

Usage:
    python gfa.py <assembly.bp.p_utg.gfa> --kmers <unique_k15.txt> [--threshold pct]

    --kmers FILE      genome-unique k-mer file (from ffind.py --output-kmers)
    --threshold PCT   override match threshold (default: 1.0)

Output: one segment ID per line → feed to gfatools: -l @segments.txt -r 0

Score = hits / len(ref_kmers) * 100  (ref_pct).
A segment passes if it covers >threshold% of the reference k-mer set (~218/21,796 for CYP2D6).
Paralogs are cleanly excluded because they share no genome-unique k-mers.

K-mer encoding (k=15, canonical):
    Each base → 2 bits: A=0 C=1 G=2 T=3
    15-mer → 30-bit unsigned int; canonical = min(forward, reverse-complement)
    Rolling hash: O(n) per sequence
"""

import sys

K = 15
MASK = (1 << (2 * K)) - 1

_ENC = [None] * 256
for _c, _v in zip('ACGTacgt', [0, 1, 2, 3, 0, 1, 2, 3]):
    _ENC[ord(_c)] = _v


def kmerize(seq):
    """Yield canonical 15-mer integers from seq using a rolling hash."""
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
            hits = 0
            for km in kmerize(seq):
                if km in ref_kmers:
                    hits += 1
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
        sys.exit(f"Usage: {sys.argv[0]} <assembly.gfa> --kmers <unique_k15.txt> [--threshold pct]")

    asm_path = positional[0]
    threshold = threshold_override if threshold_override is not None else 1.0

    sys.stderr.write(f"Loading genome-unique k-mers from {kmers_file} ...\n")
    ref_kmers = load_kmers_file(kmers_file)
    n = len(ref_kmers)
    sys.stderr.write(f"  {n:,} genome-unique {K}-mers\n")
    sys.stderr.write(f"  threshold: >{threshold:.1f}% ref_pct (hits / {n:,} unique k-mers)\n")

    score_gfa(asm_path, ref_kmers, threshold)
