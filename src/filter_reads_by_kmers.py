#!/usr/bin/env python3
"""
Retain FASTQ reads containing at least one genome-unique k-mer.

Loads both forward and RC of each reference k-mer so reads on either strand
are caught — consistent with gfa_k21.py's canonical k-mer approach.
Threshold is ≥1 hit (vs gfa_k21.py's >1% ref_pct for segments): appropriate
for short reads where even a boundary-overlapping read has few matching k-mers.

Usage:
    python3 filter_reads_by_kmers.py \
        --fastq reads.fastq.gz \
        --kmers reference/hg38CYP2D6.unique_k21.txt \
        --out   filtered.fastq.gz
"""
import argparse, gzip, sys, time

_RC = str.maketrans('ACGTacgt', 'TGCAtgca')

def rc(seq):
    return seq.translate(_RC)[::-1]

def load_kmers(path, k):
    kmers = set()
    with open(path) as fh:
        for line in fh:
            s = line.strip().upper()
            if len(s) == k:
                kmers.add(s)
                kmers.add(rc(s))
    return kmers

def read_has_kmer(seq, kmers, k):
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        if seq[i:i+k] in kmers:
            return True
    return False

def open_path(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fastq', required=True)
    ap.add_argument('--kmers', required=True)
    ap.add_argument('--out',   required=True)
    ap.add_argument('-k', type=int, default=21)
    args = ap.parse_args()

    kmers = load_kmers(args.kmers, args.k)
    print(f"Loaded {len(kmers):,} k-mers (fwd + RC)  k={args.k}", file=sys.stderr)

    total = kept = 0
    t0 = time.time()

    with open_path(args.fastq) as fin, open_path(args.out, 'wt') as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq  = fin.readline().rstrip('\n')
            plus = fin.readline()
            qual = fin.readline()
            total += 1
            if read_has_kmer(seq, kmers, args.k):
                kept += 1
                fout.write(header)
                fout.write(seq + '\n')
                fout.write(plus)
                fout.write(qual)
            if total % 50000 == 0:
                print(f"  {total:,} reads  {kept:,} kept  ({100*kept/total:.1f}%)  "
                      f"[{time.time()-t0:.0f}s]", file=sys.stderr)

    print(f"Done: {total:,} total  {kept:,} kept  ({100*kept/total:.1f}%)  "
          f"[{time.time()-t0:.0f}s]", file=sys.stderr)

if __name__ == '__main__':
    main()
