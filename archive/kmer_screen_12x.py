#!/usr/bin/env python3
"""
kmer_screen_12x.py

Build a canonical k-mer reference set from 4 non-redundant CYP2D6 regions in
the high-coverage boundary subgraph, then screen every segment in a 12x
p_utg.gfa for matches.

The 4 regions are defined by A-line read offsets (marking where CYP2D6 reads
land on each UTG) and the L-line overlap (splitting the two hap2 UTGs into
non-redundant pieces):

  R1  hap1 UTG    CYP2D6 interval  [A_min_offset, A_max_end]
  R2  from_seg    non-overlap CYP2D6 part  (before/after the junction, depending
                  on L-line orientation)
  R3  overlap     the shared junction sequence between the two hap2 UTGs
                  (taken from to_seg's coordinate space)
  R4  to_seg      non-overlap CYP2D6 part  [overlap_end, A_max_end]

Together R1–R4 cover the full CYP2D6 region for both haplotypes without
including large UTG sequence outside CYP2D6, and without double-counting the
junction overlap.

Usage example:
  python3 kmer_screen_12x.py \\
    --highcov-gfa  HG002_highcov.CYP2D6.boundary.subset.gfa \\
    --alines       HG002_highcov.CYP2D6.hap1.boundary.lines \\
                   HG002_highcov.CYP2D6.hap2.boundary.lines \\
    --gfa-12x      HG002_12x.bp.p_utg.gfa \\
    -k 61 --min-hits 2 \\
    --out-names    HG002_12x.CYP2D6.kmer61.names \\
    --out-fasta    HG002_12x.CYP2D6.kmer61.fa
"""
import argparse
import sys

# ---------------------------------------------------------------------------
# k-mer utilities (bytes throughout for speed)
# ---------------------------------------------------------------------------
_COMP = bytes.maketrans(b'ACGTacgt', b'TGCAtgca')

def _rc(b):
    return b.translate(_COMP)[::-1]

def kmerize(seq_bytes, k):
    """Yield canonical k-mers (bytes) from seq_bytes, skipping any with N/n."""
    n = len(seq_bytes) - k + 1
    for i in range(n):
        km = seq_bytes[i:i+k]
        if b'N' not in km and b'n' not in km:
            r = _rc(km)
            yield km if km <= r else r

# ---------------------------------------------------------------------------
# GFA / A-line parsing
# ---------------------------------------------------------------------------
def parse_gfa_seqs(gfa_path):
    """Return {seg_name: seq_bytes} for all S-lines with non-* sequences."""
    segs = {}
    with open(gfa_path, 'rb') as f:
        for line in f:
            if line[:2] != b'S\t':
                continue
            parts = line.rstrip(b'\n').split(b'\t')
            if parts[2] != b'*':
                segs[parts[1].decode()] = parts[2]
    return segs

def parse_gfa_link(gfa_path):
    """Return (from_seg, from_ori, to_seg, to_ori, overlap_len) from first L-line."""
    with open(gfa_path) as f:
        for line in f:
            if not line.startswith('L\t'):
                continue
            p = line.rstrip('\n').split('\t')
            cigar = p[5]
            if not cigar.endswith('M'):
                sys.exit(f"Unexpected L-line CIGAR (expected <n>M): {cigar}")
            return p[1], p[2], p[3], p[4], int(cigar[:-1])
    return None

def parse_aline_intervals(paths):
    """
    Parse A-lines files; return {seg_name: (min_offset, max_end)}.
    A-line format: A  seg  offset  strand  read_name  read_start  read_end  ...
    read span on segment = read_end - read_start
    """
    ivs = {}
    for path in paths:
        with open(path) as f:
            for line in f:
                p = line.rstrip('\n').split('\t')
                if len(p) < 7 or p[0] != 'A':
                    continue
                seg = p[1]
                offset = int(p[2])
                end = offset + int(p[6]) - int(p[5])
                if seg not in ivs:
                    ivs[seg] = (offset, end)
                else:
                    s, e = ivs[seg]
                    ivs[seg] = (min(s, offset), max(e, end))
    return ivs

# ---------------------------------------------------------------------------
# 4-region k-mer set construction
# ---------------------------------------------------------------------------
def build_kmer_set(segs, link, ivs, k):
    """
    Given 3 UTG sequences, the L-line link, and A-line CYP2D6 intervals,
    extract the 4 non-redundant regions and return their canonical k-mer set.

    L-line orientation rules:
      from_ori '-' : overlap sits at the START of from_seg (forward coords)
      from_ori '+' : overlap sits at the END   of from_seg
      to_ori   '+' : overlap sits at the START of to_seg
      to_ori   '-' : overlap sits at the END   of to_seg
    Using canonical k-mers, orientation of a region doesn't change the set.
    """
    from_seg, from_ori, to_seg, to_ori, ovlp = link

    if from_ori == '-':
        from_ovlp_start, from_ovlp_end = 0, ovlp
    else:
        fs_len = len(segs[from_seg])
        from_ovlp_start, from_ovlp_end = fs_len - ovlp, fs_len

    if to_ori == '+':
        to_ovlp_start, to_ovlp_end = 0, ovlp
    else:
        ts_len = len(segs[to_seg])
        to_ovlp_start, to_ovlp_end = ts_len - ovlp, ts_len

    # R1: hap1 UTG (the one not participating in the link)
    r1_seg = next(s for s in segs if s not in (from_seg, to_seg))
    r1_s, r1_e = ivs[r1_seg]
    r1 = segs[r1_seg][r1_s:r1_e]

    # R2: non-overlap CYP2D6 portion of from_seg
    _, from_a_max = ivs.get(from_seg, (0, len(segs[from_seg])))
    r2_s = from_ovlp_end     # start after the overlap region
    r2_e = from_a_max
    r2 = segs[from_seg][r2_s:r2_e] if r2_s < r2_e else b''

    # R3: overlap region (taken from to_seg perspective)
    r3 = segs[to_seg][to_ovlp_start:to_ovlp_end]

    # R4: non-overlap CYP2D6 portion of to_seg
    _, to_a_max = ivs.get(to_seg, (0, len(segs[to_seg])))
    r4_s = to_ovlp_end       # start after the overlap region
    r4_e = to_a_max
    r4 = segs[to_seg][r4_s:r4_e] if r4_s < r4_e else b''

    regions = [
        (f'R1 {r1_seg}[{r1_s:,}:{r1_e:,}]', r1),
        (f'R2 {from_seg}[{r2_s:,}:{r2_e:,}]', r2),
        (f'R3 {to_seg}[{to_ovlp_start:,}:{to_ovlp_end:,}] (overlap)', r3),
        (f'R4 {to_seg}[{r4_s:,}:{r4_e:,}]', r4),
    ]

    kset = set()
    total_bp = 0
    for label, seq in regions:
        before = len(kset)
        for km in kmerize(seq, k):
            kset.add(km)
        total_bp += len(seq)
        print(f'  {label}: {len(seq):,} bp → +{len(kset)-before:,} new k-mers', flush=True)

    print(f'  Total region bp: {total_bp:,}', flush=True)
    return kset

# ---------------------------------------------------------------------------
# 12x GFA screening
# ---------------------------------------------------------------------------
def screen_gfa(gfa_path, kset, k, min_hits):
    kept = []
    total = 0
    with open(gfa_path, 'rb') as f:
        for line in f:
            if line[:2] != b'S\t':
                continue
            total += 1
            parts = line.rstrip(b'\n').split(b'\t')
            name = parts[1].decode()
            seq = parts[2]
            if seq == b'*':
                continue
            hits = sum(1 for km in kmerize(seq, k) if km in kset)
            if hits >= min_hits:
                kept.append((name, seq.decode(), hits))
            if total % 10000 == 0:
                print(f'  ... {total:,} segments screened ({len(kept)} kept)', flush=True)
    return kept, total

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--highcov-gfa', required=True,
                    help='High-cov boundary subset GFA (3 segs, 1 L-line)')
    ap.add_argument('--alines', required=True, nargs='+',
                    help='A-lines files (hap1.boundary.lines hap2.boundary.lines)')
    ap.add_argument('--gfa-12x', required=True,
                    help='12x p_utg.gfa to screen')
    ap.add_argument('-k', type=int, default=61,
                    help='k-mer size (default: 61)')
    ap.add_argument('--min-hits', type=int, default=2,
                    help='Min k-mer hits to keep a segment (default: 2). '
                         'All kept segment hit counts are printed for threshold inspection.')
    ap.add_argument('--out-names', required=True,
                    help='Output: one kept segment name per line')
    ap.add_argument('--out-fasta', required=True,
                    help='Output: FASTA of kept segments (header includes hit count)')
    args = ap.parse_args()
    k = args.k

    print(f'=== CYP2D6 k-mer screen  k={k}  min_hits={args.min_hits} ===\n', flush=True)

    # Step 1: load high-cov subgraph
    print(f'[1] Loading high-cov subgraph: {args.highcov_gfa}', flush=True)
    segs = parse_gfa_seqs(args.highcov_gfa)
    for name, seq in segs.items():
        print(f'  {name}: {len(seq):,} bp', flush=True)

    link = parse_gfa_link(args.highcov_gfa)
    if link is None:
        sys.exit('No L-line found in subgraph GFA — cannot determine overlap.')
    from_seg, from_ori, to_seg, to_ori, ovlp = link
    print(f'  Link: {from_seg}({from_ori}) → {to_seg}({to_ori}),  overlap = {ovlp:,} bp\n',
          flush=True)

    # Step 2: parse A-line CYP2D6 intervals
    print(f'[2] Parsing A-line intervals from: {", ".join(args.alines)}', flush=True)
    ivs = parse_aline_intervals(args.alines)
    for seg, (s, e) in ivs.items():
        print(f'  {seg}: [{s:,}, {e:,})  =  {e-s:,} bp', flush=True)

    # Step 3: build k-mer reference set
    print(f'\n[3] Building {k}-mer reference set from 4 regions...', flush=True)
    kset = build_kmer_set(segs, link, ivs, k)
    print(f'  → {len(kset):,} unique canonical {k}-mers\n', flush=True)

    # Step 4: screen 12x p_utg.gfa
    print(f'[4] Screening {args.gfa_12x} ...', flush=True)
    kept, total = screen_gfa(args.gfa_12x, kset, k, args.min_hits)
    print(f'\n  Segments screened : {total:,}', flush=True)
    print(f'  Segments kept (>= {args.min_hits} hits): {len(kept)}', flush=True)
    print(f'  Hit count distribution (sorted by hits, descending):', flush=True)
    for name, _, hits in sorted(kept, key=lambda x: -x[2]):
        print(f'    {name}  {hits:,} hits', flush=True)

    # Step 5: write outputs
    with open(args.out_names, 'w') as f:
        for name, _, _ in kept:
            f.write(name + '\n')

    with open(args.out_fasta, 'w') as f:
        for name, seq, hits in kept:
            f.write(f'>{name} kmer_hits={hits}\n{seq}\n')

    print(f'\nNames written : {args.out_names}', flush=True)
    print(f'FASTA written : {args.out_fasta}', flush=True)

if __name__ == '__main__':
    main()
