#!/usr/bin/env python3
"""
trace_reads.py: GFA subgraph extraction by A-line read tracing.

Implements Baris's approach: given a contig-to-reference BAM, identify which
unitig segments in a p_utg.gfa cover the same locus.

Pipeline:
  1. Parse BAM primary alignments → contig name + CIGAR → contig offset range
  2. Scan p_ctg.gfa A-lines for that contig, filtered by offset range → HiFi read names
  3. Scan p_utg.gfa A-lines for those reads → unitig segment names
  4. Write utg names to stdout (one per line) for: gfatools view -l @utg_list -r 0

The BAM should be filtered to MAPQ≥60 primary alignments before running.
Contig offset range from CIGAR is expanded by --flank on both sides (default 150 kb)
to capture reads overlapping the full target locus, not just the gene body.

Usage:
  python trace_reads.py \
    --bam hg002.hap1_to_pat.CYP2D6.MAPQ60.bam \
    --ctg-gfa HG002.asm.bp.hap1.p_ctg.gfa \
    --utg-gfa HG002.asm.bp.p_utg.gfa \
    [--flank 150000] \
    > utg_list.txt

  # Then extract subgraph:
  gfatools view -l @utg_list.txt -r 0 HG002.asm.bp.p_utg.gfa > subgraph.gfa
"""

import sys
import re
import subprocess
import argparse


def parse_cigar_query_range(cigar_str):
    """Return (query_start, query_end) of the aligned (non-soft-clipped) region.

    In hifiasm contig-to-ref alignments the soft clip before the first M equals
    the contig offset of the locus, so this gives us the correct window into the
    contig's A-lines.
    """
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar_str)
    q_pos = 0
    first_m = None
    last_m_end = None
    for length_str, op in ops:
        length = int(length_str)
        if op in ('M', 'X', '=', 'I'):
            if first_m is None:
                first_m = q_pos
            q_pos += length
            last_m_end = q_pos
        elif op in ('S', 'H'):
            q_pos += length
        # D, N do not consume query bases
    return first_m, last_m_end


def parse_bam(bam_path):
    """Run samtools view and return primary alignments as (qname, qs, qe)."""
    proc = subprocess.run(
        ['samtools', 'view', '-F', '2308', bam_path],
        capture_output=True, text=True, check=True
    )
    alignments = []
    for line in proc.stdout.splitlines():
        fields = line.split('\t')
        if len(fields) < 6:
            continue
        qname = fields[0]
        cigar = fields[5]
        if cigar == '*':
            continue
        qs, qe = parse_cigar_query_range(cigar)
        if qs is not None and qe is not None:
            alignments.append((qname, qs, qe))
    return alignments


def collect_reads_from_ctg_gfa(ctg_gfa, targets):
    """Scan p_ctg.gfa A-lines; collect reads overlapping each (ctg, win_start, win_end).

    targets: list of (ctg_name, win_start, win_end)
    Returns: set of read names
    """
    target_map = {ctg: (ws, we) for ctg, ws, we in targets}
    reads = set()
    with open(ctg_gfa) as fh:
        for line in fh:
            if not line.startswith('A\t'):
                continue
            parts = line.rstrip().split('\t')
            # A <seg> <offset> <strand> <read> <rstart> <rend> ...
            if len(parts) < 7:
                continue
            seg = parts[1]
            if seg not in target_map:
                continue
            win_start, win_end = target_map[seg]
            offset = int(parts[2])
            read_span = int(parts[6]) - int(parts[5])  # read_end - read_start
            read_end_on_ctg = offset + read_span
            if offset <= win_end and read_end_on_ctg >= win_start:
                reads.add(parts[4])
    return reads


def collect_utgs_from_reads(utg_gfa, target_reads):
    """Scan p_utg.gfa A-lines; return utg names containing any read in target_reads."""
    utgs = set()
    with open(utg_gfa) as fh:
        for line in fh:
            if not line.startswith('A\t'):
                continue
            parts = line.rstrip().split('\t')
            if len(parts) < 5:
                continue
            if parts[4] in target_reads:
                utgs.add(parts[1])
    return utgs


if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description='Trace reads through hifiasm A-lines to find unitigs covering a locus.'
    )
    p.add_argument('--bam', required=True,
                   help='MAPQ-filtered BAM of contig-to-reference alignment')
    p.add_argument('--ctg-gfa', required=True,
                   help='hifiasm p_ctg GFA matching the BAM haplotype (hap1 or hap2)')
    p.add_argument('--utg-gfa', required=True,
                   help='hifiasm p_utg GFA (shared by both haplotypes)')
    p.add_argument('--flank', type=int, default=150000,
                   help='bp to extend on each side of the CIGAR-derived window (default 150000)')
    args = p.parse_args()

    sys.stderr.write(f"Step 1: parsing BAM {args.bam} ...\n")
    alignments = parse_bam(args.bam)
    if not alignments:
        sys.exit("ERROR: no primary alignments found in BAM")

    targets = []
    for qname, qs, qe in alignments:
        ws = max(0, qs - args.flank)
        we = qe + args.flank
        sys.stderr.write(
            f"  {qname}: gene-body offsets [{qs:,}, {qe:,}]  "
            f"window [{ws:,}, {we:,}]  (flank ±{args.flank:,} bp)\n"
        )
        targets.append((qname, ws, we))

    sys.stderr.write(f"\nStep 2: scanning {args.ctg_gfa} for A-lines in target windows ...\n")
    reads = collect_reads_from_ctg_gfa(args.ctg_gfa, targets)
    sys.stderr.write(f"  {len(reads):,} reads collected\n")

    if not reads:
        sys.exit("ERROR: no reads found in window — check contig name and offset range")

    sys.stderr.write(f"\nStep 3: scanning {args.utg_gfa} for those {len(reads):,} reads ...\n")
    utgs = collect_utgs_from_reads(args.utg_gfa, reads)
    sys.stderr.write(f"  {len(utgs):,} unitigs found\n\n")

    for utg in sorted(utgs):
        print(utg)
