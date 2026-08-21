#!/usr/bin/env python3
"""
ffind.py
===============

Given a reference indexed FASTA and a region (chrom, s, e),
find the region (sf, ef) such that:

    * the LEFT flank  [sf, s) contains >= m k-mers whose count <= t
    * the RIGHT flank [e, ef) contains >= m k-mers whose count <= t

Example (ensure jellyfish is in your PATH):
python ffind.py hg38.giab.fa chr22 42077656 42253758 --build -k 15 -t 1 -m 10000 --threads 96
[build] jellyfish count -m 15 -s 3G -t 96 -o ../HG008/hg38.giab.fa.k15.jf ../HG008/hg38.giab.fa
region: chr22:42077656-42253758  (len 176102)
params             : k=15 t=1 m=10000

left  flank [41945352, 42077656)  len    132304  anchors 10000/10000  OK
right flank [42253758, 42343312)  len     89554  anchors 10000/10000  OK
"""

from __future__ import annotations
import argparse
import os
import subprocess
from subprocess import Popen,PIPE
import tempfile
import sys
from typing import Dict, Iterable, Iterator, List, Optional, Tuple
from Bio.Seq import Seq

class FaidxFasta:
    def __init__(self, ref_path: str, fai_path: Optional[str] = None):
        self.ref_path = ref_path
        self.fai_path = fai_path or ref_path + ".fai"
        if not os.path.exists(self.fai_path): raise FileNotFoundError(f"FASTA index not found: {self.fai_path}")
        self.index: Dict[str, Tuple[int, int, int, int]] = {}
        with open(self.fai_path) as fh:
            for line in fh:
                name, length, offset, linebases, linewidth = line.split("\t")[:5]
                self.index[name] = (int(length), int(offset), int(linebases), int(linewidth))
        self._fh = open(ref_path, "rb")

    def length(self, chrom: str) -> int:
        return self.index[chrom][0]
    def _byte_at(self, chrom: str, pos: int) -> int:
        _, offset, linebases, linewidth = self.index[chrom]
        return offset + (pos // linebases) * linewidth + (pos % linebases)
    def fetch(self, chrom: str, start: int, end: int) -> Seq:
        if chrom not in self.index: raise KeyError(f"{chrom!r} not in {self.fai_path}")
        length = self.index[chrom][0]
        start = max(0, start)
        end = min(length, end)
        if start >= end: return Seq("")
        first = self._byte_at(chrom, start)
        last = self._byte_at(chrom, end - 1)
        self._fh.seek(first)
        raw = self._fh.read(last - first + 1)
        seq = raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii")
        return Seq(seq.upper())
    def close(self):
        self._fh.close()

class Jellyfish:
    def __init__(self, db_path: str, k: int):
        self.db_path = db_path
        self.k = k
    @classmethod
    def build(cls, ref: str, k: int, out: str, threads: int = 4, hash_size: str = "3G", verbose: bool = True) -> "Jellyfish":
        cmd = ["jellyfish", "count", "-m", str(k), "-s", hash_size, "-t", str(threads)]
        cmd += ["-o", out, ref]
        if verbose: print("[build] " + " ".join(cmd), file=sys.stderr)
        subprocess.run(cmd, check=True)
        return cls(out, k)
    def counts(self, kmers: List[str]) -> List[int]:
        if not kmers: return []
        result = []
        proc = subprocess.Popen(['jellyfish', 'query', '-i', self.db_path], stdin=PIPE, stdout=PIPE, text=True)
        for kmer in kmers:
            proc.stdin.write(kmer + "\n")
            proc.stdin.flush()
            c = proc.stdout.readline()
            result.append(int(c.strip()))
        proc.stdin.close()
        return result

def region_candidates(fa: FaidxFasta, chrom: str, s: int, e: int, k: int, chunk: int = 65536) -> Iterator[Tuple[int, str]]:
    """Yield (pos, kmer) for every non-N k-mer starting in [s, e-k]."""
    pos = s
    end_pos = e - k  # last valid k-mer start
    while pos <= end_pos:
        # Fetch chunk + k-1 extra bases so k-mers spanning a chunk boundary are complete
        block_end = min(e, pos + chunk + k - 1)
        seq = str(fa.fetch(chrom, pos, block_end))
        n_kmers = min(chunk, end_pos - pos + 1)
        for i in range(n_kmers):
            kmer = seq[i:i + k]
            if 'N' not in kmer:
                yield pos + i, kmer
        pos += chunk


def _rc(kmer: str) -> str:
    return kmer.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]


def collect_all_unique(jf: Jellyfish, candidates: Iterable[Tuple[int, str]], t: int, batch_size: int = 4096) -> List[Tuple[int, str, int]]:
    """Return all (pos, kmer, fwd_count) from candidates where fwd+rc count <= t.

    The Jellyfish DB is built without -C (non-canonical), so forward and
    reverse-complement k-mers are counted separately.  Checking only the
    forward count can include k-mers whose RC is common genome-wide; those
    would match assembly segments on the opposite strand after canonicalisation
    in gfa.py.  Querying both strands and summing gives true genome-wide count.
    """
    result = []
    buf: List[Tuple[int, str]] = []

    def _flush(buf):
        fwd_kmers = [km for _, km in buf]
        rc_kmers  = [_rc(km) for _, km in buf]
        fwd_counts = jf.counts(fwd_kmers)
        rc_counts  = jf.counts(rc_kmers)
        for (p, km), fc, rc in zip(buf, fwd_counts, rc_counts):
            if fc + rc <= t:
                result.append((p, km, fc + rc))

    for item in candidates:
        buf.append(item)
        if len(buf) >= batch_size:
            _flush(buf)
            buf = []
    if buf:
        _flush(buf)
    return result


def left_candidates(fa: FaidxFasta, chrom: str, s: int, k: int, min_pos: int, chunk: int) -> Iterator[Tuple[int, str]]:
    cur = s - k
    while cur >= min_pos:
        block_lo = max(min_pos, cur - chunk + 1)
        block = str(fa.fetch(chrom, block_lo, cur + k))
        for p in range(cur, block_lo - 1, -1):
            kmer = block[p - block_lo : p - block_lo + k]
            if "N" not in kmer: yield p, kmer
        cur = block_lo - 1

def right_candidates(fa: FaidxFasta, chrom: str, e: int, k: int, max_end: int, chunk: int) -> Iterator[Tuple[int, str]]:
    hi_limit = max_end - k
    cur = e
    while cur <= hi_limit:
        block_hi = min(hi_limit, cur + chunk - 1)
        block = str(fa.fetch(chrom, cur, block_hi + k))
        for p in range(cur, block_hi + 1):
            kmer = block[p - cur : p - cur + k]
            if "N" not in kmer: yield p, kmer
        cur = block_hi + 1

def _iter_counts(jf: Jellyfish, candidates: Iterable[Tuple[int, str]], batch_size: int) -> Iterator[Tuple[int, str, int]]:
    buf: List[Tuple[int, str]] = []
    for item in candidates:
        buf.append(item)
        if len(buf) >= batch_size:
            for (p, km), c in zip(buf, jf.counts([km for _, km in buf])): yield p, km, c
            buf = []
    if buf:
        for (p, km), c in zip(buf, jf.counts([km for _, km in buf])): yield p, km, c

def collect_anchors(jf: Jellyfish, candidates: Iterable[Tuple[int, str]], t: int, m: int, batch_size: int) -> Tuple[List[Tuple[int, str, int]], bool]:
    anchors: List[Tuple[int, str, int]] = []
    for p, km, c in _iter_counts(jf, candidates, batch_size):
        if c <= t:
            anchors.append((p, km, c))
            if len(anchors) >= m: return anchors, True
    return anchors, False

class FlankResult:
    def __init__(self, chrom, s, e, k, t, m):
        self.chrom, self.s, self.e, self.k, self.t, self.m = chrom, s, e, k, t, m
        self.sf, self.ef = s, e
        self.left_anchors: List[Tuple[int, str, int]] = []
        self.right_anchors: List[Tuple[int, str, int]] = []
        self.left_ok = self.right_ok = False
    @property
    def length(self) -> int: return self.ef - self.sf

    def report(self) -> str:
        return "\n".join(
            [
                f"region: {self.chrom}:{self.s}-{self.e}  (len {self.e - self.s})",
                f"params             : k={self.k} t={self.t} m={self.m}",
                "",
                f"left  flank [{self.sf}, {self.s})  len {self.s - self.sf:>9}  "
                f"anchors {len(self.left_anchors)}/{self.m}  "
                f"{'OK' if self.left_ok else 'INCOMPLETE (hit boundary)'}",
                f"right flank [{self.e}, {self.ef})  len {self.ef - self.e:>9}  "
                f"anchors {len(self.right_anchors)}/{self.m}  "
                f"{'OK' if self.right_ok else 'INCOMPLETE (hit boundary)'}",
                "",
                f"=> minimum region : {self.chrom}:{self.sf}-{self.ef}  (len {self.length})",
            ]
        )

def find(fa: FaidxFasta, jf: Jellyfish, chrom: str, s: int, e: int, k: int, t: int, m: int, max_extend: Optional[int] = None, chunk: int = 65536, batch_size: int = 4096) -> FlankResult:
    chrom_len = fa.length(chrom)
    if not (0 <= s < e <= chrom_len): raise ValueError(f"bad region: need 0 <= s < e <= {chrom_len} (got s={s}, e={e})")
    r = FlankResult(chrom, s, e, k, t, m)
    left_min = 0 if max_extend is None else max(0, s - max_extend)
    right_max = chrom_len if max_extend is None else min(chrom_len, e + max_extend)
    r.left_anchors, r.left_ok = collect_anchors(jf, left_candidates(fa, chrom, s, k, left_min, chunk), t, m, batch_size)
    r.sf = r.left_anchors[-1][0] if r.left_anchors else s
    r.right_anchors, r.right_ok = collect_anchors(jf, right_candidates(fa, chrom, e, k, right_max, chunk), t, m, batch_size)
    r.ef = (r.right_anchors[-1][0] + k) if r.right_anchors else e
    return r

def main(argv=None):
    ap = argparse.ArgumentParser(description="min-length region flanking (s,e) with >= m rare (count <= t) k-mers on each side")
    ap.add_argument("ref", help="reference FASTA (<ref>.fai must exist)")
    ap.add_argument("chrom")
    ap.add_argument("start", type=int, help="region start s (0-based)")
    ap.add_argument("end", type=int, help="region end e (0-based)")
    ap.add_argument("-k", type=int, required=True, help="k-mer length")
    ap.add_argument("-t", type=int, required=True, help="max genome-wide count")
    ap.add_argument("-m", type=int, required=True, help="# k-mers required")
    ap.add_argument("--db", help="jellyfish .jf database (existing, or target with --build)")
    ap.add_argument("--build", action="store_true", help="build the jellyfish DB if absent")
    ap.add_argument("--threads", type=int, default=4, help="jellyfish count threads")
    ap.add_argument("--hash-size", default="3G", help="jellyfish count hash size (-s)")
    ap.add_argument("--max-extend", type=int, default=None, help="cap distance (bp) per side (default: to chromosome end)")
    ap.add_argument("--chunk", type=int, default=65536, help="faidx fetch block size (bp)")
    ap.add_argument("--batch-size", type=int, default=4096, help="k-mers per jellyfish query")
    ap.add_argument("--fai", help="explicit path to the .fai (default: <ref>.fai)")
    ap.add_argument("--output-kmers", metavar="FILE",
                    help="after finding the flanking window, write all genome-unique k-mers "
                         "from the full window [sf, ef] to FILE (one k-mer string per line). "
                         "Used as input to gfa.py --kmers for genome-unique filtering.")
    args = ap.parse_args(argv)
    fa = FaidxFasta(args.ref, fai_path=args.fai)
    db = args.db or f"{args.ref}.k{args.k}.jf"
    if os.path.exists(db): jf = Jellyfish(db, args.k)
    elif args.build:
        jf = Jellyfish.build(args.ref, args.k, db, threads=args.threads, hash_size=args.hash_size)
    else:
        raise SystemExit(f"No jellyfish DB at {db!r}, re-run with --build to create it")

    res = find(fa, jf, args.chrom, args.start, args.end, args.k, args.t, args.m, max_extend=args.max_extend, chunk=args.chunk, batch_size=args.batch_size)
    print(res.report())

    if args.output_kmers:
        sys.stderr.write(f"\nCollecting unique k-mers (count <= {args.t}) "
                         f"from window [{res.sf}, {res.ef}) ...\n")
        unique = collect_all_unique(
            jf,
            region_candidates(fa, args.chrom, res.sf, res.ef, args.k, args.chunk),
            t=args.t,
            batch_size=args.batch_size,
        )
        with open(args.output_kmers, 'w') as out_fh:
            for _, km, _ in unique:
                out_fh.write(km + '\n')
        sys.stderr.write(f"  wrote {len(unique):,} unique k-mers to {args.output_kmers}\n")

    fa.close()
    if not (res.left_ok and res.right_ok): sys.exit(2)

if __name__ == "__main__":
    main()
