# Privacy-Preserving Genomic Co-Assembly

**PI:** Hoon Cho (Yale BIDS / Yale CS) &nbsp;·&nbsp; **Student:** Jay Tummala (Yale CS)

Many medically important genomic regions — pharmacogenes, immune receptor clusters, copy-number variable loci — are too structurally complex to assemble well from a single institution's data. This project investigates whether multiple institutions can identify shared sequence at these regions from low-coverage long reads and cooperatively assemble them better, without exposing raw genomic data.

Current work focuses on **gene-region subgraph extraction**: given a whole-genome hifiasm assembly, accurately isolating just the segments corresponding to a target gene using genome-unique k-mers. This is the foundation for co-assembly across samples.

---

## Project Status

| Stage | Status |
|-------|--------|
| Reference preparation (8-gene panel, k=21 k-mer files) | ✅ Complete |
| 12× hifiasm assembly — 58 samples (11 pedigree + 47 HPRC) | ✅ Complete |
| Gene-region subgraph extraction — k=21 method | ✅ Complete |
| HG002 validation against Q100 ground truth | ✅ Complete |
| Subgraph extraction — full cohort (58 × 8 genes) | 🔄 In progress |
| Co-assembly methods and evaluation | 📋 Planned |
| Privacy layer (differential privacy via randomized response) | 📋 Planned |

---

## Key Results

The core challenge is extracting gene-region segments from a whole-genome assembly graph containing tens of thousands of segments. We score each segment by how many genome-unique k-mers from the target gene it contains, keeping only high-scoring ones.

**The k=15 → k=21 pivot.** Initial 15-mer scoring produced 7–11 segments per assembly at CYP2D6 — too many. The cause: CYP2D6 has two near-identical pseudogenes (CYP2D7/CYP2D8) sharing >90% sequence identity, so 15-mers from CYP2D6 recurred in pseudogenes and were incorrectly labeled unique. At k=21, a single nucleotide difference anywhere in the window is enough to distinguish CYP2D6 from its paralogs.

**CYP2D6 method comparison (HG002):**

| Method | Assembly | Segments | Size | Result |
|--------|----------|----------|------|--------|
| k=15 scoring | high-cov | 11 | 1.41 Mb | ✗ paralog contamination |
| k=15 scoring | 12× | 7 | 0.91 Mb | ✗ paralog contamination |
| **k=21 scoring** | **high-cov** | **2** | **1.14 Mb** | **✓ Q100 verified — 1 segment per haplotype** |
| **k=21 scoring** | **12×** | **4** | **585 kb** | **✓ matches read-tracing exactly** |
| Read-tracing (A-lines) | high-cov | 4 | 7.55 Mb | ✓ correct reads, chromosome-scale UTGs |
| Boundary v2 (CIGAR) | high-cov | 2 | 1.14 Mb | ✓ matches k=21 exactly |

Three independent methods agree on the same segments for HG002, confirmed against the Q100 ground-truth reference. The 12× assembly fragments to 4 segments — illustrating exactly the problem co-assembly is meant to address.

**All 8 genes (HG002, 12× assembly):**

| Gene | Segments | Size | Note |
|------|----------|------|------|
| CYP2D6 | 4 | 585 kb | Q100-verified |
| TERT | 2 | 305 kb | |
| TCF3 | 3 | 190 kb | |
| SLC6A3 | 6 | 281 kb | overlaps TERT window on chr5 |
| LMF1 | 4 | 465 kb | |
| IGH | 23 | 1.6 Mb | fragmented — 1.3 Mb gene body, repeat-dense |
| KIR | 10 | 442 kb | fragmented — copy-number variable locus |
| HLA | 48 | 7.3 Mb | expected — 7.9 Mb region, extreme diversity |

---

## Future Directions

- **Co-assembly:** pool extracted subgraph reads from matched samples across institutions and evaluate whether co-assembled output improves on any single institution's 12× assembly (NGA50, QV, haplotype completeness)
- **Privacy layer:** differential privacy via randomized response to allow institutions to compare sequence without exposing raw reads

---

## Gene Panel

8 medically relevant genes, chosen for structural complexity. Coordinates in `genes.txt`.

| Gene | Chrom | Flank window | Role |
|------|-------|-------------|------|
| CYP2D6 | chr22 | 41,983,509–42,269,959 | Drug metabolism — ~25% of clinical drugs |
| TERT | chr5 | 1,186,033–1,371,812 | Telomerase reverse transcriptase |
| TCF3 | chr19 | 1,525,600–1,749,816 | Transcription factor E2A; B-cell development |
| SLC6A3 | chr5 | 1,318,105–1,511,314 | Dopamine transporter |
| LMF1 | chr16 | 794,014–1,040,403 | Lipase maturation factor |
| IGH | chr14 | 105,484,613–106,883,718 | Immunoglobulin heavy chain (1.3 Mb) |
| KIR | chr19 | 54,631,973–55,015,162 | Killer-cell immunoglobulin-like receptors |
| HLA | chr6 | 25,576,036–33,521,605 | HLA complex — most polymorphic region in genome |

Flank windows computed by `src/ffind.py` (k=15, t=1, m=10,000): each side contains ≥10,000 genome-unique anchor k-mers.

---

## Datasets

**Platinum Pedigree (CEU family 1463)** — 11 samples  
- HiFi reads: `s3://platinum-pedigree-data/data/hifi/mapped/GRCh38/<sample>.GRCh38.haplotagged.bam`
- High-coverage assemblies: `s3://platinum-pedigree-data/assemblies/`

**HPRC Year 1** — 47 samples  
- HiFi reads: `s3://human-pangenomics/working/HPRC[_PLUS]/<sample>/analysis/aligned_reads/hifi/GRCh38/<sample>_aligned_GRCh38_winnowmap.sorted.bam`
- Assemblies: `s3://human-pangenomics/working/HPRC[_PLUS]/<sample>/assemblies/year1_freeze_assembly_v2/`

---

## Prerequisites

| Tool | Version tested | Purpose |
|------|---------------|---------|
| [hifiasm](https://github.com/chhylp123/hifiasm) | 0.25.0 | Phased HiFi assembly |
| [minimap2](https://github.com/lh3/minimap2) | 2.29 | Assembly-to-reference alignment |
| [samtools](http://www.htslib.org/) | 1.21 | BAM processing |
| [jellyfish](https://github.com/gmarcais/Jellyfish) | 2.2.10 | k-mer counting |
| [gfatools](https://github.com/lh3/gfatools) | latest | GFA subgraph extraction |
| Python | 3.8+ | — |
| [BioPython](https://biopython.org/) | 1.87 | Required by `ffind.py` only |

```bash
conda create -n biotools -c bioconda -c conda-forge jellyfish biopython python=3.10
conda activate biotools
```

Set paths in `config.env` (copy from `config.env.example`):
```bash
export COASM_BASE=/path/to/project
export COASM_SCRATCH=/path/to/scratch
export COASM_TOOLS=/path/to/repo/src
```

---

## Quick Start

### 1. Build genome-unique k-mer files (one-time)
```bash
# k=21 Jellyfish database from hg38 (~22 GB)
jellyfish count -m 21 -s 10G -t 16 -o hg38.k21.jf hg38.fasta

# Genome-unique 21-mers for a gene
python src/ffind.py hg38.fasta chr22 42126499 42130810 \
    -k 21 -t 1 -m 10000 \
    --db hg38.k21.jf \
    --output-kmers hg38CYP2D6.unique_k21.txt
```

### 2. Assemble and extract subgraph
```bash
hifiasm -o sample.bp -t 16 sample.12x.fastq.gz

python src/gfa_k21.py sample.bp.p_utg.gfa \
    --kmers hg38CYP2D6.unique_k21.txt \
    > passing_segments.txt

gfatools view -l @passing_segments.txt -r 0 sample.bp.p_utg.gfa \
    > sample.CYP2D6.k21.subset.gfa
```

For cohort-scale extraction across many samples and genes, see `pipeline/subgraph/run_gfa_subset.sh`.

---

## Python Tools (`src/`)

### `ffind.py` — Flanking anchor finder

Walks outward from a gene body until each flank contains ≥m k-mers appearing at most t times genome-wide. Emits the unique k-mer list used by `gfa_k21.py`.

> **Note on Jellyfish strand handling:** a database built without `-C` stores only one strand. `ffind.py` queries both the forward k-mer and its reverse complement, requiring `fwd_count + RC_count ≤ t` — otherwise k-mers whose complement is common genome-wide appear falsely unique.

### `gfa_k21.py` — GFA subgraph scorer

Scores every segment in a hifiasm `p_utg.gfa` against the genome-unique k-mer set. A segment passes if `hits / |ref_kmers| × 100 > threshold` (default 1%). Prints passing segment IDs to stdout for piping to `gfatools view`.

### `filter_reads_by_kmers.py` — Read pre-filter

Filters a FASTQ to reads containing ≥1 genome-unique k-mer, enabling targeted hifiasm on just gene-matching reads. For HG002 CYP2D6: 358 of 2.2M reads passed (0.016%), yielding a 3-segment / 128 kb targeted assembly vs. 4 segments / 585 kb whole-genome baseline.

### `trace_reads.py` — A-line read tracer *(methodology reference)*

> **Note:** Requires a Q100 ground-truth reference (used here for HG002 validation only). For general subgraph extraction use `gfa_k21.py`.

Traces hifiasm's internal A-lines — contig-to-reference BAM → p_ctg A-lines → read names → p_utg A-lines → unitig names — as an independent check on k=21 results.

---

## Archive

[`archive/gfa_k15.py`](archive/gfa_k15.py) is the k=15 predecessor to `gfa_k21.py`. Kept to document why k=15 fails: running it on CYP2D6 produces 7–11 segments due to CYP2D7/CYP2D8 paralog contamination, motivating the switch to k=21.

---

## Citation

> Jay Tummala, Hoon Cho. Privacy-Preserving Genomic Co-Assembly at Complex Loci. *In preparation*, 2026.

## Contact

Jay Tummala — jay.tummala@yale.edu &nbsp;·&nbsp; Hoon Cho — hoon.cho@yale.edu  
Yale BIDS / Yale Computer Science
