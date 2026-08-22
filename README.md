# Privacy-Preserving Genomic Co-Assembly

Some genomic regions are medically critical but too complex to assemble well from a single institution's data. This project builds a system that lets multiple institutions identify shared sequence at these regions and cooperate to assemble them better — without any raw genomic data ever crossing institutional boundaries.

**PI:** Hoon Cho (Yale BIDS / Yale CS) &nbsp;·&nbsp; **Student:** Jay Tummala (Yale CS)

---

## What is this?

Many important genes — drug-metabolism enzymes, immune receptor clusters, disease-linked structural variants — sit in regions of the genome so repetitive and structurally variable that assembling them accurately requires more samples than any single institution has. The natural solution is collaboration across hospitals and research centers, but raw genomic reads are personally identifiable and legally protected.

This project's answer: compress each institution's reads into compact mathematical fingerprints, add carefully calibrated noise for differential privacy, and use a cryptographic protocol (MHE-PSI) to find which samples share sequence at a target gene — without either institution ever seeing the other's data. Once matched, they co-assemble the region from pooled reads, recovering assembly quality that neither could achieve alone.

The longer-term goal is a **privacy-preserving sequence index** queryable during assembly itself: rather than pooling reads, an assembler could query what other institutions observed at ambiguous positions to resolve gaps and avoid discarding real sequence as apparent errors.

This repository contains the **Python tools and primary extraction pipeline** for the assembly side: gene-specific subgraph extraction from whole-genome hifiasm assemblies using genome-unique k-mers, validated across 58 samples and 8 complex gene loci. The privacy layer (sketching → DP → MHE-PSI) builds on top of this.

---

## Project Status

| Stage | Status |
|-------|--------|
| Reference preparation (8-gene panel, k=21 k-mer files) | ✅ Complete |
| 12× hifiasm assembly — 58 samples (11 pedigree + 47 HPRC) | ✅ Complete |
| Gene-region subgraph extraction — k=21 method | ✅ Complete |
| HG002 validation against Q100 ground truth | ✅ Complete |
| Subgraph extraction — full cohort (58 × 8 genes) | 🔄 In progress |
| Minimizer sketching of subgraphs | 📋 Planned |
| Differential privacy layer (randomized response) | 📋 Planned |
| MHE-PSI inter-institutional comparison | 📋 Planned |
| Co-assembly evaluation (NGA50 / QV vs. epsilon) | 📋 Planned |

---

## Key Results

**The k=15 → k=21 pivot.** Initial 15-mer scoring produced 7–11 segments per assembly at CYP2D6 — too many. The cause: CYP2D6 has two near-identical pseudogenes (CYP2D7/CYP2D8) sharing >90% sequence identity, so 15-mers from CYP2D6 recurred in pseudogenes and were incorrectly labeled unique. At k=21, a single nucleotide difference anywhere in the 21-mer window is enough to distinguish CYP2D6 from its paralogs completely.

**CYP2D6 method comparison (HG002):**

| Method | Assembly | Segments | Size | Result |
|--------|----------|----------|------|--------|
| k=15 scoring | high-cov | 11 | 1.41 Mb | ✗ paralog contamination |
| k=15 scoring | 12× | 7 | 0.91 Mb | ✗ paralog contamination |
| **k=21 scoring** | **high-cov** | **2** | **1.14 Mb** | **✓ Q100 verified — 1 segment per haplotype** |
| **k=21 scoring** | **12×** | **4** | **585 kb** | **✓ matches read-tracing exactly** |
| Read-tracing (A-lines) | high-cov | 4 | 7.55 Mb | ✓ correct reads, chromosome-scale UTGs |
| Boundary v2 (CIGAR) | high-cov | 2 | 1.14 Mb | ✓ matches k=21 exactly |

Three independent methods agree on the same 2 segments for the high-coverage assembly, confirmed against the HG002 Q100 ground-truth reference. The 12× assembly fragments to 4 segments — illustrating the problem co-assembly is meant to solve.

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

## Gene Panel

8 medically relevant genes chosen for structural complexity. Coordinates in `genes.txt`.

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

Pipeline scripts use three environment variables (Yale HPC defaults are built in):

```bash
export COASM_BASE=/path/to/project       # assemblies, results, reference files
export COASM_SCRATCH=/path/to/scratch    # large transient outputs
export COASM_TOOLS=/path/to/repo/src     # Python tools directory
```

Copy `config.env.example` → `config.env` and fill in your paths.

---

## Quick Start

### 1. Clone and configure
```bash
git clone https://github.com/jayt456/privacy-genome-coassembly
cd privacy-genome-coassembly
cp config.env.example config.env
# edit config.env with your paths
```

### 2. Build genome-unique k-mer files (one-time, ~2–3 hrs)

Build a Jellyfish k=21 database from hg38, then run `ffind.py` for each gene:

```bash
# Build k=21 Jellyfish DB (~22 GB)
jellyfish count -m 21 -s 10G -t 16 -o hg38.k21.jf hg38.fasta

# Get genome-unique 21-mers for a gene
python src/ffind.py hg38.fasta chr22 42126499 42130810 \
    -k 21 -t 1 -m 10000 \
    --db hg38.k21.jf \
    --output-kmers hg38CYP2D6.unique_k21.txt
```

### 3. Assemble at 12× and extract subgraphs

Run hifiasm on your 12× FASTQ, then score and extract:

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

Walks outward from a gene body until each flank contains ≥m k-mers appearing at most t times genome-wide in a Jellyfish database. Emits the unique k-mers as a flat text file for use with `gfa_k21.py`.

> **Key fix:** Jellyfish built without `-C` (canonical flag) stores only one strand. Querying only the forward k-mer gives false positives — k-mers whose reverse complement is common genome-wide look unique. `ffind.py` queries both strands and requires `fwd_count + RC_count ≤ t`.

```
usage: ffind.py ref.fasta CHROM START END [-k K] [-t T] [-m M] --db JF_DB [--output-kmers FILE]
```

---

### `gfa_k21.py` — GFA subgraph scorer

Scores every segment in a hifiasm `p_utg.gfa` against a genome-unique k-mer set. A segment passes if `hits / |ref_kmers| × 100 > threshold` (default 1%). Prints passing segment IDs to stdout for piping to `gfatools view`.

```
usage: gfa_k21.py GFA --kmers KMER_FILE [--threshold FLOAT] [--k INT]
```

---

### `filter_reads_by_kmers.py` — Read pre-filter for targeted assembly

Filters a FASTQ to reads containing ≥1 genome-unique k-mer. Enables targeted hifiasm on just the gene-matching reads. At 12×, ~0.016% of whole-genome HiFi reads match CYP2D6 (358 of 2.2M for HG002); the targeted assembly produces fewer, cleaner fragments.

```
usage: filter_reads_by_kmers.py --fastq FASTQ --kmers KMER_FILE --out OUTPUT
```

---

### `trace_reads.py` — A-line read tracer *(methodology reference)*

> **Note:** Included for methodology documentation, not general use. Requires a Q100 ground-truth reference (e.g. GIAB HG002 CYP2D6 sequences) to produce the contig-to-reference BAM used as input. Used to independently validate k=21 results for HG002. For general subgraph extraction, use `gfa_k21.py`.

Traces A-lines through hifiasm's internal graph bookkeeping: contig-to-reference BAM → p_ctg A-lines → HiFi read names → p_utg A-lines → unitig names. All three independent methods (k=21 scoring, A-line tracing, CIGAR boundary parsing) agreed on the same segments for HG002.

```
usage: trace_reads.py --bam BAM --ctg-gfa P_CTG_GFA --utg-gfa P_UTG_GFA
```

---

## Cohort Extraction Script

`pipeline/subgraph/run_gfa_subset.sh` is the SLURM wrapper for running `gfa_k21.py` across all samples and all genes in `genes.txt`. It reads gene names and coordinates dynamically, skips completed outputs, and handles both pedigree and HPRC datasets.

```bash
sbatch pipeline/subgraph/run_gfa_subset.sh hprc      # 47 HPRC samples × 8 genes
sbatch pipeline/subgraph/run_gfa_subset.sh pedigree  # 11 pedigree samples × 8 genes
sbatch pipeline/subgraph/run_gfa_subset.sh both       # all 58 samples
```

Output per sample-gene: `<sample>.<GENE>.k21.subset.gfa` + `<sample>.<GENE>.k21.segments`.

---

## Archive

[`archive/gfa_k15.py`](archive/gfa_k15.py) is the k=15 predecessor to `gfa_k21.py`, kept to document why k=15 fails at pharmacogene loci. The scoring logic is identical; only the k-mer length differs. Running it on CYP2D6 produces 7–11 segments due to CYP2D7/CYP2D8 paralog contamination — the result that motivated the switch to k=21.

---

## Adding a New Gene

1. Add a row to `genes.txt` with gene body coordinates
2. Build a genome-unique k-mer file with `ffind.py` (see Quick Start step 2)
3. Run `gfa_k21.py` or `run_gfa_subset.sh` with the new gene

---

## Citation

> Jay Tummala, Hoon Cho. Privacy-Preserving Genomic Co-Assembly at Complex Loci Using Genome-Unique K-mer Subgraph Extraction. *In preparation*, 2026.

## Contact

Jay Tummala — jay.tummala@yale.edu  
PI: Hoon Cho — hoon.cho@yale.edu  
Yale BIDS / Yale Computer Science
