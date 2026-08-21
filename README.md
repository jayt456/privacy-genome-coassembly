# Privacy-Preserving Genomic Co-Assembly

Identifying shared genomic sequence across cohorts at complex loci using genome-unique k-mer sketching, enabling privacy-preserving collaborative assembly from low-coverage long reads.

**PI:** Hoon Cho (Yale BIDS / Yale CS)  
**Student:** Jay Tummala (Yale CS)

---

## Overview

Many medically important genomic regions — CYP2D6 (drug metabolism), KIR (NK cell receptors), HLA (immune response), IGH (immunoglobulin loci) — are too structurally complex for short-read assembly. Individual cohorts rarely have sufficient coverage or diversity to assemble these loci well from long reads alone. Co-assembly across institutions could help, but sharing raw reads raises privacy concerns.

**Research question:** Can multiple institutions identify which of their samples share genomic sequence at a target locus using only low-coverage reads — without exposing raw data — so they can pool reads for collaborative co-assembly that outperforms any single institution's assembly?

**Method (in development):**
1. Each institution extracts gene-specific assembly subgraphs from their individual hifiasm GFA using genome-unique k=21 k-mers
2. Subgraphs are compared across institutions using minimizer sketching + differential privacy (randomized response) + multi-party homomorphic encryption PSI
3. Matched samples pool reads for co-assembly

This repository contains the **subgraph extraction pipeline** and **validation experiments** on two public long-read datasets.

---

## Gene Panel

8 medically relevant genes chosen for assembly complexity. Coordinates in `genes.txt`.

| Gene | Chrom | Gene body | Flank window | Role |
|------|-------|-----------|--------------|------|
| CYP2D6 | chr22 | 42,126,499–42,130,810 | 41,983,509–42,269,959 | Drug metabolism (CYP450) |
| TERT | chr5 | 1,253,147–1,295,068 | 1,186,033–1,371,812 | Telomerase reverse transcriptase |
| TCF3 | chr19 | 1,609,290–1,652,615 | 1,525,600–1,749,816 | Transcription factor E2A |
| SLC6A3 | chr5 | 1,392,794–1,445,440 | 1,318,105–1,511,314 | Dopamine transporter |
| LMF1 | chr16 | 853,634–981,318 | 794,014–1,040,403 | Lipase maturation factor |
| IGH | chr14 | 105,586,437–106,879,844 | 105,484,613–106,883,718 | Immunoglobulin heavy chain |
| KIR | chr19 | 54,816,468–54,830,778 | 54,631,973–55,015,162 | Killer-cell immunoglobulin-like receptors |
| HLA | chr6 | 25,726,063–33,400,644 | 25,576,036–33,521,605 | HLA complex (7.9 Mb) |

Flank windows were computed by `src/ffind.py` (k=15, t=1, m=10,000): each side contains ≥10,000 genome-unique anchor k-mers.

---

## Datasets

Both datasets are publicly available.

**Platinum Pedigree (CEU family 1463)** — 11 samples  
NA12877, NA12878, NA12879, NA12881, NA12882, NA12885, NA12886, NA12889, NA12890, NA12891, NA12892
- High-coverage hifiasm assemblies: `s3://platinum-pedigree-data/assemblies/`
- HiFi reads (GRCh38 haplotagged BAMs): `s3://platinum-pedigree-data/data/hifi/mapped/GRCh38/<sample>.GRCh38.haplotagged.bam`

**HPRC Year 1** — 47 samples (listed in `pipeline/assembly/hprc_subsample_hifiasm.sh`)
- Assemblies: `s3://human-pangenomics/working/HPRC[_PLUS]/<sample>/assemblies/year1_freeze_assembly_v2/`
- HiFi reads (GRCh38 aligned BAMs): `s3://human-pangenomics/working/HPRC[_PLUS]/<sample>/analysis/aligned_reads/hifi/GRCh38/<sample>_aligned_GRCh38_winnowmap.sorted.bam`

---

## Prerequisites

### Tools
| Tool | Version tested | Purpose |
|------|---------------|---------|
| [hifiasm](https://github.com/chhylp123/hifiasm) | 0.25.0 | Phased HiFi assembly |
| [minimap2](https://github.com/lh3/minimap2) | 2.29 | Assembly-to-reference alignment |
| [samtools](http://www.htslib.org/) | 1.21 | BAM processing |
| [jellyfish](https://github.com/gmarcais/Jellyfish) | 2.2.10 | k-mer counting |
| [gfatools](https://github.com/lh3/gfatools) | latest | GFA subgraph extraction |
| [awscli](https://aws.amazon.com/cli/) | 2.17 | S3 streaming |
| Python | 3.8+ | — |
| [BioPython](https://biopython.org/) | 1.87 | Required by `ffind.py` only |

### Python environment
```bash
conda create -n biotools -c bioconda -c conda-forge jellyfish biopython python=3.10
conda activate biotools
```

### Environment variables
Pipeline scripts use these three variables (with Yale HPC defaults):
```bash
export COASM_BASE=/path/to/project       # assemblies, results, reference files
export COASM_SCRATCH=/path/to/scratch    # HPRC 12x assemblies (large, transient)
export COASM_TOOLS=/path/to/repo/src     # Python tools directory
```

Copy `config.env.example` → `config.env` and fill in your paths. Source it at the top of each SLURM job, or add to your `~/.bashrc`.

> **Yale Bouchet HPC:** default paths in each script point to `/nfs/roberts/project/pi_hc878/jst72/` — no changes needed if running on the same cluster.

---

## Pipeline Overview

```
hg38 reference
    │
    ├─ extract_region.sh       → per-gene reference FASTA (hg38CYP2D6.fa, etc.)
    │
    ├─ build_jellyfish_db.sh   → hg38.k21.jf  (22 GB; one-time)
    │
    └─ build_k21_all_genes.sh  → hg38<GENE>.unique_k21.txt × 8 genes
                                  (genome-unique 21-mers per gene window)

HiFi BAM (S3)
    │
    └─ hprc_subsample_hifiasm.sh → 12× FASTQ → hifiasm → p_utg.gfa + p_ctg.gfa + FASTAs

p_utg.gfa × 58 samples × 8 genes
    │
    └─ run_gfa_subset.sh        → <sample>_12x.<GENE>.k21.subset.gfa
                                   (gene-specific subgraph, ready for sketch/PSI)

High-cov assemblies (ground truth)
    │
    ├─ submit_align_{pedigree,hprc}.sh → assembly-to-gene BAMs
    │
    └─ submit_merge.sh          → cohort_bams/<dataset>_<GENE>_complete.bam
                                   (IGV visualization)
```

---

## Quick Start

### 1. Set up environment
```bash
git clone https://github.com/jayt456/privacy-genome-coassembly
cd privacy-genome-coassembly
cp config.env.example config.env
# edit config.env with your BASE and SCRATCH paths
```

### 2. Prepare reference files (one-time)
```bash
# Extract gene slices from hg38
bash pipeline/reference/extract_region.sh CYP2D6 chr22:41983509-42269959

# Build k=21 Jellyfish database (~22 GB, ~2 hrs)
sbatch pipeline/reference/build_jellyfish_db.sh  # produces hg38.k21.jf

# Compute genome-unique 21-mers for all 8 genes
sbatch pipeline/reference/build_k21_all_genes.sh  # produces hg38<GENE>.unique_k21.txt
```

### 3. Run 12× assemblies
```bash
# HPRC: stream from S3 → subsample → hifiasm (47 samples, batched 5 at a time)
bash pipeline/assembly/submit_hprc_hifiasm.sh

# Pedigree assemblies use pre-subsampled FASTQs:
sbatch pipeline/assembly/run_hifiasm_12x.sh              # NA12878
sbatch pipeline/assembly/run_hifiasm_12x_array.sh        # other 10 samples
```

### 4. Extract gene subgraphs for all samples
```bash
sbatch pipeline/subgraph/run_gfa_subset.sh hprc      # 47 HPRC samples × 8 genes
sbatch pipeline/subgraph/run_gfa_subset.sh pedigree  # 11 pedigree samples × 8 genes
```

Output: `<sample>_12x.<GENE>.k21.subset.gfa` alongside each sample's `p_utg.gfa`.

### 5. Build cohort alignment BAMs (ground truth / IGV)
```bash
bash pipeline/alignment/submit_align_pedigree.sh
bash pipeline/alignment/submit_align_hprc.sh
bash pipeline/cohort/submit_merge.sh
```

---

## Python Tools

### `src/ffind.py` — Flanking anchor finder

Walks outward from a gene body until each flank contains ≥m k-mers that appear at most t times genome-wide (Jellyfish DB). Also emits the unique k-mers from the full flanked window for use in subgraph scoring.

**Key insight:** Querying both strand orientations in the non-canonical Jellyfish DB (`fwd_count + RC_count ≤ t`) avoids false positives from the opposite-strand complement being common genome-wide.

```bash
python src/ffind.py reference/hg38.fasta chr22 42126499 42130810 \
    -k 21 -t 1 -m 10000 \
    --db reference/hg38.k21.jf \
    --output-kmers reference/hg38CYP2D6.unique_k21.txt
```

---

### `src/gfa_k21.py` — GFA subgraph scorer

Scores every unitig (S-line) in a hifiasm p_utg GFA against the genome-unique k=21 set. A segment passes if `hits / |ref_kmers| × 100 > threshold` (default 1%). Passing segment IDs pipe to `gfatools view` for subgraph extraction.

**Why k=21 and not k=15?** CYP2D6 paralogs CYP2D7/2D8 share >90% identity. At k=15 these paralogs produce thousands of false-positive k-mer matches, contaminating the subgraph with non-CYP2D6 unitigs. k=21 eliminates this completely — verified by Q100 PAF alignment showing 4312/4312 reads at MAPQ60 to the correct haplotype.

```bash
python src/gfa_k21.py assembly.bp.p_utg.gfa \
    --kmers reference/hg38CYP2D6.unique_k21.txt \
    > passing_segments.txt

gfatools view -l @passing_segments.txt -r 0 assembly.bp.p_utg.gfa \
    > CYP2D6_subgraph.gfa
```

---

### `src/filter_reads_by_kmers.py` — Targeted assembly read filter

Pre-filters a FASTQ to reads containing ≥1 genome-unique k-mer. Enables a targeted hifiasm run on just the gene-matching reads, producing a cleaner and faster assembly. At 12×, ~0.016% of whole-genome HiFi reads match CYP2D6 (358 of 2.2M reads for HG002).

```bash
python src/filter_reads_by_kmers.py \
    --fastq reads.12x.fastq.gz \
    --kmers reference/hg38CYP2D6.unique_k21.txt \
    --out filtered_CYP2D6.fastq.gz
```

---

### `src/trace_reads.py` — A-line read tracer (validation)

Alternative subgraph extraction using hifiasm's A-lines rather than k-mer scoring. Given a MAPQ-filtered contig-to-reference BAM: CIGAR offsets → p_ctg A-lines → HiFi read names → p_utg A-lines → unitig names. Used to validate k=21 results for HG002 (requires a Q100 reference, so not generalizable).

```bash
python src/trace_reads.py \
    --bam hg002.hap1_to_Q100_pat.MAPQ60.bam \
    --ctg-gfa HG002.hap1.p_ctg.gfa \
    --utg-gfa HG002.p_utg.gfa \
    > utg_names.txt

gfatools view -l @utg_names.txt -r 0 HG002.p_utg.gfa > CYP2D6_subgraph.gfa
```

---

## Pipeline Scripts

### Reference (`pipeline/reference/`)

| Script | What it does |
|--------|-------------|
| `extract_region.sh` | Cut a named gene slice from hg38 with samtools faidx |
| `build_jellyfish_db.sh` | Build genome-wide k=21 Jellyfish count DB (~22 GB; one-time) |
| `build_k21_all_genes.sh` | SLURM array (0-6): run ffind.py for all 8 genes → `hg38<GENE>.unique_k21.txt` |
| `run_ffind_new_genes.sh` | Run ffind.py for newly added genes and update `genes.txt` in-place |

### Assembly (`pipeline/assembly/`)

| Script | What it does |
|--------|-------------|
| `run_hifiasm_12x.sh` | hifiasm whole-genome phased assembly for NA12878 at 12× |
| `run_hifiasm_12x_array.sh` | SLURM array: remaining 10 pedigree samples |
| `hprc_subsample_hifiasm.sh` | SLURM worker: stream HPRC BAM from S3 → subsample 12× → hifiasm → delete FASTQ |
| `submit_hprc_hifiasm.sh` | Submit hprc_subsample_hifiasm.sh array with 5-job concurrency cap |

Coverage estimation uses 6 autosomal 1 Mb windows; subsampling uses `samtools view -s <seed.fraction>`.

### Alignment (`pipeline/alignment/`)

| Script | What it does |
|--------|-------------|
| `align_pedigree_gene.sh` | minimap2 asm5: pedigree hap1/hap2 assemblies → gene reference slice |
| `align_hprc_gene.sh` | minimap2 asm5: HPRC paternal/maternal assemblies → gene reference slice |
| `submit_align_pedigree.sh` | Submit alignment arrays for all genes in `genes.txt` |
| `submit_align_hprc.sh` | Same for HPRC |
| `verify_pedigree_alignments.sh` | Report BAM sizes + mapped-read counts for all pedigree alignments |
| `verify_hprc_alignments.sh` | Same for HPRC |

### Subgraph extraction (`pipeline/subgraph/`)

| Script | What it does |
|--------|-------------|
| `run_gfa_subset.sh` | Score all samples × all genes with gfa_k21.py; extract subset GFAs with gfatools |
| `run_HG002_all_genes_k21.sh` | SLURM array: HG002 × all 8 genes × {12×, highcov} |

`run_gfa_subset.sh` reads gene names from `genes.txt` dynamically and skips existing outputs.

### Cohort merge (`pipeline/cohort/`)

| Script | What it does |
|--------|-------------|
| `merge_cohort.sh` | MAPQ≥21 filter → per-sample haplotype merge → single cohort BAM |
| `submit_merge.sh` | Submit merge jobs for all gene × {pedigree, hprc} combinations |

Output: `cohort_bams/<dataset>_<GENE>_complete.bam` — load directly into IGV.

### Targeted assembly (`pipeline/targeted/`)

Experiment: pre-filter whole-genome 12× reads by gene k-mers, reassemble the filtered reads with hifiasm, compare to the whole-genome subgraph baseline.

| Script | What it does |
|--------|-------------|
| `run_HG002_12x_store_fastq.sh` | Stream HG002 HiFi from S3 → permanent 12× FASTQ (kept for reuse) |
| `run_HG002_CYP2D6_targeted_assembly.sh` | Filter reads → targeted hifiasm → compare vs. whole-genome baseline |
| `run_NA12878_CYP2D6_targeted.sh` | Same for NA12878 pedigree |

**Result:** For HG002 CYP2D6, targeted assembly (358 reads) yielded 3 segments / 128 kb vs. 4 segments / 585 kb from the whole-genome 12× subgraph. Pre-filtering reduces subgraph size 4.6× with one fewer fragment.

### Validation (`pipeline/validation/`)

HG002-specific scripts validating k=21 results against two independent subgraph methods:

| Script | What it does |
|--------|-------------|
| `run_HG002_CYP2D6_trace.sh` | Read-tracing for HG002 highcov: Q100 BAM → A-lines → unitigs |
| `run_HG002_CYP2D6_12x_trace.sh` | Same for 12× assembly |
| `run_HG002_CYP2D6_boundary_v2.sh` | Boundary-crossing CIGAR method (fwd/rev split) for highcov |
| `run_HG002_CYP2D6_12x_boundary_v2.sh` | Same for 12× assembly |
| `verify_k21_CYP2D6.sh` | Align k=21 subgraph segments to HG002 Q100 CYP2D6; check MAPQ60 coverage |
| `align_k21_to_hg38CYP2D6.sh` | Align subgraph segments to hg38 slice for IGV viewing |

---

## Key Results (CYP2D6, HG002)

| Method | Assembly | Segments | Size | Status |
|--------|----------|----------|------|--------|
| k=15 scoring | highcov | 11 | 1.41 Mbp | ✗ paralog contamination |
| k=15 scoring | 12× | 7 | 0.91 Mbp | ✗ paralog contamination |
| **k=21 scoring** | **highcov** | **2** | **1.14 Mbp** | **✓ Q100 verified (4312/4312 MAPQ60)** |
| **k=21 scoring** | **12×** | **4** | **585 kb** | **✓ matches read-tracing exactly** |
| Read-tracing | highcov | 4 | 7.55 Mbp | ✓ correct reads; chromosome-scale UTGs |
| Boundary v2 | highcov | 2 | 1.14 Mbp | ✓ matches k=21 exactly |

k=21 gives one clean segment per haplotype from the high-coverage assembly. The 12× subgraph is fragmented into 4 segments (585 kb total) — motivating the co-assembly approach.

### HG002 × all 8 genes (k=21, 12× assembly)

| Gene | Segs | Size | Notes |
|------|------|------|-------|
| CYP2D6 | 4 | 585 kb | Q100-verified |
| TERT | 2 | 305 kb | |
| TCF3 | 3 | 190 kb | |
| SLC6A3 | 6 | 281 kb | overlaps TERT on chr5 |
| LMF1 | 4 | 465 kb | |
| IGH | 23 | 1.6 Mb | fragmented — 1.3 Mb gene body |
| KIR | 10 | 442 kb | fragmented — copy-number variable |
| HLA | 48 | 7.3 Mb | large — 7.9 Mb region, extreme diversity |

---

## Adding a New Gene

1. Add a row to `genes.txt` with gene body coordinates (set FLANK = GENE coords for now)
2. Extract the reference slice:
   ```bash
   bash pipeline/reference/extract_region.sh MYGENE chr1:1000000-2000000
   ```
3. Run `ffind.py` to get principled flank coordinates and unique k-mers:
   ```bash
   sbatch pipeline/reference/run_ffind_new_genes.sh   # updates genes.txt
   sbatch pipeline/reference/build_k21_all_genes.sh   # add new gene to array range
   ```
4. Align existing assemblies to the new gene slice:
   ```bash
   bash pipeline/alignment/submit_align_pedigree.sh MYGENE
   bash pipeline/alignment/submit_align_hprc.sh MYGENE
   ```
5. Extract subgraphs:
   ```bash
   sbatch pipeline/subgraph/run_gfa_subset.sh both MYGENE
   ```

---

## Archive

[`archive/`](archive/) contains superseded approaches kept for methodological context:

| File | Why superseded |
|------|---------------|
| `gfa_k15.py` | k=15 scorer — paralog contamination at pharmacogene loci |
| `kmer_screen_12x.py` | 61-mer approach — paralogs still share k-mers even at k=61 |
| `score_subgraph.py` | Diagnostic inspector for existing subgraphs (not in production pipeline) |
| `build_k21_CYP2D6.sh` | CYP2D6-only k=21 build; replaced by multi-gene `build_k21_all_genes.sh` |
| `build_unique_kmers.sh` | k=15 unique k-mer generator; replaced by k=21 equivalent |
| `pedigree_cyp2d6_subsets.sh` | CYP2D6-only pedigree subgraph extraction |
| `merge_hprc_cyp2d6.sh` | CYP2D6-only cohort merge |
| `run_HG002_CYP2D6_kmer_screen.sh` | Orchestrator for the abandoned 61-mer approach |

---

## HPC Notes (Yale Bouchet)

- Login nodes: `login1.bouchet`, `login2.bouchet` — no heavy compute here
- Module load order matters: `awscli` must load before `SAMtools` (fixes OpenSSL conflict)
- Do not `module load miniconda` — conflicts with awscli. Prepend `~/.conda/envs/biotools/bin` to PATH in scripts instead.
- `samtools sort: failed to read header` = empty upstream pipe (OOM or BAM not found)

---

## Citation

> Jay Tummala, Hoon Cho. Privacy-Preserving Genomic Co-Assembly at Complex Loci Using Genome-Unique K-mer Subgraph Extraction. *In preparation*, 2026.

## Contact

Jay Tummala — jay.tummala@yale.edu  
PI: Hoon Cho — hoon.cho@yale.edu  
Yale BIDS / Yale Computer Science
