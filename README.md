# Privacy-Preserving Genomic Co-Assembly

Some genomic regions are medically critical but too complex to assemble well from a single institution's data. This project builds a system that lets multiple institutions identify shared sequence at these regions and cooperate to assemble them better — without any raw genomic data ever crossing institutional boundaries.

**PI:** Hoon Cho (Yale BIDS / Yale CS) &nbsp;·&nbsp; **Student:** Jay Tummala (Yale CS)

---

## What is this?

Many important genes — drug-metabolism enzymes, immune receptor clusters, disease-linked structural variants — sit in regions of the genome so repetitive and structurally variable that assembling them accurately requires more samples than any single institution has. The natural solution is collaboration across hospitals and research centers, but raw genomic reads are personally identifiable and legally protected. You can't just share them.

This project's answer: compress each institution's reads into compact mathematical fingerprints, add carefully calibrated noise for differential privacy, and use a cryptographic protocol (MHE-PSI) to find which samples share sequence at a target gene — without either institution ever seeing the other's data. Once matched, they co-assemble the region from pooled reads, recovering assembly quality that neither could achieve alone.

The longer-term goal is a **privacy-preserving sequence index** queryable during assembly itself: rather than pooling reads, an assembler could query what other institutions observed at ambiguous positions to resolve gaps and avoid discarding real sequence as apparent errors.

This repository contains the **assembly-side infrastructure**: genome-unique k-mer subgraph extraction, validation experiments, and the full pipeline for 58 samples × 8 complex gene loci. The privacy layer (sketching → DP → MHE-PSI) builds on top of this.

---

## Project Status

| Stage | Status |
|-------|--------|
| Reference preparation (8-gene panel, k=21 k-mer files) | ✅ Complete |
| 12× hifiasm assembly — 58 samples (11 pedigree + 47 HPRC) | ✅ Complete |
| Gene-region subgraph extraction — k=21 method | ✅ Complete |
| HG002 validation against Q100 ground truth | ✅ Complete |
| Ground-truth cohort alignment BAMs (IGV) | ✅ Complete |
| Subgraph extraction — full cohort scale (58 × 8 genes) | 🔄 In progress |
| Minimizer sketching of subgraphs | 📋 Planned |
| Differential privacy layer (randomized response) | 📋 Planned |
| MHE-PSI inter-institutional comparison | 📋 Planned |
| Co-assembly evaluation (NGA50 / QV vs. epsilon) | 📋 Planned |

---

## Key Results

The core challenge was extracting just the gene-region segments from a whole-genome assembly graph containing tens of thousands of segments. The approach: score every segment by how many genome-unique k-mers from the target gene it contains, and keep only the high-scoring ones.

**The k=15 → k=21 pivot.** Initial 15-mer scoring produced 7–11 segments per assembly at CYP2D6 — too many, and Bandage visualization showed disconnected components instead of the expected two clean haplotypes. The cause: CYP2D6 has two near-identical pseudogenes (CYP2D7/CYP2D8) sharing >90% sequence identity, so 15-mers from CYP2D6 recurred in pseudogenes and were incorrectly labeled unique. At k=21, a single nucleotide difference anywhere in the 21-mer window makes it distinct — enough to separate CYP2D6 from its paralogs completely.

**CYP2D6 method comparison (HG002):**

| Method | Assembly | Segments | Size | Result |
|--------|----------|----------|------|--------|
| k=15 scoring | high-cov | 11 | 1.41 Mb | ✗ paralog contamination |
| k=15 scoring | 12× | 7 | 0.91 Mb | ✗ paralog contamination |
| **k=21 scoring** | **high-cov** | **2** | **1.14 Mb** | **✓ Q100 verified — 1 segment per haplotype** |
| **k=21 scoring** | **12×** | **4** | **585 kb** | **✓ matches read-tracing exactly** |
| Read-tracing (A-lines) | high-cov | 4 | 7.55 Mb | ✓ correct reads, chromosome-scale UTGs |
| Boundary v2 (CIGAR) | high-cov | 2 | 1.14 Mb | ✓ matches k=21 exactly |

Three independent methods (k=21 scoring, A-line read tracing, boundary-crossing CIGAR parsing) agree on the same 2 segments for the high-coverage assembly, each confirmed against the HG002 Q100 ground-truth reference.

The 12× assembly fragments to 4 segments — illustrating exactly the problem co-assembly is meant to solve.

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

Flank windows were computed by `src/ffind.py` (k=15, t=1, m=10,000): each side contains ≥10,000 genome-unique anchor k-mers.

---

## Datasets

Both datasets are publicly available.

**Platinum Pedigree (CEU family 1463)** — 11 samples  
NA12877, NA12878, NA12879, NA12881, NA12882, NA12885, NA12886, NA12889, NA12890, NA12891, NA12892
- HiFi reads (GRCh38 haplotagged BAMs): `s3://platinum-pedigree-data/data/hifi/mapped/GRCh38/<sample>.GRCh38.haplotagged.bam`
- High-coverage assemblies: `s3://platinum-pedigree-data/assemblies/`

**HPRC Year 1** — 47 samples (listed in `pipeline/assembly/hprc_subsample_hifiasm.sh`)
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
| [awscli](https://aws.amazon.com/cli/) | 2.17 | S3 streaming |
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

## Pipeline Overview

```
hg38 reference
    ├─ extract_region.sh        → per-gene reference FASTA (hg38CYP2D6.fa, etc.)
    ├─ build_jellyfish_db.sh    → hg38.k21.jf  (22 GB; one-time)
    └─ build_k21_all_genes.sh   → hg38<GENE>.unique_k21.txt × 8 genes

HiFi BAM (S3)
    └─ hprc_subsample_hifiasm.sh → 12× FASTQ → hifiasm → p_utg.gfa + p_ctg.gfa

p_utg.gfa × 58 samples × 8 genes
    └─ run_gfa_subset.sh        → <sample>.<GENE>.k21.subset.gfa
                                   (gene-specific subgraph, ready for sketch/PSI)

High-cov assemblies (ground truth)
    ├─ submit_align_{pedigree,hprc}.sh → assembly-to-gene BAMs
    └─ submit_merge.sh          → cohort_bams/<dataset>_<GENE>_complete.bam (IGV)

[Coming] subgraph → minimizer sketch → DP noise → MHE-PSI → co-assembly
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

### 2. Prepare reference files (one-time, ~2–3 hrs)
```bash
bash pipeline/reference/extract_region.sh CYP2D6 chr22:41983509-42269959
sbatch pipeline/reference/build_jellyfish_db.sh       # produces hg38.k21.jf (~22 GB)
sbatch pipeline/reference/build_k21_all_genes.sh      # produces hg38<GENE>.unique_k21.txt × 8
```

### 3. Run 12× assemblies
```bash
bash pipeline/assembly/submit_hprc_hifiasm.sh         # 47 HPRC samples from S3
sbatch pipeline/assembly/run_hifiasm_12x.sh           # NA12878 (pedigree)
sbatch pipeline/assembly/run_hifiasm_12x_array.sh     # other 10 pedigree samples
```

### 4. Extract gene subgraphs for all samples
```bash
sbatch pipeline/subgraph/run_gfa_subset.sh hprc       # 47 HPRC × 8 genes
sbatch pipeline/subgraph/run_gfa_subset.sh pedigree   # 11 pedigree × 8 genes
```

Output: `<sample>.<GENE>.k21.subset.gfa` alongside each sample's `p_utg.gfa`.

### 5. Build cohort alignment BAMs (ground truth / IGV)
```bash
bash pipeline/alignment/submit_align_pedigree.sh
bash pipeline/alignment/submit_align_hprc.sh
bash pipeline/cohort/submit_merge.sh
```

---

## Python Tools (`src/`)

### `ffind.py` — Flanking anchor finder

Walks outward from a gene body until each flank contains ≥m k-mers that appear at most t times genome-wide. Emits the unique k-mers for subgraph scoring.

> **Key fix:** Jellyfish built without `-C` (canonical flag) stores only one strand. Querying only the forward k-mer gives false positives — k-mers whose reverse complement is common genome-wide look unique. Fix: query both strands and require `fwd_count + RC_count ≤ t`.

```bash
python src/ffind.py reference/hg38.fasta chr22 42126499 42130810 \
    -k 21 -t 1 -m 10000 \
    --db reference/hg38.k21.jf \
    --output-kmers reference/hg38CYP2D6.unique_k21.txt
```

---

### `gfa_k21.py` — GFA subgraph scorer

Scores every segment in a hifiasm `p_utg.gfa` against the genome-unique 21-mer set. A segment passes if `hits / |ref_kmers| × 100 > threshold` (default 1%). Passing IDs pipe to `gfatools` for subgraph extraction.

```bash
python src/gfa_k21.py assembly.bp.p_utg.gfa \
    --kmers reference/hg38CYP2D6.unique_k21.txt \
    > passing_segments.txt

gfatools view -l @passing_segments.txt -r 0 assembly.bp.p_utg.gfa > CYP2D6_subgraph.gfa
```

---

### `filter_reads_by_kmers.py` — Read pre-filter for targeted assembly

Keeps only reads containing ≥1 genome-unique k-mer. Enables targeted hifiasm on just the gene-matching reads. At 12×, ~0.016% of whole-genome HiFi reads match CYP2D6 (358 of 2.2M for HG002).

**Result:** targeted assembly (358 reads) → 3 segments / 128 kb vs. whole-genome baseline of 4 segments / 585 kb. Smaller, cleaner input for both co-assembly and privacy sketching.

```bash
python src/filter_reads_by_kmers.py \
    --fastq reads.12x.fastq.gz \
    --kmers reference/hg38CYP2D6.unique_k21.txt \
    --out filtered_CYP2D6.fastq.gz
```

---

### `trace_reads.py` — A-line read tracer (validation only)

Alternative extraction via hifiasm's internal A-lines: CIGAR offsets → p_ctg A-lines → HiFi read names → p_utg A-lines → unitig names. Used to independently validate k=21 results for HG002. Requires a Q100 reference, so not generalizable to arbitrary samples.

```bash
python src/trace_reads.py \
    --bam hg002.hap1_to_Q100_pat.MAPQ60.bam \
    --ctg-gfa HG002.hap1.p_ctg.gfa \
    --utg-gfa HG002.p_utg.gfa \
    > utg_names.txt
```

---

## Pipeline Scripts

### Reference (`pipeline/reference/`)
| Script | What it does |
|--------|-------------|
| `extract_region.sh` | Cut a named gene slice from hg38 with samtools faidx |
| `build_jellyfish_db.sh` | Build genome-wide k=21 Jellyfish DB (~22 GB; one-time) |
| `build_k21_all_genes.sh` | SLURM array: run ffind.py for all 8 genes |
| `run_ffind_new_genes.sh` | Run ffind.py for newly added genes, update genes.txt |

### Assembly (`pipeline/assembly/`)
| Script | What it does |
|--------|-------------|
| `run_hifiasm_12x.sh` | hifiasm whole-genome assembly for NA12878 at 12× |
| `run_hifiasm_12x_array.sh` | SLURM array: remaining 10 pedigree samples |
| `hprc_subsample_hifiasm.sh` | SLURM worker: stream HPRC BAM from S3 → subsample → hifiasm → delete FASTQ |
| `submit_hprc_hifiasm.sh` | Submit array with 5-job concurrency cap |

### Alignment (`pipeline/alignment/`)
| Script | What it does |
|--------|-------------|
| `align_{pedigree,hprc}_gene.sh` | minimap2 asm5: haplotype assemblies → gene reference slice |
| `submit_align_{pedigree,hprc}.sh` | Submit alignment arrays for all genes in genes.txt |
| `verify_{pedigree,hprc}_alignments.sh` | Check BAM sizes + mapped-read counts |

### Subgraph extraction (`pipeline/subgraph/`)
| Script | What it does |
|--------|-------------|
| `run_gfa_subset.sh` | Score all samples × all genes with gfa_k21.py; extract subset GFAs |
| `run_HG002_all_genes_k21.sh` | HG002 × all 8 genes × {12×, highcov} |

`run_gfa_subset.sh` reads gene names from `genes.txt` dynamically and skips existing outputs.

### Cohort merge (`pipeline/cohort/`)
| Script | What it does |
|--------|-------------|
| `merge_cohort.sh` | MAPQ≥21 filter → per-sample hap merge → single cohort BAM |
| `submit_merge.sh` | Submit for all gene × {pedigree, hprc} combinations |

Output: `cohort_bams/<dataset>_<GENE>_complete.bam` — load directly into IGV.

### Targeted assembly (`pipeline/targeted/`)

Experiment: pre-filter reads by gene k-mers, assemble only the filtered reads, compare to whole-genome subgraph baseline.

| Script | What it does |
|--------|-------------|
| `run_HG002_12x_store_fastq.sh` | Stream HG002 HiFi from S3 → permanent 12× FASTQ |
| `run_HG002_CYP2D6_targeted_assembly.sh` | Filter → targeted hifiasm → compare vs. baseline |
| `run_NA12878_CYP2D6_targeted.sh` | Same for NA12878 pedigree |

### Validation (`pipeline/validation/`)

HG002-specific scripts providing two independent checks on k=21 results:

| Script | What it does |
|--------|-------------|
| `run_HG002_CYP2D6_{,12x_}trace.sh` | Read-tracing via A-lines (high-cov and 12×) |
| `run_HG002_CYP2D6_{,12x_}boundary_v2.sh` | Boundary-crossing CIGAR method, fwd/rev corrected |
| `verify_k21_CYP2D6.sh` | Align k=21 subgraph segments to HG002 Q100; check MAPQ60 coverage |
| `align_k21_to_hg38CYP2D6.sh` | Align subgraph segments to hg38 slice for IGV |

---

## Adding a New Gene

1. Add a row to `genes.txt` with gene body coordinates (set FLANK = GENE coords initially)
2. Extract the reference slice:
   ```bash
   bash pipeline/reference/extract_region.sh MYGENE chr1:1000000-2000000
   ```
3. Compute principled flank coordinates and unique k-mers:
   ```bash
   sbatch pipeline/reference/run_ffind_new_genes.sh    # updates genes.txt
   sbatch pipeline/reference/build_k21_all_genes.sh    # extend array range to cover new gene
   ```
4. Align existing assemblies:
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

[`archive/`](archive/) documents superseded approaches — kept to explain the methodological evolution, not for production use.

| File | Why superseded |
|------|---------------|
| `gfa_k15.py` | k=15 scorer — paralog contamination at pharmacogene loci |
| `kmer_screen_12x.py` | 61-mer approach — paralogs still share k-mers even at k=61 |
| `score_subgraph.py` | Diagnostic inspector for existing subgraphs (not in pipeline) |
| `build_k21_CYP2D6.sh` | CYP2D6-only k=21 build; replaced by multi-gene `build_k21_all_genes.sh` |
| `build_unique_kmers.sh` | k=15 unique k-mer generator; replaced by k=21 equivalent |
| `pedigree_cyp2d6_subsets.sh` | CYP2D6-only pedigree subgraph extraction |
| `merge_hprc_cyp2d6.sh` | CYP2D6-only cohort merge |
| `run_HG002_CYP2D6_kmer_screen.sh` | Orchestrator for the abandoned 61-mer approach |

---

## HPC Notes (Yale Bouchet)

- Login nodes: `login1.bouchet`, `login2.bouchet` — no heavy compute
- Module load order matters: `awscli` must load before `SAMtools` (fixes OpenSSL conflict)
- Do not `module load miniconda` — conflicts with awscli; prepend `~/.conda/envs/biotools/bin` to PATH in scripts instead
- `samtools sort: failed to read header` = empty upstream pipe (OOM or missing input BAM)

---

## Citation

> Jay Tummala, Hoon Cho. Privacy-Preserving Genomic Co-Assembly at Complex Loci Using Genome-Unique K-mer Subgraph Extraction. *In preparation*, 2026.

## Contact

Jay Tummala — jay.tummala@yale.edu  
PI: Hoon Cho — hoon.cho@yale.edu  
Yale BIDS / Yale Computer Science
