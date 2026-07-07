<!-- ![GitHub release](https://img.shields.io/github/v/release/ManuelTgn/CRISPRitz-plus) -->
![Status: Alpha](https://img.shields.io/badge/status-alpha-orange)
![Python](https://img.shields.io/badge/python-%E2%89%A53.8-blue?logo=python&logoColor=white)
![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C?logo=cplusplus&logoColor=white)
![CMake](https://img.shields.io/badge/CMake-%E2%89%A53.19-064F8C?logo=cmake&logoColor=white)
![License: AGPL v3](https://img.shields.io/badge/license-AGPL--3.0-green)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey)


# CRISPRitz+

**High-throughput, variant-aware in-silico off-target site identification for CRISPR genome editing.**

CRISPRitz+ finds candidate off-target sites for CRISPR guide RNAs against a
reference — or a variant-enriched — genome. It combines a compact
Ternary-Search-Tree (TST) genome index with a C++/OpenMP search core (exposed to
Python through pybind11) to enumerate every site within a user-defined edit
budget of mismatches plus optional DNA/RNA bulges, then scores, annotates, and
visualises the results.

The tool is organised as **five subcommands** that form a pipeline:

```
  add-variants  ──►  index-genome  ──►  search  ──►  annotate-results  ──►  generate-report
  (optional:                              │                                        │
   build a                                │                                        │
   variant-aware                          ▼                                        ▼
   genome)                          targets.tsv / profiles          annotated.tsv + figures
```

`add-variants` is optional — start from `index-genome` if you only need the
reference genome.

---

## Installation

### Requirements

- **Python** ≥ 3.8
- A **C++17 compiler** (GCC ≥ 9 / Clang ≥ 10)
- **CMake** ≥ 3.19
- **OpenMP** (used by the parallel search core)
- On **macOS**, OpenMP is not bundled with Apple Clang — install it first:
  ```bash
  brew install libomp cmake
  ```

The compiled TST extension is built automatically from source during `pip
install` via [scikit-build-core](https://scikit-build-core.readthedocs.io/).
Python dependencies (`pysam`, `pandas`, `numpy`, `matplotlib`, `logomaker`,
`tqdm`, `colorama`) are installed for you.

### Install from source

```bash
git clone https://github.com/ManuelTgn/CRISPRitz-plus.git
cd CRISPRitz-plus
pip install .
```

For a development install (editable, with linting/test tooling):

```bash
pip install -e ".[dev]"
```

The build also compiles and runs the C++ unit tests (via CTest). Verify the
install:

```bash
crispritz-plus --version
crispritz-plus --help
```

> The console command is **`crispritz-plus`**. The examples below use it; the
> in-app help text abbreviates it to `crispritz`.

---

## Quickstart

A complete, reproducible end-to-end run on tiny toy inputs. Copy-paste the whole
block — it creates the input files, builds an index, searches, annotates, and
renders figures.

```bash
# 0) Work in a scratch directory
mkdir crispritz-demo && cd crispritz-demo

# 1) A one-record reference genome (chr1). It contains one perfect match and
#    one 1-mismatch site for the guide used below, each followed by an NGG PAM.
mkdir genome
cat > genome/chr1.fa <<'EOF'
>chr1
TTTTTTTTTTACGTTGCAAGTCACGATGCATGGTTTTTTTTTTAAGTTGCAAGTCACGATGCAAGGTTTTTTTTTT
EOF

# 2) PAM model: 20 spacer positions + an NGG PAM. The trailing integer is the
#    PAM length (NGG = 3).
echo "NNNNNNNNNNNNNNNNNNNNNGG 3" > pam.txt

# 3) One guide RNA (20 nt spacer + NNN PAM placeholder), one per line.
echo "ACGTTGCAAGTCACGATGCANNN" > guides.txt

# 4) A BED annotation track (any genomic feature you want to overlap targets with).
printf 'chr1\t0\t40\texon\n' > annotation.bed

# ── Pipeline ───────────────────────────────────────────────────────────────

# (index) Build the TST genome index
crispritz-plus index-genome \
    --genome genome \
    --genome-name demo \
    --pam pam.txt \
    --outdir index

# (search) Find off-targets: up to 3 mismatches, CFD scoring on (SpCas9 NGG PAM)
crispritz-plus search \
    --index-genome index \
    --pam pam.txt \
    --guides guides.txt \
    --mm 3 \
    --score \
    --outdir results

# (annotate) Overlap each target with the BED track → adds an annotation column
crispritz-plus annotate-results \
    --targets results/guides.targets.tsv \
    --annotations annotation.bed \
    --annotation-names exon \
    --outdir results

# (report) Off-target profile + annotation radar figures
crispritz-plus generate-report \
    --input results/guides.targets.annotated.tsv \
    --outdir results
```

**What you get in `results/`:**

| File | Produced by | Contents |
|------|-------------|----------|
| `guides.targets.tsv` | `search` | Tab-separated off-target sites table (one row per hit) |
| `guides.profile*.xls` | `search` | Per-guide mismatch/bulge profile summaries |
| `guides.targets.annotated.tsv` | `annotate-results` | The targets table with one appended column per BED track |
| `report*.png` | `generate-report` | LogoMaker off-target profile and annotation radar figures |

> **Variant-aware variant of the quickstart:** to search against a
> variant-enriched genome, run `add-variants --vcf <vcf-dir> --genome genome`
> first, then point `index-genome --genome` at the generated
> `variants_genome/` folder instead of `genome/`.

---

## Command reference

All subcommands share these conventions:

- `--threads N` — number of worker threads; `0` uses all available cores (default `1`).
- `--verbosity {0,1,2,3}` — `0` Silent, `1` Normal, `2` Verbose, `3` Debug (default `1`).
- `--debug` — propagate full Python stack traces on error.
- `--outdir DIR` — output directory (default: current working directory).
- `-h` / `--help` — per-subcommand help.

### 1. `add-variants` — build a variant-aware genome

Parses per-chromosome VCF files and integrates their variants into the matching
reference FASTA files, producing an *enriched* genome. SNPs are encoded with
IUPAC ambiguity codes (representing reference + alternative alleles); indels,
when enabled, are applied individually around each variant position. Output is
written under a `variants_genome/` folder.

| Option | Req. | Default | Description |
|--------|:----:|---------|-------------|
| `--vcf VCF-DIR` | ✔ | — | Directory of per-chromosome VCFs (e.g. `chr1.vcf.gz`, `chr2.vcf.gz`). |
| `--genome FASTA-DIR` | ✔ | — | Directory of per-chromosome reference FASTAs (e.g. `chr1.fa`). |
| `--indels` | | off | Also apply insertions/deletions (each variant individually). |
| `--keep` | | off | Keep all variants regardless of FILTER (default: `FILTER=PASS` only). |
| `--outdir OUTDIR` | | cwd | Where the `variants_genome/` folder is written. |
| `--threads N` | | 1 | Worker threads (`0` = all cores). |
| `--verbosity {0..3}` | | 1 | Output verbosity. |
| `--debug` | | off | Full error traceback. |

```bash
crispritz-plus add-variants --vcf vcfs/ --genome genome/ --indels --outdir enriched/
```

### 2. `index-genome` — build the TST index

Scans the input FASTA files, extracts every candidate target matching the PAM,
and stores them in a compact Ternary-Search-Tree index (`.bin` partitions) used
for fast, optionally bulge-aware, off-target retrieval.

| Option | Req. | Default | Description |
|--------|:----:|---------|-------------|
| `--genome FASTA-DIR` | ✔ | — | Directory of per-chromosome FASTA files (reference or enriched). |
| `--genome-name NAME` | ✔ | — | Identifier used to name the generated index folder. |
| `--pam PAM-FILE` | ✔ | — | PAM model file (see [PAM file format](#pam-file-format)). |
| `--bmax N` | | 0 | Max bulges allowed during index construction (higher = more sensitive, slower). |
| `--outdir OUTDIR` | | cwd | Output directory (defaults to `<GENOME-NAME>_<PAM>_<BMAX>`). |
| `--threads N` | | 1 | Worker threads (`0` = all cores). |
| `--verbosity {0..3}` | | 1 | Output verbosity. |
| `--debug` | | off | Full error traceback. |

```bash
crispritz-plus index-genome --genome genome/ --genome-name hg38 --pam pam.txt --bmax 2
```

### 3. `search` — off-target search

Traverses a pre-computed TST index and enumerates every site whose alignment to
a guide stays within the edit budget (mismatches plus optional DNA/RNA bulges).
The PAM and guide geometry are fixed by the index; this command supplies only
the search tolerances and output options. Results are written as a targets table
and/or per-guide profiles, ready for scoring and annotation.

| Option | Req. | Default | Description |
|--------|:----:|---------|-------------|
| `--index-genome DIR` | ✔ | — | Index directory produced by `index-genome` (contains `.bin` partitions). |
| `--pam PAM-FILE` | ✔ | — | PAM model file — must match the PAM used to build the index. |
| `--guides GUIDES-FILE` | ✔ | — | One guide RNA per line (spacer + PAM placeholder, e.g. `…NNN`). |
| `--mm MISMATCHES` | ✔ | — | Max mismatches (substitutions) per alignment. |
| `--bdna N` | | 0 | Max DNA bulges (≤ the index's bulge budget). |
| `--brna N` | | 0 | Max RNA bulges (≤ the index's bulge budget). |
| `--bulge-mode {mixed,single}` | | mixed | `mixed`: DNA+RNA bulges may co-occur in one alignment; `single`: only one bulge kind per site. |
| `--score` | | off | Compute the CFD off-target score (SpCas9/xCas9 NGG-style PAMs). |
| `--output-mode {targets,profile,both}` | | both | Which files to write: targets table, per-guide profiles, or both. |
| `--outdir OUTDIR` | | cwd | Output directory. |
| `--threads N` | | 1 | Worker threads (`0` = all cores). |
| `--verbosity {0..3}` | | 1 | Output verbosity. |
| `--debug` | | off | Full error traceback. |

```bash
crispritz-plus search --index-genome hg38_index/ --pam pam.txt \
    --guides guides.txt --mm 4 --bdna 1 --brna 1 --score --outdir results/
```

Output: `<guides-stem>.targets.tsv` (targets) and/or `<guides-stem>.profile*.xls`
(profiles).

### 4. `annotate-results` — add genomic context

Appends genomic annotations to a search targets table: for each BED track, one
column is added reporting the features overlapping every off-target site.
Overlaps use a tabix index for fast random access; tracks that are not already
indexed are sorted, bgzipped, and indexed in a temporary location (your input
files are left untouched). Accepts **only** the TSV produced by `search`.

| Option | Req. | Default | Description |
|--------|:----:|---------|-------------|
| `--targets TARGETS-FILE` | ✔ | — | Targets TSV produced by `search`. |
| `--annotations BED [BED …]` | ✔ | — | One or more BED tracks (`.bed` or block-gzipped `.bed.gz`), one column each. |
| `--annotation-names NAME [NAME …]` | | `annotation1…N` | Column names, one per track, in order. |
| `--outdir OUTDIR` | | cwd | Output directory. Result is `<targets-stem>.annotated.tsv`. |
| `--threads N` | | 1 | Threads used to prepare (sort/bgzip/index) tracks (`0` = all cores). |
| `--verbosity {0..3}` | | 1 | Output verbosity. |
| `--debug` | | off | Full error traceback. |

```bash
crispritz-plus annotate-results --targets results/guides.targets.tsv \
    --annotations genes.bed enhancers.bed.gz \
    --annotation-names genes enhancers --outdir results/
```

### 5. `generate-report` — graphical reports

Produces figures from an annotated targets TSV: a bulge-aware LogoMaker
off-target profile per guide, and a radar chart whose axes are read dynamically
from the annotation columns. Consumes **only** the annotated TSV from
`annotate-results`.

| Option | Req. | Default | Description |
|--------|:----:|---------|-------------|
| `--input ANNOTATED-TSV` | ✔ | — | Annotated targets TSV from `annotate-results`. |
| `--mm N` | | totals | Mismatch level to display on the radar (default: totals across all levels). |
| `--guide GUIDE` | | all | Restrict the report to a specific guide (gap-free gRNA+PAM). Repeatable. |
| `--prefix PREFIX` | | `report` | Filename prefix for generated figures. |
| `--outdir OUTDIR` | | cwd | Output directory. |
| `--verbosity {0..3}` | | 1 | Output verbosity. |
| `--debug` | | off | Full error traceback. |

```bash
crispritz-plus generate-report --input results/guides.targets.annotated.tsv \
    --guide ACGTTGCAAGTCACGATGCANNN --prefix demo --outdir results/
```

---

## PAM file format

A PAM file is a single line with two space-separated fields:

```
<pattern> <pam-length>
```

- `<pattern>` — the full spacer+PAM pattern: a run of `N` characters equal to
  the guide length, followed by the PAM motif.
- `<pam-length>` — the number of PAM bases (the trailing bases of the pattern).

Example (20-nt spacer, SpCas9 `NGG` PAM):

```
NNNNNNNNNNNNNNNNNNNNNGG 3
```

Guides in the `--guides` file follow the same geometry: the spacer bases
followed by a PAM-length placeholder, e.g. `ACGTTGCAAGTCACGATGCANNN`.

---

## Development

```bash
pip install -e ".[dev]"
pytest                 # Python test suite (coverage gate: 75%)
```

The C++ core additionally ships CTest unit tests, run automatically during the
build.
