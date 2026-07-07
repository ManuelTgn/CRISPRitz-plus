# Changelog

All notable changes to **CRISPRitz+** are documented in this file.


## [0.1.0] - 2026-07-03

Initial public release. CRISPRitz+ is a high-throughput, variant-aware pipeline
for in-silico CRISPR off-target site identification, built around a compact
Ternary-Search-Tree (TST) genome index and a C++/OpenMP search core exposed to
Python through pybind11.

### Added

- **Command-line interface** (`crispritz-plus`) exposing five subcommands that
  form the analysis pipeline: `add-variants`, `index-genome`, `search`,
  `annotate-results`, and `generate-report`. Consistent `--threads`,
  `--verbosity`, `--outdir`, and `--debug` options are shared across commands.
- **`add-variants`** — genome enrichment pipeline that integrates VCF variants
  into reference FASTA files. SNPs are encoded with IUPAC ambiguity codes;
  indels are optionally applied individually (`--indels`); variant filtering
  honours `FILTER=PASS` by default with an opt-in `--keep` for all variants.
- **`index-genome`** — builds a compact Ternary-Search-Tree index (`.bin`
  partitions) of all PAM-matching candidate targets, with a configurable bulge
  budget (`--bmax`) for bulge-aware search.
- **`search`** — off-target search over the TST index within a user-defined
  edit budget: mismatches (`--mm`) plus optional DNA/RNA bulges (`--bdna`,
  `--brna`), with `mixed`/`single` bulge co-occurrence modes (`--bulge-mode`).
  Supports optional CFD off-target scoring (`--score`) for SpCas9/xCas9
  NGG-style PAMs and selectable output (`--output-mode targets|profile|both`).
- **`annotate-results`** — appends genomic context to a search targets table,
  adding one column per BED track. Overlaps are resolved through tabix indexes;
  tracks are sorted, bgzipped, and indexed in a temporary location, leaving
  user input files untouched.
- **`generate-report`** — renders graphical reports from an annotated targets
  table: a bulge-aware LogoMaker off-target profile per guide and an annotation
  radar chart whose axes are derived from the annotation columns.
- **C++/OpenMP search core** exposed to Python via pybind11, with GIL-released
  parallel per-partition search, k-way sorted-shard merging, and per-guide
  profile aggregation.
- **Unified verbosity system** (`0` Silent, `1` Normal, `2` Verbose, `3` Debug)
  spanning the Python and C++ layers, with progress bars for long-running
  stages.
- **Packaging & build**: scikit-build-core build backend that compiles the
  native extension automatically on `pip install`, with a C++ CTest suite and a
  Python `pytest` suite gated at a minimum 75% coverage.
- **Documentation**: project README with installation instructions, a
  reproducible end-to-end quickstart, and a full command reference.

[0.1.0]: https://github.com/ManuelTgn/CRISPRitz-plus/releases/tag/v0.1.0
