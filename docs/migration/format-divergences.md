# Output Format Divergences: Legacy CRISPRitz - CRISPRitz+

This document is for users migrating existing pipelines, scripts, or
downstream tooling from the legacy CRISPRitz codebase to CRISPRitz+. It
records every place where an output file's **column set, column order, or
column semantics** differs between the two tools, so that parsers written
against the legacy format can be updated deliberately rather than breaking
silently on a reordered or renamed column.

This document covers *file format* only. It does not cover CLI flag
renames, algorithmic changes, or numeric changes in scores themselves.

---

## 1. Two output formats are available, by design

CRISPRitz+'s `search` command can emit results in either of two
layouts, selected via `OutputFormat`:

| Format name | Purpose |
|---|---|
| `targets` | Legacy-compatible column order, for drop-in compatibility with existing legacy-format parsers. |
| `tsv` | New canonical schema used internally and by `annotate-results` / `generate-report`. |

**If you have existing scripts that parse the legacy `*.targets.txt` /
detailedOutput layout, request `targets` format explicitly.** Do not assume
the default is legacy-compatible — see §3.

---

## 2. `targets` format (legacy-compatible layout)

This format reproduces the historical detailedOutput "targets" column
order exactly:

```
bulge_type  grna  target  chrom  pos  strand  mismatches  bulge_size  total
```

Notes for migrators:

- **Column identity is preserved**, but two columns are *derived*, not
  stored directly, and are computed at write time:
  - `bulge_size` = `bulge_dna + bulge_rna`
  - `total` = `mismatches + bulge_size` (i.e. total edit distance)
- **No CFD column.** This layout does not include a CFD score column,
  matching the legacy targets file. If you previously consumed a
  CFD-scored variant of this file (e.g. a separate `*_Scores.txt` produced
  by a legacy scoring script), that corresponds to the new `tsv` format
  below, not to `targets` — see §3 and §4.
- Field *values* are unchanged from legacy (same alignment-string
  conventions, same bulge-type vocabulary: `X`, `DNA`, `RNA`, `DNA,RNA`).

---

## 3. `tsv` format (new canonical layout) — **not legacy column order**

This is the default/internal schema and the one produced before CFD
scoring is applied:

```
chrom  pos  strand  grna  target  mismatches  bulge_dna  bulge_rna  bulge_type
```

Key divergences from the legacy layout:

- **`bulge_type` moves from first column to last.**
- **`bulge_size` and `total` are no longer emitted as combined columns.**
  They are split back into their components (`bulge_dna`, `bulge_rna`) and
  the total is not written at all — recompute it as
  `mismatches + bulge_dna + bulge_rna` if your tooling needs it.
- **Coordinate/identity columns move to the front** (`chrom, pos, strand`
  first) rather than trailing after the sequence columns.

If your migration target is "same columns, same order as legacy," use
`targets` format (§2), not `tsv`.

---

## 4. CFD score column placement

This is the most consequential divergence for scoring pipelines.

- In CRISPRitz+, the C++ search layer **never emits a CFD column**.
  `to_tsv_row()` / `TsvFormatter` deliberately stop at `bulge_type`.
- CFD scoring is a **separate pass**, applied in Python
  (`scores/shard_scoring.py`) after search, which fills in a `cfd_score`
  column that is always the **last column**, appended after `bulge_type`:

```
chrom  pos  strand  grna  spacer  mismatches  bulge_type  bulge_dna  bulge_rna  cfd_score
```

  (Note `grna`/`target` is renamed to `spacer` at this stage — see §5.)

- This 10-column, CFD-last layout is the schema shared by:
  - the per-partition shard files written by the C++ search executor,
  - `result_merger` (the k-way merge that assembles the final table — sort
    order also depends on `cfd_score`: EditDistance mode sorts CFD
    **descending**, with the `NA` sentinel treated as lowest and sorted
    last within a tie group),
  - `annotate-results`' expected input schema (`SEARCH_OUTPUT_HEADER`).

- **If your legacy pipeline expects CFD in a different position** (e.g.
  immediately after `total`, or as part of the `targets`-style bulge-first
  layout), you must re-slice/reorder columns yourself — CRISPRitz+ does
  not offer a formatter that combines the legacy `targets` column order
  with a trailing CFD column.

- **Do not rely on `crispritz_plus.offtarget.TSV_HEADER`** (the Python
  `OffTarget` class's own header constant) as the source of truth for
  scored-TSV column order. That constant currently lists `bulge_dna,
  bulge_rna, bulge_type, cfd_score` (bulge_type *before* cfd_score, but
  *after* the individual bulge columns) — this does not match the order
  actually written by the C++ formatter / shard scorer / result merger
  (`bulge_type, bulge_dna, bulge_rna, cfd_score`). This is a tracked
  internal inconsistency, not an intentional format. Until it is resolved,
  treat `annotation.SEARCH_OUTPUT_HEADER` and the shard column contract
  above as authoritative for what is actually on disk.

---

## 5. Column renames

| Legacy name | CRISPRitz+ name | Where |
|---|---|---|
| `target` / `DNA` (alignment string) | `target` in raw search TSV, renamed to `spacer` once CFD scoring is applied | `tsv` → scored TSV |

If you grep output headers for `DNA` or a legacy alignment-column name,
update to `target` (pre-scoring) or `spacer` (post-scoring, and in
`annotate-results` / `generate-report` output).

---

## 6. Annotated results (`annotate-results`) format

`annotate-results` validates its input against `SEARCH_OUTPUT_HEADER`:

```
chrom  pos  strand  grna  spacer  mismatches  bulge_type  bulge_dna  bulge_rna  cfd_score
```

i.e. it expects the **scored `tsv` layout** (§4), not the legacy `targets`
layout. Feeding it a `targets`-format file will fail validation with a
missing-column error. Annotation (BED-track overlap) columns are appended
to the right of this fixed schema, in the order the `--annotations` BED
files were supplied.

---

## 7. Profile files (`.profile.xls` family)

The profile-file *contents and layout* (per-guide mismatch/bulge tables,
`.profile.xls`, `.extended_profile.xls`, `.profile_dna.xls`,
`.profile_rna.xls`, `.profile_complete.xls`) are intentionally byte-layout
compatible with the legacy `searchOnTST` / `saveProfileGuide` output — no
divergence is tracked here. If you discover a mismatch in profile-file
layout during migration, treat it as a bug rather than an intentional
change and report it; it is not a documented format change.
