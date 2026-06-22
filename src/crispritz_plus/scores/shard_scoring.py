"""Per-shard CFD scoring for the search pipeline.

Scoring is done *from the file*: the C++ search layer writes each partition's
off-target table to its own shard file (cfd_score column left as the ``NA``
sentinel), then this module fills that column in place. Each shard is scored in
its own worker **process** via :func:`score_shards`, so the pure-Python CFD
computation is not serialised by the GIL. The CFD model pickles are loaded once
per worker process (not per shard, not per row) through the pool initializer.

Column contract (shared with the C++ shard writer and ``offtarget.OffTarget``):

    0 chrom   1 pos   2 strand   3 grna   4 spacer
    5 mismatches   6 bulge_type   7 bulge_dna   8 bulge_rna   9 cfd_score

The scorer only reads columns 3 (``grna``) and 4 (``spacer``) and writes column
9 (``cfd_score``); every other field is passed through untouched.
"""

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Optional, Tuple

from .cfd import compute_cfd, load_mismatch_pam_scores

# --- Shard column contract (keep in sync with the C++ writer) ----------------
GRNA_COL: int = 3
SPACER_COL: int = 4
CFD_COL: int = 9
N_COLS: int = 10
_SCORE_NA: str = "NA"

# --- Per-process model cache, populated by _init_worker() --------------------
_MMSCORES: Optional[dict] = None
_PAMSCORES: Optional[dict] = None


def _init_worker(debug: bool) -> None:
    """ProcessPoolExecutor initializer: load the CFD models once per process."""
    global _MMSCORES, _PAMSCORES
    _MMSCORES, _PAMSCORES = load_mismatch_pam_scores(debug)


def _ensure_models(debug: bool) -> Tuple[dict, dict]:
    """Return the (mmscores, pamscores) pair, loading them if not cached.

    Supports calling :func:`score_shard_file` directly (outside a pool), where
    the initializer has not run.
    """
    global _MMSCORES, _PAMSCORES
    if _MMSCORES is None or _PAMSCORES is None:
        _MMSCORES, _PAMSCORES = load_mismatch_pam_scores(debug)
    return _MMSCORES, _PAMSCORES


def _looks_like_header(fields: List[str]) -> bool:
    """A header line has a non-numeric ``pos`` column; a data line does not."""
    return len(fields) > 1 and not fields[1].lstrip("-").isdigit()


def score_shard_file(shard_path: str, debug: bool = False) -> int:
    """Score one shard file in place; return the number of data rows scored.

    The shard is streamed line-by-line (bounded memory) into a temporary file
    that atomically replaces the original. The header row, if present, is passed
    through unchanged.

    Parameters
    ----------
    shard_path:
        Path to a single shard table written by the C++ search layer.
    debug:
        Propagated to the CFD model loader for verbose error handling.

    Returns
    -------
    int
        Count of data rows that received a CFD score.
    """
    mmscores, pamscores = _ensure_models(debug)
    tmp_path = f"{shard_path}.scored.tmp"
    scored = 0
    with open(shard_path, "r") as fin, open(tmp_path, "w") as fout:
        first = True
        for raw in fin:
            line = raw.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if first:
                first = False
                if _looks_like_header(fields):
                    fout.write(line + "\n")
                    continue
            # A shard row may arrive without the trailing cfd_score column;
            # pad to the canonical width so column 9 always exists.
            if len(fields) < N_COLS:
                fields += [_SCORE_NA] * (N_COLS - len(fields))
            spacer = fields[SPACER_COL]
            score = compute_cfd(
                fields[GRNA_COL], spacer, spacer[-2:], mmscores, pamscores
            )
            fields[CFD_COL] = f"{score:.2f}"
            fout.write("\t".join(fields) + "\n")
            scored += 1
    os.replace(tmp_path, shard_path)
    return scored


def score_shards(
    shard_paths: List[str], threads: int, debug: bool = False
) -> int:
    """Score every shard concurrently — one process task per shard.

    Parameters
    ----------
    shard_paths:
        Shard files to score in place.
    threads:
        Maximum number of worker processes.
    debug:
        Propagated to the workers.

    Returns
    -------
    int
        Total number of data rows scored across all shards.
    """
    if not shard_paths:
        return 0
    workers = max(1, min(threads, len(shard_paths)))
    total = 0
    with ProcessPoolExecutor(
        max_workers=workers, initializer=_init_worker, initargs=(debug,)
    ) as pool:
        futures = {
            pool.submit(score_shard_file, path, debug): path
            for path in shard_paths
        }
        for future in as_completed(futures):
            total += future.result()
    return total