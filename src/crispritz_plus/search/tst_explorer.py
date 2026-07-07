"""Off-target search orchestration over Ternary Search Tree partitions.

Python driver for the ``search`` stage.  It parses the PAM and guides, builds a
:class:`~crispritz_plus.crispritz_cpp.search_configuration.SearchConfiguration`,
then runs the pipeline in phases:

1. **Parallel per-partition search** — each ``.bin`` partition is searched by
   the C++ executor on a :class:`~concurrent.futures.ThreadPoolExecutor`
   worker (the GIL is released inside the C++ call), each writing its own shard
   file.
2. **Scoring** *(optional)* — when enabled for an SpCas9/xCas9 PAM, the shard
   ``cfd_score`` columns are filled via
   :func:`~crispritz_plus.scores.shard_scoring.score_shards`.
3. **Sort + k-way merge** — the per-partition shards are merged into a single
   sorted targets table by the C++ merger, and the temporary shard directory
   is removed.
4. **Profile merge** *(optional)* — per-partition guide profiles are merged and
   written to the profile files.

Concurrency is partition-level, matching the workload shape (few guides, many
partitions); the C++ search releases the GIL so the thread pool achieves real
parallelism.

Module-level constants
----------------------
_SHARD_DIRNAME : str
    Name of the temporary per-run directory that holds the per-partition shard
    files before they are merged (``".crispritz_shards"``).
"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from time import time
from typing import List

import os


from ..crispritz_cpp import (
    PartitionResult,
    SearchConfiguration,
    make_search_configuration,
    merge_sorted_shards_cpp,
    run_search_executor_cpp,
    write_merged_profiles_cpp,
)
from ..crispritz_errors import CrispritzTstError, CrispritzSearchError
from ..exception_handlers import exception_handler
from ..guide import GuideList
from ..pam import PAM, SPCAS9, XCAS9
from ..progress import progress_bar_parallel
from ..scores import score_shards
from ..verbosity import VERBOSITY_LVL, print_verbosity


# ==========================================================================
# Shard scoring temporary folder
# ==========================================================================
_SHARD_DIRNAME = ".crispritz_shards"


# ==========================================================================
# Internal helpers
# ==========================================================================


def _query_guides(guides: GuideList, debug: bool) -> List[str]:
    """Return the spacer sequences of every parsed guide.

    Parameters
    ----------
    guides : GuideList
        The parsed guides.
    debug : bool
        When *True*, errors propagate with a full traceback.

    Returns
    -------
    List[str]
        One spacer sequence per guide, in input order.

    Raises
    ------
    CrispritzSearchError
        If the guide sequences cannot be extracted.
    """
    try:
        return [g.sequence for g in guides.guides]
    except Exception as e:
        exception_handler(
            CrispritzSearchError,
            "Invalid query guides parsed from the guide(s) file",
            os.EX_DATAERR,
            debug,
            e,
        )


def _contig_from_partition(partition_path: str, debug: bool) -> str:
    """Recover the contig name encoded in a TST partition filename.

    Partition filenames follow a ``{prefix}_{contig}_{suffix}.bin`` convention
    (the contig itself may contain underscores).  The ``.bin`` extension is
    stripped and the contig is taken as everything between the first and last
    underscore-delimited tokens.

    Parameters
    ----------
    partition_path : str
        Path to the ``.bin`` partition file.
    debug : bool
        When *True*, errors propagate with a full traceback.

    Returns
    -------
    str
        The recovered contig name.

    Raises
    ------
    CrispritzTstError
        If the filename has too few tokens to recover a contig name.
    """
    stem = os.path.basename(partition_path)
    if stem.endswith(".bin"):
        stem = stem[:-4]  # remove bin extension
    tokens = stem.split("_")
    if len(tokens) < 3:
        exception_handler(
            CrispritzTstError,
            f"Invalid TST partition name: {partition_path}. Cannot recover contig name",
            os.EX_DATAERR,
            debug,
        )
    return "_".join(tokens[1:-1])


def _search_parallel(
    config: SearchConfiguration,
    partitions: List[str],
    guides: List[str],
    pam: PAM,
    shard_dir: str,
    bulge_mode: str,
    workers: int,
    verbosity: int,
    debug: bool,
) -> List[PartitionResult]:
    """Search every partition in parallel and return the per-partition results.

    Submits one C++ search task per partition to a thread pool (the C++ call
    releases the GIL, so the searches run concurrently), assigning each a shard
    output path when targets are enabled, and aggregates the results.

    Parameters
    ----------
    config : SearchConfiguration
        Search parameters shared by all partitions.
    partitions : List[str]
        Paths to the ``.bin`` partitions to search.
    guides : List[str]
        Guide spacer sequences to search for.
    pam : PAM
        PAM geometry (motif and orientation) passed to the executor.
    shard_dir : str
        Directory for the per-partition shard files.
    bulge_mode : str
        Whether DNA and RNA bulges may be combined.
    workers : int
        Maximum number of worker threads.
    verbosity : int
        Verbosity level (see
        :data:`~crispritz_plus.verbosity.VERBOSITY_LVL`).
    debug : bool
        When *True*, errors propagate with a full traceback.

    Returns
    -------
    List[PartitionResult]
        One result per partition (order not guaranteed).

    Raises
    ------
    CrispritzSearchError
        If the search fails on any partition.
    """
    #  parallel per-partition search (GIL released in C++)
    print_verbosity(
        f"Dispatching {len(partitions)} partition(s) across {workers} worker "
        f"thread(s)",
        verbosity,
        VERBOSITY_LVL[3],
    )
    start = time()  # track parallel search phase run time
    results: List[PartitionResult] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        with progress_bar_parallel(
            len(partitions), "Searched partitions", verbosity
        ) as pbar:
            future_to_partition = {}
            for partition in partitions:
                shard = (
                    os.path.join(shard_dir, os.path.basename(partition) + ".shard.tsv")
                    if config.write_targets
                    else ""
                )
                future = pool.submit(
                    run_search_executor_cpp,
                    partition,
                    _contig_from_partition(partition, debug),
                    guides,
                    config,
                    pam.pamseq,
                    pam.upstream,
                    shard,
                    bulge_mode,
                    verbosity,
                )
                future_to_partition[future] = partition
                pbar.update(1)
        for future in as_completed(future_to_partition):
            partition = future_to_partition[future]
            try:
                results.append(future.result())
            except Exception as e:
                exception_handler(
                    CrispritzSearchError,
                    f"Failed search on partition {partition}",
                    os.EX_USAGE,
                    debug,
                    e,
                )
    print_verbosity(
        f"Partition search completed in {time() - start:.2f}s",
        verbosity,
        VERBOSITY_LVL[3],
    )
    total_hits = sum(r.total_hits for r in results)
    print_verbosity(
        f"Found {total_hits} off-target(s) total", verbosity, VERBOSITY_LVL[1]
    )
    return results


def _remove_shard_folder(shard_dir: str, debug: bool) -> None:
    """Remove the shard directory once it is empty.

    The k-way merge deletes the shard files as it consumes them; this drops the
    now-empty directory.  A non-empty directory is treated as an error.

    Parameters
    ----------
    shard_dir : str
        Path to the temporary shard directory.
    debug : bool
        When *True*, errors propagate with a full traceback.

    Returns
    -------
    None

    Raises
    ------
    CrispritzSearchError
        If the directory still contains files and cannot be removed.
    """
    try:  # merge removed the shards; drop the now-empty shard directory
        if os.path.isdir(shard_dir) and not os.listdir(shard_dir):
            os.rmdir(shard_dir)
    except OSError as e:
        exception_handler(
            CrispritzSearchError,
            "Non-empty shard folder, cannot remove",
            os.EX_OSERR,
            debug,
            e,
        )


# ==========================================================================
# Public API
# ==========================================================================


def search_offtargets_tst(
    indexes: List[str],
    pam_file: str,
    guides_file: str,
    mm: int,
    bdna: int,
    brna: int,
    outdir: str,
    threads: int,
    verbosity: int,
    debug: bool,
    bulge_mode: str = "mixed",
    output_mode: str = "both",
    sort_mode: str = "edit_distance",
    score: bool = False,
) -> None:
    """Search for off-targets across TST partitions and write the results.

    Top-level orchestrator for the ``search`` stage.  Parses the PAM and
    guides, builds the search configuration, then runs the search in phases:
    parallel per-partition search, optional CFD scoring (SpCas9/xCas9 only),
    sort + k-way merge into a single targets table, and optional profile
    merging — emitting progress messages throughout.

    Parameters
    ----------
    indexes : List[str]
        Paths to the ``.bin`` TST partitions to search.
    pam_file : str
        Path to the PAM specification file.
    guides_file : str
        Path to the guides file (one guide per line).
    mm : int
        Maximum number of mismatches permitted in a hit.
    bdna : int
        Maximum number of DNA bulges permitted in a hit.
    brna : int
        Maximum number of RNA bulges permitted in a hit.
    outdir : str
        Directory for the output table, profiles, and temporary shards.
    threads : int
        Maximum number of worker threads.
    verbosity : int
        Verbosity level (see
        :data:`~crispritz_plus.verbosity.VERBOSITY_LVL`).
    debug : bool
        When *True*, errors propagate with a full traceback.
    bulge_mode : str, optional
        Whether DNA and RNA bulges may be combined. Defaults to ``"mixed"``.
    output_mode : str, optional
        Which result files to produce. Defaults to ``"both"``.
    sort_mode : str, optional
        Ordering applied to the merged table. Defaults to
        ``"edit_distance"``.
    score : bool, optional
        When ``True`` (and the PAM is SpCas9/xCas9), fill the CFD score column
        before merging. Defaults to ``False``.

    Returns
    -------
    None

    Raises
    ------
    CrispritzSearchError
        If the search, scoring, or merge fails.
    CrispritzTstError
        If a partition filename cannot be parsed for its contig name.
    """
    pam = PAM(pam_file, debug)  # initialize pam
    guides = _query_guides(
        GuideList(guides_file, pam, debug), debug
    )  # initialize guides
    config = make_search_configuration(mm, bdna, brna, threads, output_mode=output_mode)
    print_verbosity(
        f"Loaded {len(guides)} guide(s); PAM={pam.pamseq}, edit budget: {mm} "
        f"mismatch(es) + {bdna} DNA bulge(s) + {brna} RNA bulge(s)",
        verbosity,
        VERBOSITY_LVL[2],
    )
    print_verbosity(
        f"Output mode: targets={config.write_targets}, "
        f"profiles={config.write_profile}",
        verbosity,
        VERBOSITY_LVL[2],
    )
    partitions = sorted(indexes)  # sort lexicographically
    shard_dir = os.path.join(outdir, _SHARD_DIRNAME)
    if config.write_targets:
        os.makedirs(shard_dir, exist_ok=True)
    print_verbosity(
        f"Searching {len(guides)} guide(s) across {len(partitions)} partitions "
        f"using {threads} thread(s)",
        verbosity,
        VERBOSITY_LVL[1],
    )
    # -- phase 1: parallel per-partition search
    results = _search_parallel(
        config,
        partitions,
        guides,
        pam,
        shard_dir,
        bulge_mode,
        threads,
        verbosity,
        debug,
    )
    guides_stem = os.path.splitext(os.path.basename(guides_file))[0]
    # -- phase 2-3: score shards, then sort + k-way merge (targets)
    if config.write_targets:
        shard_paths = [r.shard_path for r in results]
        if score and pam.cas_system in [SPCAS9, XCAS9]:
            print_verbosity(
                f"Scoring {len(shard_paths)} shard(s) with CFD",
                verbosity,
                VERBOSITY_LVL[2],
            )
            scored = score_shards(shard_paths, threads, debug)
            print_verbosity(
                f"Scored {scored} row(s) across {len(shard_paths)} shard(s)",
                verbosity,
                VERBOSITY_LVL[1],
            )
        final_path = os.path.join(outdir, f"{guides_stem}.targets.tsv")
        print_verbosity(
            f"Merging {len(shard_paths)} shard(s) into final table "
            f"(sort={sort_mode})",
            verbosity,
            VERBOSITY_LVL[2],
        )
        written = merge_sorted_shards_cpp(
            shard_paths, final_path, sort_mode, verbosity=verbosity
        )
        print_verbosity(
            f"Wrote {written} off-target(s) to {final_path} (sorted by {sort_mode})",
            verbosity,
            VERBOSITY_LVL[1],
        )
        _remove_shard_folder(shard_dir, debug)
    # -- phase 4: merge per-partition profiles and write the files
    if config.write_profile:
        by_partition = [list(r.profiles) for r in results]
        stem = os.path.join(outdir, guides_stem)
        print_verbosity(
            f"Merging per-partition profiles for {len(guides)} guide(s)",
            verbosity,
            VERBOSITY_LVL[2],
        )
        n_guides = write_merged_profiles_cpp(by_partition, stem, verbosity=verbosity)
        print_verbosity(
            f"Wrote profiles for {n_guides} guide(s) to {stem}.profile*.xls",
            verbosity,
            VERBOSITY_LVL[1],
        )
