""" """

from concurrent.futures import ThreadPoolExecutor, as_completed
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
from ..verbosity import VERBOSITY_LVL, print_verbosity


_SHARD_DIRNAME = ".crispritz_shards"


def _query_guides(guides: GuideList, debug: bool) -> List[str]:
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
    workers: int,
    verbosity: int,
    debug: bool,
) -> List[PartitionResult]:
    #  parallel per-partition search (GIL released in C++)
    results: List[PartitionResult] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
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
                pam.size,
                pam.upstream,
                shard,
            )
            future_to_partition[future] = partition
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
    total_hits = sum(r.total_hits for r in results)
    print_verbosity(
        f"Found {total_hits} off-target(s) total", verbosity, VERBOSITY_LVL[1]
    )
    return results


def _score_targets(
    shard_paths: List[str], config: SearchConfiguration, verbosity: int, debug: bool
) -> None:
    pass


def _remove_shard_folder(shard_dir: str, debug: bool) -> None:
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
    output_mode: str = "both",
    sort_mode: str = "edit_distance",
) -> None:
    pam = PAM(pam_file, debug)  # initialize pam
    guides = _query_guides(
        GuideList(guides_file, pam, debug), debug
    )  # initialize guides
    config = make_search_configuration(mm, bdna, brna, threads, output_mode=output_mode)
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
        config, partitions, guides, pam, shard_dir, threads, verbosity, debug
    )
    guides_stem = os.path.splitext(os.path.basename(guides_file))[0]
    # -- phase 2-3: score shards, then sort + k-way merge (targets)
    if config.write_targets:
        shard_paths = [r.shard_path for r in results]
        # if pam.cas_system in [SPCAS9, XCAS9]:
        #     print_verbosity(f"Scored {scored} row(s) across {len(shard_paths)} shard(s)", verbosity, VERBOSITY_LVL[1])
        #     pass
        final_path = os.path.join(outdir, f"{guides_stem}.targets.tsv")
        written = merge_sorted_shards_cpp(shard_paths, final_path, sort_mode)
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
        n_guides = write_merged_profiles_cpp(by_partition, stem)
        print_verbosity(
            f"Wrote profiles for {n_guides} guide(s) to {stem}.profile*.xls",
            verbosity,
            VERBOSITY_LVL[1],
        )
