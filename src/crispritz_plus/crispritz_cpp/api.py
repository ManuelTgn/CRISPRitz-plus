""" """

from typing import List, Sequence, Union


from crispritz_plus import _ternary_search_tree as tst  # type: ignore

from .guide_profile import GuideProfile
from .output_mode import OutputMode
from .partition_results import PartitionResult
from .search_configuration import SearchConfiguration
from .sort_mode import SortMode


def build_tree_cpp(
    sequence: str,
    contig: str,
    pam: str,
    pam_length: int,
    pam_size: int,
    upstream: bool,
    outdir: str,
    max_bulges: int = 0,
    threads: int = 1,
) -> None:
    """Call the C++ TST builder for a single chromosome sequence.

    Parameters
    ----------
    sequence: str
        Full genomic sequence (single chromosome, uppercase IUPAC).
    contig: str
        Chromosome / contig identifier used in the output filename(s).
    pam: str
        PAM-only string (e.g. ``"NGG"``), without guide placeholder Ns.
    pam_length: int
        Total length of the PAM + guide pattern.
    pam_size: int
        Length of the PAM portion only.
    upstream: bool
        True when the PAM precedes the guide (e.g. Cas12a).
    outdir: str
        Path to the directory where the genome index will be stored.
    max_bulges: int
        Maximum number of bulges allowed during index construction.
    threads: int
        Number of OpenMP threads for the PAM search phase.
    """
    tst.build_tree(
        sequence,
        contig,
        pam,
        pam_length,
        pam_size,
        upstream,
        outdir,
        max_bulges,
        threads,
    )


def make_search_configuration(
    max_mismatches: int,
    max_bulges_dna: int,
    max_bulges_rna: int,
    threads: int,
    output_mode: Union[str, OutputMode] = "both",
) -> SearchConfiguration:
    return SearchConfiguration(
        max_mismatches, max_bulges_dna, max_bulges_rna, threads, output_mode=output_mode
    )


def run_search_executor_cpp(
    partition_path: str,
    chrom: str,
    guides: List[str],
    config: SearchConfiguration,
    pam_len: int,
    pam_at_start: bool,
    shard_path: str,
) -> PartitionResult:
    return tst.run_search_executor(
        partition_path, chrom, guides, config.native, pam_len, pam_at_start, shard_path
    )


def merge_sorted_shards_cpp(
    shard_paths: List[str],
    final_path: str,
    sort_mode: Union[str, SortMode] = "edit_distance",
    write_header: bool = True,
    remove_inputs: bool = True,
) -> int:
    mode = SortMode.from_string(sort_mode) if isinstance(sort_mode, str) else sort_mode
    return tst.merge_sorted_shards(
        shard_paths, final_path, mode, write_header, remove_inputs
    )


def write_merged_profiles_cpp(
    profiles_by_partition: List[List[GuideProfile]], path_stem: str
) -> int:
    profiles_by_partition_native = [
        [profile.native for profile in partition] for partition in profiles_by_partition
    ]
    return tst.write_merged_profiles(profiles_by_partition, path_stem)
