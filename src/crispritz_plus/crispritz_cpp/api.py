""" """

from typing import List, Union


from crispritz_plus import _ternary_search_tree as tst  # type: ignore

from .bulge_mode import BulgeMode
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
    verbosity: int = 1,
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
    verbosity: int
        Output verbosity level forwarded to the C++ builder (0=Silent,
        1=Normal, 2=Verbose, 3=Debug).
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
        verbosity,
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
    pam: str,
    pam_at_start: bool,
    shard_path: str,
    bulge_mode: Union[str, BulgeMode] = "mixed",
    verbosity: int = 1,
) -> PartitionResult:
    mode = (
        BulgeMode.from_string(bulge_mode) if isinstance(bulge_mode, str) else bulge_mode
    )
    return tst.run_search_executor(
        partition_path,
        chrom,
        guides,
        config.native,
        pam,
        pam_at_start,
        shard_path,
        mode,
        verbosity,
    )


def merge_sorted_shards_cpp(
    shard_paths: List[str],
    final_path: str,
    sort_mode: Union[str, SortMode] = "edit_distance",
    write_header: bool = True,
    remove_inputs: bool = True,
    verbosity: int = 1,
) -> int:
    mode = SortMode.from_string(sort_mode) if isinstance(sort_mode, str) else sort_mode
    return tst.merge_sorted_shards(
        shard_paths, final_path, mode, write_header, remove_inputs, verbosity
    )


def write_merged_profiles_cpp(
    profiles_by_partition: List[List[GuideProfile]],
    path_stem: str,
    verbosity: int = 1,
) -> int:
    profiles_by_partition_native = [
        [getattr(profile, "native", profile) for profile in partition]
        for partition in profiles_by_partition
    ]
    return tst.write_merged_profiles(profiles_by_partition_native, path_stem, verbosity)
