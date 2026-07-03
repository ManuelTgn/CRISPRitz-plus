"""Free-function entry points into the C++ Ternary Search Tree core.

Each function here is a thin adapter over the compiled
``_ternary_search_tree`` extension: it accepts friendly Python arguments
(lowercase enum tokens, wrapper objects), converts them to the native types the
extension expects, and calls the corresponding C++ routine.  The four stages
covered are index construction (:func:`build_tree_cpp`), per-partition search
(:func:`run_search_executor_cpp`), shard merging
(:func:`merge_sorted_shards_cpp`), and profile writing
(:func:`write_merged_profiles_cpp`); :func:`make_search_configuration` is a
convenience constructor for the search configuration.
"""

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
) -> None:
    """Build a Ternary Search Tree index for one chromosome sequence.

    Adapter over the C++ ``build_tree`` routine.

    Parameters
    ----------
    sequence : str
        Full genomic sequence (single chromosome, uppercase IUPAC).
    contig : str
        Chromosome / contig identifier used in the output filename(s).
    pam : str
        PAM-only string (e.g. ``"NGG"``), without guide placeholder Ns.
    pam_length : int
        Total length of the PAM + guide pattern.
    pam_size : int
        Length of the PAM portion only.
    upstream : bool
        ``True`` when the PAM precedes the guide (e.g. Cas12a).
    outdir : str
        Directory where the genome index will be stored.
    max_bulges : int, optional
        Maximum number of bulges allowed during index construction.
        Defaults to ``0``.
    threads : int, optional
        Number of OpenMP threads for the PAM search phase. Defaults to ``1``.

    Returns
    -------
    None
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
    """Construct a :class:`SearchConfiguration` from search parameters.

    Convenience wrapper that forwards the edit limits, thread count, and output
    mode to the :class:`SearchConfiguration` constructor (using the default
    ``"tsv"`` output format).

    Parameters
    ----------
    max_mismatches : int
        Maximum number of mismatches permitted in a hit.
    max_bulges_dna : int
        Maximum number of DNA bulges permitted in a hit.
    max_bulges_rna : int
        Maximum number of RNA bulges permitted in a hit.
    threads : int
        Number of threads the C++ search may use.
    output_mode : Union[str, OutputMode], optional
        Which result files to produce, as a token or an :class:`OutputMode`
        value. Defaults to ``"both"``.

    Returns
    -------
    SearchConfiguration
        The constructed configuration.
    """
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
) -> PartitionResult:
    """Run the C++ search executor over a single index partition.

    Resolves *bulge_mode* (token or enum) to the native value, unwraps the
    native ``SearchConfiguration``, and invokes the C++ executor.

    Parameters
    ----------
    partition_path : str
        Path to the ``.bin`` index partition to search.
    chrom : str
        Contig name associated with the partition.
    guides : List[str]
        Guide bodies to search for.
    config : SearchConfiguration
        Search parameters; its native object is passed to C++.
    pam : str
        PAM-only string (e.g. ``"NGG"``).
    pam_at_start : bool
        ``True`` when the PAM precedes the guide.
    shard_path : str
        Output shard path for this partition's targets.
    bulge_mode : Union[str, BulgeMode], optional
        Whether DNA and RNA bulges may be combined, as a token or a
        :class:`BulgeMode` value. Defaults to ``"mixed"``.

    Returns
    -------
    PartitionResult
        The native per-partition result (counts and per-guide profiles).

    Raises
    ------
    TSTSearchError
        If *bulge_mode* is an unrecognised token, or the C++ search fails.
    """
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
    )


def merge_sorted_shards_cpp(
    shard_paths: List[str],
    final_path: str,
    sort_mode: Union[str, SortMode] = "edit_distance",
    write_header: bool = True,
    remove_inputs: bool = True,
) -> int:
    """Merge per-partition target shards into one sorted table.

    Resolves *sort_mode* (token or enum) to the native value and invokes the
    C++ merge routine.

    Parameters
    ----------
    shard_paths : List[str]
        Paths of the per-partition shard files to merge.
    final_path : str
        Destination path for the merged table.
    sort_mode : Union[str, SortMode], optional
        Ordering applied to the merged table, as a token or a
        :class:`SortMode` value. Defaults to ``"edit_distance"``.
    write_header : bool, optional
        When ``True``, write a header row to *final_path*. Defaults to ``True``.
    remove_inputs : bool, optional
        When ``True``, delete the input shards after a successful merge.
        Defaults to ``True``.

    Returns
    -------
    int
        The number of rows written to *final_path*.

    Raises
    ------
    TSTSearchError
        If *sort_mode* is an unrecognised token, or the C++ merge fails.
    """
    mode = SortMode.from_string(sort_mode) if isinstance(sort_mode, str) else sort_mode
    return tst.merge_sorted_shards(
        shard_paths, final_path, mode, write_header, remove_inputs
    )


def write_merged_profiles_cpp(
    profiles_by_partition: List[List[GuideProfile]], path_stem: str
) -> int:
    """Merge per-partition guide profiles and write them to disk.

    Unwraps every :class:`GuideProfile` to its native object (idempotently, so
    already-native entries pass through) before invoking the C++ profile
    writer.

    Parameters
    ----------
    profiles_by_partition : List[List[GuideProfile]]
        Per-partition lists of guide profiles to merge.
    path_stem : str
        Output path stem for the merged profile file(s).

    Returns
    -------
    int
        The number of merged profiles written.
    """
    profiles_by_partition_native = [
        [getattr(profile, "native", profile) for profile in partition]
        for partition in profiles_by_partition
    ]
    return tst.write_merged_profiles(profiles_by_partition_native, path_stem)
