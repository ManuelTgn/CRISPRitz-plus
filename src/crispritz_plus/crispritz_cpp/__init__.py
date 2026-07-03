"""pybind11 bridge package for the CRISPRitz-plus C++ core.
 
Thin, friendly Python wrappers over the compiled ``_ternary_search_tree``
extension.  The subpackage exposes two kinds of object:
 
Free functions (``api``)
    ``build_tree_cpp``, ``make_search_configuration``,
    ``run_search_executor_cpp``, ``merge_sorted_shards_cpp``, and
    ``write_merged_profiles_cpp`` — the entry points that drive the C++
    indexing, search, merge, and profile-writing stages.
Value wrappers
    ``SearchConfiguration``, ``PartitionResult``, ``GuideProfile``, and the
    ``BulgeMode`` / ``OutputMode`` / ``OutputFormat`` / ``SortMode`` enum
    facades — each holds a native C++ object and delegates to it, exposing the
    real type via a ``native`` property for the pybind11 boundary.
 
Design rule: enum tokens and all validation live in C++; these wrappers only
translate friendly Python arguments (lowercase tokens, wrapper objects) into
the native types the extension expects, keeping a single source of truth.
"""

from .api import (
    build_tree_cpp,
    make_search_configuration,
    merge_sorted_shards_cpp,
    run_search_executor_cpp,
    write_merged_profiles_cpp,
)
from .partition_results import PartitionResult
from .search_configuration import SearchConfiguration
