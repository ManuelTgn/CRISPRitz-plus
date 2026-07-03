#ifndef TST_VERBOSITY_HPP
#define TST_VERBOSITY_HPP

#include <iostream>
#include <string_view>

namespace crispritz {

// ---------------------------------------------------------------------------
// Verbosity levels
//
// These mirror the Python-side VERBOSITY_LVL constants so the C++ core and the
// Python orchestration layer speak the same language across the pybind11
// boundary:
//   0 = Silent  : no output except fatal errors
//   1 = Normal  : major execution stages (building an index, writing output)
//   2 = Verbose : intermediate summaries (site counts, filenames)
//   3 = Debug   : developer diagnostics (phase entry, computed values, timing)
//
// The verbosity level is passed into build_tree() from Python and stored on the
// TernarySearchTree so every diagnostic message can be gated consistently.
// ---------------------------------------------------------------------------
enum : int {
  VERBOSITY_SILENT = 0,
  VERBOSITY_NORMAL = 1,
  VERBOSITY_VERBOSE = 2,
  VERBOSITY_DEBUG = 3,
};

/**
 * @brief Write @p message to stdout when @p verbosity meets @p threshold.
 *
 * The C++ analogue of the Python ``print_verbosity`` helper. Diagnostic output
 * from the index builder must route through this gate rather than a bare
 * ``std::cout`` so that ``--verbosity 0`` (Silent) truly suppresses all
 * informational output.
 *
 * @param message   The message to print (a trailing newline is added).
 * @param verbosity The verbosity level requested by the caller.
 * @param threshold The minimum verbosity level required to print the message.
 * @complexity O(message length).
 */
inline void print_verbosity(std::string_view message, int verbosity,
                            int threshold) {
  if (verbosity >= threshold)
    std::cout << message << '\n';
}

} // namespace crispritz

#endif // TST_VERBOSITY_HPP