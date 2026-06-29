"""
DNA / IUPAC alphabet constants and sequence utilities for CRISPRitz-plus.
 
Centralises the nucleotide alphabets, the reverse-complement table, and the
IUPAC degeneracy mappings used throughout the package, along with two small
sequence-transformation helpers.
 
Module-level constants
----------------------
DNA : List[str]
    The five canonical DNA symbols (four bases plus the ambiguity wildcard
    ``'N'``).
IUPAC : List[str]
    The full IUPAC nucleotide alphabet: :data:`DNA` plus the ten degenerate
    ambiguity codes.
RC : Dict[str, str]
    Reverse-complement lookup for every IUPAC symbol, in both upper- and
    lower-case.  ``'U'`` / ``'u'`` map to ``'A'`` / ``'a'`` so RNA input is
    tolerated.
IUPACTABLE : Dict[str, str]
    Maps each IUPAC code to the sorted set of concrete bases it represents
    (e.g. ``'R' -> 'AG'``).
IUPAC_ENCODER : Dict[str, str]
    Inverse, order-insensitive mapping from a concrete base combination to
    its IUPAC code.  Every permutation of each base set maps to the same
    code, so ``'AG'`` and ``'GA'`` both resolve to ``'R'``.
"""

from itertools import permutations


# ==============================================================================
# Define DNA-related constant variables
# ==============================================================================
 
#: Canonical DNA alphabet: the four bases plus the ``'N'`` ambiguity wildcard.
DNA = ["A", "C", "G", "T", "N"]

#: Full IUPAC nucleotide alphabet: :data:`DNA` extended with the ten
#: degenerate ambiguity codes (R, Y, S, W, K, M, B, D, H, V).
IUPAC = DNA + ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"]

#: Reverse-complement lookup table covering every IUPAC symbol in both
#: upper- and lower-case.  ``'U'`` / ``'u'`` (RNA uracil) are folded onto
#: ``'A'`` / ``'a'`` so that RNA input is complemented sensibly.
RC = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
    "R": "Y",
    "Y": "R",
    "M": "K",
    "K": "M",
    "H": "D",
    "D": "H",
    "B": "V",
    "V": "B",
    "N": "N",
    "S": "S",
    "W": "W",
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "u": "a",
    "r": "y",
    "y": "r",
    "m": "k",
    "k": "m",
    "h": "d",
    "d": "h",
    "b": "v",
    "v": "b",
    "n": "n",
    "s": "s",
    "w": "w",
}

#: Maps each IUPAC code to the concrete base set it represents
#: (e.g. ``'R' -> 'AG'``, ``'N' -> 'ACGT'``).  Used to expand ambiguity
#: codes into their constituent bases.
IUPACTABLE = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "M": "AC",
    "K": "GT",
    "S": "CG",
    "W": "AT",
    "H": "ACT",
    "B": "CGT",
    "V": "ACG",
    "D": "AGT",
    "N": "ACGT",
}

#: Inverse, order-insensitive encoder built from :data:`IUPACTABLE`.  Every
#: permutation of a code's base set is a key mapping back to that code, so a
#: base combination can be collapsed to its IUPAC symbol regardless of the
#: order in which the bases are supplied (e.g. both ``'AG'`` and ``'GA'``
#: map to ``'R'``).
IUPAC_ENCODER = {
    perm: k
    for k, v in IUPACTABLE.items()
    for perm in {"".join(p) for p in permutations(v)}
}

# ==============================================================================
# Public API
# ==============================================================================

def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of an IUPAC nucleotide sequence.
 
    Reverses *sequence* and replaces each symbol with its complement from
    :data:`RC`.  Case is preserved because :data:`RC` contains both upper-
    and lower-case entries.
 
    Parameters
    ----------
    sequence : str
        Nucleotide sequence containing only symbols present in :data:`RC`.
 
    Returns
    -------
    str
        The reverse-complemented sequence.
 
    Raises
    ------
    KeyError
        If *sequence* contains a character absent from :data:`RC`.
 
    Examples
    --------
    >>> reverse_complement("ACGT")
    'ACGT'
    >>> reverse_complement("AACG")
    'CGTT'
    """
    return "".join(RC[nt] for nt in sequence[::-1])


def dna2rna(sequence: str) -> str:
    """Convert a DNA sequence to its RNA equivalent.

    Transcribes *sequence* by substituting thymine for uracil: ``'T'`` is
    replaced with ``'U'`` and ``'t'`` with ``'u'``.  Case is preserved and
    all other characters (including IUPAC ambiguity codes and gap symbols)
    are passed through unchanged.

    Parameters
    ----------
    sequence : str
        The DNA sequence to convert.

    Returns
    -------
    str
        The RNA sequence with every thymine replaced by uracil.

    Examples
    --------
    >>> dna2rna("ACGT")
    'ACGU'
    >>> dna2rna("acgt")
    'acgu'
    """
    return sequence.replace("T", "U").replace("t", "u")
