"""Utilities for repeat sequence operations."""

from typing import Dict, Tuple, Union

import pandas as pd


def extract_repeat_sequence(str_row: Union[Dict, pd.Series]) -> str:
    """Reconstruct repeat sequence from STR metadata.

    Generates the full repeat sequence by repeating the repeat unit (RU)
    the calculated number of times (COUNT).

    Parameters
    ----------
    str_row : Dict or pd.Series
        STR region data containing 'RU' (repeat unit) and 'COUNT' (number of repeats)

    Returns
    -------
    str
        Full repeat sequence

    Examples
    --------
    >>> str_row = {'RU': 'CAG', 'COUNT': 5}
    >>> extract_repeat_sequence(str_row)
    'CAGCAGCAGCAGCAG'
    """
    return str_row["RU"] * int(str_row["COUNT"])


def count_repeat_units(sequence: str, motif: str) -> int:
    """Return the longest contiguous run of `motif` in `sequence`.

    The function looks for exact, non-overlapping copies of `motif` that occur
    consecutively and returns the maximum number of such copies in any run.

    This corresponds to how STR repeat counts are typically defined: the
    length of the longest perfect contiguous block of the repeat unit.

    Parameters
    ----------
    sequence : str
        DNA sequence to search.
    motif : str
        Repeat unit motif to count (e.g. 'A', 'CAG').

    Returns
    -------
    int
        Length of the longest contiguous run of `motif` in `sequence`.

    Raises
    ------
    ValueError
        If `motif` is empty or if either argument is not a string.

    Examples
    --------
    Perfect repeats
    ~~~~~~~~~~~~~~~
    >>> count_repeat_units("CAGCAGCAG", "CAG")
    3

    Imperfect tail
    ~~~~~~~~~~~~~~
    >>> count_repeat_units("CAGCAGCA", "CAG")
    2

    No repeats
    ~~~~~~~~~~
    >>> count_repeat_units("ATCG", "CAG")
    0

    Homopolymer runs
    ~~~~~~~~~~~~~~~~
    >>> count_repeat_units("ATAAAAA", "A")
    5
    >>> count_repeat_units("AAAATAAA", "A")
    4  # longest contiguous run is 'AAAA'

    Overlapping motifs
    ~~~~~~~~~~~~~~~~~~
    'AAAA' with motif 'AA' contains two non-overlapping copies: 'AA' 'AA'
    >>> count_repeat_units("AAAA", "AA")
    2
    """
    # Basic type and value validation
    if not isinstance(sequence, str):
        raise ValueError(f"'sequence' must be a string, got {type(sequence).__name__}")
    if not isinstance(motif, str):
        raise ValueError(f"'motif' must be a string, got {type(motif).__name__}")
    if motif == "":
        raise ValueError("'motif' must be a non-empty string")

    mlen = len(motif)
    slen = len(sequence)

    if mlen > slen or slen == 0:
        return 0

    max_run = 0
    i = 0

    while i <= slen - mlen:
        # Check if a run starts at position i
        if sequence[i : i + mlen] == motif:
            run_length = 0
            j = i
            # Count how many contiguous motif copies we have from position i
            while j <= slen - mlen and sequence[j : j + mlen] == motif:
                run_length += 1
                j += mlen

            if run_length > max_run:
                max_run = run_length

            # Skip past this entire run
            i = j
        else:
            i += 1

    return max_run


def normalize_variant(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """
    Locally normalize (pos, ref, alt) by trimming shared prefix/suffix.

    - pos is 1-based VCF coordinate.
    - Trimming is case-insensitive.
    - We always keep at least 1 base in ref and alt if they differ.

    Returns
    -------
    new_pos, new_ref, new_alt
    """
    # early exit: identical → no-op
    if ref.upper() == alt.upper():
        return pos, ref, alt

    r = ref
    a = alt

    # trim common prefix
    while len(r) > 1 and len(a) > 1 and r[0].upper() == a[0].upper():
        r = r[1:]
        a = a[1:]
        pos += 1

    # trim common suffix
    while len(r) > 1 and len(a) > 1 and r[-1].upper() == a[-1].upper():
        r = r[:-1]
        a = a[:-1]

    return pos, r, a


def apply_variant_to_repeat(
    pos: int, ref: str, alt: str, repeat_start: int, repeat_seq: str
) -> str:
    """Apply a variant to the STR repeat sequence, with normalization.

    1) Normalize (pos, ref, alt) by trimming shared prefix/suffix.
    2) If the normalized variant lies fully inside the STR, apply the full ALT.
    3) If it only partially overlaps, apply only the overlapping part:
       - SNP-like (len(ref) == len(alt)): positional alignment.
       - Indel-like (len(ref) != len(alt)): align overlapping part from the
         end of ALT.

    Conceptually, we assume the genomic reference at this locus looks like
    `repeat_seq + UNKNOWN_SUFFIX`. Any differences outside the STR window
    are ignored when computing the mutated STR.

    Case handling
    -------------
    - Normalization and overlap logic are case-insensitive.
    - The output casing follows the STR sequence at the overlapping segment:
      * If the overlapping STR slice is all lowercase, ALT is lowercased.
      * If it is all uppercase (typical), ALT is uppercased.
      * Otherwise, ALT is used as-is.

    Parameters
    ----------
    pos : int
        Variant position (1-based VCF coordinate).
    ref : str
        Reference allele from VCF.
    alt : str
        Alternate allele from VCF.
    repeat_start : int
        Start position of repeat region (1-based).
    repeat_seq : str
        Reference STR sequence (panel).

    Returns
    -------
    str
        Mutated repeat sequence after applying the normalized variant
        restricted to the STR window. If there is no overlap, returns
        repeat_seq unchanged.
    """
    repeat_len = len(repeat_seq)

    # --- 1) Normalize variant locally ---
    pos, ref, alt = normalize_variant(pos, ref, alt)

    # After normalization, if ref == alt, nothing to do for the STR
    if ref.upper() == alt.upper():
        return repeat_seq

    repeat_end = repeat_start + repeat_len - 1
    var_start = pos
    var_end = pos + len(ref) - 1  # inclusive

    # --- 2) No overlap with STR ---
    if var_end < repeat_start or var_start > repeat_end:
        return repeat_seq

    # --- 3) Variant fully inside STR: apply full ALT (insertions/deletions kept) ---
    if var_start >= repeat_start and var_end <= repeat_end:
        relative_pos = var_start - repeat_start  # 0-based offset within repeat_seq
        before = repeat_seq[:relative_pos]
        after = repeat_seq[relative_pos + len(ref) :]

        panel_slice = repeat_seq[relative_pos : relative_pos + len(ref)]

        # Case-normalize ALT to match panel_slice
        alt_adj = alt
        if panel_slice.islower():
            alt_adj = alt.lower()
        elif panel_slice.isupper():
            alt_adj = alt.upper()

        return before + alt_adj + after

    # --- 4) Partial overlap with STR (starts before and/or ends after) ---
    # Compute overlap region in reference coordinates.
    overlap_start = max(repeat_start, var_start)
    overlap_end = min(repeat_end, var_end)
    overlap_len = overlap_end - overlap_start + 1
    # Position inside the STR where overlap begins (0-based)
    relative_pos = overlap_start - repeat_start

    # Extract ALT substring corresponding to overlapping bases.
    if len(ref) == len(alt):
        # SNP-like: align by position (same offset as ref)
        i0 = overlap_start - var_start
        alt_overlap_raw = alt[i0 : i0 + overlap_len]
    else:
        # Indel-like: align from the end of ALT so that the last base of ALT
        # corresponds to the last overlapping base of REF.

        # ALT shorter or equal: use entire ALT
        # Use only the suffix of ALT that covers the overlapping segment
        alt_overlap_raw = alt if overlap_len >= len(alt) else alt[-overlap_len:]

    panel_slice = repeat_seq[relative_pos : relative_pos + overlap_len]

    # Case-normalize ALT overlap to match panel_slice
    if panel_slice.islower():
        alt_overlap = alt_overlap_raw.lower()
    elif panel_slice.isupper():
        alt_overlap = alt_overlap_raw.upper()
    else:
        alt_overlap = alt_overlap_raw

    before_mut = repeat_seq[:relative_pos]
    after_mut = repeat_seq[relative_pos + overlap_len :]

    mutated = before_mut + alt_overlap + after_mut
    return mutated


def is_perfect_repeat(sequence: str, motif: str) -> bool:
    """Check if sequence is a perfect repeat of the motif.

    A perfect repeat means the sequence consists entirely of exact copies
    of the motif with no interruptions or variations.

    Parameters
    ----------
    sequence : str
        DNA sequence to check
    motif : str
        Repeat unit motif

    Returns
    -------
    bool
        True if sequence is a perfect repeat, False otherwise

    Examples
    --------
    >>> is_perfect_repeat('CAGCAGCAG', 'CAG')
    True
    >>> is_perfect_repeat('CAGCAGCA', 'CAG')
    False
    """
    if not sequence or not motif:
        return False

    count = count_repeat_units(sequence, motif)
    return sequence == motif * count
