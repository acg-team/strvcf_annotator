"""STR reference management for BED file loading and region lookups."""

from typing import Dict, Optional

import pandas as pd

from ..utils.vcf_utils import chrom_to_order


def load_str_reference(str_path: str) -> pd.DataFrame:
    """Load STR reference data from BED file.

    Loads a BED file containing STR (Short Tandem Repeat) regions and converts
    coordinates from 0-based BED format to 1-based VCF format. Calculates the
    number of repeat units for each region.

    Parameters
    ----------
    str_path : str
        Path to BED file with STR regions

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: CHROM, START, END, PERIOD, RU, COUNT
        - CHROM: Chromosome name
        - START: 1-based start position (converted from BED 0-based)
        - END: 1-based end position
        - PERIOD: Length of repeat unit
        - RU: Repeat unit sequence
        - COUNT: Number of repeat units in the region

    Notes
    -----
    BED files use 0-based coordinates, but VCF files use 1-based coordinates.
    This function converts START positions by adding 1. END positions are kept
    as-is since BED END is exclusive and VCF END is inclusive.
    """
    df = pd.read_csv(str_path, sep="\t", header=None)
    df.columns = ["CHROM", "START", "END", "PERIOD", "RU"]

    # Convert from 0-based BED to 1-based VCF coordinates
    df["START"] = df["START"] + 1

    # Calculate number of repeat units
    df["COUNT"] = (df["END"] - df["START"] + 1) / df["PERIOD"]

    # Add chromosome order column for proper sorting
    df["CHROM_ORDER"] = df["CHROM"].apply(chrom_to_order)

    # Sort by chromosome (natural order) and position for efficient lookups
    df.sort_values(by=["CHROM_ORDER", "START"], inplace=True)
    df.drop(columns="CHROM_ORDER", inplace=True)
    return df


def find_overlapping_str(str_df: pd.DataFrame, chrom: str, pos: int, end: int) -> Optional[Dict]:
    """Find STR region overlapping with variant coordinates.

    Searches for an STR region that overlaps with the given variant position.
    Uses efficient binary search on sorted DataFrame.

    Parameters
    ----------
    str_df : pd.DataFrame
        DataFrame with STR regions (from load_str_reference)
    chrom : str
        Chromosome name
    pos : int
        Variant start position (1-based)
    end : int
        Variant end position (1-based)

    Returns
    -------
    Optional[Dict]
        Dictionary with STR region data if overlap found, None otherwise
        Contains keys: CHROM, START, END, PERIOD, RU, COUNT
    """
    # Filter by chromosome
    chrom_df = str_df[str_df["CHROM"] == chrom]

    if chrom_df.empty:
        return None

    # Find overlapping regions
    # Overlap occurs when: variant_end >= str_start AND variant_start <= str_end
    overlapping = chrom_df[(chrom_df["START"] <= end) & (chrom_df["END"] >= pos)]

    if overlapping.empty:
        return None

    # Return first overlapping region
    return overlapping.iloc[0].to_dict()


def get_str_at_position(str_df: pd.DataFrame, chrom: str, pos: int) -> Optional[Dict]:
    """Get STR region containing a specific position.

    Parameters
    ----------
    str_df : pd.DataFrame
        DataFrame with STR regions (from load_str_reference)
    chrom : str
        Chromosome name
    pos : int
        Position to query (1-based)

    Returns
    -------
    Optional[Dict]
        Dictionary with STR region data if position is within an STR, None otherwise
    """
    return find_overlapping_str(str_df, chrom, pos, pos)
