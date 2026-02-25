"""STR reference management for BED file loading and region lookups."""

import os
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import pysam

from ..utils.vcf_utils import chrom_to_order


def is_valid_tabix(gz_path: str) -> bool:
    """Check that a BGZF file has a valid tabix index.

    Returns True only if:
    - .tbi exists
    - pysam can open the file
    - index is readable
    """
    tbi_path = gz_path + ".tbi"
    if not os.path.exists(gz_path) or not os.path.exists(tbi_path):
        return False

    try:
        tbx = pysam.TabixFile(gz_path)
        # Accessing contigs forces index parsing
        _ = tbx.contigs
        tbx.close()
        return True
    except Exception:
        return False


def sort_bed_file(
    bed_path: str,
    output_path: str,
    chrom_col: int = 0,
    start_col: int = 1,
) -> str:
    """Sort a BED-like file by chromosome and start coordinate.

    Parameters
    ----------
    bed_path : str
        Path to input BED file (tab-delimited).
    output_path : str
        Path to write the sorted BED file.
    chrom_col : int, optional
        Zero-based column index for chromosome. Default is 0.
    start_col : int, optional
        Zero-based column index for start coordinate. Default is 1.

    Returns
    -------
    str
        Path to the sorted BED file.

    Notes
    -----
    - This function loads the BED into memory via pandas. For extremely large
      BED files, consider an external sort.
    - Sorting is lexicographic by chromosome, then numeric by start.
    """
    bed_path = str(bed_path)
    output_path = str(output_path)

    df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        comment="#",
        dtype={chrom_col: "string"},
    )

    # Ensure start is numeric
    df[start_col] = df[start_col].astype("int64")

    # Compute genomic order
    chrom_series = df[chrom_col].astype("string")
    chrom_order = chrom_series.map(chrom_to_order)

    # If unknown chromosomes appear, push them to the end
    chrom_order = chrom_order.fillna(10**9)

    df = df.assign(_chrom_order=chrom_order)

    df.sort_values(
        by=["_chrom_order", start_col],
        kind="mergesort",  # stable sort
        inplace=True,
    )

    df.drop(columns=["_chrom_order"], inplace=True)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", header=False, index=False)

    return output_path


def load_str_reference(str_path: str) -> pd.DataFrame:
    """Ensure a BED file is BGZF-compressed and tabix-indexed.

    This function:
    - Accepts a BED path (.bed or .bed.gz).
    - If the input is already .gz and has a .tbi index, returns it.
    - Otherwise, creates a sorted BED (if needed), BGZF-compresses it, and
      creates a tabix index (preset="bed").

    Parameters
    ----------
    bed_path : str
        Path to input BED file (.bed or .bed.gz).

    Returns
    -------
    str
        Path to the BGZF-compressed, tabix-indexed BED file (*.gz).

    Notes
    -----
    - Tabix indexing requires the BED to be sorted by chromosome and start.
    - This function uses `pysam.tabix_compress` and `pysam.tabix_index`.
    """
    bed_path = Path(str_path)
    cache_dir_path = bed_path.parent

    # If input is already gz + tbi, trust and return.
    if bed_path.suffix == ".gz" and bed_path.exists() and is_valid_tabix(str(bed_path)):
        return str(bed_path)

    # Determine output names
    # Use the base name without trailing .gz if present.
    base_name = bed_path.name
    if base_name.endswith(".gz"):
        base_name = base_name[:-3]

    sorted_bed = cache_dir_path / f"{base_name}.sorted"
    gz_bed = cache_dir_path / f"{base_name}.sorted.gz"

    # If cached gz + tbi exist, reuse.
    if gz_bed.exists() and is_valid_tabix(str(gz_bed)):
        return str(gz_bed)

    # Ensure we have a sorted BED file to compress/index.
    sort_bed_file(str(bed_path), str(sorted_bed))

    # BGZF compress
    pysam.tabix_compress(str(sorted_bed), str(gz_bed), force=True)

    # Tabix index (BED preset expects chrom/start/end in the first columns)
    pysam.tabix_index(str(gz_bed), preset="bed", force=True)

    return str(gz_bed)


def find_overlapping_str(
    str_panel_gz: str,
    chrom: str,
    pos: int,
    end: int,
) -> Optional[Dict]:
    """Find STR region overlapping with variant coordinates using tabix index.

    Parameters
    ----------
    str_panel_gz : str
        Path to BGZF-compressed, tabix-indexed STR reference file.
    chrom : str
        Chromosome name.
    pos : int
        Variant start position (1-based).
    end : int
        Variant end position (1-based).

    Returns
    -------
    Optional[Dict]
        Dictionary with STR region data if overlap found, None otherwise.
        Keys: CHROM, START, END, PERIOD, RU, COUNT
    """
    # Convert to tabix coordinate system: 0-based half-open
    query_start = max(0, pos - 1)
    query_end = end

    try:
        tbx = pysam.TabixFile(str_panel_gz)

        for row in tbx.fetch(chrom, query_start, query_end):
            parts = row.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue

            try:
                str_chrom = parts[0]
                str_start = int(parts[1])
                str_end = int(parts[2])
                period = int(parts[3])
                ru = parts[4]
                count = int(parts[5]) if len(parts) > 5 else None
            except ValueError:
                continue

            # Check true overlap in 1-based coordinates
            if str_start <= end and str_end >= pos:
                result = {
                    "CHROM": str_chrom,
                    "START": str_start,
                    "END": str_end,
                    "PERIOD": period,
                    "RU": ru,
                }
                if count is not None:
                    result["COUNT"] = count
                return result

        return None

    except ValueError:
        # Chromosome not present in index
        return None


def get_str_at_position(
    str_panel_gz: str,
    chrom: str,
    pos: int,
) -> Optional[Dict]:
    """Get STR region containing a specific position using tabix index.

    Parameters
    ----------
    str_panel_gz : str
        Path to BGZF-compressed, tabix-indexed STR reference file.
    chrom : str
        Chromosome name.
    pos : int
        Position to query (1-based).

    Returns
    -------
    Optional[Dict]
        Dictionary with STR region data if position is within an STR, None otherwise.
    """
    return find_overlapping_str(str_panel_gz, chrom, pos, pos)
