"""VCF file processing and workflow management."""

import logging
from pathlib import Path
from typing import Iterator, List

import pandas as pd
import pysam

from ..parsers.base import BaseVCFParser
from ..parsers.generic import GenericParser
from ..utils.vcf_utils import chrom_to_order
from .annotation import build_new_record, make_modified_header, should_skip_genotype

logger = logging.getLogger(__name__)


def check_vcf_sorted(vcf_in: pysam.VariantFile) -> bool:
    """Validate VCF sorting by chromosome and position.

    Checks if VCF records are sorted by chromosome and position.
    Rewinds the file after checking.

    Parameters
    ----------
    vcf_in : pysam.VariantFile
        Input VCF file

    Returns
    -------
    bool
        True if VCF is sorted, False otherwise
    """
    last_chrom = None
    last_pos = -1

    for rec in vcf_in:
        chrom = chrom_to_order(rec.contig)
        pos = rec.pos

        if last_chrom is not None:  # noqa: SIM102
            if chrom < last_chrom or (chrom == last_chrom and pos < last_pos):
                logger.warning(f"Not sorted, because {chrom} < {last_chrom} or {pos} < {last_pos}")
                vcf_in.reset()
                return False

        last_chrom, last_pos = chrom, pos

    vcf_in.reset()
    return True


def reset_and_sort_vcf(vcf_in: pysam.VariantFile) -> List[pysam.VariantRecord]:
    """Sort VCF records in memory when needed.

    Loads all VCF records into memory and sorts them by chromosome
    and position according to the contig order in the header.

    Parameters
    ----------
    vcf_in : pysam.VariantFile
        Input VCF file

    Returns
    -------
    List[pysam.VariantRecord]
        Sorted list of VCF records

    Notes
    -----
    This loads the entire VCF into memory, so use with caution for large files.
    """
    header = vcf_in.header
    records = list(vcf_in)

    # Create contig order mapping
    contig_order = {c: i for i, c in enumerate(header.contigs.keys())}

    # Sort by contig order and position
    records.sort(key=lambda r: (contig_order.get(r.contig, float("inf")), r.pos))

    return records


def generate_annotated_records(
    vcf_in: pysam.VariantFile,
    str_df: pd.DataFrame,
    parser: BaseVCFParser = None,
    somatic_mode: bool = False,
) -> Iterator[pysam.VariantRecord]:
    """Generator yielding annotated VCF records.

    Processes VCF records and yields annotated records for variants that
    overlap with STR regions. Handles sorting if needed and optionally filters
    records based on genotype criteria.

    Parameters
    ----------
    vcf_in : pysam.VariantFile
        Input VCF file
    str_df : pd.DataFrame
        DataFrame with STR regions (from load_str_reference)
    parser : BaseVCFParser, optional
        Parser for genotype extraction. Uses GenericParser if None.
    somatic_mode : bool, optional
        Enable somatic filtering. When True, skips variants where both samples
        have identical genotypes. Default is False.

    Yields
    ------
    pysam.VariantRecord
        Annotated VCF records

    Notes
    -----
    - Automatically sorts VCF if not sorted
    - Skips records without STR overlap
    - If somatic_mode=True, filters records with identical genotypes
    """
    if parser is None:
        parser = GenericParser()

    header = make_modified_header(vcf_in)

    # Check if VCF is sorted
    if not check_vcf_sorted(vcf_in):
        logger.warning("Input VCF is not sorted - sorting in memory.")
        records = reset_and_sort_vcf(vcf_in)
    else:
        vcf_in.reset()
        records = vcf_in.fetch()

    # Prepare STR list for efficient lookup
    str_idx = 0
    str_list = str_df.to_dict("records")

    skipped_count = 0
    for record in records:
        # Advance STR index to current chromosome/position
        while str_idx < len(str_list) and (
            str_list[str_idx]["CHROM"] != record.chrom
            or (
                str_list[str_idx]["CHROM"] == record.chrom and str_list[str_idx]["END"] < record.pos
            )
        ):
            str_idx += 1

        if str_idx >= len(str_list):
            break

        str_row = str_list[str_idx]

        # Check for overlap - variant position should be within STR region
        if (
            str_row["CHROM"] != record.chrom
            or record.pos < str_row["START"]
            or record.pos > str_row["END"]
        ):
            continue  # No overlap

        # Skip based on genotype filtering (only if somatic_mode enabled)
        if somatic_mode and should_skip_genotype(record, parser):
            skipped_count += 1
            logger.debug(
                f"Skipped {record.contig}:{record.pos} - identical genotypes (somatic mode)"
            )
            continue

        # Build and yield annotated record
        new_record = build_new_record(record, str_row, header, parser)

        # Skip if alt == ref
        # It happens if variant was not normalized and actually happens outside of STR region
        if new_record.alleles[0] == new_record.alleles[1]:
            continue

        yield new_record

    # Log summary if records were skipped
    if skipped_count > 0:
        logger.warning(
            f"Skipped {skipped_count} records due to identical genotypes (somatic filtering enabled)"
        )


def annotate_vcf_to_file(
    vcf_path: str,
    str_df: pd.DataFrame,
    output_path: str,
    parser: BaseVCFParser = None,
    somatic_mode: bool = False,
) -> None:
    """Process VCF file and write annotated output.

    Reads a VCF file, annotates variants that overlap with STR regions,
    and writes the annotated records to an output file.

    Parameters
    ----------
    vcf_path : str
        Path to input VCF file
    str_df : pd.DataFrame
        DataFrame with STR regions (from load_str_reference)
    output_path : str
        Path to output VCF file
    parser : BaseVCFParser, optional
        Parser for genotype extraction. Uses GenericParser if None.
    somatic_mode : bool, optional
        Enable somatic filtering. Default is False.

    Notes
    -----
    Prints summary statistics after processing.
    """
    if parser is None:
        parser = GenericParser()

    vcf_in = pysam.VariantFile(vcf_path)
    new_header = make_modified_header(vcf_in)
    vcf_out = pysam.VariantFile(output_path, "w", header=new_header)

    # Process and write records
    written_count = 0
    for record in generate_annotated_records(vcf_in, str_df, parser, somatic_mode=somatic_mode):
        vcf_out.write(record)
        written_count += 1

    vcf_out.close()
    vcf_in.close()

    logger.info(f"Wrote {written_count} annotated records to {output_path}")


def process_directory(
    input_dir: str,
    str_bed_path: str,
    output_dir: str,
    parser: BaseVCFParser = None,
    somatic_mode: bool = False,
) -> None:
    """Batch process directory of VCF files.

    Processes all VCF files in a directory and writes annotated versions
    to the output directory.

    Parameters
    ----------
    input_dir : str
        Directory containing input VCF files
    str_bed_path : str
        Path to BED file with STR regions
    output_dir : str
        Directory for output VCF files
    parser : BaseVCFParser, optional
        Parser for genotype extraction. Uses GenericParser if None.
    somatic_mode : bool, optional
        Enable somatic filtering. Default is False.
    """
    from .str_reference import load_str_reference

    if parser is None:
        parser = GenericParser()

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Load STR reference
    str_df = load_str_reference(str_bed_path)

    # Process each VCF file
    input_path = Path(input_dir)
    for vcf_file in input_path.glob("*.vcf*"):
        if vcf_file.suffix in [".vcf", ".gz"]:
            # Generate output filename
            base_name = vcf_file.stem
            if base_name.endswith(".vcf"):
                base_name = base_name[:-4]
            output_file = Path(output_dir) / f"{base_name}.annotated.vcf"

            # Skip if already processed
            if output_file.exists():
                logger.info(f"Skipping {vcf_file.name} — already processed.")
                continue

            logger.info(f"Processing {vcf_file.name}...")
            annotate_vcf_to_file(
                str(vcf_file), str_df, str(output_file), parser, somatic_mode=somatic_mode
            )
            logger.info(f" → Output: {output_file}")
