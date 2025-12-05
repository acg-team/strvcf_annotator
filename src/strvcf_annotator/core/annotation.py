"""Core annotation engine for building STR-annotated VCF records."""

import contextlib
import logging
from typing import Dict, Union

import pandas as pd
import pysam

from ..parsers.base import BaseVCFParser
from .repeat_utils import (
    apply_variant_to_repeat,
    count_repeat_units,
    extract_repeat_sequence,
    is_perfect_repeat,
)

logger = logging.getLogger(__name__)


def make_modified_header(vcf_in: pysam.VariantFile) -> pysam.VariantHeader:
    """Create VCF header with STR-specific INFO and FORMAT fields.

    Creates a modified VCF header that includes all original header information
    plus STR-specific annotations. Replaces existing RU, PERIOD, REF, PERFECT
    INFO fields and REPCN FORMAT field with STR-specific definitions.

    Parameters
    ----------
    vcf_in : pysam.VariantFile
        Input VCF file object

    Returns
    -------
    pysam.VariantHeader
        New header with STR-specific fields

    Notes
    -----
    INFO fields added/replaced:
        - RU: Repeat unit
        - PERIOD: Repeat period (length of unit)
        - REF: Reference copy number
        - PERFECT: Indicates perfect repeats in REF and ALT

    FORMAT field added/replaced:
        - REPCN: Genotype as number of repeat motif copies
    """
    info_to_replace = ["RU", "PERIOD", "REF", "PERFECT"]
    format_to_replace = ["REPCN"]
    new_header = pysam.VariantHeader()

    # Copy raw header lines, except the ones we want to modify
    for rec in vcf_in.header.records:
        if rec.key == "INFO" and rec["ID"] in info_to_replace:
            continue
        if rec.key == "FORMAT" and rec["ID"] in format_to_replace:
            continue

        # Fallback: skip unrecognized or malformed lines
        with contextlib.suppress(Exception):
            new_header.add_record(rec)

    # Copy contigs explicitly
    for contig in vcf_in.header.contigs.values():
        if contig.name not in new_header.contigs:
            new_header.contigs.add(contig.name, length=contig.length)

    # Add samples
    for sample in vcf_in.header.samples:
        new_header.add_sample(sample)

    # Add modified INFO and FORMAT fields
    new_header.info.add("RU", 1, "String", "Repeat unit")
    new_header.info.add("PERIOD", 1, "Integer", "Repeat period (length of unit)")
    new_header.info.add("REF", 1, "Integer", "Reference copy number")
    new_header.info.add(
        "PERFECT",
        1,
        "String",
        "Indicates if the repeat sequence is a perfect RU repeat (both REF and ALT)",
    )
    new_header.formats.add(
        "REPCN", 2, "Integer", "Genotype given in number of copies of the repeat motif"
    )

    return new_header


def build_new_record(
    record: pysam.VariantRecord,
    str_row: Union[Dict, pd.Series],
    header: pysam.VariantHeader,
    parser: BaseVCFParser,
) -> pysam.VariantRecord:
    """Build annotated VCF record with STR alleles and metadata.

    Constructs a new VCF record where alleles represent full repeat sequences
    (before and after mutation) and adds STR-specific annotations to INFO and
    FORMAT fields.

    Parameters
    ----------
    record : pysam.VariantRecord
        Original VCF record with mutation
    str_row : Dict or pd.Series
        STR metadata (CHROM, START, END, RU, PERIOD)
    header : pysam.VariantHeader
        Modified header with STR fields
    parser : BaseVCFParser
        Parser for extracting genotype information

    Returns
    -------
    pysam.VariantRecord
        New record with STR alleles and annotations

    Notes
    -----
    - Logs warning if reference mismatch detected
    - Calculates repeat copy numbers for REF and ALT
    - Marks PERFECT=TRUE only if both alleles are perfect repeats
    - Preserves all original FORMAT fields
    """
    repeat_start = str_row["START"]
    repeat_seq = extract_repeat_sequence(str_row)
    ref_base = record.ref
    alt_base = record.alts[0] if record.alts else ref_base
    repeat_len = len(repeat_seq)
    repeat_end = repeat_start + repeat_len - 1

    # Verify reference matches
    var_start = record.pos
    var_end = record.pos + len(ref_base) - 1  # inclusive

    # --- 1) Check overlap between REF and panel sequence ---
    # Compute genomic overlap between [var_start, var_end] and [repeat_start, repeat_end]
    overlap_start = max(repeat_start, var_start)
    overlap_end = min(repeat_end, var_end)

    if overlap_start <= overlap_end:
        overlap_len = overlap_end - overlap_start + 1

        # Indices into repeat_seq (0-based)
        rep_start_idx = overlap_start - repeat_start
        rep_end_idx_excl = rep_start_idx + overlap_len

        # Indices into REF (0-based)
        ref_start_idx = overlap_start - var_start
        ref_end_idx_excl = ref_start_idx + overlap_len

        panel_sub = repeat_seq[rep_start_idx:rep_end_idx_excl]
        ref_sub = ref_base[ref_start_idx:ref_end_idx_excl]

        if panel_sub.upper() != ref_sub.upper():
            logger.warning(
                "Reference mismatch in STR overlap: VCF %s:%d %s>%s, "
                "STR panel %s:%d-%d RU=%s\n"
                "Panel overlap: %s\n"
                "VCF REF overlap: %s",
                record.contig,
                record.pos,
                record.alleles[0],
                record.alleles[1],
                str_row["CHROM"],
                str_row["START"],
                str_row["END"],
                str_row["RU"],
                panel_sub,
                ref_sub,
            )

    # Apply mutation to get alternate sequence
    mutated_seq = apply_variant_to_repeat(record.pos, ref_base, alt_base, repeat_start, repeat_seq)

    # Count repeat units
    ru = str_row["RU"]
    ref_len = count_repeat_units(repeat_seq, ru)
    alt_len = count_repeat_units(mutated_seq, ru)
    perfect = is_perfect_repeat(repeat_seq, ru) and is_perfect_repeat(mutated_seq, ru)

    # Prepare INFO fields - only copy safe, commonly used fields
    # to avoid type/Number mismatches with fields from original VCF
    info = {}

    # List of safe INFO fields to copy if they exist in original record
    safe_info_fields = ["NS", "DP", "AF", "AC", "AN", "MQ", "SB"]

    for field in safe_info_fields:
        if field in record.info:
            # Skip fields that can't be copied
            with contextlib.suppress(TypeError, ValueError, KeyError):
                info[field] = record.info[field]

    # Add our STR-specific INFO fields
    info.update(
        {
            "RU": ru,
            "PERIOD": int(str_row["PERIOD"]),
            "REF": int(ref_len),
            "PERFECT": "TRUE" if perfect else "FALSE",
        }
    )

    # Create new record
    new_record = header.new_record(
        contig=record.contig,
        start=repeat_start - 1,
        stop=str_row["END"],
        id=".",
        alleles=(repeat_seq, mutated_seq),
        info=info,
        filter=record.filter.keys(),
    )

    # Copy FORMAT fields and add REPCN
    for sample_name in record.samples:
        # Copy all FORMAT fields except GT and REPCN (we'll set those)
        old_sample = record.samples[sample_name]
        new_sample = new_record.samples[sample_name]

        for field_name in old_sample:
            if field_name not in ["GT", "REPCN"]:
                with contextlib.suppress(TypeError, ValueError):
                    new_sample[field_name] = old_sample[field_name]
                    # Skip fields that can't be copied (e.g., incompatible types)

        # Get genotype using parser and set GT + REPCN
        sample_idx = list(record.samples.keys()).index(sample_name)
        gt = parser.get_genotype(record, sample_idx)

        if gt is not None:
            alleles = [ref_len if allele == 0 else alt_len if allele == 1 else 0 for allele in gt]
            new_sample["GT"] = gt
            new_sample["REPCN"] = alleles
        else:
            new_sample["GT"] = (None, None)
            new_sample["REPCN"] = (0, 0)

    return new_record


def should_skip_genotype(record: pysam.VariantRecord, parser: BaseVCFParser) -> bool:
    """Determine if record should be skipped based on genotype filtering.

    Skips records where:
    - Not exactly 2 samples present
    - Genotypes are invalid or missing
    - Both samples have identical genotypes

    Parameters
    ----------
    record : pysam.VariantRecord
        VCF record to check
    parser : BaseVCFParser
        Parser for extracting genotypes

    Returns
    -------
    bool
        True if record should be skipped, False otherwise
    """
    samples = list(record.samples.keys())

    # Only process records with exactly 2 samples
    if len(samples) != 2:
        return False

    # Get genotypes for both samples
    gt_0 = parser.get_genotype(record, 0)
    gt_1 = parser.get_genotype(record, 1)

    # Skip if either genotype is invalid
    if gt_0 is None or gt_1 is None:
        return True

    # Skip if genotypes have more than 2 alleles
    if len(gt_0) > 2 or len(gt_1) > 2:
        return True

    # Skip if both samples have identical genotypes
    return gt_0 == gt_1
