"""Core annotation engine for building STR-annotated VCF records."""

import contextlib
import logging
from typing import Dict, Union

import pandas as pd
import pysam
from trtools.utils.utils import GetCanonicalMotif

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
    ignore_mismatch_warnings: bool = False,
    mismatch_truth: str = "panel",  # "panel" | "vcf" | "skip"
) -> Union[pysam.VariantRecord, None]:
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
    ignore_mismatch_warnings : bool
        If True, do not log mismatch warnings between VCF REF and panel repeat sequence.
        Annotation continues regardless.
    mismatch_truth : str
        Which source to treat as correct when there is a mismatch:
          - "panel": trust panel repeat sequence (default behavior)
          - "vcf": trust VCF REF. patch the panel repeat sequence overlap to match VCF REF
          - "skip": skip record with mismatch

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
    - Returns None if there is a reference mismatch and mismatch_truth is "skip"
    """
    mismatch_truth = (mismatch_truth or "panel").lower()
    if mismatch_truth not in {"panel", "vcf", "skip"}:
        raise ValueError("mismatch_truth must be one of: 'panel', 'vcf', 'skip'")

    repeat_start = int(str_row["START"])
    repeat_end = int(str_row["END"])
    repeat_seq = extract_repeat_sequence(str_row)
    repeat_seq = (repeat_seq or "").upper()

    # Variant alleles
    ref_base = record.ref.upper()
    pos = int(record.pos)

    period = int(str_row["PERIOD"])
    ru = (str_row.get("RU") or "").upper()

    # Overlap check + (optional) patching if VCF is truth
    # We compare overlap between VCF REF and panel repeat sequence.
    # If mismatch_truth == "vcf" and there is overlap, we patch repeat_seq's overlap
    # to match VCF REF overlap before applying the variant.
    repeat_len = len(repeat_seq)
    panel_sub = ""
    ref_sub = ""
    suppress = False
    raw_mismatch = False

    if repeat_len > 0 and len(ref_base) > 0:
        repeat_end_seq = repeat_start + repeat_len - 1  # inclusive
        var_start = pos
        var_end = pos + len(ref_base) - 1  # inclusive

        overlap_start = max(repeat_start, var_start)
        overlap_end = min(repeat_end_seq, var_end)

        if overlap_start <= overlap_end:
            overlap_len = overlap_end - overlap_start + 1

            rep_start_idx = overlap_start - repeat_start
            rep_end_idx_excl = rep_start_idx + overlap_len

            ref_start_idx = overlap_start - var_start
            ref_end_idx_excl = ref_start_idx + overlap_len

            panel_sub = repeat_seq[rep_start_idx:rep_end_idx_excl]
            ref_sub = ref_base[ref_start_idx:ref_end_idx_excl]

            raw_mismatch = panel_sub != ref_sub

            # Suppress if canonically it is the same RU
            if raw_mismatch and ru and period > 0 and len(ref_sub) >= period:
                ru_panel_c = GetCanonicalMotif(ru)
                ru_vcf_c = GetCanonicalMotif(ref_sub[:period])
                suppress = ru_panel_c == ru_vcf_c

            if raw_mismatch and (not suppress):
                if not ignore_mismatch_warnings:
                    logger.warning(
                        "Reference mismatch in STR overlap: VCF %s:%d %s>%s, "
                        "STR panel %s:%d-%d RU=%s\n"
                        "Panel overlap: %s\n"
                        "VCF REF overlap: %s. ",
                        record.contig,
                        record.pos,
                        record.alleles[0],
                        record.alleles[1] if record.alts else record.alleles[0],
                        str_row.get("CHROM", record.contig),
                        repeat_start,
                        repeat_end,
                        ru,
                        panel_sub,
                        ref_sub,
                    )
                if mismatch_truth == "skip":
                    # Skip this record due to mismatch
                    if not ignore_mismatch_warnings:
                        logger.warning("Skipping record...")
                    return None

            # If user says VCF is truth, patch the panel repeat sequence overlap to match VCF REF
            # (only when overlap exists; if no overlap, nothing to patch)
            if mismatch_truth == "vcf":
                # Patch even if canonical-equivalent; "vcf truth" means the literal bases matter.
                before = repeat_seq[:rep_start_idx]
                after = repeat_seq[rep_end_idx_excl:]
                repeat_seq = (before + ref_sub + after).upper()

    alt_base = (record.alts[0] if record.alts else record.ref)
    # Apply the variant to the STR sequence
    try:
        mutated_seq = apply_variant_to_repeat(
            pos=pos,
            ref=record.ref, # To keep correct case
            alt=alt_base,
            repeat_start=repeat_start,
            repeat_seq=repeat_seq,
        )
    except Exception as e:
        logger.warning(
            "Failed to apply variant to repeat sequence; keeping REF allele. "
            "VCF %s:%d %s>%s; error=%s",
            record.contig,
            pos,
            record.ref,
            alt_base,
            str(e),
        )
        mutated_seq = repeat_seq

    mutated_seq = (mutated_seq or repeat_seq).upper()

    # RU fallback if missing
    if not ru:
        ru = repeat_seq[:period].upper() if period > 0 else ""

    # Count repeat units and perfectness
    if ru:
        ref_len = count_repeat_units(repeat_seq, ru)
        alt_len = count_repeat_units(mutated_seq, ru)
        perfect = is_perfect_repeat(repeat_seq, ru) and is_perfect_repeat(mutated_seq, ru)
    else:
        ref_len = 0
        alt_len = 0
        perfect = False

    # Prepare INFO fields - only copy safe, commonly used fields
    # to avoid type/Number mismatches with fields from original VCF
    info: Dict = {}

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
            "PERIOD": int(period),
            "REF": int(ref_len),
            "PERFECT": "TRUE" if perfect else "FALSE",
        }
    )

    # Create new record ---
    new_record = header.new_record(
        contig=record.contig,
        start=repeat_start - 1,  # pysam uses 0-based start
        stop=repeat_end,  # keep as provided
        id=".",
        alleles=(repeat_seq, mutated_seq),
        info=info,
        filter=record.filter.keys(),
    )

    # Copy FORMAT fields and set GT/REPCN
    sample_names = list(record.samples.keys())
    for sample_name in sample_names:
        old_sample = record.samples[sample_name]
        new_sample = new_record.samples[sample_name]

        for field_name in old_sample:
            if field_name in ("GT", "REPCN"):
                continue
            with contextlib.suppress(TypeError, ValueError):
                new_sample[field_name] = old_sample[field_name]

        sample_idx = sample_names.index(sample_name)
        gt = parser.get_genotype(record, sample_idx)

        if gt is None:
            new_sample["GT"] = (None, None)
            new_sample["REPCN"] = (0, 0)
            continue

        repcns = []
        for allele in gt:
            if allele == 0:
                repcns.append(int(ref_len))
            elif allele == 1:
                repcns.append(int(alt_len))
            else:
                repcns.append(0)

        if len(repcns) == 1:
            repcns = repcns + [0]
        elif len(repcns) > 2:
            repcns = repcns[:2]

        new_sample["GT"] = gt
        new_sample["REPCN"] = tuple(repcns)

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
