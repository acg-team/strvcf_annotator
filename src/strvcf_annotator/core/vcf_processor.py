"""VCF file processing and workflow management."""

import logging
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import psutil
import pysam

from ..parsers.base import BaseVCFParser
from ..parsers.generic import GenericParser
from ..utils.vcf_utils import chrom_to_order
from .annotation import build_new_record, make_modified_header, should_skip_genotype

logger = logging.getLogger(__name__)
# Globals initialized once per worker process
WORKER_CONFIG = None


@dataclass(frozen=True)
class WorkerConfig:
    """Configuration container for worker processes.

    Stores parameters required by each worker to annotate VCF files.
    The configuration is passed once during worker initialization and
    reused for all tasks processed by that worker.

    Attributes
    ----------
    str_panel_gz : str
        Path to the BGZF-compressed, tabix-indexed STR reference file.
    somatic_mode : bool, optional
        Enable somatic filtering. When True, skips variants where both samples
        have identical genotypes. Default is False.
    ignore_mismatch_warnings : bool, optional
        If True, suppresses warnings about reference mismatches between the
        STR panel and VCF REF alleles. Default is False.
    mismatch_truth : str, optional
        Specifies which source to consider as ground truth for mismatches.
        Options are "panel", "vcf", or "skip". Default is "panel".

    Notes
    -----
    - The dataclass is frozen to ensure the configuration remains
      immutable once workers are initialized.
    - Instances of this class are passed to `worker_init`, which loads
      the STR reference and exposes these settings to worker tasks.
    """

    str_panel_gz: str
    somatic_mode: bool
    ignore_mismatch_warnings: bool
    mismatch_truth: str
    parser: BaseVCFParser


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
    str_panel_gz: str,
    parser: BaseVCFParser = None,
    somatic_mode: bool = False,
    ignore_mismatch_warnings: bool = False,
    mismatch_truth: str = "panel",  # "panel" | "vcf" | "skip"
) -> Iterator[pysam.VariantRecord]:
    """Generator yielding annotated VCF records.

    Processes VCF records and yields annotated records for variants that
    overlap with STR regions. Handles sorting if needed and optionally filters
    records based on genotype criteria. When multiple STR regions overlap the same POS,
    try all overlapping STR candidates and pick the first that produces a
    meaningful STR allele change.

    Parameters
    ----------
    vcf_in : pysam.VariantFile
        Input VCF file
    str_panel_gz : str
        Path to BGZF-compressed, tabix-indexed STR reference file.
    parser : BaseVCFParser, optional
        Parser for genotype extraction. Uses GenericParser if None.
    somatic_mode : bool, optional
        Enable somatic filtering. When True, skips variants where both samples
        have identical genotypes. Default is False.
    ignore_mismatch_warnings : bool, optional
        If True, suppresses warnings about reference mismatches between the
        STR panel and VCF REF alleles. Default is False.
    mismatch_truth : str, optional
        Specifies which source to consider as ground truth for mismatches.
        Options are "panel", "vcf", or "skip". Default is "panel".

    Yields
    ------
    pysam.VariantRecord
        Annotated VCF records

    Notes
    -----
    - Automatically sorts VCF if not sorted
    - Skips records without STR overlap
    - If somatic_mode=True, filters records with identical genotypes
    - ignore_mismatch_warnings controls logging of reference mismatches
    - mismatch_truth controls which source is considered ground truth for mismatches
    """
    if parser is None:
        parser = GenericParser()

    mismatch_truth = (mismatch_truth or "panel").lower()
    if mismatch_truth not in {"panel", "vcf", "skip"}:
        raise ValueError("mismatch_truth must be one of: 'panel', 'vcf', 'skip'")

    header = make_modified_header(vcf_in)

    # Check if VCF is sorted
    if not check_vcf_sorted(vcf_in):
        logger.warning("Input VCF is not sorted - sorting in memory.")
        records = reset_and_sort_vcf(vcf_in)
    else:
        vcf_in.reset()
        records = vcf_in.fetch()

    # Open tabix once for the whole generator (fast, avoids reopen per record)
    tbx = pysam.TabixFile(str_panel_gz)

    skipped_count = 0
    try:
        for record in records:
            # Skip based on genotype filtering (only if somatic_mode enabled)
            if somatic_mode and should_skip_genotype(record, parser):
                skipped_count += 1
                logger.debug(
                    f"Skipped {record.contig}:{record.pos} - identical genotypes (somatic mode)"
                )
                continue

            # Query tabix for a 1bp window at record.pos.
            # Tabix uses 0-based half-open coordinates.
            query_start = max(0, record.pos - 1)
            query_end = record.pos
            try:
                candidates = tbx.fetch(record.chrom, query_start, query_end)
            except ValueError:
                # Chromosome not present in the index
                continue

            def parse_str_line(line: str) -> Optional[Dict]:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    return None
                try:
                    row = {
                        "CHROM": parts[0],
                        "START": int(parts[1]) + 1, # Convert to 1-based
                        "END": int(parts[2]),
                        "PERIOD": int(parts[3]),
                        "RU": parts[4],
                    }
                    row["COUNT"] = (row["END"] - row["START"] + 1) / row["PERIOD"]
                    return row
                except ValueError:
                    return None

            chosen = None
            for line in candidates:
                str_row = parse_str_line(line)
                if str_row is None:
                    continue
                # True overlap check in 1-based coordinates
                if not (str_row["START"] <= record.pos <= str_row["END"]):
                    continue

                new_record = build_new_record(
                    record,
                    str_row,
                    header,
                    parser,
                    ignore_mismatch_warnings=ignore_mismatch_warnings,
                    mismatch_truth=mismatch_truth,
                )
                # Skip in case of mismatch and mismatch_truth is "skip"
                if new_record is None:
                    continue

                if new_record.alleles[0] != new_record.alleles[1]:
                    chosen = new_record
                    break
                # If variant effectively doesn't change STR allele (e.g., indel outside STR after normalization),
                # treat this STR row as not applicable and try next overlapping STR.

            if chosen is not None:
                yield chosen
    finally:
        tbx.close()

    # Log summary if records were skipped
    if skipped_count > 0:
        logger.warning(
            f"Skipped {skipped_count} records due to identical genotypes (somatic filtering enabled)"
        )


def annotate_vcf_to_file(
    vcf_path: str,
    str_panel_gz: str,
    output_path: str,
    parser: BaseVCFParser = None,
    somatic_mode: bool = False,
    ignore_mismatch_warnings: bool = False,
    mismatch_truth: str = "panel",
) -> None:
    """Process VCF file and write annotated output.

    Reads a VCF file, annotates variants that overlap with STR regions,
    and writes the annotated records to an output file.

    Parameters
    ----------
    vcf_path : str
        Path to input VCF file
    str_panel_gz : str
        Path to BGZF-compressed, tabix-indexed STR reference file.
    output_path : str
        Path to output VCF file
    parser : BaseVCFParser, optional
        Parser for genotype extraction. Uses GenericParser if None.
    somatic_mode : bool, optional
        Enable somatic filtering. When True, skips variants where both samples
        have identical genotypes. Default is False.
    ignore_mismatch_warnings : bool
        If True, do not log mismatch warnings between VCF REF and panel repeat sequence.
        Annotation continues regardless.
    mismatch_truth : str
        Which source to treat as correct when there is a mismatch:
          - "panel": trust panel repeat sequence (default behavior)
          - "vcf": trust VCF REF. patch the panel repeat sequence overlap to match VCF REF
          - "skip": skip record with mismatch

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
    for record in generate_annotated_records(
        vcf_in,
        str_panel_gz,
        parser,
        somatic_mode=somatic_mode,
        ignore_mismatch_warnings=ignore_mismatch_warnings,
        mismatch_truth=mismatch_truth,
    ):
        vcf_out.write(record)
        written_count += 1

    vcf_out.close()
    vcf_in.close()

    logger.info(f"Wrote {written_count} annotated records to {output_path}")


def annotate_one_vcf(task: Tuple[str, str]) -> str:
    """Annotate a single VCF file in a worker process.

    Runs `annotate_vcf_to_file` for one input VCF and writes the annotated VCF
    to the given output path. Created to be executed inside a process pool.

    Parameters
    ----------
    task : Tuple[str, str]
        (vcf_path, output_path) pair, where:
          - vcf_path is the input VCF (optionally gzipped)
          - output_path is the target annotated VCF path

    Returns
    -------
    str
        Path to the produced output VCF file.

    Notes
    -----
    - Expects STR reference (STR_DF) and worker configuration (WORKER_CONFIG)
      to be initialized once per worker via `worker_init`.
    - Instantiates a parser inside the worker to avoid pickling issues.
    """
    global WORKER_CONFIG

    vcf_path, output_path = task

    annotate_vcf_to_file(
        vcf_path=vcf_path,
        str_panel_gz=WORKER_CONFIG.str_panel_gz,
        output_path=output_path,
        parser=WORKER_CONFIG.parser,
        somatic_mode=WORKER_CONFIG.somatic_mode,
        ignore_mismatch_warnings=WORKER_CONFIG.ignore_mismatch_warnings,
        mismatch_truth=WORKER_CONFIG.mismatch_truth,
    )
    return output_path


def get_available_ram_bytes() -> int:
    """Get available system RAM.

    Returns
    -------
    int
        Available RAM in bytes.
    """
    return int(psutil.virtual_memory().available)


def estimate_ram_per_worker_bytes(vcf_paths: List[str]) -> int:
    """Estimate RAM usage per worker for VCF annotation.

    Provides an estimate of how much memory a single worker process
    might require while annotating one VCF. This estimate is used to cap the
    number of concurrent workers to reduce the risk of out-of-memory (OOM)
    crashes.

    Parameters
    ----------
    vcf_paths : list[str]
        List of input VCF paths that will be processed.

    Returns
    -------
    int
        Estimated RAM usage per worker in bytes.

    Notes
    -----
    - If a VCF is not sorted, the current pipeline may load all records into
      memory for sorting, which can drastically increase memory usage.
    - Even for sorted VCFs, pysam/htslib buffers plus Python object overhead
      can be substantial.
    - Each worker loads the STR reference once. The STR DataFrame and derived
      Python objects often consume several times the BED file size on disk.

    Heuristic
    ---------
    - Identify the largest input file size on disk.
    - If the largest file is gzipped, assume a higher expansion factor for the
      working set (e.g., decompression + object overhead).
    - Add a fixed overhead to account for Python/pysam allocations.
    - Add STR panel RAM estimate as: str_panel_factor * BED_size_on_disk.

    This is intentionally conservative to avoid OOM.
    """
    max_size = 0
    max_path = ""
    for p in vcf_paths:
        try:
            s = os.path.getsize(p)
        except OSError:
            s = 0
        if s > max_size:
            max_size = s
            max_path = p

    is_gz = max_path.endswith(".gz")
    expansion_factor = 5 if is_gz else 2  # VCF decompression + object overhead

    fixed_overhead = 700 * 1024**2  # ~700MB overhead per worker

    estimate = fixed_overhead + (expansion_factor * max_size)

    # Clamp to sane minimum/maximum to avoid weird estimates on tiny/huge files
    min_estimate = 1 * 1024**3  # 1 GB
    max_estimate = 120 * 1024**3  # 120 GB
    return int(min(max(estimate, min_estimate), max_estimate))


def compute_jobs_auto(n_files: int, vcf_paths: List[str]) -> int:
    """Compute an automatic number of concurrent workers.

    Chooses a default number of parallel jobs for processing a directory of VCF
    files, balancing CPU capacity and memory constraints.

    Parameters
    ----------
    n_files : int
        Number of VCF files that will be processed (after skipping outputs that
        already exist).
    vcf_paths : list[str]
        List of VCF paths used to estimate per-worker memory usage.

    Returns
    -------
    int
        Recommended number of concurrent worker processes (at least 1).

    Notes
    -----
    The selection follows:
    - jobs_auto = min(cpu_cores, n_files)
    - jobs_auto = min(jobs_auto, floor(available_ram / ram_per_worker_estimate))

    If available RAM cannot be determined, the CPU-based limit is used.
    """
    cpu_cores = multiprocessing.cpu_count()
    jobs_auto = min(cpu_cores, n_files)

    available = get_available_ram_bytes()
    ram_per_worker = estimate_ram_per_worker_bytes(vcf_paths)

    if available > 0 and ram_per_worker > 0:
        ram_cap = max(1, available // ram_per_worker)
        jobs_auto = min(jobs_auto, int(ram_cap))

    return max(1, int(jobs_auto))


def worker_init(config: WorkerConfig) -> None:
    """Initialize worker process state.

    Called once when a worker process starts. Stores configuration
    values so they can be reused for all VCF files processed by that worker.

    Parameters
    ----------
    config : WorkerConfig
        Configuration object containing:
          - str_panel_gz : path to STR panel BGZF-compressed, tabix-indexed reference file
          - somatic_mode : whether somatic filtering is enabled
          - ignore_mismatch_warnings : whether to suppress mismatch warnings
          - mismatch_truth : rule for handling panel/VCF mismatches

    Notes
    -----
    - This function is used as the `initializer` for `ProcessPoolExecutor`.
    """
    global WORKER_CONFIG
    WORKER_CONFIG = config


def process_directory(
    input_dir: str,
    str_panel_gz: str,
    output_dir: str,
    parser: BaseVCFParser = None,
    somatic_mode: bool = False,
    ignore_mismatch_warnings: bool = False,
    mismatch_truth: str = "panel",
    jobs: int = None,
) -> None:
    """Batch process directory of VCF files.

    Processes all VCF files in a directory and writes annotated versions
    to the output directory.

    Parameters
    ----------
    input_dir : str
        Directory containing input VCF files
    str_panel_gz : str
        Path to BGZF-compressed, tabix-indexed STR panel reference file
    output_dir : str
        Directory for output VCF files
    parser : BaseVCFParser, optional
        Parser for genotype extraction. Uses GenericParser if None.
    somatic_mode : bool, optional
        Enable somatic filtering. When True, skips variants where both samples
        have identical genotypes. Default is False.
    ignore_mismatch_warnings : bool
        If True, do not log mismatch warnings between VCF REF and panel repeat sequence.
        Annotation continues regardless.
    mismatch_truth : str
        Which source to treat as correct when there is a mismatch:
          - "panel": trust panel repeat sequence (default behavior)
          - "vcf": trust VCF REF. patch the panel repeat sequence overlap to match VCF REF
          - "skip": skip record with mismatch
    jobs: int, optional
        - If jobs is None: compute jobs automatically:
            jobs_auto = min(cpu_cores, n_files)
            jobs_auto = min(jobs_auto, floor(available_ram / ram_per_worker_estimate))
        - If jobs is provided: use it exactly.
    """

    if parser is None:
        parser = GenericParser()

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Process each VCF file
    input_path = Path(input_dir)

    tasks = []
    vcf_paths_for_estimate = []

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

            tasks.append((str(vcf_file), str(output_file)))
            vcf_paths_for_estimate.append(str(vcf_file))
    if not tasks:
        logger.info("No VCF files to process.")
        return

    # Process larger files first (reduces idle time at end)
    def get_size(path: str) -> int:
        try:
            return os.path.getsize(path)
        except OSError:
            return 0

    tasks.sort(key=lambda t: get_size(t[0]), reverse=True)

    n_files = len(tasks)

    if jobs is None:
        jobs_to_use = compute_jobs_auto(n_files=n_files, vcf_paths=vcf_paths_for_estimate)
        avail_gb = get_available_ram_bytes() / (1024**3) if get_available_ram_bytes() else 0
        est_gb = estimate_ram_per_worker_bytes(vcf_paths_for_estimate) / (1024**3)
        logger.info(
            f"Auto jobs={jobs_to_use} (files={n_files}, cpu={multiprocessing.cpu_count()}, "
            f"available_ram≈{avail_gb:.1f}GB, est_ram_per_worker≈{est_gb:.1f}GB)"
        )
    else:
        jobs_to_use = max(1, int(jobs))
        logger.info(f"Using fixed jobs={jobs_to_use} (files={n_files})")

    config = WorkerConfig(
        str_panel_gz=str_panel_gz,
        somatic_mode=somatic_mode,
        ignore_mismatch_warnings=ignore_mismatch_warnings,
        mismatch_truth=mismatch_truth,
        parser=parser,
    )

    with ProcessPoolExecutor(
        max_workers=jobs_to_use,
        initializer=worker_init,
        initargs=(config,),
    ) as executor:
        future_to_task = {executor.submit(annotate_one_vcf, t): t for t in tasks}

        for future in as_completed(future_to_task):
            vcf_path, out_path = future_to_task[future]
            try:
                produced = future.result()
                logger.info(f"Done: {Path(vcf_path).name} → {produced}")
            except Exception:
                logger.exception(f"Failed: {Path(vcf_path).name} → {out_path}")
