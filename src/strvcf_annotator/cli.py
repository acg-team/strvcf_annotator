"""Console script for strvcf_annotator."""

import argparse
import logging
import sys

from . import __version__
from .api import STRAnnotator
from .utils.validation import ValidationError


def setup_logging(verbose: bool = False):
    """Configure logging for CLI.

    Parameters
    ----------
    verbose : bool
        If True, set logging level to DEBUG
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def create_parser() -> argparse.ArgumentParser:
    """Create CLI argument parser.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="Annotate STR regions in VCF files using a BED file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Annotate single VCF file
  strvcf-annotator --input input.vcf --str-bed repeats.bed --output output.vcf

  # Batch process directory
  strvcf-annotator --input-dir vcf_files/ --str-bed repeats.bed --output-dir annotated/

  # Enable verbose logging
  strvcf-annotator --input input.vcf --str-bed repeats.bed --output output.vcf --verbose

  # Somatic mode (filter variants where tumor==normal genotypes)
  strvcf-annotator --input somatic.vcf --str-bed repeats.bed --output output.vcf --somatic-mode

  # Suppress mismatch warnings and trust VCF reference when mismatch happens
  strvcf-annotator --input input.vcf --str-bed repeats.bed --output output.vcf \\
                   --ignore-mismatch-warnings --mismatch-truth vcf
        """,
    )

    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--input", type=str, help="Path to input VCF file")
    input_group.add_argument("--input-dir", type=str, help="Directory containing input VCF files")

    # Required arguments
    parser.add_argument(
        "--str-bed",
        required=True,
        type=str,
        help="Path to BED file with STR regions (CHROM, START, END, PERIOD, RU)",
    )

    # Output options
    output_group = parser.add_mutually_exclusive_group(required=True)
    output_group.add_argument(
        "--output", type=str, help="Path to output VCF file (for single file mode)"
    )
    output_group.add_argument(
        "--output-dir", type=str, help="Directory for output VCF files (for batch mode)"
    )

    # Optional arguments
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")

    parser.add_argument(
        "--somatic-mode",
        action="store_true",
        help="Enable somatic filtering: skip variants where both samples have identical genotypes",
    )
    parser.add_argument(
        "--ignore-mismatch-warnings",
        action="store_true",
        help=(
            "Suppress warnings about reference mismatches between STR panel sequence "
            "and VCF REF allele."
        ),
    )
    parser.add_argument(
        "--jobs",
        type=int,
        help=(
            "Number of parallel jobs to use for processing. Each job processes one VCF file. If not specified, the number of jobs is automatically determined based on CPU cores, number of files, and available RAM."
        ),
    )
    parser.add_argument(
        "--mismatch-truth",
        choices=["panel", "vcf", "skip"],
        default="panel",
        help=(
            "Which source to treat as ground truth when a mismatch is detected: "
            "'panel' (default), 'vcf', or 'skip' (skip mismatching loci)."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    return parser


def validate_args(args: argparse.Namespace) -> None:
    """Validate CLI arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments

    Raises
    ------
    ValidationError
        If arguments are invalid or inconsistent
    """
    # Validate input/output consistency
    if args.input and not args.output:
        raise ValidationError("--input requires --output")

    if args.input_dir and not args.output_dir:
        raise ValidationError("--input-dir requires --output-dir")

    if args.output and not args.input:
        raise ValidationError("--output requires --input")

    if args.output_dir and not args.input_dir:
        raise ValidationError("--output-dir requires --input-dir")


def main():
    """CLI entry point with argument parsing and validation."""
    parser = create_parser()
    args = parser.parse_args()

    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    try:
        # Validate arguments
        validate_args(args)

        # Create annotator
        logger.info("Initializing STR annotator...")
        somatic_mode = getattr(args, "somatic_mode", False)
        ignore_mismatch_warnings = getattr(args, "ignore_mismatch_warnings", False)
        mismatch_truth = getattr(args, "mismatch_truth", "panel")
        jobs = getattr(args, "jobs", None)
        annotator = STRAnnotator(
            args.str_bed,
            somatic_mode=somatic_mode,
            ignore_mismatch_warnings=ignore_mismatch_warnings,
            mismatch_truth=mismatch_truth,
        )

        # Display statistics
        stats = annotator.get_statistics()
        logger.info(
            f"Loaded {stats['total_regions']} STR regions from {stats['chromosomes']} chromosomes"
        )

        # Process based on mode
        if args.input:
            # Single file mode
            logger.info(f"Processing single file: {args.input}")
            annotator.annotate_vcf_file(args.input, args.output)
            logger.info(f"Successfully wrote annotated VCF to {args.output}")

        elif args.input_dir:
            # Batch directory mode
            logger.info(f"Processing directory: {args.input_dir}")
            annotator.process_directory(args.input_dir, args.output_dir, jobs=jobs)
            logger.info(f"Successfully processed all VCF files to {args.output_dir}")

        logger.info("Annotation complete!")
        return 0

    except ValidationError as e:
        logger.error(f"Validation error: {e}")
        return 1

    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        return 1

    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=args.verbose)
        return 1


if __name__ == "__main__":
    sys.exit(main())
