"""Library API for programmatic access to STR annotation functionality."""

import logging
from typing import Iterator, Optional

import pysam

from .core.str_reference import load_str_reference
from .core.vcf_processor import annotate_vcf_to_file, generate_annotated_records, process_directory
from .parsers.base import BaseVCFParser
from .parsers.generic import GenericParser
from .utils.validation import validate_directory_path, validate_str_bed_file, validate_vcf_file

logger = logging.getLogger(__name__)


class STRAnnotator:
    """Main class for STR annotation functionality.

    Provides a high-level interface for annotating VCF files with STR
    (Short Tandem Repeat) information. Supports both single file and
    batch directory processing.

    Parameters
    ----------
    str_bed_path : str
        Path to BED file containing STR regions
    parser : BaseVCFParser, optional
        Custom parser for genotype extraction. Uses GenericParser if None.
    somatic_mode : bool, optional
        Enable somatic filtering mode. When True, skips variants where both
        samples (tumor/normal) have identical genotypes. Default is False.

    Attributes
    ----------
    str_bed_path : str
        Path to STR BED file
    str_df : pd.DataFrame
        Loaded STR reference data
    parser : BaseVCFParser
        Parser for genotype extraction
    somatic_mode : bool
        Whether somatic filtering is enabled

    Examples
    --------
    >>> annotator = STRAnnotator('str_regions.bed')
    >>> annotator.annotate_vcf_file('input.vcf', 'output.vcf')

    >>> # Batch process directory
    >>> annotator.process_directory('input_dir/', 'output_dir/')

    >>> # Stream processing
    >>> vcf_in = pysam.VariantFile('input.vcf')
    >>> for record in annotator.annotate_vcf_stream(vcf_in):
    ...     print(record)
    """

    def __init__(
        self, str_bed_path: str, parser: Optional[BaseVCFParser] = None, somatic_mode: bool = False
    ):
        """Initialize STR annotator with reference and parser.

        Parameters
        ----------
        str_bed_path : str
            Path to BED file with STR regions
        parser : BaseVCFParser, optional
            Custom parser for genotype extraction
        somatic_mode : bool, optional
            Enable somatic filtering (skip variants where tumor==normal genotypes).
            Default is False.

        Raises
        ------
        ValidationError
            If STR BED file is invalid
        """
        # Validate and load STR reference
        validate_str_bed_file(str_bed_path)
        self.str_bed_path = str_bed_path
        self.str_df = load_str_reference(str_bed_path)

        # Set parser
        self.parser = parser if parser is not None else GenericParser()

        # Set somatic mode
        self.somatic_mode = somatic_mode

        logger.info(f"Loaded {len(self.str_df)} STR regions from {str_bed_path}")

    def annotate_vcf_file(self, input_path: str, output_path: str) -> None:
        """Annotate single VCF file.

        Reads a VCF file, annotates variants overlapping with STR regions,
        and writes the annotated records to an output file.

        Parameters
        ----------
        input_path : str
            Path to input VCF file
        output_path : str
            Path to output VCF file

        Raises
        ------
        ValidationError
            If input VCF file is invalid

        Examples
        --------
        >>> annotator = STRAnnotator('str_regions.bed')
        >>> annotator.annotate_vcf_file('input.vcf', 'output.vcf')
        """
        # Validate input
        validate_vcf_file(input_path)

        # Annotate
        logger.info(f"Annotating {input_path}...")
        annotate_vcf_to_file(
            input_path, self.str_df, output_path, self.parser, somatic_mode=self.somatic_mode
        )
        logger.info(f"Wrote annotated VCF to {output_path}")

    def annotate_vcf_stream(self, vcf_in: pysam.VariantFile) -> Iterator[pysam.VariantRecord]:
        """Annotate VCF records from stream.

        Generator that yields annotated VCF records from an open VCF file.
        Useful for streaming processing or custom workflows.

        Parameters
        ----------
        vcf_in : pysam.VariantFile
            Open VCF file object

        Yields
        ------
        pysam.VariantRecord
            Annotated VCF records

        Examples
        --------
        >>> annotator = STRAnnotator('str_regions.bed')
        >>> vcf_in = pysam.VariantFile('input.vcf')
        >>> for record in annotator.annotate_vcf_stream(vcf_in):
        ...     # Process record
        ...     print(record.info['RU'])
        """
        yield from generate_annotated_records(
            vcf_in, self.str_df, self.parser, somatic_mode=self.somatic_mode
        )

    def process_directory(self, input_dir: str, output_dir: str) -> None:
        """Batch process directory of VCF files.

        Processes all VCF files in a directory and writes annotated versions
        to the output directory. Skips files that have already been processed.

        Parameters
        ----------
        input_dir : str
            Directory containing input VCF files
        output_dir : str
            Directory for output VCF files (created if doesn't exist)

        Raises
        ------
        ValidationError
            If input directory is invalid

        Examples
        --------
        >>> annotator = STRAnnotator('str_regions.bed')
        >>> annotator.process_directory('vcf_files/', 'annotated_vcfs/')
        """
        # Validate directories
        validate_directory_path(input_dir, must_exist=True)
        validate_directory_path(output_dir, must_exist=False, create=True)

        # Process directory
        logger.info(f"Processing VCF files in {input_dir}...")
        process_directory(
            input_dir, self.str_bed_path, output_dir, self.parser, somatic_mode=self.somatic_mode
        )
        logger.info(f"Batch processing complete. Output in {output_dir}")

    def get_str_at_position(self, chrom: str, pos: int) -> Optional[dict]:
        """Get STR region at specific genomic position.

        Parameters
        ----------
        chrom : str
            Chromosome name
        pos : int
            Genomic position (1-based)

        Returns
        -------
        Optional[dict]
            STR region data if position is within an STR, None otherwise

        Examples
        --------
        >>> annotator = STRAnnotator('str_regions.bed')
        >>> str_region = annotator.get_str_at_position('chr1', 1000000)
        >>> if str_region:
        ...     print(f"Repeat unit: {str_region['RU']}")
        """
        from .core.str_reference import get_str_at_position

        return get_str_at_position(self.str_df, chrom, pos)

    def get_statistics(self) -> dict:
        """Get statistics about loaded STR regions.

        Returns
        -------
        dict
            Statistics including total regions, chromosomes, repeat units

        Examples
        --------
        >>> annotator = STRAnnotator('str_regions.bed')
        >>> stats = annotator.get_statistics()
        >>> print(f"Total STR regions: {stats['total_regions']}")
        """
        stats = {
            "total_regions": len(self.str_df),
            "chromosomes": self.str_df["CHROM"].nunique(),
            "unique_repeat_units": self.str_df["RU"].nunique(),
            "period_distribution": self.str_df["PERIOD"].value_counts().to_dict(),
            "mean_repeat_count": self.str_df["COUNT"].mean(),
            "median_repeat_count": self.str_df["COUNT"].median(),
        }
        return stats


def annotate_vcf(
    input_vcf: str, str_bed: str, output_vcf: str, parser: Optional[BaseVCFParser] = None
) -> None:
    """Convenience function for single VCF annotation.

    Simple function interface for annotating a single VCF file.

    Parameters
    ----------
    input_vcf : str
        Path to input VCF file
    str_bed : str
        Path to STR BED file
    output_vcf : str
        Path to output VCF file
    parser : BaseVCFParser, optional
        Custom parser for genotype extraction

    Examples
    --------
    >>> from strvcf_annotator import annotate_vcf
    >>> annotate_vcf('input.vcf', 'str_regions.bed', 'output.vcf')
    """
    annotator = STRAnnotator(str_bed, parser)
    annotator.annotate_vcf_file(input_vcf, output_vcf)
