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
    """
    Main class for STR annotation functionality.

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
        Enable somatic filtering. When ``True``, variants where all samples
        have identical genotypes are skipped. If ``None``, the value set on
        the annotator instance is used.
    ignore_mismatch_warnings : bool, optional
        If ``True``, suppress warnings about reference mismatches between
        the STR panel sequence and the VCF ``REF`` allele. If ``None``,
        the value set on the annotator instance is used.
    mismatch_truth : str, optional
        Specifies which source is treated as ground truth when a mismatch
        between the STR panel and VCF ``REF`` allele is detected.

        Allowed values are:

        - ``"panel"``: trust the STR panel repeat sequence (default)
        - ``"vcf"``: trust the VCF ``REF`` allele and patch the overlapping
        panel sequence
        - ``"skip"``: skip variants with mismatches entirely

        If ``None``, the value set on the annotator instance is used.

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
    ...     print(record.info['RU'])
    """

    def __init__(
        self,
        str_bed_path: str,
        parser: Optional[BaseVCFParser] = None,
        somatic_mode: bool = False,
        ignore_mismatch_warnings: bool = False,
        mismatch_truth: str = "panel",  # "panel" | "vcf" | "skip"
    ):
        validate_str_bed_file(str_bed_path)
        self.str_bed_path = str_bed_path
        self.str_df = load_str_reference(str_bed_path)

        self.parser = parser if parser is not None else GenericParser()
        self.somatic_mode = somatic_mode

        self.ignore_mismatch_warnings = ignore_mismatch_warnings
        self.mismatch_truth = mismatch_truth

        logger.info(f"Loaded {len(self.str_df)} STR regions from {str_bed_path}")

    def annotate_vcf_file(
        self,
        input_path: str,
        output_path: str,
        *,
        somatic_mode: Optional[bool] = None,
        ignore_mismatch_warnings: Optional[bool] = None,
        mismatch_truth: Optional[str] = None,
    ) -> None:
        """
        Annotate a single VCF file with STR information.

        This method reads a VCF file, annotates variants that overlap STR regions,
        and writes the annotated records to an output VCF file.

        Parameters
        ----------
        input_path : str
            Path to the input VCF file.
        output_path : str
            Path to the output annotated VCF file.
        somatic_mode : bool, optional
            Enable somatic filtering. When ``True``, variants where all samples
            have identical genotypes are skipped. If ``None``, the value set on
            the annotator instance is used.
        ignore_mismatch_warnings : bool, optional
            If ``True``, suppress warnings about reference mismatches between
            the STR panel sequence and the VCF ``REF`` allele. If ``None``,
            the value set on the annotator instance is used.
        mismatch_truth : str, optional
            Specifies which source is treated as ground truth when a mismatch
            between the STR panel and VCF ``REF`` allele is detected.

            Allowed values are:

            - ``"panel"``: trust the STR panel repeat sequence (default)
            - ``"vcf"``: trust the VCF ``REF`` allele and patch the overlapping
            panel sequence
            - ``"skip"``: skip variants with mismatches entirely

            If ``None``, the value set on the annotator instance is used.

        Raises
        ------
        ValidationError
            If the input VCF file is invalid.

        Examples
        --------
        >>> annotator = STRAnnotator("str_regions.bed")
        >>> annotator.annotate_vcf_file("input.vcf", "output.vcf")
        """
        # Validate input
        validate_vcf_file(input_path)
        imw = (
            self.ignore_mismatch_warnings
            if ignore_mismatch_warnings is None
            else ignore_mismatch_warnings
        )
        mtruth = self.mismatch_truth if mismatch_truth is None else mismatch_truth
        smode = self.somatic_mode if somatic_mode is None else somatic_mode

        # Annotate
        logger.info(f"Annotating {input_path}...")
        annotate_vcf_to_file(
            input_path,
            self.str_df,
            output_path,
            self.parser,
            somatic_mode=smode,
            ignore_mismatch_warnings=imw,
            mismatch_truth=mtruth,
        )
        logger.info(f"Wrote annotated VCF to {output_path}")

    def annotate_vcf_stream(
        self,
        vcf_in: pysam.VariantFile,
        *,
        somatic_mode: Optional[bool] = None,
        ignore_mismatch_warnings: Optional[bool] = None,
        mismatch_truth: Optional[str] = None,
    ) -> Iterator[pysam.VariantRecord]:
        """
        Annotate VCF records from an open stream.

        This generator yields annotated VCF records from an already opened
        ``pysam.VariantFile`` object. It is useful for streaming workflows
        or custom processing pipelines.

        Parameters
        ----------
        vcf_in : pysam.VariantFile
            Open VCF file object to read variants from.
        somatic_mode : bool, optional
            Enable somatic filtering. When ``True``, variants where all samples
            have identical genotypes are skipped. If ``None``, the value set on
            the annotator instance is used.
        ignore_mismatch_warnings : bool, optional
            If ``True``, suppress warnings about reference mismatches between
            the STR panel sequence and the VCF ``REF`` allele. If ``None``,
            the value set on the annotator instance is used.
        mismatch_truth : str, optional
            Specifies which source is treated as ground truth when a mismatch
            between the STR panel and VCF ``REF`` allele is detected.

            Allowed values are:

            - ``"panel"``: trust the STR panel repeat sequence (default)
            - ``"vcf"``: trust the VCF ``REF`` allele and patch the overlapping
            panel sequence
            - ``"skip"``: skip variants with mismatches entirely

            If ``None``, the value set on the annotator instance is used.

        Yields
        ------
        pysam.VariantRecord
            Annotated VCF records.

        Examples
        --------
        >>> annotator = STRAnnotator("str_regions.bed")
        >>> vcf_in = pysam.VariantFile("input.vcf")
        >>> for record in annotator.annotate_vcf_stream(vcf_in):
        ...     print(record.info["RU"])
        """
        imw = (
            self.ignore_mismatch_warnings
            if ignore_mismatch_warnings is None
            else ignore_mismatch_warnings
        )
        mtruth = self.mismatch_truth if mismatch_truth is None else mismatch_truth
        smode = self.somatic_mode if somatic_mode is None else somatic_mode

        yield from generate_annotated_records(
            vcf_in=vcf_in,
            str_df=self.str_df,
            parser=self.parser,
            somatic_mode=smode,
            ignore_mismatch_warnings=imw,
            mismatch_truth=mtruth,
        )

    def process_directory(
        self,
        input_dir: str,
        output_dir: str,
        *,
        somatic_mode: Optional[bool] = None,
        ignore_mismatch_warnings: Optional[bool] = None,
        mismatch_truth: Optional[str] = None,
    ) -> None:
        """
        Batch process a directory of VCF files.

        This method processes all VCF files in the input directory and writes
        annotated versions to the output directory. Files that have already
        been processed are skipped automatically.

        Parameters
        ----------
        input_dir : str
            Directory containing input VCF files.
        output_dir : str
            Directory where annotated VCF files will be written. The directory
            is created if it does not already exist.
        somatic_mode : bool, optional
            Enable somatic filtering. When ``True``, variants where all samples
            have identical genotypes are skipped. If ``None``, the value set on
            the annotator instance is used.
        ignore_mismatch_warnings : bool, optional
            If ``True``, suppress warnings about reference mismatches between
            the STR panel sequence and the VCF ``REF`` allele. If ``None``,
            the value set on the annotator instance is used.
        mismatch_truth : str, optional
            Specifies which source is treated as ground truth when a mismatch
            between the STR panel and VCF ``REF`` allele is detected.

            Allowed values are:

            - ``"panel"``: trust the STR panel repeat sequence (default)
            - ``"vcf"``: trust the VCF ``REF`` allele and patch the overlapping
            panel sequence
            - ``"skip"``: skip variants with mismatches entirely

            If ``None``, the value set on the annotator instance is used.

        Raises
        ------
        ValidationError
            If the input directory is invalid.

        Examples
        --------
        >>> annotator = STRAnnotator("str_regions.bed")
        >>> annotator.process_directory("vcf_files/", "annotated_vcfs/")
        """
        # Validate directories
        validate_directory_path(input_dir, must_exist=True)
        validate_directory_path(output_dir, must_exist=False, create=True)

        imw = (
            self.ignore_mismatch_warnings
            if ignore_mismatch_warnings is None
            else ignore_mismatch_warnings
        )
        mtruth = self.mismatch_truth if mismatch_truth is None else mismatch_truth
        smode = self.somatic_mode if somatic_mode is None else somatic_mode
        # Process directory
        logger.info(f"Processing VCF files in {input_dir}...")
        process_directory(
            input_dir=input_dir,
            str_bed_path=self.str_bed_path,
            output_dir=output_dir,
            parser=self.parser,
            somatic_mode=smode,
            ignore_mismatch_warnings=imw,
            mismatch_truth=mtruth,
        )
        logger.info(f"Batch processing complete. Output in {output_dir}")

    def get_str_at_position(self, chrom: str, pos: int) -> Optional[dict]:
        """
        Get STR region at specific genomic position.

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
        """
        Get statistics about loaded STR regions.

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
    input_vcf: str,
    str_bed: str,
    output_vcf: str,
    parser: Optional[BaseVCFParser] = None,
    *,
    somatic_mode: Optional[bool] = None,
    ignore_mismatch_warnings: Optional[bool] = None,
    mismatch_truth: Optional[str] = None,
) -> None:
    """
    Convenience function for single VCF annotation.

    This is a simplified functional interface for annotating a single VCF
    file with STR information, without explicitly creating an
    :class:`STRAnnotator` instance.

    Parameters
    ----------
    input_vcf : str
        Path to the input VCF file.
    str_bed : str
        Path to the STR BED file.
    output_vcf : str
        Path to the output annotated VCF file.
    parser : BaseVCFParser, optional
        Custom parser for genotype extraction. If not provided,
        :class:`GenericParser` is used.
    somatic_mode : bool, optional
        Enable somatic filtering. When ``True``, variants where all samples
        have identical genotypes are skipped. If ``None``, the default
        behavior is used.
    ignore_mismatch_warnings : bool, optional
        If ``True``, suppress warnings about reference mismatches between
        the STR panel sequence and the VCF ``REF`` allele. If ``None``,
        the default behavior is used.
    mismatch_truth : str, optional
        Specifies which source is treated as ground truth when a mismatch
        between the STR panel and the VCF ``REF`` allele is detected.

        Allowed values are:

        - ``"panel"``: trust the STR panel repeat sequence (default)
        - ``"vcf"``: trust the VCF ``REF`` allele and patch the overlapping
          panel sequence
        - ``"skip"``: skip variants with mismatches entirely

        If ``None``, the default behavior is used.

    Examples
    --------
    >>> from strvcf_annotator import annotate_vcf
    >>> annotate_vcf("input.vcf", "str_regions.bed", "output.vcf")
    """
    annotator = STRAnnotator(str_bed, parser)
    annotator.annotate_vcf_file(
        input_vcf,
        output_vcf,
        somatic_mode=somatic_mode,
        ignore_mismatch_warnings=ignore_mismatch_warnings,
        mismatch_truth=mismatch_truth,
    )
