"""Generic parser for standard VCF format fields."""

from typing import Any, Dict, Optional, Tuple

import pysam

from .base import BaseVCFParser


class GenericParser(BaseVCFParser):
    """Generic parser for standard VCF format fields.

    Handles standard VCF FORMAT fields including GT (genotype), AD (allelic depth),
    and DP (total depth). Provides robust error handling for missing or invalid data.
    """

    def get_genotype(
        self, record: pysam.VariantRecord, sample_idx: int
    ) -> Optional[Tuple[int, int]]:
        """Extract GT field, return None for missing/invalid genotypes.

        Parameters
        ----------
        record : pysam.VariantRecord
            The VCF record to extract genotype from
        sample_idx : int
            Index of the sample in the record

        Returns
        -------
        Optional[Tuple[int, int]]
            Genotype as tuple of allele indices, or None if missing/invalid
        """
        try:
            samples = list(record.samples.values())
            if sample_idx >= len(samples):
                return None

            sample_data = samples[sample_idx]
            gt = sample_data.get("GT", (None, None))

            # Handle missing or invalid genotypes
            if gt == (None, None) or gt == (".", ".") or gt is None:
                return None

            # Handle string genotypes
            if isinstance(gt, tuple) and len(gt) == 2:
                # Convert string alleles to integers if needed
                alleles = []
                for allele in gt:
                    if allele is None or allele == ".":
                        return None
                    if isinstance(allele, str):
                        if allele == ".":
                            return None
                        allele = int(allele)
                    alleles.append(allele)
                return tuple(alleles)

            return None

        except (KeyError, ValueError, TypeError, IndexError):
            return None

    def has_variant(self, record: pysam.VariantRecord, sample_idx: int) -> bool:
        """Check variant presence using GT or alternative evidence.

        Parameters
        ----------
        record : pysam.VariantRecord
            The VCF record to check
        sample_idx : int
            Index of the sample in the record

        Returns
        -------
        bool
            True if variant is present, False otherwise
        """
        # First try genotype
        gt = self.get_genotype(record, sample_idx)
        if gt is not None:
            # Check if any allele is non-reference (not 0)
            return any(allele != 0 for allele in gt)

        # If GT is missing, check for alternative evidence (AD, DP)
        try:
            samples = list(record.samples.values())
            if sample_idx >= len(samples):
                return False

            sample_data = samples[sample_idx]

            # Check allelic depth
            ad = sample_data.get("AD", None)
            if ad is not None and len(ad) > 1:
                # If ALT allele has reads, variant is present
                return ad[1] > 0

            return False

        except (KeyError, ValueError, TypeError, IndexError):
            return False

    def extract_info(self, record: pysam.VariantRecord, sample_idx: int) -> Dict[str, Any]:
        """Extract AD, DP and other standard FORMAT fields.

        Parameters
        ----------
        record : pysam.VariantRecord
            The VCF record to extract information from
        sample_idx : int
            Index of the sample in the record

        Returns
        -------
        Dict[str, Any]
            Dictionary of FORMAT fields (AD, DP, etc.)
        """
        info = {}

        try:
            samples = list(record.samples.values())
            if sample_idx >= len(samples):
                return info

            sample_data = samples[sample_idx]

            # Extract standard fields
            if "AD" in sample_data:
                info["AD"] = sample_data["AD"]

            if "DP" in sample_data:
                info["DP"] = sample_data["DP"]

            if "GQ" in sample_data:
                info["GQ"] = sample_data["GQ"]

            if "PL" in sample_data:
                info["PL"] = sample_data["PL"]

        except (KeyError, ValueError, TypeError, IndexError):
            pass

        return info

    def validate_record(self, record: pysam.VariantRecord) -> bool:
        """Validate record has required fields for generic parsing.

        Parameters
        ----------
        record : pysam.VariantRecord
            The VCF record to validate

        Returns
        -------
        bool
            True if record is valid, False otherwise
        """
        try:
            # Check basic record structure
            if record is None:
                return False

            # Check that record has samples
            if not record.samples:
                return False

            # Check that record has alleles
            return not (not record.alleles or len(record.alleles) < 2)

        except (AttributeError, TypeError):
            return False
