"""Abstract base class for VCF parsers."""

from abc import ABC, abstractmethod
from typing import Any, Dict, Optional, Tuple

import pysam


class BaseVCFParser(ABC):
    """Abstract base class for VCF parsers.

    Defines the interface for extracting genotype and variant information
    from VCF records in a standardized way. All parser implementations must
    inherit from this class and implement all abstract methods.
    """

    @abstractmethod
    def get_genotype(
        self, record: pysam.VariantRecord, sample_idx: int
    ) -> Optional[Tuple[int, int]]:
        """Extract genotype as (allele1, allele2) or None if unknown.

        Parameters
        ----------
        record : pysam.VariantRecord
            The VCF record to extract genotype from
        sample_idx : int
            Index of the sample in the record

        Returns
        -------
        Optional[Tuple[int, int]]
            Genotype as tuple of allele indices (0=REF, 1=ALT), or None if unknown
        """
        pass

    @abstractmethod
    def has_variant(self, record: pysam.VariantRecord, sample_idx: int) -> bool:
        """Check if sample has variant even when GT is unknown.

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
        pass

    @abstractmethod
    def extract_info(self, record: pysam.VariantRecord, sample_idx: int) -> Dict[str, Any]:
        """Extract additional fields (AD, DP, etc.) as dictionary.

        Parameters
        ----------
        record : pysam.VariantRecord
            The VCF record to extract information from
        sample_idx : int
            Index of the sample in the record

        Returns
        -------
        Dict[str, Any]
            Dictionary of additional FORMAT fields
        """
        pass

    @abstractmethod
    def validate_record(self, record: pysam.VariantRecord) -> bool:
        """Validate that record is compatible with this parser.

        Parameters
        ----------
        record : pysam.VariantRecord
            The VCF record to validate

        Returns
        -------
        bool
            True if record is valid for this parser, False otherwise
        """
        pass
