"""Input validation functions."""

import contextlib
from pathlib import Path

import pandas as pd
import pysam


class ValidationError(Exception):
    """Exception raised for validation errors."""

    pass


def validate_file_path(file_path: str, must_exist: bool = True) -> Path:
    """Validate file path.

    Parameters
    ----------
    file_path : str
        Path to validate
    must_exist : bool, optional
        If True, file must exist. Default is True.

    Returns
    -------
    Path
        Validated Path object

    Raises
    ------
    ValidationError
        If file path is invalid or doesn't exist when required
    """
    if not file_path:
        raise ValidationError("File path cannot be empty")

    path = Path(file_path)

    if must_exist and not path.exists():
        raise ValidationError(f"File does not exist: {file_path}")

    if must_exist and not path.is_file():
        raise ValidationError(f"Path is not a file: {file_path}")

    return path


def validate_directory_path(dir_path: str, must_exist: bool = True, create: bool = False) -> Path:
    """Validate directory path.

    Parameters
    ----------
    dir_path : str
        Directory path to validate
    must_exist : bool, optional
        If True, directory must exist. Default is True.
    create : bool, optional
        If True, create directory if it doesn't exist. Default is False.

    Returns
    -------
    Path
        Validated Path object

    Raises
    ------
    ValidationError
        If directory path is invalid
    """
    if not dir_path:
        raise ValidationError("Directory path cannot be empty")

    path = Path(dir_path)

    if create and not path.exists():
        path.mkdir(parents=True, exist_ok=True)

    if must_exist and not path.exists():
        raise ValidationError(f"Directory does not exist: {dir_path}")

    if path.exists() and not path.is_dir():
        raise ValidationError(f"Path is not a directory: {dir_path}")

    return path


def validate_vcf_file(vcf_path: str) -> bool:
    """Validate VCF file format.

    Parameters
    ----------
    vcf_path : str
        Path to VCF file

    Returns
    -------
    bool
        True if VCF is valid

    Raises
    ------
    ValidationError
        If VCF file is invalid or cannot be opened
    """
    path = validate_file_path(vcf_path, must_exist=True)

    try:
        vcf = pysam.VariantFile(str(path))

        # Check that VCF has samples
        if not vcf.header.samples:
            raise ValidationError(f"VCF file has no samples: {vcf_path}")

        # Try to read first record

        # Empty VCF is technically valid
        with contextlib.suppress(StopIteration):
            next(vcf)

        vcf.close()
        return True

    except Exception as e:
        raise ValidationError(f"Invalid VCF file {vcf_path}: {str(e)}") from e


def validate_bed_file(bed_path: str) -> bool:
    """Validate BED file format.

    Parameters
    ----------
    bed_path : str
        Path to BED file

    Returns
    -------
    bool
        True if BED is valid

    Raises
    ------
    ValidationError
        If BED file is invalid or cannot be opened
    """
    path = validate_file_path(bed_path, must_exist=True)

    try:
        df = pd.read_csv(str(path), sep="\t", header=None, nrows=1)

        # Check minimum number of columns (at least 3 for BED)
        if df.shape[1] < 3:
            raise ValidationError(f"BED file has fewer than 3 columns: {bed_path}")

        # Check that first 3 columns are chrom, start, end
        # Start and end should be numeric
        if not pd.api.types.is_numeric_dtype(df.iloc[:, 1]):
            raise ValidationError(f"BED file column 2 (start) is not numeric: {bed_path}")

        if not pd.api.types.is_numeric_dtype(df.iloc[:, 2]):
            raise ValidationError(f"BED file column 3 (end) is not numeric: {bed_path}")

        return True

    except pd.errors.EmptyDataError as e:
        raise ValidationError(f"BED file is empty: {bed_path}") from e
    except Exception as e:
        raise ValidationError(f"Invalid BED file {bed_path}: {str(e)}") from e


def validate_str_bed_file(bed_path: str) -> bool:
    """Validate STR BED file format with required columns.

    Parameters
    ----------
    bed_path : str
        Path to STR BED file

    Returns
    -------
    bool
        True if STR BED is valid

    Raises
    ------
    ValidationError
        If STR BED file is invalid or missing required columns
    """
    path = validate_file_path(bed_path, must_exist=True)

    try:
        df = pd.read_csv(str(path), sep="\t", header=None, nrows=1)

        # Check for required 5 columns: CHROM, START, END, PERIOD, RU
        if df.shape[1] < 5:
            raise ValidationError(
                f"STR BED file must have at least 5 columns "
                f"(CHROM, START, END, PERIOD, RU): {bed_path}"
            )

        # Validate column types
        if not pd.api.types.is_numeric_dtype(df.iloc[:, 1]):
            raise ValidationError(f"START column (2) must be numeric: {bed_path}")

        if not pd.api.types.is_numeric_dtype(df.iloc[:, 2]):
            raise ValidationError(f"END column (3) must be numeric: {bed_path}")

        if not pd.api.types.is_numeric_dtype(df.iloc[:, 3]):
            raise ValidationError(f"PERIOD column (4) must be numeric: {bed_path}")

        return True

    except pd.errors.EmptyDataError as e:
        raise ValidationError(f"STR BED file is empty: {bed_path}") from e
    except Exception as e:
        raise ValidationError(f"Invalid STR BED file {bed_path}: {str(e)}") from e
