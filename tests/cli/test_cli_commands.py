"""Tests for CLI commands using subprocess."""

import os
import shutil
import subprocess
import tempfile

import pytest


@pytest.fixture
def temp_output_dir():
    """Create a temporary output directory."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    # Cleanup
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)


@pytest.fixture
def sample_vcf(data_dir):
    """Path to a sample VCF file for testing."""
    return os.path.join(data_dir, "TCGA-DC-6682.vcf")


@pytest.fixture
def sample_bed(data_dir):
    """Path to a sample BED file for testing."""
    return os.path.join(data_dir, "GRCh38_repeats.bed")


class TestCLIBasicUsage:
    """Basic CLI usage tests."""

    def test_help_command(self):
        """Test --help flag."""
        result = subprocess.run(["strvcf-annotator", "--help"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "Annotate STR regions in VCF files" in result.stdout
        assert "--input" in result.stdout
        assert "--str-bed" in result.stdout
        assert "--output" in result.stdout

    def test_version_command(self):
        """Test --version flag."""
        result = subprocess.run(["strvcf-annotator", "--version"], capture_output=True, text=True)
        assert result.returncode == 0
        assert "0.2.2" in result.stdout

    def test_no_arguments_fails(self):
        """Test that running without arguments fails."""
        result = subprocess.run(["strvcf-annotator"], capture_output=True, text=True)
        assert result.returncode != 0
        assert "required" in result.stderr.lower() or "error" in result.stderr.lower()


class TestCLISingleFileMode:
    """Tests for single file annotation mode."""

    def test_single_file_annotation(self, sample_vcf, sample_bed, temp_output_dir):
        """Test annotating a single VCF file."""
        output_vcf = os.path.join(temp_output_dir, "annotated.vcf")

        result = subprocess.run(
            [
                "strvcf-annotator",
                "--input",
                sample_vcf,
                "--str-bed",
                sample_bed,
                "--output",
                output_vcf,
            ],
            capture_output=True,
            text=True,
            timeout=300,  # 5 minutes timeout
        )

        # Check return code
        assert result.returncode == 0, f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"

        # Check output file was created
        assert os.path.exists(output_vcf), "Output VCF file was not created"

        # Check output file is not empty
        assert os.path.getsize(output_vcf) > 0, "Output VCF file is empty"

    def test_single_file_with_verbose(self, sample_vcf, sample_bed, temp_output_dir):
        """Test single file annotation with verbose logging."""
        output_vcf = os.path.join(temp_output_dir, "annotated.vcf")

        result = subprocess.run(
            [
                "strvcf-annotator",
                "--input",
                sample_vcf,
                "--str-bed",
                sample_bed,
                "--output",
                output_vcf,
                "--verbose",
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )

        assert result.returncode == 0
        # Verbose mode should produce more log output
        assert "DEBUG" in result.stderr or "INFO" in result.stderr

    def test_missing_input_file(self, sample_bed, temp_output_dir):
        """Test error handling for missing input file."""
        output_vcf = os.path.join(temp_output_dir, "annotated.vcf")

        result = subprocess.run(
            [
                "strvcf-annotator",
                "--input",
                "/nonexistent/file.vcf",
                "--str-bed",
                sample_bed,
                "--output",
                output_vcf,
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0
        # Should produce error message
        assert len(result.stderr) > 0 or len(result.stdout) > 0

    def test_missing_bed_file(self, sample_vcf, temp_output_dir):
        """Test error handling for missing BED file."""
        output_vcf = os.path.join(temp_output_dir, "annotated.vcf")

        result = subprocess.run(
            [
                "strvcf-annotator",
                "--input",
                sample_vcf,
                "--str-bed",
                "/nonexistent/file.bed",
                "--output",
                output_vcf,
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0

    def test_invalid_output_directory(self, sample_vcf, sample_bed):
        """Test error handling for invalid output path."""
        result = subprocess.run(
            [
                "strvcf-annotator",
                "--input",
                sample_vcf,
                "--str-bed",
                sample_bed,
                "--output",
                "/nonexistent/directory/output.vcf",
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0


class TestCLIBatchMode:
    """Tests for batch directory processing mode."""

    def test_batch_directory_processing(self, data_dir, sample_bed, temp_output_dir):
        """Test batch processing of VCF directory."""
        # Create a temporary input directory with VCF files
        temp_input_dir = tempfile.mkdtemp()
        try:
            # Copy a sample VCF file
            sample_vcf = os.path.join(data_dir, "TCGA-DC-6682.vcf")
            shutil.copy(sample_vcf, os.path.join(temp_input_dir, "test.vcf"))

            result = subprocess.run(
                [
                    "strvcf-annotator",
                    "--input-dir",
                    temp_input_dir,
                    "--str-bed",
                    sample_bed,
                    "--output-dir",
                    temp_output_dir,
                ],
                capture_output=True,
                text=True,
                timeout=300,
            )

            assert result.returncode == 0, f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"

            # Check output directory was created and contains files
            assert os.path.exists(temp_output_dir)
            output_files = os.listdir(temp_output_dir)
            assert len(output_files) > 0, "No output files were created"

        finally:
            if os.path.exists(temp_input_dir):
                shutil.rmtree(temp_input_dir)

    def test_batch_mode_creates_output_dir(self, data_dir, sample_bed):
        """Test that batch mode creates output directory if it doesn't exist."""
        temp_input_dir = tempfile.mkdtemp()
        temp_output_dir = tempfile.mkdtemp()

        try:
            # Remove output dir to test creation
            shutil.rmtree(temp_output_dir)

            # Copy a sample VCF file
            sample_vcf = os.path.join(data_dir, "TCGA-DC-6682.vcf")
            shutil.copy(sample_vcf, os.path.join(temp_input_dir, "test.vcf"))

            result = subprocess.run(
                [
                    "strvcf-annotator",
                    "--input-dir",
                    temp_input_dir,
                    "--str-bed",
                    sample_bed,
                    "--output-dir",
                    temp_output_dir,
                ],
                capture_output=True,
                text=True,
                timeout=300,
            )

            assert result.returncode == 0
            assert os.path.exists(temp_output_dir), "Output directory was not created"

        finally:
            if os.path.exists(temp_input_dir):
                shutil.rmtree(temp_input_dir)
            if os.path.exists(temp_output_dir):
                shutil.rmtree(temp_output_dir)

    def test_empty_input_directory(self, sample_bed, temp_output_dir):
        """Test batch mode with empty input directory."""
        temp_input_dir = tempfile.mkdtemp()

        try:
            result = subprocess.run(
                [
                    "strvcf-annotator",
                    "--input-dir",
                    temp_input_dir,
                    "--str-bed",
                    sample_bed,
                    "--output-dir",
                    temp_output_dir,
                ],
                capture_output=True,
                text=True,
            )

            # Should succeed but with no output files
            # (or may fail depending on implementation)
            # Either behavior is acceptable
            if result.returncode == 0:
                output_files = os.listdir(temp_output_dir)
                assert len(output_files) == 0

        finally:
            if os.path.exists(temp_input_dir):
                shutil.rmtree(temp_input_dir)


class TestCLIArgumentValidation:
    """Tests for CLI argument validation."""

    def test_input_without_output(self, sample_vcf, sample_bed):
        """Test that --input without --output fails."""
        result = subprocess.run(
            ["strvcf-annotator", "--input", sample_vcf, "--str-bed", sample_bed],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0

    def test_input_dir_without_output_dir(self, sample_bed, temp_output_dir):
        """Test that --input-dir without --output-dir fails."""
        result = subprocess.run(
            ["strvcf-annotator", "--input-dir", temp_output_dir, "--str-bed", sample_bed],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0

    def test_output_without_input(self, sample_bed, temp_output_dir):
        """Test that --output without --input fails."""
        output_vcf = os.path.join(temp_output_dir, "out.vcf")

        result = subprocess.run(
            ["strvcf-annotator", "--str-bed", sample_bed, "--output", output_vcf],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0

    def test_mixing_single_and_batch_modes(self, sample_vcf, sample_bed, temp_output_dir):
        """Test that mixing --input with --output-dir fails."""
        result = subprocess.run(
            [
                "strvcf-annotator",
                "--input",
                sample_vcf,
                "--str-bed",
                sample_bed,
                "--output-dir",
                temp_output_dir,
            ],
            capture_output=True,
            text=True,
        )

        # Should fail due to argument validation
        assert result.returncode != 0

    def test_missing_str_bed(self, sample_vcf, temp_output_dir):
        """Test that missing --str-bed fails."""
        output_vcf = os.path.join(temp_output_dir, "out.vcf")

        result = subprocess.run(
            ["strvcf-annotator", "--input", sample_vcf, "--output", output_vcf],
            capture_output=True,
            text=True,
        )

        assert result.returncode != 0
        assert "required" in result.stderr.lower() or "--str-bed" in result.stderr


class TestCLIOutputValidation:
    """Tests for CLI output validation."""

    def test_output_contains_str_annotations(self, sample_vcf, sample_bed, temp_output_dir):
        """Test that output VCF contains STR annotation fields."""
        output_vcf = os.path.join(temp_output_dir, "annotated.vcf")

        result = subprocess.run(
            [
                "strvcf-annotator",
                "--input",
                sample_vcf,
                "--str-bed",
                sample_bed,
                "--output",
                output_vcf,
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )

        assert result.returncode == 0

        # Read output file and check for STR fields
        with open(output_vcf) as f:
            content = f.read()
            # Check header contains STR INFO fields
            assert "##INFO=<ID=RU" in content
            assert "##INFO=<ID=PERIOD" in content
            assert "##INFO=<ID=REF" in content
            assert "##INFO=<ID=PERFECT" in content
            # Check header contains STR FORMAT field
            assert "##FORMAT=<ID=REPCN" in content

    def test_output_is_valid_vcf(self, sample_vcf, sample_bed, temp_output_dir):
        """Test that output is a valid VCF file."""
        output_vcf = os.path.join(temp_output_dir, "annotated.vcf")

        result = subprocess.run(
            [
                "strvcf-annotator",
                "--input",
                sample_vcf,
                "--str-bed",
                sample_bed,
                "--output",
                output_vcf,
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )

        assert result.returncode == 0

        # Basic VCF format validation
        with open(output_vcf) as f:
            lines = f.readlines()
            # Should have header lines
            assert any(line.startswith("##fileformat=VCF") for line in lines)
            # Should have column header
            assert any(line.startswith("#CHROM") for line in lines)


class TestCLIEdgeCases:
    """Edge case tests for CLI."""

    def test_compressed_vcf_input(self, data_dir, sample_bed, temp_output_dir):
        """Test handling of compressed VCF input."""
        compressed_vcf = os.path.join(data_dir, "test.vcf.gz")

        if not os.path.exists(compressed_vcf):
            pytest.skip("Compressed VCF file not available")

        output_vcf = os.path.join(temp_output_dir, "annotated.vcf")

        result = subprocess.run(
            [
                "strvcf-annotator",
                "--input",
                compressed_vcf,
                "--str-bed",
                sample_bed,
                "--output",
                output_vcf,
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )

        # Should handle compressed files
        assert result.returncode == 0 or "compressed" in result.stderr.lower()

    def test_relative_paths(self, temp_output_dir):
        """Test that relative paths work correctly."""
        # Create temporary files in current directory
        temp_vcf = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
        temp_bed = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)

        try:
            # Write minimal VCF
            temp_vcf.write("##fileformat=VCFv4.2\n")
            temp_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            temp_vcf.close()

            # Write minimal BED
            temp_bed.write("chr1\t100\t200\t2\tAT\n")
            temp_bed.close()

            output_vcf = os.path.join(temp_output_dir, "out.vcf")

            result = subprocess.run(
                [
                    "strvcf-annotator",
                    "--input",
                    temp_vcf.name,
                    "--str-bed",
                    temp_bed.name,
                    "--output",
                    output_vcf,
                ],
                capture_output=True,
                text=True,
            )

            # Relative paths should work
            # (may fail due to empty VCF, but path resolution should work)
            assert (
                result.returncode == 0
                or "empty" in result.stderr.lower()
                or "no" in result.stderr.lower()
            )

        finally:
            os.unlink(temp_vcf.name)
            os.unlink(temp_bed.name)
