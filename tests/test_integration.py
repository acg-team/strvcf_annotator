"""Integration tests for end-to-end workflows."""

import tempfile
from pathlib import Path

import pysam
import pytest

from strvcf_annotator import STRAnnotator, annotate_vcf
from strvcf_annotator.parsers.generic import GenericParser


class TestEndToEndAnnotation:
    """Test complete annotation pipeline."""

    @pytest.fixture
    def str_bed_file(self):
        """Create temporary STR BED file."""
        content = """chr1\t100\t115\t3\tCAG
chr1\t200\t212\t4\tATCG
chr2\t300\t318\t3\tGAT"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(content)
            temp_path = f.name

        yield temp_path
        Path(temp_path).unlink()

    @pytest.fixture
    def vcf_file(self):
        """Create temporary VCF file."""
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr2>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2
chr1\t105\t.\tA\tT\t.\t.\t.\tGT:AD:DP\t0/1:10,5:15\t1/1:0,20:20
chr1\t205\t.\tC\tG\t.\t.\t.\tGT:AD:DP\t0/0:15,0:15\t0/1:8,7:15
chr2\t305\t.\tG\tA\t.\t.\t.\tGT:AD:DP\t0/1:12,8:20\t1/1:0,18:18
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write(vcf_content)
            temp_path = f.name

        yield temp_path
        Path(temp_path).unlink()

    def test_annotate_single_file(self, str_bed_file, vcf_file):
        """Test annotating a single VCF file."""
        with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as output_file:
            output_path = output_file.name

        try:
            annotator = STRAnnotator(str_bed_file)
            annotator.annotate_vcf_file(vcf_file, output_path)

            # Verify output file exists and is valid
            assert Path(output_path).exists()

            # Read and verify output
            vcf_out = pysam.VariantFile(output_path)
            records = list(vcf_out)

            # Should have annotated records
            assert len(records) > 0

            # Check that STR-specific fields are present
            for record in records:
                assert "RU" in record.info
                assert "PERIOD" in record.info
                assert "REF" in record.info
                assert "PERFECT" in record.info

                # Check FORMAT fields
                for sample in record.samples.values():
                    assert "REPCN" in sample

            vcf_out.close()

        finally:
            if Path(output_path).exists():
                Path(output_path).unlink()

    def test_convenience_function(self, str_bed_file, vcf_file):
        """Test convenience function for annotation."""
        with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as output_file:
            output_path = output_file.name

        try:
            annotate_vcf(vcf_file, str_bed_file, output_path)

            # Verify output
            assert Path(output_path).exists()
            vcf_out = pysam.VariantFile(output_path)
            records = list(vcf_out)
            assert len(records) > 0
            vcf_out.close()

        finally:
            if Path(output_path).exists():
                Path(output_path).unlink()

    def test_stream_processing(self, str_bed_file, vcf_file):
        """Test stream processing of VCF records."""
        annotator = STRAnnotator(str_bed_file)
        vcf_in = pysam.VariantFile(vcf_file)

        records = list(annotator.annotate_vcf_stream(vcf_in))

        # Should have some annotated records
        assert len(records) > 0

        # Verify annotations
        for record in records:
            assert "RU" in record.info
            assert "REPCN" in record.samples[list(record.samples.keys())[0]]

        vcf_in.close()

    def test_batch_directory_processing(self, str_bed_file, vcf_file):
        """Test batch processing of directory."""
        with tempfile.TemporaryDirectory() as input_dir:
            with tempfile.TemporaryDirectory() as output_dir:
                # Copy VCF to input directory
                input_vcf = Path(input_dir) / "test.vcf"
                input_vcf.write_text(Path(vcf_file).read_text())

                # Process directory
                annotator = STRAnnotator(str_bed_file)
                annotator.process_directory(
                    input_dir, output_dir, ignore_mismatch_warnings=True, mismatch_truth="vcf"
                )

                # Verify output
                output_files = list(Path(output_dir).glob("*.vcf"))
                assert len(output_files) > 0

                # Verify content
                vcf_out = pysam.VariantFile(str(output_files[0]))
                records = list(vcf_out)
                assert len(records) > 0
                vcf_out.close()

    def test_get_statistics(self, str_bed_file):
        """Test getting STR statistics."""
        annotator = STRAnnotator(str_bed_file)
        stats = annotator.get_statistics()

        assert "total_regions" in stats
        assert stats["total_regions"] == 3
        assert "chromosomes" in stats
        assert stats["chromosomes"] == 2
        assert "unique_repeat_units" in stats


class TestErrorHandling:
    """Test error handling in integration scenarios."""

    def test_invalid_vcf_file(self):
        """Test handling of invalid VCF file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as bed_file:
            bed_file.write("chr1\t100\t115\t3\tCAG\n")
            bed_path = bed_file.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as vcf_file:
            vcf_file.write("This is not a VCF file")
            vcf_path = vcf_file.name

        try:
            annotator = STRAnnotator(bed_path)

            with pytest.raises(Exception):
                annotator.annotate_vcf_file(vcf_path, "output.vcf")

        finally:
            Path(bed_path).unlink()
            Path(vcf_path).unlink()

    def test_invalid_bed_file(self):
        """Test handling of invalid BED file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as bed_file:
            bed_file.write("Invalid BED content")
            bed_path = bed_file.name

        try:
            with pytest.raises(Exception):
                STRAnnotator(bed_path)

        finally:
            Path(bed_path).unlink()

    def test_nonexistent_files(self):
        """Test handling of nonexistent files."""
        with pytest.raises(Exception):
            STRAnnotator("/nonexistent/file.bed")


class TestCustomParser:
    """Test using custom parser."""

    @pytest.fixture
    def str_bed_file(self):
        """Create temporary STR BED file."""
        content = "chr1\t100\t115\t3\tCAG\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(content)
            temp_path = f.name

        yield temp_path
        Path(temp_path).unlink()

    def test_with_generic_parser(self, str_bed_file):
        """Test annotation with explicit GenericParser."""
        parser = GenericParser()
        annotator = STRAnnotator(str_bed_file, parser=parser)

        stats = annotator.get_statistics()
        assert stats["total_regions"] == 1
