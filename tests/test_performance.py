"""Performance validation tests."""

import tempfile
import time
from pathlib import Path

import pysam
import pytest

from strvcf_annotator import STRAnnotator
from strvcf_annotator.core.str_reference import is_valid_tabix


class TestPerformance:
    """Performance validation tests."""

    @pytest.fixture
    def large_str_bed(self):
        """Create a larger STR BED file for performance testing."""
        content = []
        for i in range(1000):
            chrom = f"chr{(i % 22) + 1}"
            start = i * 1000
            end = start + 15
            content.append(f"{chrom}\t{start}\t{end}\t3\tCAG")

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write("\n".join(content))
            temp_path = f.name

        yield temp_path
        Path(temp_path).unlink()

    @pytest.fixture
    def large_vcf(self):
        """Create a larger VCF file for performance testing."""
        vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depth">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2
"""
        # Add 100 variant records
        for i in range(100):
            pos = i * 1000 + 5
            vcf_content += f"chr1\t{pos}\t.\tA\tT\t.\t.\t.\tGT:AD:DP\t0/1:10,5:15\t1/1:0,20:20\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write(vcf_content)
            temp_path = f.name

        yield temp_path
        Path(temp_path).unlink()

    def test_loading_performance(self, large_str_bed):
        """Test STR reference loading performance."""
        start_time = time.time()
        annotator = STRAnnotator(large_str_bed)
        load_time = time.time() - start_time

        # Should load 1000 regions in under 1 second
        assert load_time < 1.0, f"Loading took {load_time:.2f}s, expected < 1.0s"
        assert Path(annotator.str_panel_gz).exists() and is_valid_tabix(annotator.str_panel_gz), (
            "STR BED file should be cached and indexed"
        )

    def test_annotation_performance(self, large_str_bed, large_vcf):
        """Test annotation performance."""
        annotator = STRAnnotator(large_str_bed)

        with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as output_file:
            output_path = output_file.name

        try:
            start_time = time.time()
            annotator.annotate_vcf_file(large_vcf, output_path)
            annotation_time = time.time() - start_time

            # Should annotate 100 records in under 5 seconds
            assert annotation_time < 5.0, f"Annotation took {annotation_time:.2f}s, expected < 5.0s"

            # Verify output
            vcf_out = pysam.VariantFile(output_path)
            records = list(vcf_out)
            assert len(records) > 0
            vcf_out.close()

        finally:
            if Path(output_path).exists():
                Path(output_path).unlink()

    def test_streaming_memory_efficiency(self, large_str_bed, large_vcf):
        """Test that streaming doesn't load entire file into memory."""
        annotator = STRAnnotator(large_str_bed)
        vcf_in = pysam.VariantFile(large_vcf)

        # Process records one at a time (streaming)
        record_count = 0
        for record in annotator.annotate_vcf_stream(vcf_in):
            record_count += 1
            # Verify record has annotations
            assert "RU" in record.info

        assert record_count > 0
        vcf_in.close()

    def test_statistics_performance(self, large_str_bed):
        """Test statistics calculation performance."""
        annotator = STRAnnotator(large_str_bed)

        start_time = time.time()
        stats = annotator.get_statistics()
        stats_time = time.time() - start_time

        # Should calculate stats in under 0.5 seconds
        assert stats_time < 0.5, f"Stats calculation took {stats_time:.2f}s, expected < 0.5s"
        assert stats["total_regions"] == 1000
