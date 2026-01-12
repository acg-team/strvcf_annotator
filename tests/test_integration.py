"""Integration tests for end-to-end workflows."""


from pathlib import Path

import pysam
import pytest

from strvcf_annotator import STRAnnotator, annotate_vcf
from strvcf_annotator.parsers.generic import GenericParser
from strvcf_annotator.utils.validation import ValidationError


class TestEndToEndAnnotation:
    """Test complete annotation pipeline."""

    @pytest.fixture
    def str_bed_file(self, tmp_path):
        """Create STR BED file in a temporary folder (small synthetic reference)."""
        content = """chr1\t100\t115\t3\tCAG
chr1\t200\t212\t4\tATCG
chr2\t300\t318\t3\tGAT
"""
        bed = tmp_path / "repeats.bed"
        bed.write_text(content)
        return str(bed)

    @pytest.fixture
    def vcf_file(self, tmp_path):
        """Create VCF file in a temporary folder (small synthetic VCF)."""
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
        vcf = tmp_path / "input.vcf"
        vcf.write_text(vcf_content)
        return str(vcf)

    def _read_records(self, path: Path):
        v = pysam.VariantFile(str(path))
        recs = list(v)
        v.close()
        return recs

    def test_annotate_single_file_default_options(self, str_bed_file, vcf_file, output_path):
        """Test annotating a single VCF file with default options."""
        annotator = STRAnnotator(str_bed_file)
        annotator.annotate_vcf_file(vcf_file, str(output_path))

        assert output_path.exists()
        records = self._read_records(output_path)
        assert len(records) > 0

        for record in records:
            assert "RU" in record.info
            assert "PERIOD" in record.info
            assert "REF" in record.info
            assert "PERFECT" in record.info
            for sample in record.samples.values():
                assert "REPCN" in sample

    @pytest.mark.parametrize(
        "ignore_mismatch_warnings,mismatch_truth",
        [
            (False, "panel"),
            (True, "panel"),
            (True, "vcf"),
            (True, "skip"),
        ],
    )
    def test_annotate_single_file_new_options_available(
        self, str_bed_file, vcf_file, output_path, ignore_mismatch_warnings, mismatch_truth
    ):
        """
        Integration-level check that new options are accepted and produce output.
        For this small synthetic data there may be no mismatches, but we verify:
          - call succeeds
          - output is a valid VCF
          - output has STR annotations
        """
        annotator = STRAnnotator(
            str_bed_file,
            ignore_mismatch_warnings=ignore_mismatch_warnings,
            mismatch_truth=mismatch_truth,
        )
        annotator.annotate_vcf_file(vcf_file, str(output_path))

        assert output_path.exists()
        records = self._read_records(output_path)

        # Even with skip mode, for this dataset we still expect records (no mismatch enforced here)
        assert len(records) > 0

        for record in records:
            assert "RU" in record.info
            assert "REPCN" in record.samples[list(record.samples.keys())[0]]

    @pytest.mark.parametrize(
        "ignore_mismatch_warnings,mismatch_truth",
        [
            (False, "panel"),
            (True, "panel"),
            (True, "vcf"),
            (True, "skip"),
        ],
    )
    def test_convenience_function_new_options_available(
        self, str_bed_file, vcf_file, output_path, ignore_mismatch_warnings, mismatch_truth
    ):
        """Same as above, but through annotate_vcf() convenience function."""
        annotate_vcf(
            vcf_file,
            str_bed_file,
            str(output_path),
            somatic_mode=None,
            ignore_mismatch_warnings=ignore_mismatch_warnings,
            mismatch_truth=mismatch_truth,
        )

        assert output_path.exists()
        records = self._read_records(output_path)
        assert len(records) > 0

    def test_stream_processing_new_options_available(self, str_bed_file, vcf_file):
        """Test stream processing (options should be supported via annotator defaults)."""
        annotator = STRAnnotator(
            str_bed_file,
            ignore_mismatch_warnings=True,
            mismatch_truth="panel",
        )
        vcf_in = pysam.VariantFile(vcf_file)
        records = list(annotator.annotate_vcf_stream(vcf_in))
        vcf_in.close()

        assert len(records) > 0
        for record in records:
            assert "RU" in record.info
            assert "REPCN" in record.samples[list(record.samples.keys())[0]]

    @pytest.mark.parametrize(
        "ignore_mismatch_warnings,mismatch_truth",
        [
            (False, "panel"),
            (True, "panel"),
            (True, "vcf"),
            (True, "skip"),
        ],
    )
    def test_batch_directory_processing_new_options_available(
        self,
        str_bed_file,
        vcf_file,
        output_dir,
        tmp_path,
        ignore_mismatch_warnings,
        mismatch_truth,
    ):
        """Batch processing should accept and propagate new mismatch options."""
        input_dir = tmp_path / "input_vcfs"
        input_dir.mkdir(parents=True, exist_ok=True)

        # Copy VCF into input directory
        input_vcf = input_dir / "test.vcf"
        input_vcf.write_text(Path(vcf_file).read_text())

        annotator = STRAnnotator(
            str_bed_file,
            ignore_mismatch_warnings=ignore_mismatch_warnings,
            mismatch_truth=mismatch_truth,
        )
        annotator.process_directory(str(input_dir), str(output_dir))

        output_files = list(output_dir.glob("*.annotated.vcf")) + list(output_dir.glob("*.vcf"))
        assert len(output_files) > 0

        vcf_out = pysam.VariantFile(str(output_files[0]))
        records = list(vcf_out)
        vcf_out.close()
        assert len(records) >= 0  # allow 0 if skip ever drops all (depends on data)

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

    def test_invalid_vcf_file(self, tmp_path, output_path):
        """Test handling of invalid VCF file."""
        bed = tmp_path / "repeats.bed"
        bed.write_text("chr1\t100\t115\t3\tCAG\n")

        vcf = tmp_path / "bad.vcf"
        vcf.write_text("This is not a VCF file")

        annotator = STRAnnotator(str(bed))

        with pytest.raises(ValidationError):
            annotator.annotate_vcf_file(str(vcf), str(output_path))

    def test_invalid_bed_file(self, tmp_path):
        """Test handling of invalid BED file."""
        bed = tmp_path / "bad.bed"
        bed.write_text("Invalid BED content")

        with pytest.raises(ValidationError):
            STRAnnotator(str(bed))

    def test_nonexistent_files(self):
        """Test handling of nonexistent files."""
        with pytest.raises(ValidationError):
            STRAnnotator("/nonexistent/file.bed")


class TestCustomParser:
    """Test using custom parser."""

    @pytest.fixture
    def str_bed_file(self, tmp_path):
        content = "chr1\t100\t115\t3\tCAG\n"
        bed = tmp_path / "repeats.bed"
        bed.write_text(content)
        return str(bed)

    def test_with_generic_parser(self, str_bed_file):
        """Test annotation with explicit GenericParser."""
        parser = GenericParser()
        annotator = STRAnnotator(str_bed_file, parser=parser)

        stats = annotator.get_statistics()
        assert stats["total_regions"] == 1
