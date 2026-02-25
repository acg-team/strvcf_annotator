"""Unit tests for VCF utility functions."""

import os
import tempfile

import pysam
import pytest

from strvcf_annotator.utils.vcf_utils import (
    chrom_to_order,
    get_sample_by_index,
    get_sample_by_name,
    has_format_field,
    normalize_info_fields,
)


@pytest.fixture
def basic_header():
    """Create a basic VCF header for tests."""
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")
    header.contigs.add("chr1", length=1000000)

    header.add_sample("TUMOR")
    header.add_sample("NORMAL")

    header.info.add("FLAG1", 0, "Flag", "Flag field")
    header.info.add("S1", 1, "String", "Single string")
    header.info.add("I1", 1, "Integer", "Single integer")
    header.info.add("R1", "R", "Integer", "REF+ALT values")

    header.formats.add("GT", 1, "String", "Genotype")
    header.formats.add("AD", "R", "Integer", "Allelic depths")

    return header


@pytest.fixture
def basic_record(basic_header):
    """Create a pysam VariantRecord for tests."""
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
    tmp.close()

    try:
        with pysam.VariantFile(tmp.name, "w", header=basic_header) as vcf_out:
            rec = vcf_out.new_record(
                contig="chr1",
                start=99,
                alleles=("A", "T"),
                qual=30,
                filter=["PASS"],
            )

            # INFO: include various types
            rec.info["FLAG1"] = True
            rec.info["S1"] = ("a", "b")  # should become "a|b"
            rec.info["I1"] = (10,)  # should become 10
            rec.info["R1"] = (1, 2, 3, 4)  # should be clipped to [1,2]
            rec.info["UNKNOWN"] = "skip_me"  # not in header -> skipped

            # FORMAT fields
            rec.samples["TUMOR"]["GT"] = (0, 1)
            rec.samples["TUMOR"]["AD"] = (10, 5)
            rec.samples["NORMAL"]["GT"] = (0, 0)
            rec.samples["NORMAL"]["AD"] = (12, 0)

            vcf_out.write(rec)

        with pysam.VariantFile(tmp.name) as vcf_in:
            records = list(vcf_in)
            assert len(records) == 1
            return records[0]
    finally:
        os.unlink(tmp.name)


class TestChromToOrder:
    """Tests for chrom_to_order."""

    def test_none(self):
        assert chrom_to_order(None) == 1_000_000

    def test_autosomes_with_chr_prefix(self):
        assert chrom_to_order("chr1") == 1
        assert chrom_to_order("chr2") == 2
        assert chrom_to_order("chr10") == 10

    def test_autosomes_without_chr_prefix(self):
        assert chrom_to_order("1") == 1
        assert chrom_to_order("22") == 22

    def test_sex_chromosomes(self):
        assert chrom_to_order("chrX") == 23
        assert chrom_to_order("X") == 23
        assert chrom_to_order("chrY") == 24
        assert chrom_to_order("Y") == 24

    def test_mitochondrial(self):
        assert chrom_to_order("chrM") == 25
        assert chrom_to_order("M") == 25
        assert chrom_to_order("chrMT") == 25
        assert chrom_to_order("MT") == 25

    def test_other_contigs_go_last(self):
        assert chrom_to_order("chrUn_gl000220") == 1_000_000
        assert chrom_to_order("GL000220.1") == 1_000_000


class TestNormalizeInfoFields:
    """Tests for normalize_info_fields."""

    def test_skips_unknown_info_fields(self, basic_record, basic_header):
        fixed = normalize_info_fields(basic_record, basic_header)
        assert "UNKNOWN" not in fixed

    def test_flag_field_included_only_if_true(self, basic_record, basic_header):
        fixed = normalize_info_fields(basic_record, basic_header)
        assert fixed["FLAG1"] is True

    def test_string_number_one_tuple_joined(self, basic_record, basic_header):
        fixed = normalize_info_fields(basic_record, basic_header)
        assert fixed["S1"] == "a|b"

    def test_scalar_number_one_tuple_clipped(self, basic_record, basic_header):
        fixed = normalize_info_fields(basic_record, basic_header)
        assert fixed["I1"] == 10

    def test_r_field_is_clipped_to_two(self, basic_record, basic_header):
        fixed = normalize_info_fields(basic_record, basic_header)
        assert fixed["R1"] == [1, 2]


class TestGetSampleByName:
    """Tests for get_sample_by_name."""

    def test_returns_sample(self, basic_record):
        tumor = get_sample_by_name(basic_record, "TUMOR")
        assert tumor is not None
        assert tumor["GT"] == (0, 1)

    def test_raises_keyerror(self, basic_record):
        with pytest.raises(KeyError):
            get_sample_by_name(basic_record, "NOT_A_SAMPLE")


class TestGetSampleByIndex:
    """Tests for get_sample_by_index."""

    def test_returns_by_index(self, basic_record):
        s0 = get_sample_by_index(basic_record, 0)
        s1 = get_sample_by_index(basic_record, 1)

        assert s0 is not None
        assert s1 is not None
        # Order is the header sample order: TUMOR then NORMAL
        assert s0["GT"] == (0, 1)
        assert s1["GT"] == (0, 0)

    def test_raises_indexerror(self, basic_record):
        with pytest.raises(IndexError):
            get_sample_by_index(basic_record, 999)


class TestHasFormatField:
    """Tests for has_format_field."""

    def test_true_when_present(self, basic_record):
        assert has_format_field(basic_record, "GT") is True
        assert has_format_field(basic_record, "AD") is True

    def test_false_when_absent(self, basic_record):
        assert has_format_field(basic_record, "DP") is False
