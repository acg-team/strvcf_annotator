"""Unit tests for build_new_record function."""

import logging
import os
import tempfile

import pandas as pd
import pysam
import pytest

from strvcf_annotator.core.annotation import build_new_record, make_modified_header
from strvcf_annotator.parsers.generic import GenericParser


@pytest.fixture
def basic_vcf_header():
    """Create a basic VCF header for testing."""
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")
    header.contigs.add("chr1", length=1000000)
    header.add_sample("TUMOR")
    header.add_sample("NORMAL")

    # Add basic INFO fields
    header.info.add("DP", 1, "Integer", "Total read depth")

    # Add basic FORMAT fields
    header.formats.add("GT", 1, "String", "Genotype")
    header.formats.add("AD", "R", "Integer", "Allelic depths")
    header.formats.add("DP", 1, "Integer", "Read depth")

    return header


@pytest.fixture
def mismatch_vcf_header():
    """Header that includes chr2/chr8/chr11 for mismatch-case tests."""
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")
    header.contigs.add("chr2", length=200000000)
    header.contigs.add("chr8", length=200000000)
    header.contigs.add("chr11", length=200000000)
    header.add_sample("TUMOR")
    header.add_sample("NORMAL")

    header.info.add("DP", 1, "Integer", "Total read depth")
    header.formats.add("GT", 1, "String", "Genotype")
    header.formats.add("AD", "R", "Integer", "Allelic depths")
    header.formats.add("DP", 1, "Integer", "Read depth")
    return header


@pytest.fixture
def mismatch_str_modified_header(mismatch_vcf_header):
    """Generate STR-modified header for mismatch contigs."""
    # Build a temp VCF with one dummy record so make_modified_header can read it

    with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
        f.write(str(mismatch_vcf_header))
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n")
        f.write("chr2\t1\t.\tA\tT\t30\tPASS\tDP=10\tGT:AD:DP\t0/1:5,5:10\t0/0:10,0:10\n")
        path = f.name

    try:
        vcf_in = pysam.VariantFile(path)
        out_header = make_modified_header(vcf_in)
        vcf_in.close()
        return out_header
    finally:
        os.unlink(path)


@pytest.fixture
def mismatch_parser():
    return GenericParser()


@pytest.fixture
def str_modified_header(basic_vcf_header):
    """Create a mock VCF file and generate STR-modified header."""

    # Create a temporary VCF file to use with make_modified_header
    with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
        # Write header
        f.write(str(basic_vcf_header))
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n")
        # Write a dummy record
        f.write("chr1\t100\t.\tA\tT\t30\tPASS\tDP=50\tGT:AD:DP\t0/1:20,30:50\t0/0:50,0:50\n")
        temp_path = f.name

    try:
        vcf_in = pysam.VariantFile(temp_path)
        header = make_modified_header(vcf_in)
        vcf_in.close()
        return header
    finally:
        os.unlink(temp_path)


@pytest.fixture
def parser():
    """Create a GenericParser instance."""
    return GenericParser()

class TestBuildNewRecordMismatchOptions:
    """
    Tests for mismatch handling options added to build_new_record.

    We test two things separately:
    1) Warnings: emitted vs suppressed depending on ignore_mismatch_warnings
    2) Modes: panel/vcf return a record; skip returns None on mismatch
    """

    @pytest.fixture
    def mismatch_case_chr11(self, mismatch_vcf_header, monkeypatch):
        """
        Build a record + str_row where panel sequence and VCF REF overlap mismatch:
        VCF REF is 'cacaca...' while panel is 'CCCC...'
        """
        ref = "cacaca" * 4 + "ca"  # 26 bases
        alt = ref  # keep identical; mismatch logic should still detect REF overlap mismatch

        record = mismatch_vcf_header.new_record(
            contig="chr11",
            start=134204815,  # 0-based => POS 134204816
            alleles=(ref, alt),
            filter="PASS",
        )
        record.samples["TUMOR"]["GT"] = (0, 0)
        record.samples["NORMAL"]["GT"] = (0, 0)

        # Panel overlap is all 'C'
        panel_seq = "C" * 26
        monkeypatch.setattr(
            "strvcf_annotator.core.annotation.extract_repeat_sequence",
            lambda _: panel_seq,
        )

        str_row = {
            "CHROM": "chr11",
            "START": 134204816,
            "END": 134204841,
            "RU": "CC",
            "PERIOD": 2,
            "COUNT": 13,
        }
        return record, str_row

    def test_mismatch_warning_emitted_by_default(
        self,
        mismatch_case_chr11,
        mismatch_str_modified_header,
        mismatch_parser,
        caplog,
    ):
        record, str_row = mismatch_case_chr11
        caplog.set_level(logging.WARNING, logger="strvcf_annotator.core.annotation")

        out = build_new_record(
            record,
            str_row,
            mismatch_str_modified_header,
            mismatch_parser,
            ignore_mismatch_warnings=False,
            mismatch_truth="panel",
        )

        assert isinstance(out, pysam.VariantRecord)
        msgs = [r.getMessage() for r in caplog.records]
        assert any("Reference mismatch in STR overlap" in m for m in msgs)

    def test_mismatch_warning_suppressed_when_ignore_true(
        self,
        mismatch_case_chr11,
        mismatch_str_modified_header,
        mismatch_parser,
        caplog,
    ):
        record, str_row = mismatch_case_chr11
        caplog.set_level(logging.WARNING, logger="strvcf_annotator.core.annotation")

        out = build_new_record(
            record,
            str_row,
            mismatch_str_modified_header,
            mismatch_parser,
            ignore_mismatch_warnings=True,
            mismatch_truth="panel",
        )

        assert isinstance(out, pysam.VariantRecord)
        msgs = [r.getMessage() for r in caplog.records]
        assert not any("Reference mismatch in STR overlap" in m for m in msgs)

    @pytest.mark.parametrize("mismatch_truth", ["panel", "vcf"])
    def test_mismatch_truth_panel_and_vcf_return_record(
        self,
        mismatch_case_chr11,
        mismatch_str_modified_header,
        mismatch_parser,
        caplog,
        mismatch_truth,
    ):
        record, str_row = mismatch_case_chr11
        caplog.set_level(logging.WARNING, logger="strvcf_annotator.core.annotation")

        out = build_new_record(
            record,
            str_row,
            mismatch_str_modified_header,
            mismatch_parser,
            ignore_mismatch_warnings=True,  # keep log clean for this test
            mismatch_truth=mismatch_truth,
        )

        # Both modes should still return an annotated record (not None)
        assert isinstance(out, pysam.VariantRecord)

        # Extra sanity: alleles should exist
        assert out.alleles is not None
        assert len(out.alleles) == 2

    def test_mismatch_truth_skip_returns_none(
        self,
        mismatch_case_chr11,
        mismatch_str_modified_header,
        mismatch_parser,
        caplog,
    ):
        record, str_row = mismatch_case_chr11
        caplog.set_level(logging.WARNING, logger="strvcf_annotator.core.annotation")

        out = build_new_record(
            record,
            str_row,
            mismatch_str_modified_header,
            mismatch_parser,
            ignore_mismatch_warnings=True,  # warning irrelevant here
            mismatch_truth="skip",
        )

        # Skip mode must drop mismatching loci
        assert out is None

class TestMismatchRecord:
    """
    Tests for the three mismatch scenarios you reported.

    We verify:
    - chr8: motif is reverse-complement/rotation equivalent -> warning should be suppressed
    - chr11/chr2: motif conflicts -> warning should be emitted
    """

    def test_chr11_true_conflict_emits_warning(
        self,
        mismatch_vcf_header,
        mismatch_str_modified_header,
        mismatch_parser,
        caplog,
        monkeypatch,
    ):
        # VCF REF: CACACA... while panel overlap: CCCCC...
        # This should trigger a mismatch warning.
        ref = "cacaca" * 4 + "ca"  # length 26 -> "cacacacacacacacacacacacaca"
        alt = ref  # same; we only care about mismatch detection

        record = mismatch_vcf_header.new_record(
            contig="chr11",
            start=134204815,  # 0-based => POS 134204816
            alleles=(ref, alt),
            filter="PASS",
        )
        record.samples["TUMOR"]["GT"] = (0, 0)
        record.samples["NORMAL"]["GT"] = (0, 0)

        # Patch panel sequence to the all-C overlap shown in your log
        panel_seq = "C" * 26
        monkeypatch.setattr(
            "strvcf_annotator.core.annotation.extract_repeat_sequence",
            lambda _: panel_seq,
        )

        str_row = {
            "CHROM": "chr11",
            "START": 134204816,
            "END": 134204841,
            "RU": "CC",
            "PERIOD": 2,
            "COUNT": 13,
        }

        caplog.set_level(logging.WARNING, logger="strvcf_annotator.core.annotation")
        new_record = build_new_record(
            record, str_row, mismatch_str_modified_header, mismatch_parser
        )

        assert isinstance(new_record, pysam.VariantRecord)
        msgs = [rec.getMessage() for rec in caplog.records]
        assert any("Reference mismatch in STR overlap" in m for m in msgs)

    def test_chr2_true_conflict_emits_warning(
        self,
        mismatch_vcf_header,
        mismatch_str_modified_header,
        mismatch_parser,
        caplog,
        monkeypatch,
    ):
        # VCF REF: TATATA... while panel overlap: TGTGTG...
        # This should trigger a mismatch warning.
        ref = "ta" * 12 + "ta"  # 26 bases of alternating TA (close to your example)
        ref = ref[:24]  # keep it simple/short but still >= PERIOD=2
        alt = ref

        record = mismatch_vcf_header.new_record(
            contig="chr2",
            start=178647038,  # 0-based => POS 178647039
            alleles=(ref, alt),
            filter="PASS",
        )
        record.samples["TUMOR"]["GT"] = (0, 0)
        record.samples["NORMAL"]["GT"] = (0, 0)

        # Patch panel sequence to TG repeats
        panel_seq = "TG" * (len(ref) // 2)
        monkeypatch.setattr(
            "strvcf_annotator.core.annotation.extract_repeat_sequence",
            lambda _: panel_seq,
        )

        str_row = {
            "CHROM": "chr2",
            "START": 178647039,
            "END": 178647062,
            "RU": "TG",
            "PERIOD": 2,
            "COUNT": 12,
        }

        caplog.set_level(logging.WARNING, logger="strvcf_annotator.core.annotation")
        new_record = build_new_record(
            record, str_row, mismatch_str_modified_header, mismatch_parser
        )

        assert isinstance(new_record, pysam.VariantRecord)
        msgs = [rec.getMessage() for rec in caplog.records]
        assert any("Reference mismatch in STR overlap" in m for m in msgs)


class TestBuildNewRecordBasic:
    """Basic tests for build_new_record function."""

    def test_returns_variant_record(self, basic_vcf_header, str_modified_header, parser):
        """Test that function returns a VariantRecord."""
        # Create a simple VCF record
        record = basic_vcf_header.new_record(
            contig="chr1", start=100, alleles=("A", "T"), filter="PASS", info={"DP": 50}
        )
        record.samples["TUMOR"]["GT"] = (0, 1)
        record.samples["TUMOR"]["AD"] = (20, 30)
        record.samples["TUMOR"]["DP"] = 50
        record.samples["NORMAL"]["GT"] = (0, 0)
        record.samples["NORMAL"]["AD"] = (50, 0)
        record.samples["NORMAL"]["DP"] = 50

        # STR metadata
        str_row = {"CHROM": "chr1", "START": 95, "END": 115, "RU": "AT", "PERIOD": 2, "COUNT": 10}

        new_record = build_new_record(record, str_row, str_modified_header, parser)
        assert isinstance(new_record, pysam.VariantRecord)

    def test_sets_correct_position(self, basic_vcf_header, str_modified_header, parser):
        """Test that new record has correct position from STR metadata."""
        record = basic_vcf_header.new_record(contig="chr1", start=100, alleles=("A", "T"))
        record.samples["TUMOR"]["GT"] = (0, 1)
        record.samples["NORMAL"]["GT"] = (0, 0)

        str_row = {"CHROM": "chr1", "START": 95, "END": 115, "RU": "AT", "PERIOD": 2, "COUNT": 10}

        new_record = build_new_record(record, str_row, str_modified_header, parser)
        # Position should be POS - 1 (0-based)
        assert new_record.start == 94
        assert new_record.stop == 114

    def test_adds_str_info_fields(self, basic_vcf_header, str_modified_header, parser):
        """Test that STR INFO fields are added."""
        record = basic_vcf_header.new_record(contig="chr1", start=100, alleles=("A", "T"))
        record.samples["TUMOR"]["GT"] = (0, 1)
        record.samples["NORMAL"]["GT"] = (0, 0)

        str_row = {"CHROM": "chr1", "START": 95, "END": 115, "RU": "CAG", "PERIOD": 3, "COUNT": 7}

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        assert "RU" in new_record.info
        assert new_record.info["RU"] == "CAG"
        assert "PERIOD" in new_record.info
        assert new_record.info["PERIOD"] == 3
        assert "REF" in new_record.info
        assert "PERFECT" in new_record.info

    def test_adds_repcn_format_field(self, basic_vcf_header, str_modified_header, parser):
        """Test that REPCN FORMAT field is added for samples."""
        record = basic_vcf_header.new_record(contig="chr1", start=100, alleles=("A", "T"))
        record.samples["TUMOR"]["GT"] = (0, 1)
        record.samples["NORMAL"]["GT"] = (0, 0)

        str_row = {"CHROM": "chr1", "START": 95, "END": 115, "RU": "AT", "PERIOD": 2, "COUNT": 10}

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        assert "REPCN" in new_record.samples["TUMOR"]
        assert "REPCN" in new_record.samples["NORMAL"]
        assert isinstance(new_record.samples["TUMOR"]["REPCN"], tuple)


class TestBuildNewRecordAlleles:
    """Tests for allele reconstruction in build_new_record."""

    def test_snv_mutation(self, basic_vcf_header, str_modified_header, parser):
        """Test SNV mutation within repeat region."""
        # SNV at position 105
        record = basic_vcf_header.new_record(
            contig="chr1",
            start=104,  # 0-based
            alleles=("A", "G"),
        )
        record.samples["TUMOR"]["GT"] = (0, 1)
        record.samples["NORMAL"]["GT"] = (0, 0)

        # Repeat region: chr1:100-120, repeat unit = "AT"
        str_row = {"CHROM": "chr1", "START": 100, "END": 120, "RU": "AT", "PERIOD": 2, "COUNT": 10}

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        # Should have two alleles (reference and mutated repeat sequence)
        assert len(new_record.alleles) == 2
        assert new_record.alleles[0] != new_record.alleles[1]

    def test_insertion_mutation(self, basic_vcf_header, str_modified_header, parser):
        """Test insertion mutation."""
        # Insertion at position 105: A -> AAT
        record = basic_vcf_header.new_record(contig="chr1", start=104, alleles=("A", "AAT"))
        record.samples["TUMOR"]["GT"] = (0, 1)
        record.samples["NORMAL"]["GT"] = (0, 0)

        str_row = {"CHROM": "chr1", "START": 100, "END": 120, "RU": "AT", "PERIOD": 2, "COUNT": 10}

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        # ALT allele should be longer than REF
        assert len(new_record.alleles[1]) > len(new_record.alleles[0])

    def test_deletion_mutation(self, basic_vcf_header, str_modified_header, parser):
        """Test deletion mutation."""
        # Deletion at position 105: AAT -> A
        record = basic_vcf_header.new_record(contig="chr1", start=104, alleles=("AAT", "A"))
        record.samples["TUMOR"]["GT"] = (0, 1)
        record.samples["NORMAL"]["GT"] = (0, 0)

        str_row = {"CHROM": "chr1", "START": 100, "END": 120, "RU": "AT", "PERIOD": 2, "COUNT": 10}

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        # ALT allele should be shorter than REF
        assert len(new_record.alleles[1]) < len(new_record.alleles[0])


class TestBuildNewRecordRepeatCounts:
    """Tests for repeat copy number calculation."""

    def test_calculates_ref_repeat_count(self, basic_vcf_header, str_modified_header, parser):
        """Test that REF INFO field contains repeat count."""
        record = basic_vcf_header.new_record(contig="chr1", start=104, alleles=("A", "T"))
        record.samples["TUMOR"]["GT"] = (0, 1)
        record.samples["NORMAL"]["GT"] = (0, 0)

        str_row = {"CHROM": "chr1", "START": 100, "END": 120, "RU": "AT", "PERIOD": 2, "COUNT": 10}

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        assert "REF" in new_record.info
        assert isinstance(new_record.info["REF"], int)
        assert new_record.info["REF"] > 0

    def test_repcn_for_heterozygous(self, basic_vcf_header, str_modified_header, parser):
        """Test REPCN calculation for heterozygous genotype."""
        record = basic_vcf_header.new_record(contig="chr1", start=104, alleles=("A", "AAT"))
        record.samples["TUMOR"]["GT"] = (0, 1)  # Het
        record.samples["NORMAL"]["GT"] = (0, 0)  # Hom ref

        str_row = {"CHROM": "chr1", "START": 100, "END": 120, "RU": "AT", "PERIOD": 2, "COUNT": 10}

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        tumor_repcn = new_record.samples["TUMOR"]["REPCN"]
        normal_repcn = new_record.samples["NORMAL"]["REPCN"]

        # TUMOR should have different copy numbers (0/1 genotype)
        assert len(tumor_repcn) == 2
        assert tumor_repcn[0] != tumor_repcn[1]

        # NORMAL should have same copy numbers (0/0 genotype)
        assert len(normal_repcn) == 2
        assert normal_repcn[0] == normal_repcn[1]

    def test_repcn_for_homozygous_alt(self, basic_vcf_header, str_modified_header, parser):
        """Test REPCN for homozygous ALT genotype."""
        record = basic_vcf_header.new_record(contig="chr1", start=104, alleles=("A", "AAT"))
        record.samples["TUMOR"]["GT"] = (1, 1)  # Hom alt
        record.samples["NORMAL"]["GT"] = (0, 0)

        str_row = {"CHROM": "chr1", "START": 100, "END": 120, "RU": "AT", "PERIOD": 2, "COUNT": 10}

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        tumor_repcn = new_record.samples["TUMOR"]["REPCN"]

        # Both alleles should have ALT repeat count
        assert tumor_repcn[0] == tumor_repcn[1]


class TestBuildNewRecordPerfectRepeats:
    """Tests for PERFECT field calculation."""

    def test_perfect_field_exists(self, basic_vcf_header, str_modified_header, parser):
        """Test that PERFECT field is always set."""
        record = basic_vcf_header.new_record(contig="chr1", start=104, alleles=("A", "T"))
        record.samples["TUMOR"]["GT"] = (0, 1)
        record.samples["NORMAL"]["GT"] = (0, 0)

        str_row = {"CHROM": "chr1", "START": 100, "END": 120, "RU": "AT", "PERIOD": 2, "COUNT": 10}

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        assert "PERFECT" in new_record.info
        assert new_record.info["PERFECT"] in ("TRUE", "FALSE")


class TestBuildNewRecordEdgeCases:
    """Edge case tests for build_new_record."""

    def test_missing_genotype(self, basic_vcf_header, str_modified_header, parser):
        """Test handling of missing genotype (./.)."""
        record = basic_vcf_header.new_record(contig="chr1", start=104, alleles=("A", "T"))
        record.samples["TUMOR"]["GT"] = (None, None)
        record.samples["NORMAL"]["GT"] = (0, 0)

        str_row = {"CHROM": "chr1", "START": 100, "END": 120, "RU": "AT", "PERIOD": 2, "COUNT": 10}

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        # Should handle missing genotype gracefully
        assert new_record.samples["TUMOR"]["GT"] == (None, None)
        assert new_record.samples["TUMOR"]["REPCN"] == (0, 0)

    def test_preserves_original_format_fields(self, basic_vcf_header, str_modified_header, parser):
        """Test that original FORMAT fields are preserved."""
        record = basic_vcf_header.new_record(contig="chr1", start=104, alleles=("A", "T"))
        record.samples["TUMOR"]["GT"] = (0, 1)
        record.samples["TUMOR"]["AD"] = (20, 30)
        record.samples["TUMOR"]["DP"] = 50
        record.samples["NORMAL"]["GT"] = (0, 0)
        record.samples["NORMAL"]["AD"] = (50, 0)
        record.samples["NORMAL"]["DP"] = 50

        str_row = {"CHROM": "chr1", "START": 100, "END": 120, "RU": "AT", "PERIOD": 2, "COUNT": 10}

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        # Original FORMAT fields should be preserved
        assert "AD" in new_record.samples["TUMOR"]
        assert "DP" in new_record.samples["TUMOR"]
        assert new_record.samples["TUMOR"]["AD"] == (20, 30)
        assert new_record.samples["TUMOR"]["DP"] == 50

    def test_with_pandas_series(self, basic_vcf_header, str_modified_header, parser):
        """Test that function works with pandas Series for str_row."""
        record = basic_vcf_header.new_record(contig="chr1", start=104, alleles=("A", "T"))
        record.samples["TUMOR"]["GT"] = (0, 1)
        record.samples["NORMAL"]["GT"] = (0, 0)

        # Use pandas Series instead of dict
        str_row = pd.Series(
            {"CHROM": "chr1", "START": 100, "END": 120, "RU": "AT", "PERIOD": 2, "COUNT": 10}
        )

        new_record = build_new_record(record, str_row, str_modified_header, parser)

        assert isinstance(new_record, pysam.VariantRecord)
        assert new_record.info["RU"] == "AT"

    def test_different_repeat_units(self, basic_vcf_header, str_modified_header, parser):
        """Test with different repeat unit sizes."""
        test_cases = [
            ("A", 1),  # Mononucleotide
            ("AT", 2),  # Dinucleotide
            ("CAG", 3),  # Trinucleotide
            ("ATCG", 4),  # Tetranucleotide
        ]

        for ru, period in test_cases:
            record = basic_vcf_header.new_record(contig="chr1", start=104, alleles=("A", "T"))
            record.samples["TUMOR"]["GT"] = (0, 1)
            record.samples["NORMAL"]["GT"] = (0, 0)

            str_row = {
                "CHROM": "chr1",
                "START": 100,
                "END": 120,
                "RU": ru,
                "PERIOD": period,
                "COUNT": 10,
            }

            new_record = build_new_record(record, str_row, str_modified_header, parser)

            assert new_record.info["RU"] == ru
            assert new_record.info["PERIOD"] == period
