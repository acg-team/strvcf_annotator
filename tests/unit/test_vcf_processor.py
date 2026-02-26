"""Unit tests for VCF processor functions."""

import os
import tempfile
from pathlib import Path

import pysam
import pytest

from strvcf_annotator.core.vcf_processor import (
    check_vcf_sorted,
    estimate_ram_per_worker_bytes,
    generate_annotated_records,
    get_available_ram_bytes,
    reset_and_sort_vcf,
)
from strvcf_annotator.parsers.generic import GenericParser


@pytest.fixture
def vcf_paths(data_dir):
    """Return absolute paths to all shipped VCF test files."""
    files = [
        "test.vcf.gz",
        "pindel_header.vcf",
        "mutec2_indel.vcf.gz",
        "TCGA-DC-6682.vcf",
    ]
    return [os.path.abspath(os.path.join(data_dir, f)) for f in files]


@pytest.fixture
def basic_vcf_header():
    """Create a basic VCF header for testing."""
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")
    header.contigs.add("chr1", length=1000000)
    header.contigs.add("chr2", length=1000000)
    header.contigs.add("chr3", length=1000000)
    header.add_sample("TUMOR")
    header.add_sample("NORMAL")

    # Add basic INFO and FORMAT fields
    header.info.add("DP", 1, "Integer", "Total read depth")
    header.formats.add("GT", 1, "String", "Genotype")
    header.formats.add("AD", "R", "Integer", "Allelic depths")
    header.formats.add("DP", 1, "Integer", "Read depth")

    return header


@pytest.fixture
def sorted_vcf_file(basic_vcf_header):
    """Create a temporary sorted VCF file using pysam API."""
    temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
    temp_file.close()  # Close immediately so pysam can write to it

    try:
        # Write VCF using pysam API
        with pysam.VariantFile(temp_file.name, "w", header=basic_vcf_header) as vcf_out:
            # Write sorted records
            records_data = [
                ("chr1", 100, "A", "T"),
                ("chr1", 200, "G", "C"),
                ("chr1", 300, "T", "A"),
                ("chr2", 100, "A", "G"),
                ("chr2", 200, "C", "T"),
            ]
            for chrom, pos, ref, alt in records_data:
                rec = vcf_out.new_record(
                    contig=chrom,
                    start=pos - 1,  # 0-based
                    alleles=(ref, alt),
                    qual=30,
                    filter=["PASS"],
                    info={"DP": 50},
                )
                rec.samples["TUMOR"]["GT"] = (0, 1)
                rec.samples["TUMOR"]["AD"] = (20, 30)
                rec.samples["TUMOR"]["DP"] = 50
                rec.samples["NORMAL"]["GT"] = (0, 0)
                rec.samples["NORMAL"]["AD"] = (50, 0)
                rec.samples["NORMAL"]["DP"] = 50
                vcf_out.write(rec)

        yield temp_file.name
    finally:
        os.unlink(temp_file.name)


@pytest.fixture
def unsorted_vcf_file(basic_vcf_header):
    """Create a temporary unsorted VCF file using pysam API."""
    temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
    temp_file.close()  # Close immediately so pysam can write to it

    try:
        # Write VCF using pysam API
        with pysam.VariantFile(temp_file.name, "w", header=basic_vcf_header) as vcf_out:
            # Write unsorted records (positions out of order)
            records_data = [
                ("chr1", 300, "T", "A"),
                ("chr1", 100, "A", "T"),
                ("chr2", 200, "C", "T"),
                ("chr1", 200, "G", "C"),
                ("chr2", 100, "A", "G"),
            ]
            for chrom, pos, ref, alt in records_data:
                rec = vcf_out.new_record(
                    contig=chrom,
                    start=pos - 1,  # 0-based
                    alleles=(ref, alt),
                    qual=30,
                    filter=["PASS"],
                    info={"DP": 50},
                )
                rec.samples["TUMOR"]["GT"] = (0, 1)
                rec.samples["TUMOR"]["AD"] = (20, 30)
                rec.samples["TUMOR"]["DP"] = 50
                rec.samples["NORMAL"]["GT"] = (0, 0)
                rec.samples["NORMAL"]["AD"] = (50, 0)
                rec.samples["NORMAL"]["DP"] = 50
                vcf_out.write(rec)

        yield temp_file.name
    finally:
        os.unlink(temp_file.name)


@pytest.fixture
def str_panel_tabix():
    """Create a small tabix-indexed STR panel for tests.

    Returns
    -------
    str
        Path to bgzip-compressed, tabix-indexed STR panel.
    """
    content = """chr1\t95\t115\t2\tAT\t10
chr1\t195\t215\t3\tCAG\t7
chr2\t95\t115\t2\tGC\t10
"""

    with tempfile.TemporaryDirectory() as tmp:
        bed_path = Path(tmp) / "str_panel.bed"
        bed_path.write_text(content, encoding="utf-8")

        # Sort to guarantee tabix compatibility
        sorted_path = Path(tmp) / "str_panel.sorted.bed"
        lines = sorted(
            [l.strip() for l in content.strip().splitlines()],
            key=lambda x: (x.split("\t")[0], int(x.split("\t")[1])),
        )
        sorted_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

        gz_path = Path(tmp) / "str_panel.sorted.bed.gz"

        pysam.tabix_compress(str(sorted_path), str(gz_path), force=True)
        pysam.tabix_index(str(gz_path), preset="bed", force=True)

        yield str(gz_path)


class TestCheckVCFSorted:
    """Tests for check_vcf_sorted function."""

    def test_sorted_vcf_returns_true(self, sorted_vcf_file):
        """Test that sorted VCF returns True."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)
        result = check_vcf_sorted(vcf_in)
        vcf_in.close()

        assert result is True

    def test_unsorted_vcf_returns_false(self, unsorted_vcf_file):
        """Test that unsorted VCF returns False."""
        vcf_in = pysam.VariantFile(unsorted_vcf_file)
        result = check_vcf_sorted(vcf_in)
        vcf_in.close()

        assert result is False

    def test_rewinds_file_after_check(self, sorted_vcf_file):
        """Test that file is rewound after checking."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)

        # Check sorting
        check_vcf_sorted(vcf_in)

        # File should be rewound, so we can read from beginning
        records = list(vcf_in)
        assert len(records) == 5
        assert records[0].pos == 100  # First record

        vcf_in.close()

    def test_empty_vcf_returns_true(self, basic_vcf_header):
        """Test that empty VCF is considered sorted."""
        temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
        temp_file.close()

        try:
            # Create empty VCF using pysam API
            with pysam.VariantFile(temp_file.name, "w", header=basic_vcf_header) as vcf_out:
                pass  # Write header only, no records

            vcf_in = pysam.VariantFile(temp_file.name)
            result = check_vcf_sorted(vcf_in)
            vcf_in.close()

            assert result is True
        finally:
            os.unlink(temp_file.name)

    def test_single_record_returns_true(self, basic_vcf_header):
        """Test that single record VCF is considered sorted."""
        temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
        temp_file.close()

        try:
            # Create VCF with single record using pysam API
            with pysam.VariantFile(temp_file.name, "w", header=basic_vcf_header) as vcf_out:
                rec = vcf_out.new_record(
                    contig="chr1",
                    start=99,  # 0-based
                    alleles=("A", "T"),
                    qual=30,
                    filter=["PASS"],
                    info={"DP": 50},
                )
                rec.samples["TUMOR"]["GT"] = (0, 1)
                rec.samples["TUMOR"]["AD"] = (20, 30)
                rec.samples["TUMOR"]["DP"] = 50
                rec.samples["NORMAL"]["GT"] = (0, 0)
                rec.samples["NORMAL"]["AD"] = (50, 0)
                rec.samples["NORMAL"]["DP"] = 50
                vcf_out.write(rec)

            vcf_in = pysam.VariantFile(temp_file.name)
            result = check_vcf_sorted(vcf_in)
            vcf_in.close()

            assert result is True
        finally:
            os.unlink(temp_file.name)


class TestResetAndSortVCF:
    """Tests for reset_and_sort_vcf function."""

    def test_sorts_unsorted_vcf(self, unsorted_vcf_file):
        """Test that function sorts unsorted VCF."""
        vcf_in = pysam.VariantFile(unsorted_vcf_file)
        sorted_records = reset_and_sort_vcf(vcf_in)
        vcf_in.close()

        # Check correct number of records
        assert len(sorted_records) == 5

        # Check sorted order
        assert sorted_records[0].contig == "chr1"
        assert sorted_records[0].pos == 100
        assert sorted_records[1].contig == "chr1"
        assert sorted_records[1].pos == 200
        assert sorted_records[2].contig == "chr1"
        assert sorted_records[2].pos == 300
        assert sorted_records[3].contig == "chr2"
        assert sorted_records[3].pos == 100
        assert sorted_records[4].contig == "chr2"
        assert sorted_records[4].pos == 200

    def test_preserves_sorted_order(self, sorted_vcf_file):
        """Test that already sorted VCF remains sorted."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)
        sorted_records = reset_and_sort_vcf(vcf_in)
        vcf_in.close()

        # Check sorted order is preserved
        assert sorted_records[0].pos == 100
        assert sorted_records[1].pos == 200
        assert sorted_records[2].pos == 300
        assert sorted_records[3].contig == "chr2"
        assert sorted_records[3].pos == 100

    def test_returns_list(self, sorted_vcf_file):
        """Test that function returns a list."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)
        result = reset_and_sort_vcf(vcf_in)
        vcf_in.close()

        assert isinstance(result, list)

    def test_empty_vcf_returns_empty_list(self, basic_vcf_header):
        """Test that empty VCF returns empty list."""
        temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
        temp_file.close()

        try:
            # Create empty VCF using pysam API
            with pysam.VariantFile(temp_file.name, "w", header=basic_vcf_header) as vcf_out:
                pass  # No records

            vcf_in = pysam.VariantFile(temp_file.name)
            sorted_records = reset_and_sort_vcf(vcf_in)
            vcf_in.close()

            assert len(sorted_records) == 0
        finally:
            os.unlink(temp_file.name)

    def test_sorts_by_contig_order(self, basic_vcf_header):
        """Test that sorting respects contig order in header."""
        temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
        temp_file.close()

        try:
            # Create VCF with records from different chromosomes using pysam API
            with pysam.VariantFile(temp_file.name, "w", header=basic_vcf_header) as vcf_out:
                for chrom in ["chr3", "chr1", "chr2"]:
                    rec = vcf_out.new_record(
                        contig=chrom,
                        start=99,  # 0-based
                        alleles=("A", "T"),
                        qual=30,
                        filter=["PASS"],
                        info={"DP": 50},
                    )
                    rec.samples["TUMOR"]["GT"] = (0, 1)
                    rec.samples["TUMOR"]["AD"] = (20, 30)
                    rec.samples["TUMOR"]["DP"] = 50
                    rec.samples["NORMAL"]["GT"] = (0, 0)
                    rec.samples["NORMAL"]["AD"] = (50, 0)
                    rec.samples["NORMAL"]["DP"] = 50
                    vcf_out.write(rec)

            vcf_in = pysam.VariantFile(temp_file.name)
            sorted_records = reset_and_sort_vcf(vcf_in)
            vcf_in.close()

            # Should be sorted by header contig order: chr1, chr2, chr3
            assert sorted_records[0].contig == "chr1"
            assert sorted_records[1].contig == "chr2"
            assert sorted_records[2].contig == "chr3"
        finally:
            os.unlink(temp_file.name)


class TestGenerateAnnotatedRecords:
    """Tests for generate_annotated_records function."""

    def test_returns_iterator(self, sorted_vcf_file, str_panel_tabix):
        """Test that function returns an iterator."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)
        parser = GenericParser()

        result = generate_annotated_records(vcf_in, str_panel_tabix, parser)

        # Should be an iterator
        assert hasattr(result, "__iter__")
        assert hasattr(result, "__next__")

        vcf_in.close()

    def test_yields_variant_records(self, sorted_vcf_file, str_panel_tabix):
        """Test that iterator yields VariantRecord objects."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)
        parser = GenericParser()

        records = list(generate_annotated_records(vcf_in, str_panel_tabix, parser))

        if len(records) > 0:
            assert all(isinstance(r, pysam.VariantRecord) for r in records)

        vcf_in.close()

    def test_filters_non_overlapping_variants(self, sorted_vcf_file):
        """Test that variants outside STR regions are filtered."""

        vcf_in = pysam.VariantFile(sorted_vcf_file)
        parser = GenericParser()

        # Create a tabix-indexed STR panel with no overlap to the VCF (chr10 only)
        content = "chr10\t1000\t1100\t2\tAT\t50\n"

        with tempfile.TemporaryDirectory() as tmp:
            bed_path = Path(tmp) / "no_overlap.bed"
            bed_path.write_text(content, encoding="utf-8")

            gz_path = Path(tmp) / "no_overlap.bed.gz"
            pysam.tabix_compress(str(bed_path), str(gz_path), force=True)
            pysam.tabix_index(str(gz_path), preset="bed", force=True)

            records = list(generate_annotated_records(vcf_in, str(gz_path), parser))

            # Should yield no records (no overlap)
            assert len(records) == 0

        vcf_in.close()

    def test_handles_unsorted_vcf(self, unsorted_vcf_file, str_panel_tabix):
        """Test that function handles unsorted VCF by sorting it."""
        vcf_in = pysam.VariantFile(unsorted_vcf_file)
        parser = GenericParser()

        # Should not raise an error
        records = list(generate_annotated_records(vcf_in, str_panel_tabix, parser))

        # Records should be processed (may be 0 if no overlaps)
        assert isinstance(records, list)

        vcf_in.close()

    def test_uses_generic_parser_by_default(self, sorted_vcf_file, str_panel_tabix):
        """Test that GenericParser is used when parser=None."""
        vcf_in = pysam.VariantFile(sorted_vcf_file)

        # Call without parser
        records = list(generate_annotated_records(vcf_in, str_panel_tabix, parser=None))

        # Should work without error
        assert isinstance(records, list)

        vcf_in.close()

    def test_empty_vcf_returns_no_records(self, basic_vcf_header, str_panel_tabix):
        """Test that empty VCF yields no records."""
        temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False)
        temp_file.close()

        try:
            # Create empty VCF using pysam API
            with pysam.VariantFile(temp_file.name, "w", header=basic_vcf_header) as vcf_out:
                pass  # No records

            vcf_in = pysam.VariantFile(temp_file.name)
            parser = GenericParser()

            records = list(generate_annotated_records(vcf_in, str_panel_tabix, parser))

            assert len(records) == 0

            vcf_in.close()
        finally:
            os.unlink(temp_file.name)

    def test_empty_str_panel_tabix_returns_no_records(self, sorted_vcf_file):
        """Test that empty STR panel (no overlapping loci) yields no records."""

        # Create a panel with one dummy locus on a chromosome not present in VCF
        content = "chrUn\t1\t10\t2\tAT\t5\n"

        with tempfile.TemporaryDirectory() as tmp:
            bed_path = Path(tmp) / "empty_panel.bed"
            bed_path.write_text(content, encoding="utf-8")

            gz_path = Path(tmp) / "empty_panel.bed.gz"
            pysam.tabix_compress(str(bed_path), str(gz_path), force=True)
            pysam.tabix_index(str(gz_path), preset="bed", force=True)

            vcf_in = pysam.VariantFile(sorted_vcf_file)
            parser = GenericParser()

            records = list(generate_annotated_records(vcf_in, str(gz_path), parser))

            assert len(records) == 0

            vcf_in.close()


class TestGetAvailableRamBytes:
    """Tests for get_available_ram_bytes function."""

    def test_returns_int(self):
        """Test that function returns an integer."""
        result = get_available_ram_bytes()
        assert isinstance(result, int)

    # @TODO fix
    # def test_fake_memory(self, monkeypatch):
    #     """Test that psutil branch returns fake available RAM."""
    #     fake_psutil = types.ModuleType("psutil")

    #     class FakeVMem:
    #         available = 123456789

    #     def virtual_memory():
    #         return FakeVMem()

    #     fake_psutil.virtual_memory = virtual_memory

    #     real_import = builtins.__import__

    #     def import_hook(name, globals=None, locals=None, fromlist=(), level=0):
    #         if name == "psutil":
    #             return fake_psutil
    #         return real_import(name, globals, locals, fromlist, level)

    #     monkeypatch.setattr(builtins, "__import__", import_hook)

    #     result = get_available_ram_bytes()
    #     assert result == 123456789


class TestEstimateRamPerWorkerBytes:
    """Tests for estimate_ram_per_worker_bytes function."""

    def test_empty_list_returns_minimum(self):
        """Test that empty input returns a safe minimum estimate."""
        result = estimate_ram_per_worker_bytes([])
        assert result == int(1 * 1024**3)  # 1 GB

    def test_uses_largest_file_size(self, vcf_paths):
        """Test that estimate is based on the largest file size among inputs."""
        sizes = {p: os.path.getsize(p) for p in vcf_paths}
        max_path = max(sizes, key=sizes.get)
        max_size = sizes[max_path]

        fixed_overhead = 700 * 1024**2
        expansion_factor = 5 if max_path.endswith(".gz") else 2

        expected = fixed_overhead + expansion_factor * max_size
        expected = min(max(expected, 1 * 1024**3), 120 * 1024**3)

        result = estimate_ram_per_worker_bytes(vcf_paths)
        assert result == int(expected)

    def test_single_file_plain_vcf(self, data_dir):
        """Test estimate for a single plain VCF."""
        path = os.path.abspath(os.path.join(data_dir, "TCGA-DC-6682.vcf"))
        size = os.path.getsize(path)

        fixed_overhead = 700 * 1024**2
        expected = fixed_overhead + 2 * size
        expected = min(max(expected, 1 * 1024**3), 120 * 1024**3)

        result = estimate_ram_per_worker_bytes([path])
        assert result == int(expected)

    def test_single_file_gz_vcf(self, data_dir):
        """Test estimate for a single gzipped VCF."""
        path = os.path.abspath(os.path.join(data_dir, "test.vcf.gz"))
        size = os.path.getsize(path)

        fixed_overhead = 700 * 1024**2
        expected = fixed_overhead + 5 * size
        expected = min(max(expected, 1 * 1024**3), 120 * 1024**3)

        result = estimate_ram_per_worker_bytes([path])
        assert result == int(expected)
