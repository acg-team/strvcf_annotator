"""Unit tests for STR reference management."""

import os
import tempfile
from pathlib import Path

import pandas as pd
import pysam
import pytest

from strvcf_annotator.core.str_reference import (
    find_overlapping_str,
    get_str_at_position,
    is_valid_tabix,
    load_str_reference,
)


class TestLoadSTRReference:
    """Test suite for load_str_reference."""

    @pytest.fixture
    def temp_bed_file(self):
        """Create temporary BED file."""
        content = """chr1\t100\t115\t3\tCAG
chr2\t200\t212\t4\tATCG
chr1\t300\t318\t3\tGAT"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(content)
            temp_path = f.name

        yield temp_path

        # Cleanup
        Path(temp_path).unlink()

    def test_creates_gz_and_tbi(self, temp_bed_file):
        """Test that load_str_reference creates .gz and .tbi."""
        gz_path = load_str_reference(temp_bed_file)

        assert gz_path.endswith(".gz")
        assert Path(gz_path).exists() and is_valid_tabix(gz_path), (
            "BGZF-compressed and indexed file should exist"
        )

    def test_reuses_existing_indexed_panel(self, temp_bed_file):
        """Test that load_str_reference reuses already created gz+tbi."""
        gz_path_1 = load_str_reference(temp_bed_file)
        mtime_gz_1 = os.path.getmtime(gz_path_1)
        mtime_tbi_1 = os.path.getmtime(gz_path_1 + ".tbi")

        gz_path_2 = load_str_reference(temp_bed_file)
        mtime_gz_2 = os.path.getmtime(gz_path_2)
        mtime_tbi_2 = os.path.getmtime(gz_path_2 + ".tbi")

        assert gz_path_1 == gz_path_2
        assert mtime_gz_1 == mtime_gz_2
        assert mtime_tbi_1 == mtime_tbi_2

    def test_output_is_sorted_for_tabix(self, temp_bed_file):
        """Test that output is sorted by chromosome order and start."""
        gz_path = load_str_reference(temp_bed_file)

        tbx = pysam.TabixFile(gz_path)
        rows = list(tbx.fetch())  # full file iteration
        tbx.close()

        parsed = []
        for line in rows:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            start = int(parts[1])
            parsed.append((chrom, start))

        # chr1 entries first
        chr1_starts = [s for c, s in parsed if c == "chr1"]
        chr2_starts = [s for c, s in parsed if c == "chr2"]

        assert len(chr1_starts) == 2
        assert len(chr2_starts) == 1
        assert chr1_starts == sorted(chr1_starts)

        # Ensure overall order doesn't put chr2 before chr1
        first_chr = parsed[0][0]
        assert first_chr == "chr1"

    def test_recreates_index_if_tbi_is_invalid(self, temp_bed_file):
        """Test that an invalid .tbi triggers rebuilding a valid tabix index."""
        gz_path = load_str_reference(temp_bed_file)
        tbi_path = gz_path + ".tbi"

        assert Path(gz_path).exists() and is_valid_tabix(gz_path), (
            "BGZF-compressed and indexed file should exist"
        )

        # Record times so we can detect replacement if the same path is reused
        old_gz = gz_path
        old_tbi_mtime = os.path.getmtime(tbi_path)

        # Corrupt the index file
        with open(tbi_path, "wb") as f:
            f.write(b"NOT_A_REAL_TABIX_INDEX")

        # Ensure filesystem mtime changes (some FS have 1s resolution)
        # time.sleep(1.1)

        # Now call load_str_reference on the gz itself
        rebuilt_gz_path = load_str_reference(old_gz)
        rebuilt_tbi_path = rebuilt_gz_path + ".tbi"

        assert Path(rebuilt_gz_path).exists() and is_valid_tabix(rebuilt_gz_path), (
            "BGZF-compressed and indexed file should exist"
        )

        new_tbi_mtime = os.path.getmtime(rebuilt_tbi_path)
        assert new_tbi_mtime > old_tbi_mtime, "Expected .tbi to be recreated"


class TestFindOverlappingSTR:
    """Test suite for find_overlapping_str (tabix-backed)."""

    @pytest.fixture
    def tabix_panel(self):
        """Create a small bgzip+tabix STR panel for overlap tests."""
        content = """chr1\t101\t115\t3\tCAG\t5
chr1\t201\t212\t4\tATCG\t3
chr2\t301\t318\t3\tGAT\t6
"""
        with tempfile.TemporaryDirectory() as tmp:
            bed_path = Path(tmp) / "panel.bed"
            bed_path.write_text(content, encoding="utf-8")

            gz_path = load_str_reference(str(bed_path))
            yield gz_path

    def test_exact_overlap(self, tabix_panel):
        """Test exact position overlap."""
        result = find_overlapping_str(tabix_panel, "chr1", 101, 115)

        assert result is not None
        assert result["START"] == 101
        assert result["END"] == 115
        assert result["RU"] == "CAG"

    def test_partial_overlap(self, tabix_panel):
        """Test partial overlap."""
        result = find_overlapping_str(tabix_panel, "chr1", 105, 110)

        assert result is not None
        assert result["START"] == 101
        assert result["RU"] == "CAG"

    def test_no_overlap(self, tabix_panel):
        """Test no overlap."""
        result = find_overlapping_str(tabix_panel, "chr1", 150, 160)
        assert result is None

    def test_wrong_chromosome(self, tabix_panel):
        """Test wrong chromosome."""
        result = find_overlapping_str(tabix_panel, "chr3", 101, 115)
        assert result is None

    def test_variant_extends_beyond(self, tabix_panel):
        """Test variant extending beyond STR region."""
        result = find_overlapping_str(tabix_panel, "chr1", 110, 120)

        assert result is not None
        assert result["START"] == 101
        assert result["END"] == 115


class TestGetSTRAtPosition:
    """Test suite for get_str_at_position (tabix-backed)."""

    @pytest.fixture
    def tabix_panel(self):
        """Create a small bgzip+tabix STR panel for position tests."""
        content = """chr1\t101\t115\t3\tCAG\t5
chr1\t201\t212\t4\tATCG\t3
chr2\t301\t318\t3\tGAT\t6
"""
        with tempfile.TemporaryDirectory() as tmp:
            bed_path = Path(tmp) / "panel.bed"
            bed_path.write_text(content, encoding="utf-8")
            gz_path = load_str_reference(str(bed_path))
            yield gz_path

    def test_position_in_str(self, tabix_panel):
        """Test position within STR."""
        result = get_str_at_position(tabix_panel, "chr1", 105)

        assert result is not None
        assert result["RU"] == "CAG"

    def test_position_outside_str(self, tabix_panel):
        """Test position outside STR."""
        result = get_str_at_position(tabix_panel, "chr1", 101)
        assert result is None

    def test_position_at_boundary(self, tabix_panel):
        """Test position at STR boundary."""
        result = get_str_at_position(tabix_panel, "chr1", 102)

        assert result is not None
        assert result["RU"] == "CAG"
