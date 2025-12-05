"""Unit tests for repeat utilities."""

from strvcf_annotator.core.repeat_utils import (
    apply_variant_to_repeat,
    count_repeat_units,
    extract_repeat_sequence,
    is_perfect_repeat,
)


class TestExtractRepeatSequence:
    """Test suite for extract_repeat_sequence."""

    def test_simple_repeat(self):
        """Test simple repeat extraction."""
        str_row = {"RU": "CAG", "COUNT": 3}
        result = extract_repeat_sequence(str_row)
        assert result == "CAGCAGCAG"

    def test_single_repeat(self):
        """Test single repeat unit."""
        str_row = {"RU": "AT", "COUNT": 1}
        result = extract_repeat_sequence(str_row)
        assert result == "AT"

    def test_long_repeat(self):
        """Test longer repeat unit."""
        str_row = {"RU": "ATCG", "COUNT": 5}
        result = extract_repeat_sequence(str_row)
        assert result == "ATCGATCGATCGATCGATCG"


class TestCountRepeatUnits:
    """Test suite for count_repeat_units."""

    def test_perfect_repeat(self):
        """Test counting in perfect repeat."""
        assert count_repeat_units("CAGCAGCAG", "CAG") == 3

    def test_imperfect_repeat(self):
        """Test counting in imperfect repeat."""
        assert count_repeat_units("CAGCAGCA", "CAG") == 2

    def test_no_repeats(self):
        """Test sequence with no repeats."""
        assert count_repeat_units("ATCG", "CAG") == 0

    def test_single_repeat(self):
        """Test single repeat unit."""
        assert count_repeat_units("CAG", "CAG") == 1

    def test_overlapping_pattern(self):
        """Test non-overlapping counting."""
        # 'AAA' contains 'AA' once (non-overlapping)
        assert count_repeat_units("AAAA", "AA") == 2

    def test_empty_sequence(self):
        """Test empty sequence."""
        assert count_repeat_units("", "CAG") == 0

    def test_motif_longer_than_sequence(self):
        """Test motif longer than sequence."""
        assert count_repeat_units("CA", "CAG") == 0


class TestApplyVariantToRepeat:
    """Test suite for apply_variant_to_repeat."""

    def test_simple_substitution(self):
        """Test simple substitution."""
        result = apply_variant_to_repeat(100, "A", "T", 100, "AAAA")
        assert result == "TAAA"

    def test_deletion(self):
        """Test deletion."""
        result = apply_variant_to_repeat(100, "AA", "A", 100, "AAAA")
        assert result == "AAA"

    def test_insertion(self):
        """Test insertion."""
        result = apply_variant_to_repeat(100, "A", "AT", 100, "AAAA")
        assert result == "ATAAA"

    def test_variant_in_middle(self):
        """Test variant in middle of repeat."""
        result = apply_variant_to_repeat(102, "A", "T", 100, "AAAA")
        assert result == "AATA"

    def test_variant_before_repeat(self):
        """Test variant starting before repeat and overlapping its start."""
        result = apply_variant_to_repeat(98, "CCA", "T", 100, "AAAA")
        # Overlap is only at position 100 (the last 'A' of "CCA")
        # so effectively ALT='T' is applied at repeat position 0:
        # repeat_seq = "AAAA" -> "TAAA"
        assert result == "TAAA"

    def test_variant_at_end(self):
        """Test variant at end of repeat."""
        result = apply_variant_to_repeat(103, "A", "T", 100, "AAAA")
        assert result == "AAAT"

    def test_homopolymer_run_at_end(self):
        """AT AAAAA with unit 'A' should yield 5 (longest contiguous run)."""
        assert count_repeat_units("ATAAAAA", "A") == 5

    def test_homopolymer_snp_interrupts_run(self):
        """AAAATAAA with unit 'A' should yield 4 (run is broken by SNP T)."""
        assert count_repeat_units("AAAATAAA", "A") == 4

    def test_homopolymer_snp_outside_run(self):
        """
        SNP outside the main run should not change the run length.
        Reference: TTTAAAACTT (AAAA = 4)
        Alt example: CTTAAAACTT (SNP at position 1).
        """
        assert count_repeat_units("CTTAAAACTT", "A") == 4

    def test_homopolymer_multiple_runs_take_longest(self):
        """
        Longest contiguous run should be reported, not the sum of runs.
        AA T AAAAA T AA  -> runs: 2, 5, 2 -> expect 5.
        """
        assert count_repeat_units("AATAAAAATAA", "A") == 5

    def test_multibase_repeat_internal_snp_interrupts(self):
        """
        SNP inside one CAG unit breaks the contiguous run.

        Reference:  CAG CAG CAG  (3 repeats)
        Alt:        CAG CTG CAG  -> longest contiguous CAG run is 1.
        """
        assert count_repeat_units("CAGCTGCAG", "CAG") == 1

    def test_multibase_repeat_snp_between_units(self):
        """
        SNP between repeat units breaks contiguity of the run.

        Reference: CAG CAG CAG (3)
        Alt:       CAG T CAG CAG -> longest contiguous CAG run is 2 (last two).
        """
        assert count_repeat_units("CAGTCAGCAG", "CAG") == 2

    def test_multibase_repeat_complex_indel_like_pattern(self):
        """
        Complex allele mimicking indel + SNP from reference perspective.

        Sequence: CAG CAG T CAG
        Longest contiguous CAG run is 2 (first two CAGs).
        """
        assert count_repeat_units("CAGCAGTCAG", "CAG") == 2

    def test_homopolymer_deletion_shortens_run(self):
        """
        Deletion inside a homopolymer (from AAAA to AAA) corresponds
        to a shorter run: longest run should be 3.
        """
        assert count_repeat_units("AAA", "A") == 3

    def test_motif_not_aligned_to_frame(self):
        """
        Pattern that could confuse naive implementations that don't
        respect motif boundaries.

        Sequence:  ACAGCAGCA
        Motif:     CAG
        Valid repeats: CAG CAG (positions 2-7), then 'CA' is incomplete.
        Expect 2.
        """
        assert count_repeat_units("ACAGCAGCA", "CAG") == 2

    def test_alternating_motif_never_forms_repeat(self):
        """
        Alternating bases without contiguity.
        Sequence: ATATAT, motif A -> longest run of 'A' is 1.
        """
        assert count_repeat_units("ATATAT", "A") == 1


    def test_variant_entirely_before_repeat_no_overlap(self):
        """
        Variant ends before the repeat starts: STR sequence must be unchanged.
        repeat: positions 100..103 (len=4), variant at 95 with len(ref)=2 -> ends at 96.
        """
        repeat_seq = "AAAA"
        result_ref = apply_variant_to_repeat(
            pos=95,
            ref="AC",
            alt="TG",
            repeat_start=100,
            repeat_seq=repeat_seq,
        )
        assert result_ref == repeat_seq

    def test_variant_entirely_after_repeat_no_overlap(self):
        """
        Variant starts after the repeat ends: STR sequence must be unchanged.
        repeat: 100..103, variant starts at 105.
        """
        repeat_seq = "AAAA"
        result_ref = apply_variant_to_repeat(
            pos=105,
            ref="AC",
            alt="TG",
            repeat_start=100,
            repeat_seq=repeat_seq,
        )
        assert result_ref == repeat_seq

    def test_variant_extends_beyond_end_ref_and_alt_consistent(self):
        """
        Variant extends beyond the end of the STR:
        - applying REF as ALT must leave repeat_seq unchanged
        - applying real ALT must only modify overlapping part.

        repeat_start = 100
        repeat_seq   = TTTTTTTTTT (10 Ts, positions 100..109)
        variant pos  = 104 (inside STR), ref length 11 -> ends at 114 (> 109)
        """
        repeat_start = 100
        repeat_seq = "TTTTTTTTTT"  # len = 10

        pos = 104
        ref = "TTTTTTTTTTT"       # 11 Ts
        alt_same = ref            # same as ref, should NOT change STR
        alt_mut = "AAAAAAAAAAA"   # 11 As

        # 1) Applying REF as ALT: no change
        seq_ref = apply_variant_to_repeat(
            pos=pos,
            ref=ref,
            alt=alt_same,
            repeat_start=repeat_start,
            repeat_seq=repeat_seq,
        )
        assert seq_ref == repeat_seq

        # 2) Applying ALT: only overlapping part is changed
        seq_alt = apply_variant_to_repeat(
            pos=pos,
            ref=ref,
            alt=alt_mut,
            repeat_start=repeat_start,
            repeat_seq=repeat_seq,
        )

        # We only replace the overlapping part:
        # relative_pos = 104-100 = 4
        # overlap_len  = min(len(ref)=11, repeat_len-relative_pos=6) = 6
        # repeat_seq[0:4] = "TTTT"
        # alt_overlap[0:6] = "AAAAAA"
        # result = "TTTT" + "AAAAAA" = "TTTTAAAAAA"
        assert seq_alt == "TTTTAAAAAA"

    def test_lowercase_ref_alt_with_uppercase_repeat_identity(self):
        """
        Lowercase REF/ALT from VCF on an uppercase STR panel.

        Applying REF as ALT (identity) must leave the STR sequence unchanged,
        even though REF/ALT are lowercase.
        """
        repeat_start = 100
        repeat_seq = "AAAAAAAAAA"       # STR panel (uppercase)
        pos = 100
        ref = "aaaaaaaaaa"             # VCF gives lowercase
        alt = "aaaaaaaaaa"             # same as ref

        mutated = apply_variant_to_repeat(
            pos=pos,
            ref=ref,
            alt=alt,
            repeat_start=repeat_start,
            repeat_seq=repeat_seq,
        )
        assert mutated == repeat_seq   # no mismatch due to case

    def test_lowercase_repeat_and_alt_preserve_lowercase(self):
        """
        When the STR sequence is lowercase, ALT should also be applied in lowercase.
        """
        repeat_start = 100
        repeat_seq = "aaaaaaaaaa"      # lowercase STR
        pos = 100
        ref = "aaaaaaaaaa"
        alt = "tttttttttt"

        mutated = apply_variant_to_repeat(
            pos=pos,
            ref=ref,
            alt=alt,
            repeat_start=repeat_start,
            repeat_seq=repeat_seq,
        )
        assert mutated == "tttttttttt"  # all lowercase, matching STR case

    def test_apply_real_variant(self):
        pos = 45512562
        ref = "cacacacacacacacacacacacaca"
        alt = "cacacacacacacacacacacaca"
        mutated = apply_variant_to_repeat(
            pos=pos,
            ref=ref,
            alt=alt,
            repeat_start=pos,
            repeat_seq=ref,
        )
        assert mutated == alt
        # chr15	78938368	78938381	1	A
        # chr15	78939176	.	aaaaaaaaaaaaaaaa	aaaaaaaaaaaaaaaaa
        # It is a mismatch of panels
        repeat_start = 78938369
        repeat_seq = "AAAAAAAAAAAAA"
        pos = 78939176
        ref = "aaaaaaaaaaaaaaaa"
        alt = "aaaaaaaaaaaaaaaaa"
        mutated = apply_variant_to_repeat(
            pos=pos,
            ref=ref,
            alt=alt,
            repeat_start=repeat_start,
            repeat_seq=repeat_seq,
        )
        assert mutated == "AAAAAAAAAAAAA"

    def test_overlap_snp(self):
        pos = 98
        ref = "CCA"
        alt = "TGG"
        repeat_start = 100
        mutated = apply_variant_to_repeat(
            pos=pos,
            ref=ref,
            alt=alt,
            repeat_start=repeat_start,
            repeat_seq="AA",
        )
        assert mutated == "GA"

class TestIsPerfectRepeat:
    """Test suite for is_perfect_repeat."""

    def test_perfect_repeat(self):
        """Test perfect repeat detection."""
        assert is_perfect_repeat("CAGCAGCAG", "CAG") is True

    def test_imperfect_repeat(self):
        """Test imperfect repeat detection."""
        assert is_perfect_repeat("CAGCAGCA", "CAG") is False

    def test_single_unit(self):
        """Test single repeat unit."""
        assert is_perfect_repeat("CAG", "CAG") is True

    def test_no_repeat(self):
        """Test non-repeat sequence."""
        assert is_perfect_repeat("ATCG", "CAG") is False

    def test_empty_sequence(self):
        """Test empty sequence."""
        assert is_perfect_repeat("", "CAG") is False

    def test_empty_motif(self):
        """Test empty motif."""
        assert is_perfect_repeat("CAG", "") is False
