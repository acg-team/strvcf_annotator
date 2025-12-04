import hashlib
import os

import pytest

from strvcf_annotator import STRAnnotator

# Get base directory for test data
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def file_hash(path):
    """Calculate MD5 hash of a file."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


@pytest.fixture(
    params=[
        ("test.vcf.gz", "d5354b559173a69e7045a68bb3e1b6f3"),
        ("pindel_header.vcf", "9bd195a201d6b3317645ce5d44d40a2e"),
        ("mutec2_indel.vcf.gz", "6bd0bdf9d034b63e89ef35725a07b1cc"),
        ("TCGA-DC-6682.vcf", "b7889a2db9ce89adfc49cf8c11c4b7ba"),
    ]
)
def vcf_case(request, data_dir):
    input_vcf = os.path.abspath(os.path.join(data_dir, request.param[0]))
    expected_hash = request.param[1]
    return input_vcf, expected_hash, request


class TestProcessVcf:
    def test_process_vcf(self, vcf_case, output_dir):
        input_vcf, expected_hash, request = vcf_case
        update_hashes = request.config.getoption("--update-vcf-hashes")

        str_bed = os.path.abspath(os.path.join(base_dir, "data", "GRCh38_repeats.bed"))
        output_filename = (
            os.path.basename(input_vcf).replace(".vcf.gz", "").replace(".vcf", "")
            + ".processed.vcf"
        )
        output_path = os.path.abspath(os.path.join(base_dir, "output", output_filename))
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # Use new API
        annotator = STRAnnotator(str_bed)
        annotator.annotate_vcf_file(input_vcf, output_path)

        actual_hash = file_hash(output_path)

        if update_hashes:
            print(f"New hash for {input_vcf}: {actual_hash}")
            # optionally save to a JSON/yaml file here
        else:
            assert actual_hash == expected_hash, f"{input_vcf} hash mismatch: {actual_hash}"

    def test_reannotating_test_vcf_is_idempotent(self, data_dir, output_dir):
        """Re-annotating an already annotated VCF should produce the same file."""
        # Use only test.vcf.gz for idempotence check
        input_vcf = os.path.abspath(os.path.join(data_dir, "test.vcf.gz"))
        str_bed = os.path.abspath(os.path.join(base_dir, "data", "GRCh38_repeats.bed"))

        # First annotation
        first_output = os.path.abspath(os.path.join(output_dir, "test.processed.vcf"))
        # Second annotation (annotate result again)
        second_output = os.path.abspath(os.path.join(output_dir, "test.reannotated.vcf"))

        os.makedirs(os.path.dirname(first_output), exist_ok=True)

        annotator = STRAnnotator(str_bed)
        annotator.annotate_vcf_file(input_vcf, first_output)
        annotator.annotate_vcf_file(first_output, second_output)

        first_hash = file_hash(first_output)
        second_hash = file_hash(second_output)

        assert first_hash == second_hash, (
            f"Re-annotating {input_vcf} is not idempotent: {first_hash} != {second_hash}"
        )
