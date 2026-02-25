"""Core modules for STR annotation functionality."""

from .annotation import build_new_record, make_modified_header, should_skip_genotype
from .repeat_utils import apply_variant_to_repeat, count_repeat_units, extract_repeat_sequence
from .str_reference import load_str_reference
from .vcf_processor import (
    annotate_vcf_to_file,
    check_vcf_sorted,
    generate_annotated_records,
    reset_and_sort_vcf,
)

__all__ = [
    'load_str_reference',
    'extract_repeat_sequence',
    'count_repeat_units',
    'apply_variant_to_repeat',
    'make_modified_header',
    'build_new_record',
    'should_skip_genotype',
    'check_vcf_sorted',
    'reset_and_sort_vcf',
    'generate_annotated_records',
    'annotate_vcf_to_file'
]
