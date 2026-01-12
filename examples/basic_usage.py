"""Basic usage examples for strvcf_annotator."""

from strvcf_annotator import STRAnnotator, annotate_vcf


def example_1_simple_annotation():
    """Example 1: Simple single file annotation using convenience function."""
    print("Example 1: Simple annotation")
    print("-" * 50)

    # Simplest way to annotate a VCF file
    annotate_vcf(input_vcf="input.vcf", str_bed="str_regions.bed", output_vcf="output.vcf")

    print("✓ Annotated input.vcf -> output.vcf")
    print()


def example_2_class_based():
    """Example 2: Class-based approach for multiple operations."""
    print("Example 2: Class-based approach")
    print("-" * 50)

    # Create annotator instance (loads STR reference once)
    annotator = STRAnnotator("str_regions.bed")

    # Get statistics
    stats = annotator.get_statistics()
    print(f"Loaded {stats['total_regions']} STR regions")
    print(f"Covering {stats['chromosomes']} chromosomes")
    print()

    # Annotate multiple files efficiently
    annotator.annotate_vcf_file("sample1.vcf", "sample1.annotated.vcf")
    annotator.annotate_vcf_file("sample2.vcf", "sample2.annotated.vcf")

    print("✓ Annotated multiple files")
    print()


def example_3_batch_processing():
    """Example 3: Batch process entire directory."""
    print("Example 3: Batch processing")
    print("-" * 50)

    annotator = STRAnnotator("str_regions.bed")

    # Process all VCF files in directory
    annotator.process_directory(input_dir="vcf_files/", output_dir="annotated_vcfs/")

    print("✓ Processed all VCF files in directory")
    print()


def example_4_streaming():
    """Example 4: Stream processing for large files."""
    print("Example 4: Streaming processing")
    print("-" * 50)

    import pysam

    annotator = STRAnnotator("str_regions.bed")
    vcf_in = pysam.VariantFile("large_file.vcf")

    # Process records one at a time (memory efficient)
    for i, record in enumerate(annotator.annotate_vcf_stream(vcf_in), 1):
        if i <= 3:  # Show first 3 records
            print(f"Record {i}: {record.contig}:{record.pos}")
            print(f"  Repeat unit: {record.info['RU']}")
            print(f"  Copy number: {record.info['REF']}")

    vcf_in.close()
    print("✓ Streamed large file efficiently")
    print()


def example_5_query_regions():
    """Example 5: Query STR regions."""
    print("Example 5: Query STR regions")
    print("-" * 50)

    annotator = STRAnnotator("str_regions.bed")

    # Check if position is in STR
    str_region = annotator.get_str_at_position("chr1", 1000000)

    if str_region:
        print("Position chr1:1000000 is in STR:")
        print(f"  Repeat unit: {str_region['RU']}")
        print(f"  Period: {str_region['PERIOD']}")
        print(f"  Region: {str_region['START']}-{str_region['END']}")
    else:
        print("Position chr1:1000000 is not in an STR region")

    print()


def example_6_error_handling():
    """Example 6: Error handling."""
    print("Example 6: Error handling")
    print("-" * 50)

    from strvcf_annotator import ValidationError

    try:
        # Try to load non-existent file
        STRAnnotator("nonexistent.bed")
    except ValidationError as e:
        print(f"✓ Caught validation error: {e}")

    try:
        # Try to annotate invalid VCF
        annotate_vcf("invalid.vcf", "str_regions.bed", "output.vcf")
    except Exception as e:
        print(f"✓ Caught error: {type(e).__name__}")

    print()


def example_7_custom_parser():
    """Example 7: Using custom parser."""
    print("Example 7: Custom parser")
    print("-" * 50)

    from strvcf_annotator.parsers.generic import GenericParser

    # Use explicit parser
    parser = GenericParser()
    annotator = STRAnnotator("str_regions.bed", parser=parser)

    annotator.annotate_vcf_file("input.vcf", "output.vcf")

    print("✓ Used custom parser")
    print()


if __name__ == "__main__":
    print("=" * 50)
    print("STR VCF Annotator - Usage Examples")
    print("=" * 50)
    print()

    # Note: These examples assume files exist
    # In practice, you would use actual file paths

    print("Run individual examples by uncommenting them:")
    print()

    # Uncomment to run examples:
    # example_1_simple_annotation()
    # example_2_class_based()
    # example_3_batch_processing()
    # example_4_streaming()
    # example_5_query_regions()
    # example_6_error_handling()
    # example_7_custom_parser()

    print("See the source code for detailed examples!")
