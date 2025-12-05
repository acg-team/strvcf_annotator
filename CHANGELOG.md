# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- Add support for multi-allelic variants
- Implement custom parser for additional VCF callers
- Add performance optimizations for large VCF files
- Integrate with cloud storage (S3, GCS)

## [0.2.2] - 2025-12-04
- Change `apply_variant_to_repeat` function for applying VCF variants to STR panel sequences.
  * Locally normalizes variants (`pos`, `ref`, `alt`) by trimming shared prefix/suffix before applying.
  * Applies the full **ALT** if the normalized variant starts inside the STR; only clips when the variant starts before the STR and overlaps it.
  * Handles REF/ALT in a **case-insensitive** way and matches the output case to the STR panel sequence.
  * Includes unit tests for:
    * variants overlapping STR boundaries,
    * long homopolymer variants extending beyond the STR,
    * lowercase REF/ALT vs uppercase/lowercase STR panel.
- Filter records with same alt and ref (happens if mutation is not normalized and actually happens outside of STR region)


## [0.2.1] - 2025-11-26
- Reformatted code for consistency and style.
- Added tox-based test matrix and GitHub Actions CI.
- Updated repeat unit counting logic to use the maximum length of uninterrupted motif runs, and added tests for complex cases.


## [0.2.0] - 2025-11-25
- Fixed error in detection of unsorted VCF files
- Fixed bug that caused annotation to stop at chromosome 9; annotation now runs through all chromosomes

## [0.1.0] - 2025-11-04

### Added
- Initial release of strvcf-annotator
- Complete modular architecture with clean separation of concerns
- **Core functionality:**
  - STR annotation for VCF files using BED reference
  - Support for indels and SNVs in repeat regions
  - Calculation of repeat copy numbers (REPCN)
  - Perfect repeat detection
  - Reconstruction of full repeat sequences

- **Parser system:**
  - BaseVCFParser abstract interface for extensibility
  - GenericParser implementation supporting standard VCF fields (GT, AD, DP)
  - Easy to extend with custom parsers

- **API (Library usage):**
  - STRAnnotator class for programmatic access
  - annotate_vcf() convenience function
  - Support for both single file and batch directory processing
  - Statistics and metadata extraction

- **CLI (Command-line usage):**
  - strvcf-annotator command-line tool
  - Single file mode: `--input` and `--output`
  - Batch directory mode: `--input-dir` and `--output-dir`
  - Verbose logging support with `--verbose`
  - Comprehensive argument validation

- **Core modules:**
  - `core/annotation.py` - VCF record annotation engine
  - `core/vcf_processor.py` - VCF file processing and workflow
  - `core/str_reference.py` - STR BED file loading and lookup
  - `core/repeat_utils.py` - Repeat sequence manipulation

- **Utilities:**
  - `utils/vcf_utils.py` - VCF helper functions
  - `utils/validation.py` - Input validation utilities

- **Testing:**
  - Comprehensive unit tests for all core modules
  - Integration tests for end-to-end workflows
  - CLI tests using subprocess
  - Performance benchmarks
  - Test coverage configuration

- **Documentation:**
  - Complete README with examples
  - API documentation
  - Migration guide from previous versions
  - Inline code documentation with docstrings

- **Packaging:**
  - pyproject.toml configuration
  - MIT License
  - PyPI-ready setup
  - Support for Python 3.8-3.12

### Changed
- **Breaking changes from v0.0.x:**
  - Removed hardcoded Pindel/constrain_dumb parser logic
  - Removed monolithic strvcf_annotator.py file
  - Refactored into modular architecture
  - New CLI argument structure (--input instead of positional args)
  - Library API changed to use STRAnnotator class

### Deprecated
- Old monolithic module structure (removed in this version)

### Removed
- Legacy Pindel-specific parser
- constrain_dumb filtering logic
- Hardcoded execution code in main module

### Fixed
- Reference mismatch handling with proper logging
- VCF sorting issues with automatic detection and sorting
- Genotype filtering for non-standard GT formats
- Header preservation during annotation
- FORMAT field copying in annotated records

### Security
- Input validation for all file paths
- Proper error handling to prevent crashes
- No execution of untrusted code

## [0.0.1] - 2025-10-25 (Legacy)

### Added
- Initial prototype implementation
- Basic STR annotation for Pindel VCF files

