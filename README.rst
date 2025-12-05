================
strvcf_annotator
================


STR (Short Tandem Repeat) annotation tool for VCF files.

`strvcf_annotator` is a Python library and CLI tool for annotating 
variants in VCF files that overlap short tandem repeat (STR) regions. 
The tool converts SNPs, SNVs and indels into full repeat sequences and adds STR metadata.

Installation
------------
Package is available through `PyPI <https://pypi.org/project/strvcf-annotator/>`_. To install, type:

.. code-block:: shell

    pip install strvcf-annotator
    
.. code-block:: shell

    # Install from source
    git clone https://github.com/acg-team/strvcf_annotator.git
    cd strvcf_annotator
    pip install -e .

.. code-block:: shell

    # Dev dependencies
    pip install -r requirements_dev.txt


Quick Start
-----------

Command Line
~~~~~~~~~~~~

.. code-block:: shell

    # Annotate a single VCF
    strvcf-annotator --input input.vcf --str-bed repeats.bed --output output.vcf

    # Batch-process a directory
    strvcf-annotator --input-dir vcf_files/ --str-bed repeats.bed --output-dir annotated/

    # With verbose logging
    strvcf-annotator --input input.vcf --str-bed repeats.bed --output output.vcf --verbose

Library Usage
~~~~~~~~~~~~~

.. code-block:: python

    from strvcf_annotator import STRAnnotator

    # Create the annotator
    annotator = STRAnnotator('repeats.bed')

    # Annotate a single file
    annotator.annotate_vcf_file('input.vcf', 'output.vcf')

    # Batch processing
    annotator.process_directory('vcf_files/', 'annotated/')

    # Streaming processing
    import pysam
    vcf_in = pysam.VariantFile('input.vcf')
    for record in annotator.annotate_vcf_stream(vcf_in):
        print(f"Repeat unit: {record.info['RU']}")

Input format
-------------

BED file with STR regions
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    CHROM   START   END     PERIOD  RU
    chr1    100     115     3       CAG
    chr1    200     212     4       ATCG
    chr2    300     318     3       GAT

- **CHROM**: Chromosome name
- **START**: Start position (0-based, BED format)
- **END**: End position (0-based, exclusive)
- **PERIOD**: Repeat unit length
- **RU**: Repeat unit sequence

VCF file
~~~~~~~~~~~~

A standard VCF with variants. Must contain:

- FORMAT field **GT** (genotype)
- Optional: **AD** (allelic depth), **DP** (total depth)

Output format
-------------

The annotated VCF contains additional fields:

INFO fields
~~~~~~~~~~~~

- **RU**: Repeat unit
- **PERIOD**: Repeat period (unit length)
- **REF**: Reference copy number
- **PERFECT**: TRUE if both alleles are perfect repeats

FORMAT fields
~~~~~~~~~~~~~

- **REPCN**: Genotype expressed as repeat copy numbers

Example
~~~~~~~~~~~~

.. code-block:: shell

    ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit">
    ##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Repeat period">
    ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference copy number">
    ##INFO=<ID=PERFECT,Number=1,Type=String,Description="Perfect repeat indicator">
    ##FORMAT=<ID=REPCN,Number=2,Type=Integer,Description="Repeat copy number">

    #CHROM  POS  ID  REF         ALT             QUAL  FILTER  INFO                              FORMAT      Sample1
    chr1    101  .   CAGCAGCAG   CAGCAGCAGCAG    .     .       RU=CAG;PERIOD=3;REF=3;PERFECT=TRUE  GT:REPCN    0/1:3,4

Extending functionality
-----------------------

Creating a custom parser
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    from strvcf_annotator.parsers.base import BaseVCFParser

    class CustomParser(BaseVCFParser):
        def get_genotype(self, record, sample_idx):
            # Your logic for extracting the genotype
            pass
        
        def has_variant(self, record, sample_idx):
            # Your logic for determining if there is a variant
            pass
        
        def extract_info(self, record, sample_idx):
            # Your logic for extracting additional fields
            pass
        
        def validate_record(self, record):
            # Your logic for validating the record
            pass

    # Usage
    annotator = STRAnnotator('repeats.bed', parser=CustomParser())

Troubleshooting
---------------

Issue: ModuleNotFoundError
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    # Install the package in editable (dev) mode
    pip install -e .

Issue: Unnormalized VCF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tool **only accepts normalized VCFs**. Please normalize with `bcftools` before running.
Example (produces a normalized, indexed VCF):

.. code-block:: shell

   # Replace reference.fa with the exact reference used for the VCF
   bcftools norm -f reference.fa -m input.vcf


Issue: Unsorted VCF
~~~~~~~~~~~~~~~~~~~

The tool automatically sorts the VCF in memory, but for large files pre-sorting is recommended:

.. code-block:: shell

    bcftools sort input.vcf -o sorted.vcf


Issue: Reference mismatch
~~~~~~~~~~~~~~~~~~~~~~~~~

If you see warnings about a reference mismatch, check:

- The correctness of the STR BED file
- Matching reference genome versions

Documentation
-------------

* `API Documentation <docs/API.md>`_
* `Examples <examples/>`_

Contributing
------------

Contributions are welcome! 
For major changes, please open an issue first 
to discuss what you’d like to change.
Please ensure:

1. All tests pass
2. Code follows existing style
3. New features include tests
4. Documentation is updated

License
-------

MIT License

Credits
-------

Test bed files were taken from `ConSTRain` repository `https://github.com/acg-team/ConSTRain`. 
