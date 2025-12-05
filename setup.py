"""Setup script for strvcf_annotator."""

import re
from pathlib import Path

from setuptools import find_packages, setup

with open("README.rst") as readme_file:
    readme = readme_file.read()

# Read version from __init__.py
init_file = Path(__file__).parent / "src" / "strvcf_annotator" / "__init__.py"
with open(init_file) as f:
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", f.read(), re.M)
    if version_match:
        version = version_match.group(1)
    else:
        raise RuntimeError("Unable to find version string.")

requirements = [
    "pysam>=0.22.0",
    "pandas>=2.0.0",
]

test_requirements = [
    "pytest>=8.0.0",
    "pytest-cov>=4.1.0",
]

setup(
    author="Olesia Kondrateva",
    author_email="xkdnoa@gmail.com",
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
    description="STR annotation tool for VCF files",
    entry_points={
        "console_scripts": [
            "strvcf-annotator=strvcf_annotator.cli:main",
        ],
    },
    install_requires=requirements,
    license="MIT",
    long_description=readme,
    include_package_data=True,
    keywords="strvcf_annotator",
    name="strvcf_annotator",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/acg-team/strvcf_annotator",
    version=version,
    zip_safe=False,
)
