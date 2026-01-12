import os
import shutil
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def data_dir():
    """Provides absolute path to the test data directory."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "data"))

def pytest_addoption(parser):
    parser.addoption(
        "--update-vcf-hashes",
        action="store_true",
        default=False,
        help="Recalculate expected hashes for VCF outputs",
    )


# Write all outputs here (committed folder, but files are NOT committed)
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "output"


@pytest.fixture(scope="session", autouse=True)
def _ensure_output_dir():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


@pytest.fixture
def output_path(request):
    """
    Provide a unique output path per test function under tests/output/.
    Ensures cleanup after the test.
    """
    out = OUTPUT_DIR / f"{request.node.name}.vcf"
    if out.exists():
        out.unlink()
    yield out
    if out.exists():
        out.unlink()


@pytest.fixture
def output_dir(request):
    """
    Provide a unique output directory per test under tests/output/.
    Ensures cleanup after the test.
    """
    outdir = OUTPUT_DIR / f"{request.node.name}_dir"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    yield outdir
    if outdir.exists():
        shutil.rmtree(outdir)
