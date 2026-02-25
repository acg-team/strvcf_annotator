import os
import shutil
import zipfile
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
OUTPUT_DIR = Path(__file__).resolve().parents[0] / "output"


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


@pytest.fixture(scope="session")
def output_dir(request):
    """
    Provide a unique output directory per test under tests/output/.
    Ensures cleanup after the test.
    """
    test_name = os.environ.get('PYTEST_CURRENT_TEST').split(':')[-1].split(' ')[0]
    outdir = OUTPUT_DIR / f"{test_name}_dir"
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    yield outdir
    if outdir.exists():
        shutil.rmtree(outdir)


@pytest.fixture(scope="session")
def vcf_dir(data_dir: str, output_dir: str) -> str:
    """
    Unpack tests/data/vcfs/test_input.zip into tests/output/vcfs
    and return the directory containing VCF files.

    If the directory already exists (from a previous run), it is reused.
    """
    data_path = Path(data_dir)
    zip_path = data_path / "vcfs" / "test_input.zip"
    assert zip_path.is_file(), f"Missing test input zip: {zip_path}"

    vcf_root = Path(output_dir) / "vcfs"
    vcf_root.mkdir(parents=True, exist_ok=True)

    # If directory is empty, extract; otherwise assume it's already populated
    if not any(vcf_root.iterdir()):
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(vcf_root)

    inner = list(vcf_root.iterdir())
    if len(inner) == 1 and inner[0].is_dir():
        # Zip contained a single directory; use that
        return str(inner[0])

    # Otherwise, use the top-level extraction dir
    return str(vcf_root)
