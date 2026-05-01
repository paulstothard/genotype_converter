from __future__ import annotations

from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent / "data"
MANIFEST = DATA_DIR / "manifest.csv"
REFERENCE = DATA_DIR / "reference.fa"


@pytest.fixture(scope="session")
def manifest_path():
    return str(MANIFEST)


@pytest.fixture(scope="session")
def reference_path():
    return str(REFERENCE)


@pytest.fixture(scope="session")
def pipeline_output(tmp_path_factory, manifest_path, reference_path):
    """Run the full pipeline once and return the output directory."""
    from genotype_converter.pipeline import run
    out = tmp_path_factory.mktemp("pipeline_out")
    run(
        manifest_path=manifest_path,
        reference_path=reference_path,
        outdir=str(out),
        species="bos_taurus",
        workers=1,
        save_alignment=True,
    )
    return out / "bos_taurus" / "reference"
