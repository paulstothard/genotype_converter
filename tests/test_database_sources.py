from __future__ import annotations

import csv
from pathlib import Path


def test_database_source_fixture_layout_is_available():
    root = Path("tests/data/database_sources/bos_taurus/ARS_UCD_v2_0")

    assert (root / "manifests/tiny_bovine_manifest.csv").is_file()
    assert (root / "references/tiny_reference.fa").is_file()
    assert (root / "genotypes/tiny_top_wide.csv").is_file()
    assert (root / "genotypes/tiny_top_long.csv").is_file()
    assert (root / "expected/tiny_lookup_markers.txt").is_file()

    expected_markers = {
        line.strip()
        for line in (root / "expected/tiny_lookup_markers.txt").read_text().splitlines()
        if line.strip()
    }
    with open(root / "manifests/tiny_bovine_manifest.csv", newline="") as handle:
        manifest_markers = {row["Name"] for row in csv.DictReader(handle)}

    assert expected_markers == manifest_markers
