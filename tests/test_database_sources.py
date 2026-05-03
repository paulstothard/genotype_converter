from __future__ import annotations

import csv
from pathlib import Path


def test_database_source_fixture_layout_is_available():
    root = Path("tests/data/database_sources/bos_taurus")

    assert (root / "manifests/tiny_bovine_manifest.csv").is_file()
    assert (root / "references/ARS_UCD_v2_0/tiny_reference.fa").is_file()
    assert (root / "references/ARS_UCD1_2/tiny_reference.fa").is_file()
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


def test_mixed_manifest_example_files_are_available():
    root = Path("examples/mixed_manifests")

    assert (root / "README.md").is_file()
    assert (root / "panel_a.lookup.csv").is_file()
    assert (root / "panel_b.lookup.csv").is_file()
    assert (root / "mixed_genotypes.csv").is_file()

    with (root / "panel_a.lookup.csv").open(newline="") as handle:
        panel_a_markers = {row["marker_name"] for row in csv.DictReader(handle)}
    with (root / "panel_b.lookup.csv").open(newline="") as handle:
        panel_b_markers = {row["marker_name"] for row in csv.DictReader(handle)}

    assert "SNP_SHARED" in panel_a_markers & panel_b_markers
