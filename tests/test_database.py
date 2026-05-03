from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from genotype_converter.database import (
    discover_source_folders,
    import_lookup,
    init_database,
    list_assemblies,
    list_manifests,
    list_species,
    load_lookup_table_from_database,
    query_marker,
)


def test_import_lookup_tracks_context_and_queries_markers(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"

    init_database(str(db_path))
    stats = import_lookup(
        database_path=str(db_path),
        lookup_path=str(lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="tiny_manifest",
    )

    assert stats.rows_imported == 8
    assert list_species(str(db_path)) == ["bos_taurus"]
    assert list_assemblies(str(db_path), "bos_taurus") == ["ARS_UCD_v2_0"]
    manifests = list_manifests(str(db_path), species="bos_taurus")
    assert manifests[0]["manifest_name"] == "tiny_manifest"

    snp_rows = query_marker(
        database_path=str(db_path),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        marker_name="SNP2",
    )
    assert len(snp_rows) == 1
    assert snp_rows[0]["manifest_name"] == "tiny_manifest"
    assert snp_rows[0]["A_in_TOP"] == "A"
    assert snp_rows[0]["A_in_PLUS"] == "T"

    indel_rows = query_marker(
        database_path=str(db_path),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        marker_name="INDEL1",
    )
    assert len(indel_rows) == 1
    assert indel_rows[0]["A_in_TOP"] == "I"
    assert indel_rows[0]["B_in_TOP"] == "D"


def test_same_marker_name_can_exist_in_multiple_manifests(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"

    import_lookup(
        database_path=str(db_path),
        lookup_path=str(lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="manifest_a",
    )
    import_lookup(
        database_path=str(db_path),
        lookup_path=str(lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="manifest_b",
    )

    rows = query_marker(
        database_path=str(db_path),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        marker_name="SNP2",
    )
    assert [row["manifest_name"] for row in rows] == ["manifest_a", "manifest_b"]

    filtered = query_marker(
        database_path=str(db_path),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        marker_name="SNP2",
        manifest_name="manifest_b",
    )
    assert len(filtered) == 1
    assert filtered[0]["manifest_name"] == "manifest_b"


def test_load_lookup_table_from_database_matches_converter_shape(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    import_lookup(
        database_path=str(db_path),
        lookup_path=str(lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="tiny_manifest",
    )

    table = load_lookup_table_from_database(
        database_path=str(db_path),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="tiny_manifest",
    )

    assert table["SNP2"][0]["TOP"] == "A"
    assert table["SNP2"][0]["PLUS"] == "T"
    assert table["SNP2"][1]["TOP"] == "C"
    assert table["SNP2"][1]["PLUS"] == "G"
    assert table["INDEL1"][0]["TOP"] == "I"
    assert table["INDEL1"][1]["TOP"] == "D"


def test_load_lookup_table_from_database_rejects_ambiguous_manifest_context(
    pipeline_output, tmp_path
):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    import_lookup(
        database_path=str(db_path),
        lookup_path=str(lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="tiny_manifest",
    )

    modified_lookup = tmp_path / "modified.lookup.csv"
    text = lookup.read_text()
    modified_lookup.write_text(text.replace("SNP1-0_T_F_1511658221", "SNP1-alt", 1))
    import_lookup(
        database_path=str(db_path),
        lookup_path=str(modified_lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="tiny_manifest",
    )

    with pytest.raises(ValueError, match="Multiple lookup sources"):
        load_lookup_table_from_database(
            database_path=str(db_path),
            species="bos_taurus",
            assembly="ARS_UCD_v2_0",
            manifest_name="tiny_manifest",
        )


def test_import_lookup_requires_replace_for_same_source(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    kwargs = {
        "database_path": str(db_path),
        "lookup_path": str(lookup),
        "species": "bos_taurus",
        "assembly": "ARS_UCD_v2_0",
        "manifest_name": "tiny_manifest",
    }

    import_lookup(**kwargs)
    with pytest.raises(ValueError, match="already exists"):
        import_lookup(**kwargs)

    stats = import_lookup(**kwargs, replace=True)
    assert stats.rows_replaced == 1
    assert stats.rows_imported == 8

    with sqlite3.connect(db_path) as conn:
        source_count = conn.execute("SELECT COUNT(*) FROM lookup_sources").fetchone()[0]
        rule_count = conn.execute("SELECT COUNT(*) FROM marker_rules").fetchone()[0]
    assert source_count == 1
    assert rule_count == 8


def test_database_source_fixture_contains_expected_layout():
    root = Path("tests/data/database_sources/bos_taurus/ARS_UCD_v2_0")
    assert (root / "manifests/tiny_bovine_manifest.csv").is_file()
    assert (root / "references/tiny_reference.fa").is_file()
    assert (root / "genotypes/tiny_top_wide.csv").is_file()


def test_discover_source_folders_reads_tiny_fixture():
    folders = discover_source_folders("tests/data/database_sources")

    assert len(folders) == 1
    folder = folders[0]
    assert folder.species == "bos_taurus"
    assert folder.assembly == "ARS_UCD_v2_0"
    assert [Path(path).name for path in folder.manifest_paths] == [
        "tiny_bovine_manifest.csv"
    ]
    assert [Path(path).name for path in folder.reference_paths] == [
        "tiny_reference.fa"
    ]
