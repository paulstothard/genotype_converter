from __future__ import annotations

import csv
import sqlite3
from pathlib import Path

import pytest

from genotype_converter.database import (
    LOOKUP_COLUMNS,
    database_stats,
    discover_source_folders,
    duplicate_marker_report,
    infer_lookup_source_for_markers,
    import_lookup,
    init_database,
    list_import_warnings,
    list_assemblies,
    list_manifests,
    list_species,
    load_lookup_table_from_database,
    query_marker,
    remove_source,
    source_details,
    unresolved_marker_report,
    validate_database,
)


def _write_lookup_subset(source: Path, dest: Path, marker_names: set[str]) -> None:
    with source.open(newline="") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    reader = csv.DictReader(lines)
    assert reader.fieldnames is not None
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            if row["marker_name"] in marker_names:
                writer.writerow(row)


def _write_lookup_rows(dest: Path, rows: list[dict[str, str]]) -> None:
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LOOKUP_COLUMNS)
        writer.writeheader()
        for row in rows:
            complete_row = {column: "" for column in LOOKUP_COLUMNS}
            complete_row.update(row)
            writer.writerow(complete_row)


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


def test_infer_lookup_source_prefers_manifest_with_most_input_markers(
    pipeline_output, tmp_path
):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    small_lookup = tmp_path / "small.lookup.csv"
    large_lookup = tmp_path / "large.lookup.csv"
    _write_lookup_subset(lookup, small_lookup, {"SNP1", "SNP2"})
    _write_lookup_subset(lookup, large_lookup, {"SNP1", "SNP2", "SNP3", "INDEL1"})

    import_lookup(
        database_path=str(db_path),
        lookup_path=str(small_lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="small_panel",
    )
    import_lookup(
        database_path=str(db_path),
        lookup_path=str(large_lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="large_panel",
    )

    inferred = infer_lookup_source_for_markers(
        database_path=str(db_path),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        marker_names=["SNP1", "SNP2", "SNP3"],
    )

    assert inferred.manifest_name == "large_panel"
    assert inferred.markers_matched == 3
    assert inferred.markers_requested == 3


def test_infer_lookup_source_rejects_tied_manifest_matches(pipeline_output, tmp_path):
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

    with pytest.raises(ValueError, match="Could not infer a unique manifest"):
        infer_lookup_source_for_markers(
            database_path=str(db_path),
            species="bos_taurus",
            assembly="ARS_UCD_v2_0",
            marker_names=["SNP1", "SNP2"],
        )


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


def test_import_lookup_replace_context_removes_stale_sources(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    modified_lookup = tmp_path / "modified.lookup.csv"
    modified_lookup.write_text(
        lookup.read_text().replace("SNP1-0_T_F_1511658221", "SNP1-alt", 1)
    )
    kwargs = {
        "database_path": str(db_path),
        "species": "bos_taurus",
        "assembly": "ARS_UCD_v2_0",
        "manifest_name": "tiny_manifest",
    }

    import_lookup(lookup_path=str(lookup), **kwargs)
    import_lookup(lookup_path=str(modified_lookup), **kwargs)
    stats = import_lookup(
        lookup_path=str(lookup),
        replace_context=True,
        **kwargs,
    )

    assert stats.rows_replaced == 2
    assert len(list_manifests(str(db_path))) == 1
    assert database_stats(str(db_path))["marker_rules"] == 8


def test_import_lookup_collapses_identical_duplicate_marker_rows(tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = tmp_path / "duplicate_identical.lookup.csv"
    row = {
        "marker_name": "DUP1",
        "chromosome": "1",
        "position": "10",
        "A_in_TOP": "A",
        "B_in_TOP": "C",
    }
    _write_lookup_rows(lookup, [row, dict(row)])

    stats = import_lookup(
        database_path=str(db_path),
        lookup_path=str(lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="duplicate_panel",
    )

    assert stats.rows_imported == 1
    assert stats.rows_skipped_duplicate_marker == 0
    assert stats.duplicate_marker_names == ()
    rows = query_marker(
        database_path=str(db_path),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        marker_name="DUP1",
    )
    assert len(rows) == 1


def test_import_lookup_skips_conflicting_duplicate_marker_rows(tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = tmp_path / "duplicate_conflict.lookup.csv"
    _write_lookup_rows(
        lookup,
        [
            {
                "marker_name": "DUP1",
                "chromosome": "1",
                "position": "10",
                "A_in_TOP": "A",
                "B_in_TOP": "C",
            },
            {
                "marker_name": "DUP1",
                "chromosome": "1",
                "position": "20",
                "A_in_TOP": "G",
                "B_in_TOP": "T",
            },
            {
                "marker_name": "UNIQUE1",
                "chromosome": "2",
                "position": "30",
                "A_in_TOP": "A",
                "B_in_TOP": "G",
            },
        ],
    )

    stats = import_lookup(
        database_path=str(db_path),
        lookup_path=str(lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="duplicate_panel",
    )

    assert stats.rows_imported == 1
    assert stats.rows_skipped_duplicate_marker == 2
    assert stats.duplicate_marker_names == ("DUP1",)
    assert query_marker(
        database_path=str(db_path),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        marker_name="DUP1",
    ) == []
    assert len(query_marker(
        database_path=str(db_path),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        marker_name="UNIQUE1",
    )) == 1

    warnings = list_import_warnings(str(db_path))
    assert len(warnings) == 1
    assert warnings[0]["warning_type"] == "conflicting_duplicate_marker"
    assert warnings[0]["marker_name"] == "DUP1"


def test_database_stats_source_details_and_validation(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    stats = import_lookup(
        database_path=str(db_path),
        lookup_path=str(lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="tiny_manifest",
    )

    db_stats = database_stats(str(db_path))
    assert db_stats["lookup_sources"] == 1
    assert db_stats["marker_rules"] == 8
    assert db_stats["species"] == 1
    assert db_stats["assemblies"] == 1

    details = source_details(str(db_path), stats.source_id)
    assert details["manifest_name"] == "tiny_manifest"
    assert details["marker_rules"] == 8

    validation = validate_database(str(db_path))
    assert validation.status == "pass"
    assert {row["check"] for row in validation.checks} >= {
        "sqlite_integrity",
        "foreign_keys",
        "duplicate_source_contexts",
        "duplicate_rules_within_source",
    }


def test_remove_source_deletes_rules_and_warnings(tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = tmp_path / "duplicate_conflict.lookup.csv"
    _write_lookup_rows(
        lookup,
        [
            {"marker_name": "DUP1", "position": "10"},
            {"marker_name": "DUP1", "position": "20"},
            {"marker_name": "UNIQUE1", "position": "30"},
        ],
    )
    stats = import_lookup(
        database_path=str(db_path),
        lookup_path=str(lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="duplicate_panel",
    )

    removed = remove_source(str(db_path), source_id=stats.source_id)

    assert removed.sources_removed == 1
    assert removed.marker_rules_removed == 1
    assert database_stats(str(db_path))["lookup_sources"] == 0
    assert list_import_warnings(str(db_path)) == []


def test_remove_source_by_context_removes_multiple_stale_sources(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    modified_lookup = tmp_path / "modified.lookup.csv"
    modified_lookup.write_text(
        lookup.read_text().replace("SNP1-0_T_F_1511658221", "SNP1-alt", 1)
    )
    import_lookup(
        database_path=str(db_path),
        lookup_path=str(lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="tiny_manifest",
    )
    import_lookup(
        database_path=str(db_path),
        lookup_path=str(modified_lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="tiny_manifest",
    )

    removed = remove_source(
        str(db_path),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="tiny_manifest",
    )

    assert removed.sources_removed == 2
    assert removed.marker_rules_removed == 16
    assert list_manifests(str(db_path)) == []


def test_duplicate_and_unresolved_reports(pipeline_output, tmp_path):
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

    duplicates = duplicate_marker_report(str(db_path), limit=1)
    assert duplicates[0]["source_count"] == 2
    assert "manifest_a" in duplicates[0]["manifest_names"]
    assert "manifest_b" in duplicates[0]["manifest_names"]

    unresolved_lookup = tmp_path / "unresolved.lookup.csv"
    _write_lookup_rows(
        unresolved_lookup,
        [{"marker_name": "NO_POS", "A_in_TOP": "A", "B_in_TOP": "C"}],
    )
    import_lookup(
        database_path=str(db_path),
        lookup_path=str(unresolved_lookup),
        species="bos_taurus",
        assembly="ARS_UCD_v2_0",
        manifest_name="unresolved_panel",
    )
    unresolved = unresolved_marker_report(str(db_path), manifest_name="unresolved_panel")
    assert [row["marker_name"] for row in unresolved] == ["NO_POS"]


def test_database_source_fixture_contains_expected_layout():
    root = Path("tests/data/database_sources/bos_taurus")
    assert (root / "manifests/tiny_bovine_manifest.csv").is_file()
    assert (root / "references/ARS_UCD_v2_0/tiny_reference.fa").is_file()
    assert (root / "references/ARS_UCD1_2/tiny_reference.fa").is_file()
    assert (root / "genotypes/tiny_top_wide.csv").is_file()


def test_discover_source_folders_reads_tiny_fixture():
    folders = discover_source_folders("tests/data/database_sources")

    assert [folder.assembly for folder in folders] == ["ARS_UCD1_2", "ARS_UCD_v2_0"]
    for folder in folders:
        assert folder.species == "bos_taurus"
        assert Path(folder.root_path).name == "bos_taurus"
        assert Path(folder.manifest_root_path).name == "manifests"
        assert Path(folder.reference_root_path).name == folder.assembly
        assert [Path(path).name for path in folder.manifest_paths] == [
            "tiny_bovine_manifest.csv"
        ]
        assert [Path(path).name for path in folder.reference_paths] == [
            "tiny_reference.fa"
        ]
