from __future__ import annotations

import csv
import json
from pathlib import Path

from click.testing import CliRunner

from genotype_converter.cli import main


def test_batch_options_are_shown_in_help():
    runner = CliRunner()
    expected_flags = {
        "convert": [
            "--genotypes-dir", "--outdir", "--suffix", "--overwrite",
            "--resolve-mixed-manifests", "--on-ambiguous-marker",
        ],
        "convert-plink": [
            "--bfile-dir", "--outdir", "--suffix", "--overwrite",
            "--resolve-mixed-manifests", "--on-ambiguous-marker",
        ],
        "convert-pfile": [
            "--pfile-dir", "--outdir", "--suffix", "--overwrite",
            "--resolve-mixed-manifests", "--on-ambiguous-marker",
        ],
    }

    for command, flags in expected_flags.items():
        result = runner.invoke(main, [command, "--help"])
        assert result.exit_code == 0, result.output
        for flag in flags:
            assert flag in result.output


def test_db_options_are_shown_in_help():
    runner = CliRunner()
    result = runner.invoke(main, ["db", "--help"])
    assert result.exit_code == 0, result.output
    for command in [
        "init",
        "build",
        "import-lookup",
        "list-species",
        "list-assemblies",
        "marker",
        "discover-sources",
    ]:
        assert command in result.output


def _write_plink_files(prefix, bim_rows):
    prefix.with_suffix(".bed").write_bytes(bytes([0x6C, 0x1B, 0x01, 0x00]))
    prefix.with_suffix(".fam").write_text("F1 I1 0 0 1 -9\n")
    prefix.with_suffix(".bim").write_text(
        "".join("\t".join(row) + "\n" for row in bim_rows)
    )


def _write_pfile(prefix, pvar_lines):
    prefix.with_suffix(".pgen").write_bytes(b"fake-pgen")
    prefix.with_suffix(".psam").write_text("#IID\tSEX\nI1\t1\n")
    prefix.with_suffix(".pvar").write_text("".join(line + "\n" for line in pvar_lines))


def _write_lookup_subset(source, dest, marker_names, overrides=None):
    overrides = overrides or {}
    with source.open(newline="") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    reader = csv.DictReader(lines)
    assert reader.fieldnames is not None
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            if row["marker_name"] in marker_names:
                row.update(overrides.get(row["marker_name"], {}))
                writer.writerow(row)


def test_convert_command_accepts_genotypes_dir(pipeline_output, tmp_path):
    input_dir = tmp_path / "genotypes"
    input_dir.mkdir()
    (input_dir / "herd.csv").write_text("sample_id,SNP1,SNP2\nS1,A/G,A/C\n")
    out_dir = tmp_path / "converted"

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "convert",
            "--genotypes-dir",
            str(input_dir),
            "--lookup",
            str(pipeline_output / "manifest.reference.lookup.csv"),
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--layout",
            "wide",
            "--outdir",
            str(out_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Converted 1/1 genotype files" in result.output
    assert "Batch summary:" in result.output
    rows = list(csv.DictReader((out_dir / "herd.converted.csv").read_text().splitlines()))
    assert rows[0]["SNP2"] == "T/G"
    summary_rows = list(csv.DictReader((out_dir / "conversion_summary.csv").read_text().splitlines()))
    assert summary_rows[0]["genotypes_changed"] == "1"


def test_convert_command_reports_single_file_unknown_alleles(pipeline_output, tmp_path):
    input_csv = tmp_path / "genotypes.csv"
    output_csv = tmp_path / "converted.csv"
    input_csv.write_text("sample_id,SNP2\nS1,Z/C\n")
    runner = CliRunner()

    result = runner.invoke(
        main,
        [
            "convert",
            "--genotypes",
            str(input_csv),
            "--lookup",
            str(pipeline_output / "manifest.reference.lookup.csv"),
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--output",
            str(output_csv),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Unknown allele labels left unchanged: 1" in result.output
    rows = list(csv.DictReader(output_csv.read_text().splitlines()))
    assert rows[0]["SNP2"] == "Z/G"


def test_convert_command_accepts_database_for_single_csv(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    input_csv = tmp_path / "genotypes.csv"
    output_csv = tmp_path / "converted.csv"
    input_csv.write_text("sample_id,SNP1,SNP2\nS1,A/G,A/C\n")
    runner = CliRunner()

    import_result = runner.invoke(
        main,
        [
            "db",
            "import-lookup",
            "--database",
            str(db_path),
            "--lookup",
            str(lookup),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--manifest-name",
            "tiny_manifest",
        ],
    )
    assert import_result.exit_code == 0, import_result.output

    result = runner.invoke(
        main,
        [
            "convert",
            "--genotypes",
            str(input_csv),
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--manifest-name",
            "tiny_manifest",
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--output",
            str(output_csv),
        ],
    )

    assert result.exit_code == 0, result.output
    rows = list(csv.DictReader(output_csv.read_text().splitlines()))
    assert rows[0]["SNP2"] == "T/G"


def test_convert_command_accepts_database_for_batch_csv(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    input_dir = tmp_path / "genotypes"
    input_dir.mkdir()
    (input_dir / "herd.csv").write_text("sample_id,SNP1,SNP2\nS1,A/G,A/C\n")
    out_dir = tmp_path / "converted"
    runner = CliRunner()

    import_result = runner.invoke(
        main,
        [
            "db",
            "import-lookup",
            "--database",
            str(db_path),
            "--lookup",
            str(lookup),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--manifest-name",
            "tiny_manifest",
        ],
    )
    assert import_result.exit_code == 0, import_result.output

    result = runner.invoke(
        main,
        [
            "convert",
            "--genotypes-dir",
            str(input_dir),
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--manifest-name",
            "tiny_manifest",
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--outdir",
            str(out_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    rows = list(csv.DictReader((out_dir / "herd.converted.csv").read_text().splitlines()))
    assert rows[0]["SNP2"] == "T/G"
    assert (out_dir / "conversion_summary.csv").is_file()


def test_convert_database_without_manifest_name_reports_failed_inference(tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    input_csv = tmp_path / "genotypes.csv"
    output_csv = tmp_path / "converted.csv"
    input_csv.write_text("sample_id,SNP1\nS1,A/G\n")
    runner = CliRunner()
    init_result = runner.invoke(main, ["db", "init", "--database", str(db_path)])
    assert init_result.exit_code == 0, init_result.output

    result = runner.invoke(
        main,
        [
            "convert",
            "--genotypes",
            str(input_csv),
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--output",
            str(output_csv),
        ],
    )

    assert result.exit_code != 0
    assert "Could not infer a manifest" in result.output


def test_convert_database_infers_manifest_from_csv_markers(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    small_lookup = tmp_path / "small.lookup.csv"
    large_lookup = tmp_path / "large.lookup.csv"
    _write_lookup_subset(lookup, small_lookup, {"SNP1", "SNP2"})
    _write_lookup_subset(lookup, large_lookup, {"SNP1", "SNP2", "SNP3"})
    input_csv = tmp_path / "genotypes.csv"
    output_csv = tmp_path / "converted.csv"
    input_csv.write_text("sample_id,SNP1,SNP2,SNP3\nS1,A/G,A/C,A/C\n")
    runner = CliRunner()

    for manifest_name, path in [("small_panel", small_lookup), ("large_panel", large_lookup)]:
        import_result = runner.invoke(
            main,
            [
                "db",
                "import-lookup",
                "--database",
                str(db_path),
                "--lookup",
                str(path),
                "--species",
                "bos_taurus",
                "--assembly",
                "ARS_UCD_v2_0",
                "--manifest-name",
                manifest_name,
            ],
        )
        assert import_result.exit_code == 0, import_result.output

    result = runner.invoke(
        main,
        [
            "convert",
            "--genotypes",
            str(input_csv),
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--output",
            str(output_csv),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Inferred manifest 'large_panel'" in result.output
    rows = list(csv.DictReader(output_csv.read_text().splitlines()))
    assert rows[0]["SNP2"] == "T/G"


def test_convert_database_resolves_mixed_manifests_per_marker(
    pipeline_output, tmp_path
):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    panel_a = tmp_path / "panel_a.lookup.csv"
    panel_b = tmp_path / "panel_b.lookup.csv"
    report = tmp_path / "resolution.csv"
    _write_lookup_subset(lookup, panel_a, {"SNP1", "SNP2", "SNP4"})
    _write_lookup_subset(
        lookup,
        panel_b,
        {"SNP2", "SNP3"},
        overrides={"SNP2": {"A_in_PLUS": "C", "B_in_PLUS": "A"}},
    )
    input_csv = tmp_path / "genotypes.csv"
    output_csv = tmp_path / "converted.csv"
    input_csv.write_text("sample_id,SNP1,SNP2,SNP4\nS1,A/G,A/C,A/G\n")
    runner = CliRunner()

    for manifest_name, path in [("panel_a", panel_a), ("panel_b", panel_b)]:
        import_result = runner.invoke(
            main,
            [
                "db",
                "import-lookup",
                "--database",
                str(db_path),
                "--lookup",
                str(path),
                "--species",
                "bos_taurus",
                "--assembly",
                "ARS_UCD_v2_0",
                "--manifest-name",
                manifest_name,
            ],
        )
        assert import_result.exit_code == 0, import_result.output

    result = runner.invoke(
        main,
        [
            "convert",
            "--genotypes",
            str(input_csv),
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--resolve-mixed-manifests",
            "--resolution-report",
            str(report),
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--output",
            str(output_csv),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Resolved mixed manifests" in result.output
    rows = list(csv.DictReader(output_csv.read_text().splitlines()))
    assert rows[0]["SNP2"] == "T/G"
    report_rows = list(csv.DictReader(report.read_text().splitlines()))
    snp2 = next(row for row in report_rows if row["marker_name"] == "SNP2")
    assert snp2["selected_manifest_name"] == "panel_a"
    assert snp2["candidate_count"] == "2"


def test_convert_database_mixed_manifest_mode_fails_on_ambiguous_marker(
    pipeline_output, tmp_path
):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    panel_a = tmp_path / "panel_a.lookup.csv"
    panel_b = tmp_path / "panel_b.lookup.csv"
    _write_lookup_subset(lookup, panel_a, {"SNP1", "SNP2"})
    _write_lookup_subset(
        lookup,
        panel_b,
        {"SNP2", "SNP3"},
        overrides={"SNP2": {"A_in_PLUS": "C", "B_in_PLUS": "A"}},
    )
    input_csv = tmp_path / "genotypes.csv"
    output_csv = tmp_path / "converted.csv"
    input_csv.write_text("sample_id,SNP1,SNP2,SNP3\nS1,A/G,A/C,A/C\n")
    runner = CliRunner()

    for manifest_name, path in [("panel_a", panel_a), ("panel_b", panel_b)]:
        import_result = runner.invoke(
            main,
            [
                "db",
                "import-lookup",
                "--database",
                str(db_path),
                "--lookup",
                str(path),
                "--species",
                "bos_taurus",
                "--assembly",
                "ARS_UCD_v2_0",
                "--manifest-name",
                manifest_name,
            ],
        )
        assert import_result.exit_code == 0, import_result.output

    result = runner.invoke(
        main,
        [
            "convert",
            "--genotypes",
            str(input_csv),
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--resolve-mixed-manifests",
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--output",
            str(output_csv),
        ],
    )

    assert result.exit_code != 0
    assert "Ambiguous marker rule" in result.output


def test_convert_rejects_lookup_and_database_together(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    input_csv = tmp_path / "genotypes.csv"
    output_csv = tmp_path / "converted.csv"
    input_csv.write_text("sample_id,SNP1\nS1,A/G\n")
    runner = CliRunner()
    init_result = runner.invoke(main, ["db", "init", "--database", str(db_path)])
    assert init_result.exit_code == 0, init_result.output

    result = runner.invoke(
        main,
        [
            "convert",
            "--genotypes",
            str(input_csv),
            "--lookup",
            str(pipeline_output / "manifest.reference.lookup.csv"),
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--manifest-name",
            "tiny_manifest",
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--output",
            str(output_csv),
        ],
    )

    assert result.exit_code != 0
    assert "exactly one of --lookup or --database" in result.output


def test_db_cli_imports_and_queries_lookup(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    runner = CliRunner()

    init_result = runner.invoke(main, ["db", "init", "--database", str(db_path)])
    assert init_result.exit_code == 0, init_result.output

    import_result = runner.invoke(
        main,
        [
            "db",
            "import-lookup",
            "--database",
            str(db_path),
            "--lookup",
            str(lookup),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--manifest-name",
            "tiny_manifest",
        ],
    )
    assert import_result.exit_code == 0, import_result.output
    assert "Imported 8 marker rules" in import_result.output

    species_result = runner.invoke(main, ["db", "list-species", "--database", str(db_path)])
    assert species_result.exit_code == 0, species_result.output
    assert species_result.output.strip() == "bos_taurus"

    marker_result = runner.invoke(
        main,
        [
            "db",
            "marker",
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--marker",
            "SNP2",
        ],
    )
    assert marker_result.exit_code == 0, marker_result.output
    assert "tiny_manifest" in marker_result.output
    assert "\tSNP2\t" in marker_result.output

    json_marker_result = runner.invoke(
        main,
        [
            "db",
            "marker",
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--marker",
            "SNP2",
            "--format",
            "json",
        ],
    )
    assert json_marker_result.exit_code == 0, json_marker_result.output
    marker_rows = json.loads(json_marker_result.output)
    assert marker_rows[0]["marker_name"] == "SNP2"
    assert marker_rows[0]["manifest_name"] == "tiny_manifest"


def test_db_cli_discovers_source_folders():
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "db",
            "discover-sources",
            "--source-root",
            "tests/data/database_sources",
            "--format",
            "csv",
        ],
    )

    assert result.exit_code == 0, result.output
    rows = list(csv.DictReader(result.output.splitlines()))
    assert rows[0]["species"] == "bos_taurus"
    assert rows[0]["assembly"] == "ARS_UCD_v2_0"
    assert rows[0]["manifest_count"] == "1"
    assert rows[0]["reference_count"] == "1"


def test_db_cli_builds_from_source_root(tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    build_outdir = tmp_path / "build"
    runner = CliRunner()

    result = runner.invoke(
        main,
        [
            "db",
            "build",
            "--source-root",
            "tests/data/database_sources",
            "--database",
            str(db_path),
            "--build-outdir",
            str(build_outdir),
            "--workers",
            "1",
            "--no-progress",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Database build complete: 1 lookup source(s), 8 marker rules imported." in result.output

    marker_result = runner.invoke(
        main,
        [
            "db",
            "marker",
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--manifest-name",
            "tiny_bovine_manifest",
            "--marker",
            "SNP2",
            "--format",
            "json",
        ],
    )
    assert marker_result.exit_code == 0, marker_result.output
    marker_rows = json.loads(marker_result.output)
    assert marker_rows[0]["A_in_PLUS"] == "T"


def test_db_cli_build_requires_one_reference(tmp_path):
    source_root = tmp_path / "sources"
    source_dir = source_root / "bos_taurus" / "ARS_UCD_v2_0"
    manifests_dir = source_dir / "manifests"
    references_dir = source_dir / "references"
    manifests_dir.mkdir(parents=True)
    references_dir.mkdir(parents=True)
    fixture_root = Path("tests/data/database_sources/bos_taurus/ARS_UCD_v2_0")
    (manifests_dir / "tiny.csv").write_text(
        (fixture_root / "manifests/tiny_bovine_manifest.csv").read_text()
    )
    reference_text = (fixture_root / "references/tiny_reference.fa").read_text()
    (references_dir / "ref_a.fa").write_text(reference_text)
    (references_dir / "ref_b.fa").write_text(reference_text)
    runner = CliRunner()

    result = runner.invoke(
        main,
        [
            "db",
            "build",
            "--source-root",
            str(source_root),
            "--database",
            str(tmp_path / "conversion.sqlite"),
            "--build-outdir",
            str(tmp_path / "build"),
            "--no-progress",
        ],
    )

    assert result.exit_code != 0
    assert "Expected exactly one reference file" in result.output


def test_convert_plink_command_accepts_bfile_dir(pipeline_output, tmp_path):
    input_dir = tmp_path / "plink"
    input_dir.mkdir()
    _write_plink_files(
        input_dir / "herd",
        [["1", "SNP2", "0", "700", "A", "C"]],
    )
    out_dir = tmp_path / "converted"

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "convert-plink",
            "--bfile-dir",
            str(input_dir),
            "--lookup",
            str(pipeline_output / "manifest.reference.lookup.csv"),
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--outdir",
            str(out_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Converted 1/1 PLINK filesets" in result.output
    assert "Batch summary:" in result.output
    assert (out_dir / "herd.converted.bim").read_text().strip() == "1\tSNP2\t0\t700\tT\tG"
    summary_rows = list(csv.DictReader((out_dir / "conversion_summary.csv").read_text().splitlines()))
    assert summary_rows[0]["alleles_changed"] == "2"


def test_convert_pfile_command_accepts_pfile_dir(pipeline_output, tmp_path):
    input_dir = tmp_path / "pfiles"
    input_dir.mkdir()
    _write_pfile(
        input_dir / "herd",
        ["#CHROM\tPOS\tID\tREF\tALT", "1\t700\tSNP2\tA\tC"],
    )
    out_dir = tmp_path / "converted"

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "convert-pfile",
            "--pfile-dir",
            str(input_dir),
            "--lookup",
            str(pipeline_output / "manifest.reference.lookup.csv"),
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--outdir",
            str(out_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Converted 1/1 PLINK 2 filesets" in result.output
    assert "Batch summary:" in result.output
    assert (out_dir / "herd.converted.pvar").read_text().splitlines() == [
        "#CHROM\tPOS\tID\tREF\tALT",
        "1\t700\tSNP2\tT\tG",
    ]
    summary_rows = list(csv.DictReader((out_dir / "conversion_summary.csv").read_text().splitlines()))
    assert summary_rows[0]["alleles_changed"] == "2"


def test_convert_plink_command_accepts_database_for_single_bfile(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    input_prefix = tmp_path / "plink" / "herd"
    input_prefix.parent.mkdir()
    output_prefix = tmp_path / "converted" / "herd_plus"
    _write_plink_files(
        input_prefix,
        [["1", "SNP2", "0", "700", "A", "C"]],
    )
    runner = CliRunner()

    import_result = runner.invoke(
        main,
        [
            "db",
            "import-lookup",
            "--database",
            str(db_path),
            "--lookup",
            str(lookup),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--manifest-name",
            "tiny_manifest",
        ],
    )
    assert import_result.exit_code == 0, import_result.output

    result = runner.invoke(
        main,
        [
            "convert-plink",
            "--bfile",
            str(input_prefix),
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--manifest-name",
            "tiny_manifest",
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--out",
            str(output_prefix),
        ],
    )

    assert result.exit_code == 0, result.output
    assert output_prefix.with_suffix(".bim").read_text().strip() == "1\tSNP2\t0\t700\tT\tG"


def test_convert_plink_command_accepts_database_for_bfile_dir(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    input_dir = tmp_path / "plink"
    input_dir.mkdir()
    _write_plink_files(
        input_dir / "herd",
        [["1", "SNP2", "0", "700", "A", "C"]],
    )
    out_dir = tmp_path / "converted"
    runner = CliRunner()

    import_result = runner.invoke(
        main,
        [
            "db",
            "import-lookup",
            "--database",
            str(db_path),
            "--lookup",
            str(lookup),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--manifest-name",
            "tiny_manifest",
        ],
    )
    assert import_result.exit_code == 0, import_result.output

    result = runner.invoke(
        main,
        [
            "convert-plink",
            "--bfile-dir",
            str(input_dir),
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--manifest-name",
            "tiny_manifest",
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--outdir",
            str(out_dir),
        ],
    )

    assert result.exit_code == 0, result.output
    assert (out_dir / "herd.converted.bim").read_text().strip() == "1\tSNP2\t0\t700\tT\tG"
    assert (out_dir / "conversion_summary.csv").is_file()


def test_convert_pfile_command_accepts_database_for_single_pfile(pipeline_output, tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    input_prefix = tmp_path / "pfiles" / "herd"
    input_prefix.parent.mkdir()
    output_prefix = tmp_path / "converted" / "herd_plus"
    _write_pfile(
        input_prefix,
        ["#CHROM\tPOS\tID\tREF\tALT", "1\t700\tSNP2\tA\tC"],
    )
    runner = CliRunner()

    import_result = runner.invoke(
        main,
        [
            "db",
            "import-lookup",
            "--database",
            str(db_path),
            "--lookup",
            str(lookup),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--manifest-name",
            "tiny_manifest",
        ],
    )
    assert import_result.exit_code == 0, import_result.output

    result = runner.invoke(
        main,
        [
            "convert-pfile",
            "--pfile",
            str(input_prefix),
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--manifest-name",
            "tiny_manifest",
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--out",
            str(output_prefix),
        ],
    )

    assert result.exit_code == 0, result.output
    assert output_prefix.with_suffix(".pvar").read_text().splitlines() == [
        "#CHROM\tPOS\tID\tREF\tALT",
        "1\t700\tSNP2\tT\tG",
    ]


def test_convert_plink_database_without_manifest_name_reports_failed_inference(tmp_path):
    db_path = tmp_path / "conversion.sqlite"
    input_prefix = tmp_path / "plink" / "herd"
    input_prefix.parent.mkdir()
    output_prefix = tmp_path / "converted" / "herd_plus"
    _write_plink_files(
        input_prefix,
        [["1", "SNP2", "0", "700", "A", "C"]],
    )
    runner = CliRunner()
    init_result = runner.invoke(main, ["db", "init", "--database", str(db_path)])
    assert init_result.exit_code == 0, init_result.output

    result = runner.invoke(
        main,
        [
            "convert-plink",
            "--bfile",
            str(input_prefix),
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--out",
            str(output_prefix),
        ],
    )

    assert result.exit_code != 0
    assert "Could not infer a manifest" in result.output


def test_convert_plink_database_infers_manifest_from_bim_markers(
    pipeline_output, tmp_path
):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    small_lookup = tmp_path / "small.lookup.csv"
    large_lookup = tmp_path / "large.lookup.csv"
    _write_lookup_subset(lookup, small_lookup, {"SNP1", "SNP2"})
    _write_lookup_subset(lookup, large_lookup, {"SNP1", "SNP2", "SNP3"})
    input_prefix = tmp_path / "plink" / "herd"
    input_prefix.parent.mkdir()
    output_prefix = tmp_path / "converted" / "herd_plus"
    _write_plink_files(
        input_prefix,
        [
            ["1", "SNP1", "0", "300", "A", "G"],
            ["1", "SNP2", "0", "700", "A", "C"],
            ["1", "SNP3", "0", "900", "A", "C"],
        ],
    )
    runner = CliRunner()

    for manifest_name, path in [("small_panel", small_lookup), ("large_panel", large_lookup)]:
        import_result = runner.invoke(
            main,
            [
                "db",
                "import-lookup",
                "--database",
                str(db_path),
                "--lookup",
                str(path),
                "--species",
                "bos_taurus",
                "--assembly",
                "ARS_UCD_v2_0",
                "--manifest-name",
                manifest_name,
            ],
        )
        assert import_result.exit_code == 0, import_result.output

    result = runner.invoke(
        main,
        [
            "convert-plink",
            "--bfile",
            str(input_prefix),
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--out",
            str(output_prefix),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Inferred manifest 'large_panel'" in result.output
    assert "SNP2\t0\t700\tT\tG" in output_prefix.with_suffix(".bim").read_text()


def test_convert_pfile_database_infers_manifest_from_pvar_markers(
    pipeline_output, tmp_path
):
    db_path = tmp_path / "conversion.sqlite"
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    small_lookup = tmp_path / "small.lookup.csv"
    large_lookup = tmp_path / "large.lookup.csv"
    _write_lookup_subset(lookup, small_lookup, {"SNP1", "SNP2"})
    _write_lookup_subset(lookup, large_lookup, {"SNP1", "SNP2", "SNP3"})
    input_prefix = tmp_path / "pfiles" / "herd"
    input_prefix.parent.mkdir()
    output_prefix = tmp_path / "converted" / "herd_plus"
    _write_pfile(
        input_prefix,
        [
            "#CHROM\tPOS\tID\tREF\tALT",
            "1\t300\tSNP1\tA\tG",
            "1\t700\tSNP2\tA\tC",
            "1\t900\tSNP3\tA\tC",
        ],
    )
    runner = CliRunner()

    for manifest_name, path in [("small_panel", small_lookup), ("large_panel", large_lookup)]:
        import_result = runner.invoke(
            main,
            [
                "db",
                "import-lookup",
                "--database",
                str(db_path),
                "--lookup",
                str(path),
                "--species",
                "bos_taurus",
                "--assembly",
                "ARS_UCD_v2_0",
                "--manifest-name",
                manifest_name,
            ],
        )
        assert import_result.exit_code == 0, import_result.output

    result = runner.invoke(
        main,
        [
            "convert-pfile",
            "--pfile",
            str(input_prefix),
            "--database",
            str(db_path),
            "--species",
            "bos_taurus",
            "--assembly",
            "ARS_UCD_v2_0",
            "--from-format",
            "TOP",
            "--to-format",
            "PLUS",
            "--out",
            str(output_prefix),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Inferred manifest 'large_panel'" in result.output
    assert "1\t700\tSNP2\tT\tG" in output_prefix.with_suffix(".pvar").read_text()
