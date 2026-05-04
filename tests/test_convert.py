from __future__ import annotations

import csv
from pathlib import Path

import pytest

from genotype_converter.convert_genotypes import (
    _split_genotype,
    convert_affymetrix_matrix,
    convert_batch,
    convert_long,
    convert_illumina_long,
    convert_illumina_matrix,
    convert_wide,
    load_lookup_table,
)


@pytest.fixture(scope="session")
def lookup_path(pipeline_output):
    return str(pipeline_output / "manifest.reference.lookup.csv")


@pytest.fixture(scope="session")
def table(lookup_path):
    return load_lookup_table(lookup_path)


def test_lookup_table_loaded(table):
    assert "SNP1" in table
    assert "SNP5" in table
    assert len(table["SNP1"]) == 2


def test_convert_allele_snp1_top_to_plus(table):
    from genotype_converter.convert_genotypes import _convert_allele
    # SNP1 plus strand, both TOP and PLUS are the same
    assert _convert_allele("A", "SNP1", "TOP", "PLUS", table) == "A"
    assert _convert_allele("G", "SNP1", "TOP", "PLUS", table) == "G"


def test_convert_allele_snp2_top_to_forward(table):
    # SNP2: TOP_A=A → FORWARD_A=T; TOP_B=C → FORWARD_B=G
    from genotype_converter.convert_genotypes import _convert_allele
    assert _convert_allele("A", "SNP2", "TOP", "FORWARD", table) == "T"
    assert _convert_allele("C", "SNP2", "TOP", "FORWARD", table) == "G"


def test_convert_allele_unknown_passthrough(table):
    from genotype_converter.convert_genotypes import _convert_allele
    assert _convert_allele("X", "SNP1", "TOP", "PLUS", table) == "X"


def test_convert_allele_missing_target_value_passthrough():
    from genotype_converter.convert_genotypes import _convert_allele, _convert_pair
    table = {
        "SNP_NO_PLUS": [
            {"TOP": "A", "PLUS": ""},
            {"TOP": "C", "PLUS": ""},
        ]
    }

    assert _convert_allele("A", "SNP_NO_PLUS", "TOP", "PLUS", table) == "A"
    converted, changed, unknown = _convert_pair(
        ("A", "C"), "SNP_NO_PLUS", "TOP", "PLUS", table
    )
    assert converted == ("A", "C")
    assert changed == 0
    assert unknown == 2


def test_convert_wide_top_to_plus(table, tmp_path):
    input_data = "sample_id,SNP1,SNP2\nS1,A/G,A/C\nS2,A/A,C/C\n"
    in_file = tmp_path / "in.csv"
    in_file.write_text(input_data)
    out_file = tmp_path / "out.csv"

    stats = convert_wide(
        input_path=str(in_file),
        output_path=str(out_file),
        table=table,
        from_fmt="TOP",
        to_fmt="PLUS",
        in_sep="/",
        out_sep="/",
        sample_col="sample_id",
    )

    rows = list(csv.DictReader(out_file.read_text().splitlines()))
    # SNP1: same in TOP and PLUS
    assert rows[0]["SNP1"] == "A/G"
    # SNP2: TOP_A=A→PLUS_A=T, TOP_B=C→PLUS_B=G
    assert rows[0]["SNP2"] == "T/G"
    assert rows[1]["SNP2"] == "G/G"
    assert stats.genotype_cells_total == 4
    assert stats.genotypes_changed == 2
    assert stats.alleles_changed == 4
    assert stats.unknown_alleles == 0


def test_convert_wide_missing_passthrough(table, tmp_path):
    input_data = "sample_id,SNP1\nS1,0/0\nS2,NA/NA\n"
    in_file = tmp_path / "in.csv"
    in_file.write_text(input_data)
    out_file = tmp_path / "out.csv"

    convert_wide(
        input_path=str(in_file),
        output_path=str(out_file),
        table=table,
        from_fmt="TOP",
        to_fmt="PLUS",
        in_sep="/",
        out_sep="/",
        sample_col="sample_id",
    )
    rows = list(csv.DictReader(out_file.read_text().splitlines()))
    assert rows[0]["SNP1"] == "0/0"
    assert rows[1]["SNP1"] == "NA/NA"


def test_convert_wide_excludes_unconvertible_marker_by_default(table, tmp_path):
    input_data = "sample_id,SNP1,SNP2,INDEL1\nS1,AB,AB,ID\nS2,BB,AA,DD\n"
    in_file = tmp_path / "ab_in.csv"
    in_file.write_text(input_data)
    out_file = tmp_path / "plus_out.csv"

    stats = convert_wide(
        input_path=str(in_file),
        output_path=str(out_file),
        table=table,
        from_fmt="AB",
        to_fmt="PLUS",
        in_sep=None,
        out_sep="/",
        sample_col="sample_id",
    )

    rows = list(csv.DictReader(out_file.read_text().splitlines()))
    assert rows[0]["SNP1"] == "A/G"
    assert rows[0]["SNP2"] == "T/G"
    assert rows[1]["SNP2"] == "T/T"
    assert "INDEL1" not in rows[0]
    assert stats.markers_excluded == 1
    assert stats.genotypes_excluded == 2
    assert stats.markers_incomplete_mapping == 1


def test_convert_wide_can_keep_unconvertible_marker_for_audit(table, tmp_path):
    input_data = "sample_id,SNP1,SNP2,INDEL1\nS1,AB,AB,ID\nS2,BB,AA,DD\n"
    in_file = tmp_path / "ab_in.csv"
    in_file.write_text(input_data)
    out_file = tmp_path / "plus_out.csv"

    stats = convert_wide(
        input_path=str(in_file),
        output_path=str(out_file),
        table=table,
        from_fmt="AB",
        to_fmt="PLUS",
        in_sep=None,
        out_sep="/",
        sample_col="sample_id",
        on_unconvertible_marker="keep",
    )

    rows = list(csv.DictReader(out_file.read_text().splitlines()))
    assert rows[0]["INDEL1"] == "I/D"
    assert rows[1]["INDEL1"] == "D/D"
    assert stats.markers_excluded == 0
    assert stats.markers_incomplete_mapping == 1
    assert stats.unknown_alleles == 4


def test_split_genotype_auto_detects_space_and_strips_missing():
    assert _split_genotype("A G", None) == ("A", "G")
    assert _split_genotype(" NA ", None) is None


def test_convert_long(table, tmp_path):
    input_data = "sample_id,marker_name,genotype\nS1,SNP2,A/C\nS1,SNP3,A/C\n"
    in_file = tmp_path / "long_in.csv"
    in_file.write_text(input_data)
    out_file = tmp_path / "long_out.csv"

    stats = convert_long(
        input_path=str(in_file),
        output_path=str(out_file),
        table=table,
        from_fmt="TOP",
        to_fmt="FORWARD",
        in_sep="/",
        out_sep="/",
        sample_col="sample_id",
        marker_col="marker_name",
        genotype_col="genotype",
    )
    rows = list(csv.DictReader(out_file.read_text().splitlines()))
    # SNP2: TOP_A=A→FORWARD_A=T; TOP_B=C→FORWARD_B=G
    assert rows[0]["genotype"] == "T/G"
    # SNP3: same logic
    assert rows[1]["genotype"] == "T/G"
    assert stats.genotype_cells_total == 2
    assert stats.genotypes_changed == 2
    assert stats.alleles_changed == 4


def test_convert_long_excludes_unconvertible_marker_and_reports(table, tmp_path):
    input_data = (
        "sample_id,marker_name,genotype\n"
        "S1,SNP2,A/C\n"
        "S2,SNP2,A/C\n"
        "S1,NOT_IN_LOOKUP,A/G\n"
    )
    in_file = tmp_path / "long_in.csv"
    in_file.write_text(input_data)
    out_file = tmp_path / "long_out.csv"

    stats = convert_long(
        input_path=str(in_file),
        output_path=str(out_file),
        table=table,
        from_fmt="TOP",
        to_fmt="PLUS",
        in_sep="/",
        out_sep="/",
        sample_col="sample_id",
        marker_col="marker_name",
        genotype_col="genotype",
    )

    rows = list(csv.DictReader(out_file.read_text().splitlines()))
    assert [row["marker_name"] for row in rows] == ["SNP2", "SNP2"]
    assert stats.markers_missing_lookup == 1
    assert stats.markers_excluded == 1
    assert stats.genotypes_excluded == 1
    report_rows = list(csv.DictReader(Path(stats.marker_report_path).read_text().splitlines()))
    assert report_rows[0]["marker_name"] == "NOT_IN_LOOKUP"
    assert report_rows[0]["reason"] == "missing_lookup"
    assert report_rows[0]["action"] == "exclude"
    assert Path(stats.exclude_marker_path).read_text().strip() == "NOT_IN_LOOKUP"


def test_convert_long_strictly_rejects_incomplete_target_mapping(tmp_path):
    input_data = "sample_id,marker_name,genotype\nS1,SNP_NO_PLUS,A/C\n"
    in_file = tmp_path / "long_in.csv"
    in_file.write_text(input_data)
    out_file = tmp_path / "long_out.csv"
    table = {
        "SNP_NO_PLUS": [
            {"TOP": "A", "PLUS": ""},
            {"TOP": "C", "PLUS": ""},
        ]
    }

    with pytest.raises(ValueError, match="SNP_NO_PLUS"):
        convert_long(
            input_path=str(in_file),
            output_path=str(out_file),
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
            in_sep="/",
            out_sep="/",
            sample_col="sample_id",
            marker_col="marker_name",
            genotype_col="genotype",
            on_unconvertible_marker="fail",
        )


def test_convert_long_vcf_output(table, tmp_path):
    input_data = "sample_id,marker_name,genotype\nS1,SNP2,A/C\nS1,INDEL1,I/D\n"
    in_file = tmp_path / "long_in.csv"
    in_file.write_text(input_data)
    out_file = tmp_path / "long_vcf.csv"

    convert_long(
        input_path=str(in_file),
        output_path=str(out_file),
        table=table,
        from_fmt="TOP",
        to_fmt="VCF",
        in_sep="/",
        out_sep="/",
        sample_col="sample_id",
        marker_col="marker_name",
        genotype_col="genotype",
    )
    rows = list(csv.DictReader(out_file.read_text().splitlines()))
    assert rows[0]["genotype"] == "REF/ALT"
    assert rows[1]["genotype"] == "ALT/REF"


def test_convert_illumina_matrix_ab_to_plus(table, tmp_path):
    input_data = "\n".join(
        [
            "[Header]",
            "GSGT Version\t2.0.4",
            "[Data]",
            "\tS1\tS2",
            "SNP1\tAB\tBB",
            "SNP2\tAA\tAB",
        ]
    ) + "\n"
    in_file = tmp_path / "gsgt_matrix.txt"
    out_file = tmp_path / "gsgt_matrix_plus.txt"
    in_file.write_text(input_data)

    stats = convert_illumina_matrix(
        input_path=str(in_file),
        output_path=str(out_file),
        table=table,
        from_fmt="AB",
        to_fmt="PLUS",
        in_sep=None,
        out_sep="",
    )

    lines = out_file.read_text().splitlines()
    assert lines[:4] == ["[Header]", "GSGT Version\t2.0.4", "[Data]", "\tS1\tS2"]
    assert lines[4] == "SNP1\tAG\tGG"
    assert lines[5] == "SNP2\tTT\tTG"
    assert stats.markers_total == 2
    assert stats.genotype_cells_total == 4


def test_convert_illumina_long_fills_target_columns(table, tmp_path):
    input_data = "\n".join(
        [
            "[Header]",
            "GSGT Version\t2.0.4",
            "[Data]",
            "SNP Name\tSample ID\tAllele1 - Top\tAllele2 - Top\tAllele1 - Forward\tAllele2 - Forward\tAllele1 - AB\tAllele2 - AB\tAllele1 - Design\tAllele2 - Design\tAllele1 - Plus\tAllele2 - Plus\tGC Score",
            "SNP2\tS1\tA\tC\t-\t-\tA\tB\t-\t-\t-\t-\t0.99",
            "SNP1\tS1\tA\tG\t-\t-\tA\tB\t-\t-\t-\t-\t0.99",
        ]
    ) + "\n"
    in_file = tmp_path / "gsgt_long.txt"
    out_file = tmp_path / "gsgt_long_plus.txt"
    in_file.write_text(input_data)

    stats = convert_illumina_long(
        input_path=str(in_file),
        output_path=str(out_file),
        table=table,
        from_fmt="TOP",
        to_fmt="PLUS",
    )

    lines = out_file.read_text().splitlines()
    assert lines[4].split("\t")[10:12] == ["T", "G"]
    assert lines[5].split("\t")[10:12] == ["A", "G"]
    assert stats.genotypes_parsed == 2
    assert stats.alleles_changed == 2


def test_convert_affymetrix_matrix_ab_to_plus(table, tmp_path):
    input_data = "\n".join(
        [
            "probeset_id\tS1\tS1\tS2\tS2",
            "SNP1\tAB\t--\tBB\t--",
            "SNP2\tAA\t--\tNoCall\t---",
        ]
    ) + "\n"
    in_file = tmp_path / "affy.txt"
    out_file = tmp_path / "affy_plus.txt"
    in_file.write_text(input_data)

    stats = convert_affymetrix_matrix(
        input_path=str(in_file),
        output_path=str(out_file),
        table=table,
        from_fmt="AB",
        to_fmt="PLUS",
    )

    lines = out_file.read_text().splitlines()
    assert lines[1] == "SNP1\tAB\tAG\tBB\tGG"
    assert lines[2] == "SNP2\tAA\tTT\tNoCall\t---"
    assert stats.genotype_cells_total == 4
    assert stats.genotypes_parsed == 3


def test_convert_batch_wide_top_to_plus(table, tmp_path):
    input_dir = tmp_path / "genotypes"
    input_dir.mkdir()
    (input_dir / "herd_a.csv").write_text("sample_id,SNP1,SNP2\nS1,A/G,A/C\n")
    (input_dir / "herd_b.csv").write_text("sample_id,SNP1,SNP2\nS2,G/G,C/C\n")
    (input_dir / "notes.txt").write_text("not a genotype file\n")
    out_dir = tmp_path / "converted"

    stats = convert_batch(
        input_dir=str(input_dir),
        output_dir=str(out_dir),
        pattern="*.csv",
        suffix=".plus.csv",
        overwrite=False,
        table=table,
        from_fmt="TOP",
        to_fmt="PLUS",
        layout="wide",
        in_sep="/",
        out_sep="/",
        sample_col="sample_id",
        marker_col="marker_name",
        genotype_col="genotype",
    )

    assert stats.files_total == 2
    assert stats.files_converted == 2
    assert [Path(path).name for path in stats.output_paths] == [
        "herd_a.plus.csv",
        "herd_b.plus.csv",
    ]
    assert Path(stats.summary_path).name == "conversion_summary.csv"
    rows_a = list(csv.DictReader((out_dir / "herd_a.plus.csv").read_text().splitlines()))
    rows_b = list(csv.DictReader((out_dir / "herd_b.plus.csv").read_text().splitlines()))
    summary_rows = list(csv.DictReader((out_dir / "conversion_summary.csv").read_text().splitlines()))
    assert rows_a[0]["SNP1"] == "A/G"
    assert rows_a[0]["SNP2"] == "T/G"
    assert rows_b[0]["SNP2"] == "G/G"
    assert [Path(row["output_path"]).name for row in summary_rows] == [
        "herd_a.plus.csv",
        "herd_b.plus.csv",
    ]
    assert summary_rows[0]["genotype_cells_total"] == "2"
    assert summary_rows[0]["genotypes_changed"] == "1"


def test_convert_batch_refuses_overwrite(table, tmp_path):
    input_dir = tmp_path / "genotypes"
    input_dir.mkdir()
    (input_dir / "herd.csv").write_text("sample_id,SNP1\nS1,A/G\n")
    out_dir = tmp_path / "converted"
    out_dir.mkdir()
    (out_dir / "herd.converted.csv").write_text("existing\n")

    with pytest.raises(FileExistsError, match="already exists"):
        convert_batch(
            input_dir=str(input_dir),
            output_dir=str(out_dir),
            pattern="*.csv",
            suffix=".converted.csv",
            overwrite=False,
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
            layout="wide",
            in_sep="/",
            out_sep="/",
            sample_col="sample_id",
            marker_col="marker_name",
            genotype_col="genotype",
        )


def test_convert_batch_requires_matching_files(table, tmp_path):
    input_dir = tmp_path / "genotypes"
    input_dir.mkdir()

    with pytest.raises(ValueError, match="No genotype files matched"):
        convert_batch(
            input_dir=str(input_dir),
            output_dir=str(tmp_path / "converted"),
            pattern="*.csv",
            suffix=".converted.csv",
            overwrite=False,
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
            layout="wide",
            in_sep="/",
            out_sep="/",
            sample_col="sample_id",
            marker_col="marker_name",
            genotype_col="genotype",
        )


def test_convert_wide_requires_sample_column(table, tmp_path):
    in_file = tmp_path / "in.csv"
    in_file.write_text("id,SNP1\nS1,A/G\n")
    out_file = tmp_path / "out.csv"

    with pytest.raises(ValueError, match="sample_id"):
        convert_wide(
            input_path=str(in_file),
            output_path=str(out_file),
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
            in_sep="/",
            out_sep="/",
            sample_col="sample_id",
        )


def test_convert_wide_strictly_rejects_unknown_marker(table, tmp_path):
    in_file = tmp_path / "in.csv"
    in_file.write_text("sample_id,SNP1,NOT_IN_LOOKUP\nS1,A/G,A/G\n")
    out_file = tmp_path / "out.csv"

    with pytest.raises(ValueError, match="NOT_IN_LOOKUP"):
        convert_wide(
            input_path=str(in_file),
            output_path=str(out_file),
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
            in_sep="/",
            out_sep="/",
            sample_col="sample_id",
            on_unconvertible_marker="fail",
        )


def test_convert_wide_excludes_unknown_marker_by_default(table, tmp_path):
    in_file = tmp_path / "in.csv"
    in_file.write_text("sample_id,SNP1,NOT_IN_LOOKUP\nS1,A/G,A/G\n")
    out_file = tmp_path / "out.csv"

    stats = convert_wide(
        input_path=str(in_file),
        output_path=str(out_file),
        table=table,
        from_fmt="TOP",
        to_fmt="PLUS",
        in_sep="/",
        out_sep="/",
        sample_col="sample_id",
    )

    rows = list(csv.DictReader(out_file.read_text().splitlines()))
    assert list(rows[0]) == ["sample_id", "SNP1"]
    assert stats.markers_missing_lookup == 1
    assert stats.markers_excluded == 1


def test_convert_long_requires_marker_and_genotype_columns(table, tmp_path):
    in_file = tmp_path / "long_in.csv"
    in_file.write_text("sample_id,marker,call\nS1,SNP1,A/G\n")
    out_file = tmp_path / "long_out.csv"

    with pytest.raises(ValueError, match="marker_name"):
        convert_long(
            input_path=str(in_file),
            output_path=str(out_file),
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
            in_sep="/",
            out_sep="/",
            sample_col="sample_id",
            marker_col="marker_name",
            genotype_col="genotype",
        )


def test_lookup_csv_has_correct_columns(pipeline_output):
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    lines = [l for l in lookup.read_text().splitlines() if not l.startswith("#")]
    reader = csv.DictReader(lines)
    row = next(reader)
    for col in ["marker_name", "chromosome", "position", "ref_allele",
                "determination_type",
                "A_in_TOP", "B_in_TOP", "A_in_PLUS", "B_in_PLUS"]:
        assert col in row, f"Missing column: {col}"


def test_lookup_snp1_values(pipeline_output):
    lookup = pipeline_output / "manifest.reference.lookup.csv"
    lines = [l for l in lookup.read_text().splitlines() if not l.startswith("#")]
    rows = {r["marker_name"]: r for r in csv.DictReader(lines)}
    assert rows["SNP1"]["chromosome"] == "1"
    assert rows["SNP1"]["position"] == "300"
    assert rows["SNP1"]["ref_allele"] == "A"
    assert rows["SNP1"]["determination_type"]
    assert rows["SNP1"]["A_in_AB"] == "A"
    assert rows["SNP1"]["B_in_AB"] == "B"
    assert rows["SNP1"]["A_in_TOP"] == "A"
    assert rows["SNP1"]["B_in_TOP"] == "G"
