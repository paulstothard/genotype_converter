from __future__ import annotations

import csv
from pathlib import Path

import pytest

from genotype_converter.convert_genotypes import (
    _split_genotype,
    convert_long,
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


def test_convert_wide_top_to_plus(table, tmp_path):
    input_data = "sample_id,SNP1,SNP2\nS1,A/G,A/C\nS2,A/A,C/C\n"
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
    # SNP1: same in TOP and PLUS
    assert rows[0]["SNP1"] == "A/G"
    # SNP2: TOP_A=A→PLUS_A=T, TOP_B=C→PLUS_B=G
    assert rows[0]["SNP2"] == "T/G"
    assert rows[1]["SNP2"] == "G/G"


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


def test_convert_wide_ab_adjacent_to_plus(table, tmp_path):
    input_data = "sample_id,SNP1,SNP2,INDEL1\nS1,AB,AB,ID\nS2,BB,AA,DD\n"
    in_file = tmp_path / "ab_in.csv"
    in_file.write_text(input_data)
    out_file = tmp_path / "plus_out.csv"

    convert_wide(
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
    assert rows[0]["INDEL1"] == "I/D"
    assert rows[1]["SNP2"] == "T/T"
    assert rows[1]["INDEL1"] == "D/D"


def test_split_genotype_auto_detects_space_and_strips_missing():
    assert _split_genotype("A G", None) == ("A", "G")
    assert _split_genotype(" NA ", None) is None


def test_convert_long(table, tmp_path):
    input_data = "sample_id,marker_name,genotype\nS1,SNP2,A/C\nS1,SNP3,A/C\n"
    in_file = tmp_path / "long_in.csv"
    in_file.write_text(input_data)
    out_file = tmp_path / "long_out.csv"

    convert_long(
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


def test_convert_wide_rejects_unknown_marker(table, tmp_path):
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
        )


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
