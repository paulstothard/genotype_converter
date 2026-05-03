from __future__ import annotations

import csv

import pytest

from genotype_converter.convert_genotypes import load_lookup_table
from genotype_converter.convert_plink import (
    convert_plink_bfile,
    convert_plink_bfile_batch,
    convert_plink_pfile,
    convert_plink_pfile_batch,
)


@pytest.fixture(scope="session")
def table(pipeline_output):
    return load_lookup_table(str(pipeline_output / "manifest.reference.lookup.csv"))


def _write_plink_files(prefix, bim_rows):
    prefix.with_suffix(".bed").write_bytes(bytes([0x6C, 0x1B, 0x01, 0x00]))
    prefix.with_suffix(".fam").write_text("F1 I1 0 0 1 -9\n")
    prefix.with_suffix(".bim").write_text(
        "".join("\t".join(row) + "\n" for row in bim_rows)
    )


def test_convert_plink_bfile_rewrites_bim_and_copies_bed_fam(table, tmp_path):
    input_prefix = tmp_path / "input"
    output_prefix = tmp_path / "output"
    _write_plink_files(
        input_prefix,
        [
            ["1", "SNP1", "0", "300", "A", "G"],
            ["1", "SNP2", "0", "700", "A", "C"],
        ],
    )

    stats = convert_plink_bfile(
        input_prefix=str(input_prefix),
        output_prefix=str(output_prefix),
        table=table,
        from_fmt="TOP",
        to_fmt="PLUS",
    )

    assert stats.variants_total == 2
    assert stats.variants_converted == 2
    assert stats.alleles_changed == 2
    assert output_prefix.with_suffix(".bed").read_bytes() == input_prefix.with_suffix(".bed").read_bytes()
    assert output_prefix.with_suffix(".fam").read_text() == input_prefix.with_suffix(".fam").read_text()
    assert output_prefix.with_suffix(".bim").read_text().splitlines() == [
        "1\tSNP1\t0\t300\tA\tG",
        "1\tSNP2\t0\t700\tT\tG",
    ]


def test_convert_plink_bfile_can_update_positions(table, tmp_path):
    input_prefix = tmp_path / "input"
    output_prefix = tmp_path / "output"
    _write_plink_files(
        input_prefix,
        [["0", "SNP2", "0", "0", "A", "C"]],
    )

    convert_plink_bfile(
        input_prefix=str(input_prefix),
        output_prefix=str(output_prefix),
        table=table,
        from_fmt="TOP",
        to_fmt="PLUS",
        update_position=True,
    )

    assert output_prefix.with_suffix(".bim").read_text().strip() == "1\tSNP2\t0\t700\tT\tG"


def test_convert_plink_bfile_requires_all_three_files(table, tmp_path):
    input_prefix = tmp_path / "input"
    input_prefix.with_suffix(".bim").write_text("1\tSNP1\t0\t300\tA\tG\n")

    with pytest.raises(ValueError, match=".bed, .bim, and .fam"):
        convert_plink_bfile(
            input_prefix=str(input_prefix),
            output_prefix=str(tmp_path / "output"),
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
        )


def test_convert_plink_bfile_rejects_missing_lookup_marker(table, tmp_path):
    input_prefix = tmp_path / "input"
    _write_plink_files(
        input_prefix,
        [["1", "NOT_IN_LOOKUP", "0", "100", "A", "G"]],
    )

    with pytest.raises(ValueError, match="NOT_IN_LOOKUP"):
        convert_plink_bfile(
            input_prefix=str(input_prefix),
            output_prefix=str(tmp_path / "output"),
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
        )


def test_convert_plink_bfile_rejects_vcf_target(table, tmp_path):
    input_prefix = tmp_path / "input"
    _write_plink_files(
        input_prefix,
        [["1", "SNP1", "0", "300", "A", "G"]],
    )

    with pytest.raises(ValueError, match="allele-label formats"):
        convert_plink_bfile(
            input_prefix=str(input_prefix),
            output_prefix=str(tmp_path / "output"),
            table=table,
            from_fmt="TOP",
            to_fmt="VCF",
        )


def test_convert_plink_bfile_batch_rewrites_filesets(table, tmp_path):
    input_dir = tmp_path / "plink"
    input_dir.mkdir()
    _write_plink_files(
        input_dir / "herd_a",
        [["1", "SNP1", "0", "300", "A", "G"]],
    )
    _write_plink_files(
        input_dir / "herd_b",
        [["1", "SNP2", "0", "700", "A", "C"]],
    )
    (input_dir / "notes.txt").write_text("not a fileset\n")
    out_dir = tmp_path / "converted"

    stats = convert_plink_bfile_batch(
        input_dir=str(input_dir),
        output_dir=str(out_dir),
        pattern="*.bed",
        suffix=".plus",
        overwrite=False,
        table=table,
        from_fmt="TOP",
        to_fmt="PLUS",
    )

    assert stats.filesets_total == 2
    assert stats.filesets_converted == 2
    assert stats.summary_path == str(out_dir / "conversion_summary.csv")
    assert (out_dir / "herd_a.plus.bim").read_text().strip() == "1\tSNP1\t0\t300\tA\tG"
    assert (out_dir / "herd_b.plus.bim").read_text().strip() == "1\tSNP2\t0\t700\tT\tG"
    assert (out_dir / "herd_a.plus.bed").read_bytes() == (input_dir / "herd_a.bed").read_bytes()
    assert (out_dir / "herd_b.plus.fam").read_text() == (input_dir / "herd_b.fam").read_text()
    summary_rows = list(csv.DictReader((out_dir / "conversion_summary.csv").read_text().splitlines()))
    assert [row["variants_total"] for row in summary_rows] == ["1", "1"]
    assert [row["alleles_changed"] for row in summary_rows] == ["0", "2"]


def test_convert_plink_bfile_batch_refuses_overwrite(table, tmp_path):
    input_dir = tmp_path / "plink"
    input_dir.mkdir()
    _write_plink_files(
        input_dir / "herd",
        [["1", "SNP1", "0", "300", "A", "G"]],
    )
    out_dir = tmp_path / "converted"
    out_dir.mkdir()
    (out_dir / "herd.converted.bim").write_text("existing\n")

    with pytest.raises(FileExistsError, match="already exist"):
        convert_plink_bfile_batch(
            input_dir=str(input_dir),
            output_dir=str(out_dir),
            pattern="*.bed",
            suffix=".converted",
            overwrite=False,
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
        )


def test_convert_plink_bfile_batch_reports_allowed_missing_markers(table, tmp_path):
    input_dir = tmp_path / "plink"
    input_dir.mkdir()
    _write_plink_files(
        input_dir / "herd",
        [
            ["1", "SNP2", "0", "700", "A", "C"],
            ["1", "NOT_IN_LOOKUP", "0", "100", "A", "G"],
        ],
    )
    out_dir = tmp_path / "converted"

    stats = convert_plink_bfile_batch(
        input_dir=str(input_dir),
        output_dir=str(out_dir),
        pattern="*.bed",
        suffix=".plus",
        overwrite=False,
        table=table,
        from_fmt="TOP",
        to_fmt="PLUS",
        require_all_markers=False,
    )

    assert stats.filesets_converted == 1
    assert stats.stats[0].variants_total == 2
    assert stats.stats[0].variants_converted == 1
    assert stats.stats[0].variants_missing_lookup == 1
    assert (out_dir / "herd.plus.bim").read_text().splitlines() == [
        "1\tSNP2\t0\t700\tT\tG",
        "1\tNOT_IN_LOOKUP\t0\t100\tA\tG",
    ]
    summary_rows = list(csv.DictReader((out_dir / "conversion_summary.csv").read_text().splitlines()))
    assert summary_rows[0]["variants_total"] == "2"
    assert summary_rows[0]["variants_converted"] == "1"
    assert summary_rows[0]["variants_missing_lookup"] == "1"


def _write_pfile(prefix, pvar_lines):
    prefix.with_suffix(".pgen").write_bytes(b"fake-pgen")
    prefix.with_suffix(".psam").write_text("#IID\tSEX\nI1\t1\n")
    prefix.with_suffix(".pvar").write_text("".join(line + "\n" for line in pvar_lines))


def test_convert_plink_pfile_rewrites_pvar_and_copies_pgen_psam(table, tmp_path):
    input_prefix = tmp_path / "input2"
    output_prefix = tmp_path / "output2"
    _write_pfile(
        input_prefix,
        [
            "##fileformat=PLINKv2.0",
            "#CHROM\tPOS\tID\tREF\tALT",
            "1\t300\tSNP1\tA\tG",
            "1\t700\tSNP2\tA\tC",
        ],
    )

    stats = convert_plink_pfile(
        input_prefix=str(input_prefix),
        output_prefix=str(output_prefix),
        table=table,
        from_fmt="TOP",
        to_fmt="PLUS",
    )

    assert stats.variants_total == 2
    assert stats.variants_converted == 2
    assert stats.alleles_changed == 2
    assert output_prefix.with_suffix(".pgen").read_bytes() == input_prefix.with_suffix(".pgen").read_bytes()
    assert output_prefix.with_suffix(".psam").read_text() == input_prefix.with_suffix(".psam").read_text()
    assert output_prefix.with_suffix(".pvar").read_text().splitlines() == [
        "##fileformat=PLINKv2.0",
        "#CHROM\tPOS\tID\tREF\tALT",
        "1\t300\tSNP1\tA\tG",
        "1\t700\tSNP2\tT\tG",
    ]


def test_convert_plink_pfile_can_update_positions(table, tmp_path):
    input_prefix = tmp_path / "input2"
    output_prefix = tmp_path / "output2"
    _write_pfile(
        input_prefix,
        [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL",
            "0\t0\tSNP2\tA\tC\t.",
        ],
    )

    convert_plink_pfile(
        input_prefix=str(input_prefix),
        output_prefix=str(output_prefix),
        table=table,
        from_fmt="TOP",
        to_fmt="PLUS",
        update_position=True,
    )

    assert output_prefix.with_suffix(".pvar").read_text().splitlines() == [
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL",
        "1\t700\tSNP2\tT\tG\t.",
    ]


def test_convert_plink_pfile_requires_all_three_files(table, tmp_path):
    input_prefix = tmp_path / "input2"
    input_prefix.with_suffix(".pvar").write_text("#CHROM\tPOS\tID\tREF\tALT\n")

    with pytest.raises(ValueError, match=".pgen, .pvar, and .psam"):
        convert_plink_pfile(
            input_prefix=str(input_prefix),
            output_prefix=str(tmp_path / "output2"),
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
        )


def test_convert_plink_pfile_rejects_missing_lookup_marker(table, tmp_path):
    input_prefix = tmp_path / "input2"
    _write_pfile(
        input_prefix,
        [
            "#CHROM\tPOS\tID\tREF\tALT",
            "1\t100\tNOT_IN_LOOKUP\tA\tG",
        ],
    )

    with pytest.raises(ValueError, match="NOT_IN_LOOKUP"):
        convert_plink_pfile(
            input_prefix=str(input_prefix),
            output_prefix=str(tmp_path / "output2"),
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
        )


def test_convert_plink_pfile_rejects_multiallelic_marker(table, tmp_path):
    input_prefix = tmp_path / "input2"
    _write_pfile(
        input_prefix,
        [
            "#CHROM\tPOS\tID\tREF\tALT",
            "1\t300\tSNP1\tA\tG,T",
        ],
    )

    with pytest.raises(ValueError, match="multiallelic"):
        convert_plink_pfile(
            input_prefix=str(input_prefix),
            output_prefix=str(tmp_path / "output2"),
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
        )


def test_convert_plink_pfile_rejects_headerless_pvar(table, tmp_path):
    input_prefix = tmp_path / "input2"
    _write_pfile(input_prefix, ["1\t300\tSNP1\tA\tG"])

    with pytest.raises(ValueError, match="#CHROM header"):
        convert_plink_pfile(
            input_prefix=str(input_prefix),
            output_prefix=str(tmp_path / "output2"),
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
        )


def test_convert_plink_pfile_batch_rewrites_filesets(table, tmp_path):
    input_dir = tmp_path / "pfiles"
    input_dir.mkdir()
    _write_pfile(
        input_dir / "herd_a",
        ["#CHROM\tPOS\tID\tREF\tALT", "1\t300\tSNP1\tA\tG"],
    )
    _write_pfile(
        input_dir / "herd_b",
        ["#CHROM\tPOS\tID\tREF\tALT", "1\t700\tSNP2\tA\tC"],
    )
    out_dir = tmp_path / "converted"

    stats = convert_plink_pfile_batch(
        input_dir=str(input_dir),
        output_dir=str(out_dir),
        pattern="*.pgen",
        suffix=".plus",
        overwrite=False,
        table=table,
        from_fmt="TOP",
        to_fmt="PLUS",
    )

    assert stats.filesets_total == 2
    assert stats.filesets_converted == 2
    assert stats.summary_path == str(out_dir / "conversion_summary.csv")
    assert (out_dir / "herd_a.plus.pvar").read_text().splitlines() == [
        "#CHROM\tPOS\tID\tREF\tALT",
        "1\t300\tSNP1\tA\tG",
    ]
    assert (out_dir / "herd_b.plus.pvar").read_text().splitlines() == [
        "#CHROM\tPOS\tID\tREF\tALT",
        "1\t700\tSNP2\tT\tG",
    ]
    assert (out_dir / "herd_a.plus.pgen").read_bytes() == (input_dir / "herd_a.pgen").read_bytes()
    assert (out_dir / "herd_b.plus.psam").read_text() == (input_dir / "herd_b.psam").read_text()
    summary_rows = list(csv.DictReader((out_dir / "conversion_summary.csv").read_text().splitlines()))
    assert [row["variants_converted"] for row in summary_rows] == ["1", "1"]
    assert [row["alleles_changed"] for row in summary_rows] == ["0", "2"]


def test_convert_plink_pfile_batch_requires_matching_pgen(table, tmp_path):
    input_dir = tmp_path / "pfiles"
    input_dir.mkdir()
    _write_pfile(
        input_dir / "herd",
        ["#CHROM\tPOS\tID\tREF\tALT", "1\t300\tSNP1\tA\tG"],
    )

    with pytest.raises(ValueError, match="pattern should match .pgen"):
        convert_plink_pfile_batch(
            input_dir=str(input_dir),
            output_dir=str(tmp_path / "converted"),
            pattern="*.pvar",
            suffix=".converted",
            overwrite=False,
            table=table,
            from_fmt="TOP",
            to_fmt="PLUS",
        )
