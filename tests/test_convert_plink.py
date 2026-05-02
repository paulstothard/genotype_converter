from __future__ import annotations

import pytest

from genotype_converter.convert_genotypes import load_lookup_table
from genotype_converter.convert_plink import convert_plink_bfile, convert_plink_pfile


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
