from __future__ import annotations

import pytest

from genotype_converter.manifest import IllmnStrand, parse_manifest


def test_parse_returns_eight_records(manifest_path):
    records = parse_manifest(manifest_path)
    assert len(records) == 8  # 5 SNPs + 3 indels


def test_indel_records_parsed(manifest_path):
    records = parse_manifest(manifest_path)
    indels = [r for r in records if r.is_indel]
    assert len(indels) == 3
    by_name = {r.name: r for r in indels}

    # INDEL1: SNP=[I/D], flanking [-/ACGT]
    assert by_name["INDEL1"].ab_a == "I"
    assert by_name["INDEL1"].ab_b == "D"
    assert by_name["INDEL1"].first_allele == "-"
    assert by_name["INDEL1"].second_allele == "ACGT"

    # INDEL2: SNP=[I/D], flanking [TCGA/-]
    assert by_name["INDEL2"].ab_a == "I"
    assert by_name["INDEL2"].first_allele == "TCGA"
    assert by_name["INDEL2"].second_allele == "-"

    # INDEL3: SNP=[D/I] — ab_a=D, ab_b=I
    assert by_name["INDEL3"].ab_a == "D"
    assert by_name["INDEL3"].ab_b == "I"
    assert by_name["INDEL3"].first_allele == "-"
    assert by_name["INDEL3"].second_allele == "TTCC"


def test_snp1_identity(manifest_path):
    records = parse_manifest(manifest_path)
    snp1 = next(r for r in records if r.name == "SNP1")
    assert snp1.alt_name == "SNP1-0_T_F_1511658221"
    assert snp1.ab_a == "A"
    assert snp1.ab_b == "G"
    assert snp1.ilmn_strand == IllmnStrand.TOP
    assert snp1.source_strand == IllmnStrand.TOP
    assert snp1.ilmn_id_tbpm_fru == "T_F"
    assert snp1.is_snp
    assert not snp1.is_indel


def test_snp2_bot_strand(manifest_path):
    records = parse_manifest(manifest_path)
    snp2 = next(r for r in records if r.name == "SNP2")
    assert snp2.ilmn_strand == IllmnStrand.BOT
    assert snp2.ilmn_id_tbpm_fru == "B_F"
    assert snp2.ab_a == "T"
    assert snp2.ab_b == "G"


def test_snp4_b_f_probe(manifest_path):
    records = parse_manifest(manifest_path)
    snp4 = next(r for r in records if r.name == "SNP4")
    assert snp4.ilmn_id_tbpm_fru == "B_F"


def test_snp5_t_r_probe(manifest_path):
    records = parse_manifest(manifest_path)
    snp5 = next(r for r in records if r.name == "SNP5")
    assert snp5.ilmn_id_tbpm_fru == "T_R"
    assert snp5.source_strand == IllmnStrand.BOT
    assert snp5.ilmn_strand == IllmnStrand.TOP


def test_unrecognised_manifest_raises(tmp_path):
    manifest = tmp_path / "not_a_manifest.csv"
    manifest.write_text("marker,allele_a,allele_b\nSNP1,A,G\n")

    with pytest.raises(ValueError, match="Manifest format not recognised"):
        parse_manifest(str(manifest))

