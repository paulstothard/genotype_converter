from __future__ import annotations

import pytest
from genotype_converter.aligner import AlignmentResult, GenomicStrand
from genotype_converter.conversion import compute_conversion
from genotype_converter.manifest import IllmnStrand, ManifestRecord, ManifestType, parse_manifest


def _get(records, name):
    return next(r for r in records if r.name == name)


def test_snp1_all_formats(manifest_path):
    records = parse_manifest(manifest_path)
    snp1 = _get(records, "SNP1")
    aln = AlignmentResult(
        chromosome="1", position=300, strand=GenomicStrand.PLUS,
        vcf_ref="A", vcf_alt="G"
    )
    r = compute_conversion(snp1, aln)
    assert r.top_a == "A"
    assert r.top_b == "G"
    assert r.forward_a == "A"
    assert r.design_a == "A"
    assert r.plus_a == "A"
    assert r.vcf_a == "REF"
    assert r.vcf_b == "ALT"


def test_snp2_bot_bf(manifest_path):
    records = parse_manifest(manifest_path)
    snp2 = _get(records, "SNP2")
    aln = AlignmentResult(
        chromosome="1", position=700, strand=GenomicStrand.PLUS,
        vcf_ref="T", vcf_alt="G"
    )
    r = compute_conversion(snp2, aln)
    assert r.top_a == "A"
    assert r.top_b == "C"
    assert r.forward_a == "T"
    assert r.design_a == "T"
    assert r.plus_a == "T"
    assert r.vcf_a == "REF"
    assert r.vcf_b == "ALT"


def test_snp3_b_f_ref_is_b(manifest_path):
    records = parse_manifest(manifest_path)
    snp3 = _get(records, "SNP3")
    aln = AlignmentResult(
        chromosome="1", position=900, strand=GenomicStrand.PLUS,
        vcf_ref="G", vcf_alt="T"
    )
    r = compute_conversion(snp3, aln)
    assert r.plus_b == "G"
    assert r.vcf_b == "REF"
    assert r.vcf_a == "ALT"


def test_snp4_minus_strand(manifest_path):
    records = parse_manifest(manifest_path)
    snp4 = _get(records, "SNP4")
    aln = AlignmentResult(
        chromosome="1", position=200, strand=GenomicStrand.MINUS,
        vcf_ref="A", vcf_alt="G"
    )
    r = compute_conversion(snp4, aln)
    assert r.top_a == "A"
    assert r.top_b == "G"
    assert r.forward_a == "T"
    assert r.design_a == "T"
    assert r.plus_a == "A"
    assert r.vcf_a == "REF"


def test_snp5_top_tr(manifest_path):
    records = parse_manifest(manifest_path)
    snp5 = _get(records, "SNP5")
    aln = AlignmentResult(
        chromosome="1", position=800, strand=GenomicStrand.PLUS,
        vcf_ref="T", vcf_alt="G"
    )
    r = compute_conversion(snp5, aln)
    assert r.top_a == "A"
    assert r.forward_a == "T"
    assert r.design_a == "A"
    assert r.plus_a == "T"
    assert r.vcf_a == "REF"


def _illumina_indel_record(ab_a="I", ab_b="D", first_allele="-", second_allele="ATGC"):
    return ManifestRecord(
        name="INDEL1", alt_name="INDEL1",
        ab_a=ab_a, ab_b=ab_b,
        flanking=f"AAAA[{first_allele}/{second_allele}]TTTT",
        first_allele=first_allele, second_allele=second_allele,
        manifest_type=ManifestType.ILLUMINA, is_indel=True,
        ilmn_strand=IllmnStrand.TOP, source_strand=IllmnStrand.TOP,
        ilmn_id_tbpm_fru="T_F",
    )


def test_indel_vcf_ab_insertion():
    record = _illumina_indel_record(ab_a="I", ab_b="D", first_allele="-", second_allele="ATGC")
    aln = AlignmentResult(
        chromosome="1", position=4, strand=GenomicStrand.PLUS,
        vcf_ref="T", vcf_alt="TATGC",
        indel_ref_is_del=True,
    )
    r = compute_conversion(record, aln)
    assert r.vcf_a == "ALT"
    assert r.vcf_b == "REF"


def test_indel_vcf_ab_deletion():
    record = _illumina_indel_record(ab_a="I", ab_b="D", first_allele="ATGC", second_allele="-")
    aln = AlignmentResult(
        chromosome="1", position=4, strand=GenomicStrand.PLUS,
        vcf_ref="TATGC", vcf_alt="T",
        indel_ref_is_del=False,
    )
    r = compute_conversion(record, aln)
    assert r.vcf_a == "REF"
    assert r.vcf_b == "ALT"
