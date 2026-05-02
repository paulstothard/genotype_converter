from __future__ import annotations

import pytest
from genotype_converter.aligner import (
    DeterminationType,
    GenomicStrand,
    VariantAligner,
    _resolve_indel_vcf,
    reverse_complement,
)
from genotype_converter.manifest import ManifestRecord, ManifestType, parse_manifest


def test_reverse_complement_basic():
    assert reverse_complement("A") == "T"
    assert reverse_complement("T") == "A"
    assert reverse_complement("G") == "C"
    assert reverse_complement("C") == "G"
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AACCGGTT") == "AACCGGTT"
    assert reverse_complement("AACC") == "GGTT"


def test_snp1_aligns_plus_strand(manifest_path, reference_path):
    records = parse_manifest(manifest_path)
    snp1 = next(r for r in records if r.name == "SNP1")
    aligner = VariantAligner(reference_path)
    result = aligner.align(snp1)
    assert result.chromosome is not None
    assert result.position == 300
    assert result.strand == GenomicStrand.PLUS
    assert result.vcf_ref == "A"
    assert result.vcf_alt == "G"
    assert result.alignment_text is None


def test_alignment_text_is_optional(manifest_path, reference_path):
    records = parse_manifest(manifest_path)
    snp1 = next(r for r in records if r.name == "SNP1")
    aligner = VariantAligner(reference_path, save_alignment=True)
    result = aligner.align(snp1)
    assert result.alignment_text is not None
    assert "Determination type:" in result.alignment_text


def test_snp2_position(manifest_path, reference_path):
    records = parse_manifest(manifest_path)
    snp2 = next(r for r in records if r.name == "SNP2")
    aligner = VariantAligner(reference_path)
    result = aligner.align(snp2)
    assert result.position == 700
    assert result.vcf_ref == "T"
    assert result.vcf_alt == "G"


def test_snp3_position(manifest_path, reference_path):
    records = parse_manifest(manifest_path)
    snp3 = next(r for r in records if r.name == "SNP3")
    aligner = VariantAligner(reference_path)
    result = aligner.align(snp3)
    assert result.position == 900
    assert result.vcf_ref == "G"
    assert result.vcf_alt == "T"


def test_snp4_minus_strand(manifest_path, reference_path):
    records = parse_manifest(manifest_path)
    snp4 = next(r for r in records if r.name == "SNP4")
    aligner = VariantAligner(reference_path)
    result = aligner.align(snp4)
    assert result.position == 200
    assert result.vcf_ref == "A"
    assert result.vcf_alt == "G"


def test_snp5_position(manifest_path, reference_path):
    records = parse_manifest(manifest_path)
    snp5 = next(r for r in records if r.name == "SNP5")
    aligner = VariantAligner(reference_path)
    result = aligner.align(snp5)
    assert result.position == 800
    assert result.vcf_ref == "T"
    assert result.vcf_alt == "G"


def _indel_record(first_allele, second_allele, ab_a="I", ab_b="D"):
    """Build a minimal ManifestRecord for indel unit tests."""
    return ManifestRecord(
        name="INDEL1",
        alt_name="INDEL1",
        ab_a=ab_a,
        ab_b=ab_b,
        flanking=f"AAAA[{first_allele}/{second_allele}]TTTT",
        first_allele=first_allele,
        second_allele=second_allele,
        manifest_type=ManifestType.ILLUMINA,
        is_indel=True,
    )


def test_indel_insertion_plus_strand():
    ref_seq = "ACGTACGT"
    record = _indel_record("-", "ATGC")
    vcf_ref, vcf_alt, pos_0, det, ref_is_del = _resolve_indel_vcf(
        record, GenomicStrand.PLUS, ref_seq, 4
    )
    assert vcf_ref == ref_seq[3]
    assert vcf_alt == ref_seq[3] + "ATGC"
    assert pos_0 == 3
    assert det == DeterminationType.INDEL_INSERTION
    assert ref_is_del is True


def test_indel_deletion_plus_strand():
    ref_seq = "ACGTATGCTT"
    record = _indel_record("ATGC", "-")
    vcf_ref, vcf_alt, pos_0, det, ref_is_del = _resolve_indel_vcf(
        record, GenomicStrand.PLUS, ref_seq, 7
    )
    assert vcf_ref == ref_seq[3] + "ATGC"
    assert vcf_alt == ref_seq[3]
    assert pos_0 == 3
    assert det == DeterminationType.INDEL_DELETION
    assert ref_is_del is False


def test_indel_insertion_minus_strand():
    ref_seq = "ACGTACGT"
    record = _indel_record("-", "ATGC")
    vcf_ref, vcf_alt, pos_0, det, ref_is_del = _resolve_indel_vcf(
        record, GenomicStrand.MINUS, ref_seq, 4
    )
    assert vcf_ref == ref_seq[3]
    assert vcf_alt == ref_seq[3] + reverse_complement("ATGC")
    assert ref_is_del is True


def test_indel_symbolic_fallback():
    ref_seq = "ACGTACGT"
    record = _indel_record("I", "D", ab_a="I", ab_b="D")
    vcf_ref, vcf_alt, pos_0, det, ref_is_del = _resolve_indel_vcf(
        record, GenomicStrand.PLUS, ref_seq, 4
    )
    assert vcf_ref == vcf_alt
    assert det == DeterminationType.INDEL_SITE
    assert ref_is_del is None


def test_indel_no_anchor():
    ref_seq = "ACGT"
    record = _indel_record("-", "ATGC")
    vcf_ref, vcf_alt, pos_0, det, ref_is_del = _resolve_indel_vcf(
        record, GenomicStrand.PLUS, ref_seq, 0
    )
    assert det == DeterminationType.INDEL_NO_ANCHOR
    assert ref_is_del is None
