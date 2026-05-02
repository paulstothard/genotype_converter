from __future__ import annotations

import pytest
from genotype_converter.aligner import (
    DeterminationType,
    GenomicStrand,
    ProbeSide,
    VariantAligner,
    _aligned_probe_placements,
    _choose_probe_ref_pos,
    _cigar_query_to_ref,
    _format_probe_alignment_line,
    _hit_has_gap_near_variant,
    _prepare_query,
    _probe_adjacent_ref_pos,
    _probe_placements,
    _reference_probe_adjacent_ref_pos,
    _refine_snp_with_local_alignment,
    _probe_side,
    _resolve_indel_vcf,
    _semiglobal_align_query_to_ref,
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


def test_prepare_query_preserves_real_flanking_ns():
    query, n_pos = _prepare_query("ACNT[A/G]TGNA")
    assert query == "ACNTNTGNA"
    assert n_pos == 4


def test_cigar_query_to_ref_reports_left_adjacent_gap_position():
    # 60M 1D 60M: if the variant is the first matched base after a deletion
    # in the query, report the deleted reference base adjacent to the probe.
    assert _cigar_query_to_ref(
        60, 1000, [(60, 0), (1, 2), (60, 0)], probe_side=ProbeSide.LEFT
    ) == 1060
    assert _cigar_query_to_ref(
        60, 1000, [(59, 0), (1, 2), (62, 0)], probe_side=ProbeSide.LEFT
    ) == 1060
    assert _cigar_query_to_ref(
        60, 1000, [(60, 0), (1, 2), (60, 0)], probe_side=ProbeSide.RIGHT
    ) == 1061

    # 60M 1I 60M: if the variant is inside an insertion in the query, report
    # the left-adjacent reference base, not the right-adjacent base.
    assert _cigar_query_to_ref(
        60, 1000, [(60, 0), (1, 1), (60, 0)], probe_side=ProbeSide.RIGHT
    ) == 1059

    # Indel anchoring uses the direct CIGAR mapping, preserving existing
    # insertion/deletion anchor behavior.
    assert _cigar_query_to_ref(60, 1000, [(60, 0), (1, 1), (60, 0)]) == 1060


def test_gap_near_variant_detection():
    class Hit:
        strand = 1
        q_st = 0
        q_en = 121
        cigar = [(57, 0), (1, 2), (64, 0)]

    assert _hit_has_gap_near_variant(60, Hit())
    assert not _hit_has_gap_near_variant(100, Hit())


def test_probe_side_detects_manifest_probe_orientation(manifest_path):
    records = parse_manifest(manifest_path)
    snp1 = next(r for r in records if r.name == "SNP1")
    snp5 = next(r for r in records if r.name == "SNP5")
    assert _probe_side(snp1) == ProbeSide.LEFT
    assert _probe_side(snp5) == ProbeSide.RIGHT


def test_probe_placement_detects_allele_including_left_probe():
    record = ManifestRecord(
        name="SNP",
        alt_name=None,
        ab_a="A",
        ab_b="G",
        flanking="ACCC[A/G]TTT",
        first_allele="A",
        second_allele="G",
        manifest_type=ManifestType.ILLUMINA,
        is_indel=False,
        allele_a_probe_seq="ACCCA",
    )

    placements = _probe_placements(record)
    assert [(p.side, p.includes_allele) for p in placements] == [
        (ProbeSide.LEFT, True)
    ]
    assert _probe_side(record) == ProbeSide.LEFT


def test_probe_placement_reports_disagreeing_probe_sides_as_ambiguous():
    record = ManifestRecord(
        name="SNP",
        alt_name=None,
        ab_a="A",
        ab_b="G",
        flanking="ACCC[A/G]TTT",
        first_allele="A",
        second_allele="G",
        manifest_type=ManifestType.ILLUMINA,
        is_indel=False,
        allele_a_probe_seq="ACCC",
        allele_b_probe_seq="TTT",
    )

    assert {p.side for p in _probe_placements(record)} == {
        ProbeSide.LEFT,
        ProbeSide.RIGHT,
    }
    assert _probe_side(record) is None


def test_aligned_probe_placements_flip_on_minus_strand():
    class Hit:
        strand = -1

    record = ManifestRecord(
        name="SNP",
        alt_name=None,
        ab_a="A",
        ab_b="G",
        flanking="ACCC[A/G]TTT",
        first_allele="A",
        second_allele="G",
        manifest_type=ManifestType.ILLUMINA,
        is_indel=False,
        allele_a_probe_seq="ACCCA",
    )

    placements = _aligned_probe_placements(record, Hit())
    assert [(p.side, p.includes_allele) for p in placements] == [
        (ProbeSide.RIGHT, True)
    ]
    assert placements[0].sequence == "TGGGT"


def test_probe_adjacent_ref_pos_uses_reference_base_next_to_probe():
    record = ManifestRecord(
        name="SNP",
        alt_name=None,
        ab_a="A",
        ab_b="G",
        flanking="ACCC[A/G]TTT",
        first_allele="A",
        second_allele="G",
        manifest_type=ManifestType.ILLUMINA,
        is_indel=False,
        allele_a_probe_seq="ACCC",
    )

    assert _probe_adjacent_ref_pos(
        record,
        q_aln="ACCCNTTT",
        r_aln="ACCC-ATT",
        r_start=100,
        probe_side=ProbeSide.LEFT,
        n_index=4,
    ) == 104

    right_record = ManifestRecord(
        name="SNP",
        alt_name=None,
        ab_a="A",
        ab_b="G",
        flanking="AAAC[A/G]TTT",
        first_allele="A",
        second_allele="G",
        manifest_type=ManifestType.ILLUMINA,
        is_indel=False,
        allele_a_probe_seq="TTT",
    )
    assert _probe_adjacent_ref_pos(
        right_record,
        q_aln="AAACNTTT",
        r_aln="AAAC-ATT",
        r_start=100,
        probe_side=ProbeSide.RIGHT,
        n_index=4,
    ) == 103


def test_probe_adjacent_ref_pos_handles_allele_including_probe():
    record = ManifestRecord(
        name="SNP",
        alt_name=None,
        ab_a="A",
        ab_b="G",
        flanking="ACCC[A/G]TTT",
        first_allele="A",
        second_allele="G",
        manifest_type=ManifestType.ILLUMINA,
        is_indel=False,
        allele_a_probe_seq="ACCCA",
    )

    assert _probe_adjacent_ref_pos(
        record,
        q_aln="ACCCNTTT",
        r_aln="ACCCATTT",
        r_start=100,
        probe_side=ProbeSide.LEFT,
        n_index=4,
    ) == 104


def test_probe_adjacent_ref_pos_does_not_confuse_real_n_with_variant_n():
    record = ManifestRecord(
        name="SNP",
        alt_name=None,
        ab_a="A",
        ab_b="G",
        flanking="NACC[A/G]TTT",
        first_allele="A",
        second_allele="G",
        manifest_type=ManifestType.ILLUMINA,
        is_indel=False,
        allele_a_probe_seq="ACC",
    )

    assert _probe_adjacent_ref_pos(
        record,
        q_aln="NACCNTTT",
        r_aln="NACC-ATT",
        r_start=100,
        probe_side=ProbeSide.LEFT,
        n_index=4,
    ) == 104


def test_reference_probe_adjacent_ref_pos_uses_subject_side_probe_endpoint():
    record = ManifestRecord(
        name="SNP",
        alt_name=None,
        ab_a="A",
        ab_b="G",
        flanking="GCTAAT[A/G]AAA",
        first_allele="A",
        second_allele="G",
        manifest_type=ManifestType.ILLUMINA,
        is_indel=False,
        allele_a_probe_seq="GCTAAT",
    )

    assert _reference_probe_adjacent_ref_pos(
        record,
        q_aln="TTGCTAA-TNGATT",
        r_aln="TTGCTAATTGGATT",
        r_start=100,
        probe_side=ProbeSide.LEFT,
        n_index=8,
        candidate_ref_pos=108,
    ) == 108


def test_reference_probe_adjacent_ref_pos_handles_allele_including_probe():
    record = ManifestRecord(
        name="SNP",
        alt_name=None,
        ab_a="A",
        ab_b="G",
        flanking="ACCC[A/G]TTT",
        first_allele="A",
        second_allele="G",
        manifest_type=ManifestType.ILLUMINA,
        is_indel=False,
        allele_a_probe_seq="ACCCA",
    )

    assert _reference_probe_adjacent_ref_pos(
        record,
        q_aln="ACCCNTTT",
        r_aln="ACCCATTT",
        r_start=100,
        probe_side=ProbeSide.LEFT,
        n_index=4,
        candidate_ref_pos=104,
    ) == 104


def test_reference_probe_position_uses_reverse_complemented_minus_strand_probe():
    class Hit:
        strand = -1

    record = ManifestRecord(
        name="ARS-BFGL-NGS-101112",
        alt_name=None,
        ab_a="A",
        ab_b="G",
        flanking=(
            "ACCAGACCACGCAGGTCAGAGGCCTAGCGGCAGGCTGGCCAGGGAGTTAACTCGAGGTTC"
            "[T/C]"
            "TTTCACTGAATTGATAACTGCTGCTTTCTACTTGACTGAGAATCAAACACACGCACGGCT"
        ),
        first_allele="T",
        second_allele="C",
        manifest_type=ManifestType.ILLUMINA,
        is_indel=False,
        allele_a_probe_seq="GTGTTTGATTCTCAGTCAAGTAGAAAGCAGCAGTTATCAATTCAGTGAAA",
    )

    placements = _aligned_probe_placements(record, Hit())
    assert placements[0].side == ProbeSide.LEFT
    assert placements[0].sequence == record.allele_a_probe_seq
    assert _reference_probe_adjacent_ref_pos(
        record,
        q_aln=(
            "AGCCGTGCGTGTGTTTGATTCTCAGTCAAGTAGAAAGCAGCAGTTATCAATTCAGTG-"
            "AAANGAACCTCGAGTTAACTCCCTGGCCAGCCTGCCGCTAGGCCTCTGACCTGCGTGGTCTGGT"
        ),
        r_aln=(
            "AGCCGTGCGTGTGTTTGATTCTCAGTCAAGTAGAAAGCAGCAGTTATCAATTCAGTG"
            "AAAAGGAACCTCGAGTTAACTCCCTGGCCAGCCTGCCGCTAGGCCTCTGACCTGCGTGGTCTGGT"
        ),
        r_start=115336747,
        probe_side=ProbeSide.LEFT,
        n_index=60,
        candidate_ref_pos=115336808,
        placements=placements,
    ) == 115336807


def test_reference_probe_position_overrides_direct_cigar_when_probe_is_unique():
    class Hit:
        r_st = 115336797
        r_en = 115336819

    record = ManifestRecord(
        name="ARS-BFGL-NGS-101112",
        alt_name=None,
        ab_a="A",
        ab_b="G",
        flanking=(
            "ACCAGACCACGCAGGTCAGAGGCCTAGCGGCAGGCTGGCCAGGGAGTTAACTCGAGGTTC"
            "[T/C]"
            "TTTCACTGAATTGATAACTGCTGCTTTCTACTTGACTGAGAATCAAACACACGCACGGCT"
        ),
        first_allele="T",
        second_allele="C",
        manifest_type=ManifestType.ILLUMINA,
        is_indel=False,
        allele_a_probe_seq="GTGTTTGATTCTCAGTCAAGTAGAAAGCAGCAGTTATCAATTCAGTGAAA",
    )
    query = (
        "ACCAGACCACGCAGGTCAGAGGCCTAGCGGCAGGCTGGCCAGGGAGTTAACTCGAGGTTC"
        "N"
        "TTTCACTGAATTGATAACTGCTGCTTTCTACTTGACTGAGAATCAAACACACGCACGGCT"
    )
    ref_window = (
        "AGCCGTGCGTGTGTTTGATTCTCAGTCAAGTAGAAAGCAGCAGTTATCAATTCAGTG"
        "AAAAGGAACCTCGAGTTAACTCCCTGGCCAGCCTGCCGCTAGGCCTCTGACCTGCGTGGTCTGGT"
    )
    placements = [
        p for p in _aligned_probe_placements(record, type("Hit", (), {"strand": -1})())
    ]

    pos, det, _q_aln, _r_aln, _start = _refine_snp_with_local_alignment(
        record=record,
        query=query,
        n_pos=60,
        chrom="3",
        hit=Hit(),
        strand=GenomicStrand.MINUS,
        probe_side=ProbeSide.LEFT,
        direct_ref_pos=115336808,
        allele1="G",
        allele2="A",
        fetch_ref=lambda _chrom, _start, _end: ref_window,
        placements=placements,
    )

    assert pos == 115336807
    assert det == DeterminationType.SNP_PROBE_ADJACENT


def test_probe_alignment_line_shows_oriented_probe_on_reference_alignment():
    record = ManifestRecord(
        name="SNP",
        alt_name=None,
        ab_a="A",
        ab_b="G",
        flanking="ACCC[A/G]TTT",
        first_allele="A",
        second_allele="G",
        manifest_type=ManifestType.ILLUMINA,
        is_indel=False,
        allele_a_probe_seq="ACCC",
    )
    line = _format_probe_alignment_line(
        q_aln="AACCCN-TTT",
        r_aln="AACCCGGTTT",
        placements=_probe_placements(record),
    )

    assert line == " ACCC"


def test_semiglobal_align_query_to_ref_allows_reference_window_overhangs():
    q_aln, r_aln, offset = _semiglobal_align_query_to_ref(
        "AACCGG",
        "TTTTAACCGGAAAA",
    )
    assert q_aln == "AACCGG"
    assert r_aln == "AACCGG"
    assert offset == 4


def test_allele_assisted_probe_resolution_uses_unique_compatible_base():
    bases = {
        100: "C",
        101: "A",
    }
    pos, det = _choose_probe_ref_pos(
        100,
        101,
        "A",
        "G",
        lambda p: bases[p],
    )
    assert pos == 101
    assert det == DeterminationType.SNP_ALLELE_ASSISTED


def test_allele_assisted_probe_resolution_keeps_probe_when_ambiguous():
    bases = {
        100: "A",
        101: "G",
    }
    pos, det = _choose_probe_ref_pos(
        100,
        101,
        "A",
        "G",
        lambda p: bases[p],
    )
    assert pos == 100
    assert det == DeterminationType.SNP_AMBIGUOUS


def test_allele_assisted_probe_resolution_keeps_probe_when_uninformative():
    bases = {
        100: "C",
        101: "T",
    }
    pos, det = _choose_probe_ref_pos(
        100,
        101,
        "A",
        "G",
        lambda p: bases[p],
    )
    assert pos == 100
    assert det == DeterminationType.SNP_AMBIGUOUS


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
