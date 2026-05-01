from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from .aligner import AlignmentResult, GenomicStrand, reverse_complement
from .manifest import IllmnStrand, ManifestRecord, ManifestType


@dataclass
class VariantResult:
    # Identity
    name: str
    alt_name: Optional[str]
    # Raw
    ab_a: str
    ab_b: str
    # Position
    chromosome: Optional[str]
    position: Optional[int]
    strand: Optional[GenomicStrand]
    # Conversions
    top_a: Optional[str]
    top_b: Optional[str]
    forward_a: Optional[str]
    forward_b: Optional[str]
    design_a: Optional[str]
    design_b: Optional[str]
    plus_a: Optional[str]
    plus_b: Optional[str]
    # VCF
    vcf_ref: Optional[str]
    vcf_alt: Optional[str]
    vcf_a: Optional[str]
    vcf_b: Optional[str]
    # Debug
    alignment_text: Optional[str]
    is_indel: bool


def _rc(s: Optional[str]) -> Optional[str]:
    if s is None:
        return None
    return reverse_complement(s)


def _set_vcf_ab(
    plus_a: Optional[str],
    plus_b: Optional[str],
    vcf_ref: Optional[str],
    is_indel: bool,
    indel_ref_is_del: Optional[bool] = None,
) -> tuple[Optional[str], Optional[str]]:
    if vcf_ref is None:
        return None, None
    if is_indel:
        if indel_ref_is_del is None:
            return None, None
        # indel_ref_is_del=True → deletion allele matches reference.
        ref_ab = "D" if indel_ref_is_del else "I"
        vcf_a = "REF" if plus_a == ref_ab else "ALT"
        vcf_b = "REF" if plus_b == ref_ab else "ALT"
        return vcf_a, vcf_b
    else:
        if plus_a == vcf_ref:
            return "REF", "ALT"
        if plus_b == vcf_ref:
            return "ALT", "REF"
        return "ALT", "ALT"


def _compute_illumina_snp(record: ManifestRecord, aln: AlignmentResult) -> dict:
    ab_a = record.ab_a.upper()
    ab_b = record.ab_b.upper()

    # Step 1: TOP/BOT
    if record.ilmn_strand == IllmnStrand.TOP:
        top_a, top_b = ab_a, ab_b
        bot_a, bot_b = _rc(ab_a), _rc(ab_b)
    else:  # BOT
        top_a, top_b = _rc(ab_a), _rc(ab_b)
        bot_a, bot_b = ab_a, ab_b

    # Step 2: SourceSeq alleles
    if record.source_strand == IllmnStrand.TOP:
        ss_a, ss_b = top_a, top_b
    else:
        ss_a, ss_b = bot_a, bot_b

    # Step 3: DESIGN
    if record.source_strand == record.ilmn_strand:
        design_a, design_b = ss_a, ss_b
    else:
        design_a, design_b = _rc(ss_a), _rc(ss_b)

    # Step 4: FORWARD (driven by the T/B _ F/R/U code from IlmnID)
    code = record.ilmn_id_tbpm_fru or ""
    if record.ilmn_strand == IllmnStrand.TOP and code == "T_F":
        forward_a, forward_b = top_a, top_b
    elif record.ilmn_strand == IllmnStrand.TOP and code in ("T_R", "T_U"):
        forward_a, forward_b = bot_a, bot_b
    elif record.ilmn_strand == IllmnStrand.BOT and code in ("B_F", "B_U"):
        forward_a, forward_b = bot_a, bot_b
    elif record.ilmn_strand == IllmnStrand.BOT and code == "B_R":
        forward_a, forward_b = top_a, top_b
    else:
        forward_a, forward_b = None, None

    # Step 5: PLUS (genomic forward strand)
    if aln.strand == GenomicStrand.PLUS:
        plus_a, plus_b = ss_a, ss_b
    elif aln.strand == GenomicStrand.MINUS:
        plus_a, plus_b = _rc(ss_a), _rc(ss_b)
    else:
        plus_a, plus_b = None, None

    return dict(
        top_a=top_a, top_b=top_b,
        forward_a=forward_a, forward_b=forward_b,
        design_a=design_a, design_b=design_b,
        plus_a=plus_a, plus_b=plus_b,
    )


def _compute_illumina_indel(record: ManifestRecord) -> dict:
    # Indels are strand-independent: all formats carry the AB designation
    a, b = record.ab_a, record.ab_b
    return dict(
        top_a=a, top_b=b,
        forward_a=a, forward_b=b,
        design_a=a, design_b=b,
        plus_a=a, plus_b=b,
    )


def _compute_affymetrix_snp(record: ManifestRecord, aln: AlignmentResult) -> dict:
    # Affy alleles are on the forward strand by definition
    forward_a = record.ab_a.upper()
    forward_b = record.ab_b.upper()

    if aln.strand == GenomicStrand.PLUS:
        plus_a, plus_b = forward_a, forward_b
    elif aln.strand == GenomicStrand.MINUS:
        plus_a, plus_b = _rc(forward_a), _rc(forward_b)
    else:
        plus_a, plus_b = None, None

    return dict(
        top_a=None, top_b=None,
        forward_a=forward_a, forward_b=forward_b,
        design_a=None, design_b=None,
        plus_a=plus_a, plus_b=plus_b,
    )


def _compute_affymetrix_indel(record: ManifestRecord) -> dict:
    a, b = record.ab_a, record.ab_b
    return dict(
        top_a=a, top_b=b,
        forward_a=a, forward_b=b,
        design_a=a, design_b=b,
        plus_a=a, plus_b=b,
    )


def compute_conversion(record: ManifestRecord, aln: AlignmentResult) -> VariantResult:
    if record.manifest_type == ManifestType.ILLUMINA:
        fields = _compute_illumina_indel(record) if record.is_indel else _compute_illumina_snp(record, aln)
    else:  # AFFYMETRIX
        fields = _compute_affymetrix_indel(record) if record.is_indel else _compute_affymetrix_snp(record, aln)

    vcf_a, vcf_b = _set_vcf_ab(
        fields["plus_a"],
        fields["plus_b"],
        aln.vcf_ref,
        record.is_indel,
        aln.indel_ref_is_del,
    )

    return VariantResult(
        name=record.name,
        alt_name=record.alt_name,
        ab_a=record.ab_a,
        ab_b=record.ab_b,
        chromosome=aln.chromosome,
        position=aln.position,
        strand=aln.strand,
        vcf_ref=aln.vcf_ref,
        vcf_alt=aln.vcf_alt,
        vcf_a=vcf_a,
        vcf_b=vcf_b,
        alignment_text=aln.alignment_text,
        is_indel=record.is_indel,
        **fields,
    )
