from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import Optional

import mappy

from .manifest import ManifestRecord


class GenomicStrand(Enum):
    PLUS = auto()
    MINUS = auto()


class DeterminationType(Enum):
    SNP_ALIGNED = auto()       # ref base matched one flanking allele → REF/ALT assigned
    SNP_NO_ALLELES = auto()    # flanking had no [X/Y] bracket; only ref base available
    INDEL_INSERTION = auto()   # insertion [-/SEQ] resolved with anchor base
    INDEL_DELETION = auto()    # deletion [SEQ/-] resolved with anchor base
    INDEL_SITE = auto()        # symbolic I/D or unrecognised; VCF alleles unavailable
    INDEL_NO_ANCHOR = auto()   # variant at position 0; no anchor base possible


@dataclass
class AlignmentResult:
    chromosome: Optional[str] = None
    position: Optional[int] = None
    strand: Optional[GenomicStrand] = None
    vcf_ref: Optional[str] = None
    vcf_alt: Optional[str] = None
    alignment_text: Optional[str] = None
    determination_type: Optional[DeterminationType] = None
    indel_ref_is_del: Optional[bool] = None


def reverse_complement(seq: str) -> str:
    table = str.maketrans(
        "gatcryswkmbdhvnGATCRYSWKMBDHVN",
        "ctagyrswmkvhdbnCTAGYRSWMKVHDBN",
    )
    return seq.translate(table)[::-1]


def _load_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks).upper()
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def _prepare_query(flanking: str) -> tuple[str, int]:
    """
    Build alignment query from flanking sequence.
    Steps (matching the original pipeline):
      1. Remove existing N/n
      2. Replace [X/Y] variant notation with N
      3. Remove non-GATCN characters
    Returns (query_uppercase, n_position_0based) or (query, -1) if no N found.
    """
    seq = re.sub(r"[Nn]", "", flanking)
    seq = re.sub(r"\[.*?\]", "N", seq)
    seq = re.sub(r"[^GATCNgatcn]", "", seq)
    seq = seq.upper()
    n_pos = seq.find("N")
    return seq, n_pos


def _cigar_query_to_ref(q_offset: int, r_start: int, cigar) -> Optional[int]:
    """Walk a mappy cigar list to map a 0-based query offset to a 0-based reference position."""
    q = 0
    r = r_start
    for length, op in cigar:
        if op in (4, 5):  # S, H
            continue
        if op in (0, 7, 8):  # M, =, X
            if q + length > q_offset:
                return r + (q_offset - q)
            q += length
            r += length
        elif op == 1:  # I – insertion in query
            if q + length > q_offset:
                return r  # N is inside an insertion; report left-adjacent ref pos
            q += length
        elif op == 2:  # D – deletion in query
            r += length
    return None


def _query_pos_to_ref_pos(n_pos: int, hit) -> Optional[int]:
    """Map a 0-based query position to a 0-based reference position via CIGAR."""
    if hit.strand == 1:
        if not (hit.q_st <= n_pos < hit.q_en):
            return None
        return _cigar_query_to_ref(n_pos - hit.q_st, hit.r_st, hit.cigar)
    else:
        rc_offset = hit.q_en - 1 - n_pos
        if not (0 <= rc_offset < hit.q_en - hit.q_st):
            return None
        return _cigar_query_to_ref(rc_offset, hit.r_st, hit.cigar)


def _build_alignment_strings(query: str, ref_seq: str, hit) -> tuple[str, str]:
    """Reconstruct gapped alignment strings from mappy hit + reference sequence."""
    if hit.strand == 1:
        q_bases = query[hit.q_st : hit.q_en]
    else:
        q_bases = reverse_complement(query[hit.q_st : hit.q_en])

    ref_sub = ref_seq[hit.r_st : hit.r_en]
    q_aln: list[str] = []
    r_aln: list[str] = []
    q_pos = 0
    r_pos = 0

    for length, op in hit.cigar:
        if op in (4, 5):
            continue
        if op in (0, 7, 8):
            q_aln.append(q_bases[q_pos : q_pos + length])
            r_aln.append(ref_sub[r_pos : r_pos + length])
            q_pos += length
            r_pos += length
        elif op == 1:
            q_aln.append(q_bases[q_pos : q_pos + length])
            r_aln.append("-" * length)
            q_pos += length
        elif op == 2:
            q_aln.append("-" * length)
            r_aln.append(ref_sub[r_pos : r_pos + length])
            r_pos += length

    return "".join(q_aln), "".join(r_aln)


def _build_ruler(ref_aln: str, r_start: int) -> str:
    """Build a position ruler string aligned to ref_aln (gaps don't advance position)."""
    ruler = []
    pos = r_start
    for ch in ref_aln:
        if ch != "-":
            if pos % 10 == 0:
                ruler.append("|")
            elif pos % 5 == 0:
                ruler.append(".")
            else:
                ruler.append(" ")
            pos += 1
        else:
            ruler.append(" ")
    return "".join(ruler)


def _format_alignment(
    record: ManifestRecord,
    query: str,
    q_aln: str,
    r_aln: str,
    r_start: int,
    ref_pos_0: int,
    ref_base: str,
    allele1: str,
    allele2: str,
    vcf_ref: str,
    vcf_alt: str,
    variant_type: str,
    determination_type: DeterminationType,
) -> str:
    pad = 12
    ref_col = 0
    ref_count = r_start
    for i, ch in enumerate(r_aln):
        if ch != "-":
            if ref_count == ref_pos_0:
                ref_col = i + 1
                break
            ref_count += 1

    lines = []
    lines.append(f"{record.name}\n")
    lines.append(f"Type: {variant_type}\n")
    lines.append(f"{'QUERY':>{pad}} {q_aln}\n")
    lines.append(f"{'SUBJECT':>{pad}} {r_aln}\n")
    ruler = _build_ruler(r_aln, r_start)
    lines.append(f"{r_start + 1:>{pad}} {ruler}\n")
    lines.append(f"{'ALLELE1':>{pad}} {allele1:>{ref_col}}\n")
    lines.append(f"{'ALLELE2':>{pad}} {allele2:>{ref_col}}\n")
    lines.append(f"{'POSITION':>{pad}} {str(ref_pos_0 + 1) + '|':>{ref_col}}\n")
    lines.append(f"{'REF':>{pad}} {ref_base:>{ref_col}}\n")
    lines.append(f"{'VCF_REF':>{pad}} {vcf_ref:>{ref_col}}\n")
    lines.append(f"{'VCF_ALT':>{pad}} {vcf_alt:>{ref_col}}\n")
    lines.append(f"Determination type: {determination_type.name}\n")
    return "".join(lines)


def _resolve_vcf(ref_base: str, allele1: str, allele2: str) -> tuple[str, str, DeterminationType]:
    """Determine VCF REF/ALT for a SNP."""
    if ref_base == allele1:
        return allele1, allele2, DeterminationType.SNP_ALIGNED
    if ref_base == allele2:
        return allele2, allele1, DeterminationType.SNP_ALIGNED
    return ref_base, f"{allele1}/{allele2}", DeterminationType.SNP_ALIGNED


def _is_indel_seq(s: Optional[str]) -> bool:
    """True if s looks like an actual DNA sequence, not a symbolic placeholder."""
    return s is not None and s not in ("-", "I", "D")


def _resolve_indel_vcf(
    record: ManifestRecord,
    strand: GenomicStrand,
    ref_seq: str,
    ref_pos_0: int,
) -> tuple[str, str, int, DeterminationType, Optional[bool]]:
    """
    Build proper VCF REF/ALT for an indel using the flanking sequence and reference.

    Returns (vcf_ref, vcf_alt, anchor_pos_0, determination_type, indel_ref_is_del).
    indel_ref_is_del: True = deletion allele matches reference (no inserted bases).
    Falls back to ref_base/ref_base when sequence information is unavailable.
    """
    ref_base = ref_seq[ref_pos_0] if ref_pos_0 < len(ref_seq) else "N"

    if ref_pos_0 == 0:
        return ref_base, ref_base, ref_pos_0, DeterminationType.INDEL_NO_ANCHOR, None

    anchor_pos_0 = ref_pos_0 - 1
    anchor = ref_seq[anchor_pos_0]

    fa = record.first_allele
    sa = record.second_allele

    if fa == "-" and _is_indel_seq(sa):
        ins_seq = sa.upper() if strand == GenomicStrand.PLUS else reverse_complement(sa.upper())
        return anchor, anchor + ins_seq, anchor_pos_0, DeterminationType.INDEL_INSERTION, True

    if sa == "-" and _is_indel_seq(fa):
        del_len = len(fa)
        true_anchor_pos_0 = ref_pos_0 - del_len
        if true_anchor_pos_0 < 0:
            return ref_base, ref_base, ref_pos_0, DeterminationType.INDEL_NO_ANCHOR, None
        true_anchor = ref_seq[true_anchor_pos_0]
        del_seq_ref = ref_seq[true_anchor_pos_0 + 1 : true_anchor_pos_0 + 1 + del_len]
        return true_anchor + del_seq_ref, true_anchor, true_anchor_pos_0, DeterminationType.INDEL_DELETION, False

    return ref_base, ref_base, ref_pos_0, DeterminationType.INDEL_SITE, None


class VariantAligner:
    """Aligns variant flanking sequences to a reference genome using minimap2."""

    def __init__(self, reference_path: str):
        self._ref_path = reference_path
        self._aligner = mappy.Aligner(reference_path, preset="sr", best_n=5)
        if not self._aligner:
            raise RuntimeError(f"Failed to load reference: {reference_path}")
        self._ref_seqs = _load_fasta(reference_path)

    def align(self, record: ManifestRecord) -> AlignmentResult:
        query, n_pos = _prepare_query(record.flanking)
        if n_pos < 0:
            return AlignmentResult()

        hits = list(self._aligner.map(query))
        if not hits:
            return AlignmentResult()

        hit = max(hits, key=lambda h: h.mapq)
        chrom = hit.ctg
        strand = GenomicStrand.PLUS if hit.strand == 1 else GenomicStrand.MINUS

        ref_pos_0 = _query_pos_to_ref_pos(n_pos, hit)
        if ref_pos_0 is None:
            return AlignmentResult(chromosome=chrom, strand=strand)

        ref_seq = self._ref_seqs.get(chrom, "")
        if ref_pos_0 >= len(ref_seq):
            return AlignmentResult(chromosome=chrom, strand=strand)

        ref_base = ref_seq[ref_pos_0]

        # Alleles on the forward strand (for VCF and PLUS encoding)
        allele1: Optional[str] = None
        allele2: Optional[str] = None
        if record.first_allele is not None and record.second_allele is not None:
            if strand == GenomicStrand.PLUS:
                allele1 = record.first_allele.upper()
                allele2 = record.second_allele.upper()
            else:
                allele1 = reverse_complement(record.second_allele.upper())
                allele2 = reverse_complement(record.first_allele.upper())

        if record.is_indel:
            vcf_ref, vcf_alt, pos_0, determination, indel_ref_is_del = _resolve_indel_vcf(
                record, strand, ref_seq, ref_pos_0
            )
            variant_type = "INDEL"
        else:
            if allele1 is not None and allele2 is not None:
                vcf_ref, vcf_alt, determination = _resolve_vcf(ref_base, allele1, allele2)
            else:
                vcf_ref, vcf_alt, determination = ref_base, None, DeterminationType.SNP_NO_ALLELES
            pos_0 = ref_pos_0
            indel_ref_is_del = None
            variant_type = "SNP"

        q_aln, r_aln = _build_alignment_strings(query, ref_seq, hit)
        aln_text = _format_alignment(
            record=record,
            query=query,
            q_aln=q_aln,
            r_aln=r_aln,
            r_start=hit.r_st,
            ref_pos_0=ref_pos_0,
            ref_base=ref_base,
            allele1=allele1 or "",
            allele2=allele2 or "",
            vcf_ref=vcf_ref,
            vcf_alt=vcf_alt or "",
            variant_type=variant_type,
            determination_type=determination,
        )

        return AlignmentResult(
            chromosome=chrom,
            position=pos_0 + 1,
            strand=strand,
            vcf_ref=vcf_ref,
            vcf_alt=vcf_alt,
            alignment_text=aln_text,
            determination_type=determination,
            indel_ref_is_del=indel_ref_is_del,
        )


# Module-level state for multiprocessing workers
_worker_aligner: Optional[VariantAligner] = None


def _init_worker(reference_path: str) -> None:
    global _worker_aligner
    _worker_aligner = VariantAligner(reference_path)


def _align_record(record: ManifestRecord) -> AlignmentResult:
    assert _worker_aligner is not None
    return _worker_aligner.align(record)
