from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from enum import Enum, auto
from typing import Callable, Optional

import mappy

from .manifest import ManifestRecord


class GenomicStrand(Enum):
    PLUS = auto()
    MINUS = auto()


class DeterminationType(Enum):
    SNP_ALIGNED = auto()       # ref base matched one flanking allele → REF/ALT assigned
    SNP_NO_ALLELES = auto()    # flanking had no [X/Y] bracket; only ref base available
    SNP_PROBE_ADJACENT = auto()
    SNP_ALLELE_ASSISTED = auto()
    SNP_AMBIGUOUS = auto()
    INDEL_INSERTION = auto()   # insertion [-/SEQ] resolved with anchor base
    INDEL_DELETION = auto()    # deletion [SEQ/-] resolved with anchor base
    INDEL_SITE = auto()        # symbolic I/D or unrecognised; VCF alleles unavailable
    INDEL_NO_ANCHOR = auto()   # variant at position 0; no anchor base possible


class ProbeSide(Enum):
    LEFT = auto()
    RIGHT = auto()


@dataclass(frozen=True)
class ProbePlacement:
    sequence: str
    side: ProbeSide
    includes_allele: bool


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


def _prepare_query(flanking: str) -> tuple[str, int]:
    """
    Build alignment query from flanking sequence.
    The bracketed variant notation is replaced with a synthetic N and its
    position is tracked explicitly. Other N/n bases in the flanking sequence
    are preserved because they are part of the manifest's sequence context.

    Returns (query_uppercase, n_position_0based) or (query, -1) if no N found.
    """
    def clean(seq: str) -> str:
        return re.sub(r"[^GATCNgatcn]", "", seq).upper()

    m = re.search(r"\[.*?\]", flanking)
    if m:
        left = clean(flanking[: m.start()])
        right = clean(flanking[m.end() :])
        return left + "N" + right, len(left)

    seq = clean(flanking)
    n_pos = seq.find("N")
    return seq, n_pos


def _clean_probe_context(seq: str) -> str:
    seq = re.sub(r"[Nn]", "", seq)
    seq = re.sub(r"[^GATCgatc]", "", seq)
    return seq.upper()


def _flanking_sides(flanking: str) -> tuple[str, str]:
    m = re.search(r"\[.*?\]", flanking)
    if not m:
        return _clean_probe_context(flanking), ""
    return (
        _clean_probe_context(flanking[: m.start()]),
        _clean_probe_context(flanking[m.end() :]),
    )


def _matches_left_probe_side(probe: str, left: str) -> bool:
    return left.endswith(probe) or probe.endswith(left)


def _matches_right_probe_side(probe: str, right: str) -> bool:
    return right.startswith(probe) or probe.startswith(right)


def _probe_placements(record: ManifestRecord) -> list[ProbePlacement]:
    left, right = _flanking_sides(record.flanking)
    alleles = {
        allele.upper()
        for allele in (record.first_allele, record.second_allele)
        if allele and len(allele) == 1
    }
    placements: list[ProbePlacement] = []
    seen: set[tuple[str, ProbeSide, bool]] = set()

    def add(sequence: str, side: ProbeSide, includes_allele: bool) -> None:
        key = (sequence, side, includes_allele)
        if key not in seen:
            seen.add(key)
            placements.append(
                ProbePlacement(
                    sequence=sequence,
                    side=side,
                    includes_allele=includes_allele,
                )
            )

    for probe in (record.allele_a_probe_seq, record.allele_b_probe_seq):
        if not probe:
            continue
        cleaned = _clean_probe_context(probe)
        if not cleaned:
            continue
        candidates = (cleaned, reverse_complement(cleaned))
        for candidate in candidates:
            if _matches_left_probe_side(candidate, left):
                add(candidate, ProbeSide.LEFT, False)
            if _matches_right_probe_side(candidate, right):
                add(candidate, ProbeSide.RIGHT, False)
            if (
                alleles
                and len(candidate) > 1
                and candidate[-1] in alleles
                and _matches_left_probe_side(candidate[:-1], left)
            ):
                add(candidate, ProbeSide.LEFT, True)
            if (
                alleles
                and len(candidate) > 1
                and candidate[0] in alleles
                and _matches_right_probe_side(candidate[1:], right)
            ):
                add(candidate, ProbeSide.RIGHT, True)
    return placements


def _probe_side(record: ManifestRecord) -> Optional[ProbeSide]:
    sides = {placement.side for placement in _probe_placements(record)}
    if len(sides) == 1:
        return next(iter(sides))
    return None


def _flip_probe_side(side: ProbeSide) -> ProbeSide:
    return ProbeSide.RIGHT if side == ProbeSide.LEFT else ProbeSide.LEFT


def _aligned_probe_placements(record: ManifestRecord, hit) -> list[ProbePlacement]:
    placements = _probe_placements(record)
    if hit.strand == 1:
        return placements
    return [
        ProbePlacement(
            sequence=reverse_complement(placement.sequence),
            side=_flip_probe_side(placement.side),
            includes_allele=placement.includes_allele,
        )
        for placement in placements
    ]


def _aligned_probe_side(
    record: ManifestRecord,
    hit,
) -> Optional[ProbeSide]:
    sides = {placement.side for placement in _aligned_probe_placements(record, hit)}
    if len(sides) == 1:
        return next(iter(sides))
    return None


def _query_probe_matches(
    placement: ProbePlacement,
    q_ungapped: str,
    n_index: int,
) -> list[tuple[int, int, int]]:
    """Return (start, end, assayed_query_index) matches for a probe placement."""
    seq = placement.sequence
    matches: list[tuple[int, int, int]] = []

    if placement.includes_allele:
        if placement.side == ProbeSide.LEFT:
            context = seq[:-1]
            start = q_ungapped.find(context)
            while start != -1:
                end = start + len(context)
                if end == n_index:
                    matches.append((start, end, n_index))
                start = q_ungapped.find(context, start + 1)
        else:
            context = seq[1:]
            start = q_ungapped.find(context)
            while start != -1:
                end = start + len(context)
                if start == n_index + 1:
                    matches.append((start, end, n_index))
                start = q_ungapped.find(context, start + 1)
        return matches

    start = q_ungapped.find(seq)
    while start != -1:
        end = start + len(seq)
        if placement.side == ProbeSide.LEFT and end <= n_index:
            matches.append((start, end, end))
        elif placement.side == ProbeSide.RIGHT and start > n_index:
            matches.append((start, end, start - 1))
        start = q_ungapped.find(seq, start + 1)
    return matches


def _alignment_col_to_ref_pos(
    ref_aln: str,
    r_start: int,
    col: int,
) -> Optional[int]:
    if col < 0 or col >= len(ref_aln) or ref_aln[col] == "-":
        return None
    return r_start + sum(1 for ch in ref_aln[:col] if ch != "-")


def _ref_index_to_alignment_col(ref_aln: str, ref_index: int) -> Optional[int]:
    ref_count = 0
    for col, base in enumerate(ref_aln):
        if base == "-":
            continue
        if ref_count == ref_index:
            return col
        ref_count += 1
    return None


def _query_index_to_alignment_col(q_aln: str, q_index: int) -> Optional[int]:
    query_count = 0
    for col, base in enumerate(q_aln):
        if base == "-":
            continue
        if query_count == q_index:
            return col
        query_count += 1
    return None


def _probe_adjacent_ref_pos(
    record: ManifestRecord,
    q_aln: str,
    r_aln: str,
    r_start: int,
    probe_side: ProbeSide,
    n_index: int,
    placements: Optional[list[ProbePlacement]] = None,
) -> Optional[int]:
    q_ungapped = q_aln.replace("-", "")
    if n_index < 0 or n_index >= len(q_ungapped):
        return None

    ungapped_to_col = [
        col for col, base in enumerate(q_aln)
        if base != "-"
    ]
    if len(ungapped_to_col) != len(q_ungapped):
        return None

    probe_matches: list[tuple[int, int, int]] = []
    for placement in placements if placements is not None else _probe_placements(record):
        if placement.side != probe_side:
            continue
        probe_matches.extend(_query_probe_matches(placement, q_ungapped, n_index))

    if not probe_matches:
        return None

    if probe_side == ProbeSide.LEFT:
        _start, end, assayed_index = max(probe_matches, key=lambda match: match[1])
    else:
        start, _end, assayed_index = min(probe_matches, key=lambda match: match[0])

    assayed_col = ungapped_to_col[assayed_index]
    ref_pos = _alignment_col_to_ref_pos(r_aln, r_start, assayed_col)
    if ref_pos is not None:
        return ref_pos

    if probe_side == ProbeSide.LEFT:
        search_range = range(assayed_col + 1, len(r_aln))
    else:
        search_range = range(assayed_col - 1, -1, -1)

    for col in search_range:
        ref_pos = _alignment_col_to_ref_pos(r_aln, r_start, col)
        if ref_pos is not None:
            return ref_pos
    return None


def _reference_probe_adjacent_ref_pos(
    record: ManifestRecord,
    q_aln: str,
    r_aln: str,
    r_start: int,
    probe_side: ProbeSide,
    n_index: int,
    candidate_ref_pos: Optional[int],
    placements: Optional[list[ProbePlacement]] = None,
) -> Optional[int]:
    q_n_col = _query_index_to_alignment_col(q_aln, n_index)
    candidate_ref_index = (
        candidate_ref_pos - r_start if candidate_ref_pos is not None else None
    )
    if q_n_col is None and candidate_ref_index is None:
        return None

    r_ungapped = r_aln.replace("-", "")
    scored_matches: list[tuple[int, int]] = []
    for placement in placements if placements is not None else _probe_placements(record):
        if placement.side != probe_side:
            continue
        probe = placement.sequence
        start = r_ungapped.find(probe)
        while start != -1:
            end = start + len(probe)
            if placement.includes_allele:
                assayed_ref_index = end - 1 if probe_side == ProbeSide.LEFT else start
                probe_edge_index = assayed_ref_index
            elif probe_side == ProbeSide.LEFT:
                assayed_ref_index = end
                probe_edge_index = end - 1
            else:
                assayed_ref_index = start - 1
                probe_edge_index = start

            probe_edge_col = _ref_index_to_alignment_col(r_aln, probe_edge_index)

            if 0 <= assayed_ref_index < len(r_ungapped):
                distances = []
                if candidate_ref_index is not None:
                    distances.append(abs(assayed_ref_index - candidate_ref_index))
                if q_n_col is not None and probe_edge_col is not None:
                    distances.append(abs(probe_edge_col - q_n_col))
                if distances:
                    scored_matches.append((min(distances), assayed_ref_index))
            start = r_ungapped.find(probe, start + 1)

    if not scored_matches:
        return None

    best_distance = min(distance for distance, _idx in scored_matches)
    best_indices = {
        ref_index for distance, ref_index in scored_matches
        if distance == best_distance
    }
    if len(best_indices) == 1:
        return r_start + next(iter(best_indices))
    return None


def _aligned_query_variant_index(n_pos: int, hit) -> Optional[int]:
    if hit.strand == 1:
        if not (hit.q_st <= n_pos < hit.q_en):
            return None
        return n_pos - hit.q_st
    rc_offset = hit.q_en - 1 - n_pos
    if not (0 <= rc_offset < hit.q_en - hit.q_st):
        return None
    return rc_offset


def _score_bases(query_base: str, ref_base: str) -> int:
    if query_base == "N" or ref_base == "N":
        return 0
    return 2 if query_base == ref_base else -1


def _semiglobal_align_query_to_ref(query: str, ref: str) -> tuple[str, str, int]:
    """
    Align all query bases to the best substring of ref.

    Reference overhangs are free; query overhangs are not. Returns
    (query_alignment, ref_alignment, ref_start_offset).
    """
    gap = -2
    n = len(query)
    m = len(ref)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[""] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap
        trace[i][0] = "U"
    for j in range(1, m + 1):
        score[0][j] = 0
        trace[0][j] = "L"

    for i in range(1, n + 1):
        qb = query[i - 1]
        for j in range(1, m + 1):
            rb = ref[j - 1]
            diag = score[i - 1][j - 1] + _score_bases(qb, rb)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap
            best = max(diag, up, left)
            score[i][j] = best
            if best == diag:
                trace[i][j] = "D"
            elif best == up:
                trace[i][j] = "U"
            else:
                trace[i][j] = "L"

    j = max(range(m + 1), key=lambda col: score[n][col])
    i = n
    q_out: list[str] = []
    r_out: list[str] = []
    while i > 0:
        move = trace[i][j]
        if move == "D":
            q_out.append(query[i - 1])
            r_out.append(ref[j - 1])
            i -= 1
            j -= 1
        elif move == "U":
            q_out.append(query[i - 1])
            r_out.append("-")
            i -= 1
        else:
            q_out.append("-")
            r_out.append(ref[j - 1])
            j -= 1

    return "".join(reversed(q_out)), "".join(reversed(r_out)), j


def _choose_probe_ref_pos(
    probe_ref_pos: int,
    direct_ref_pos: Optional[int],
    allele1: Optional[str],
    allele2: Optional[str],
    fetch_ref_base: Callable[[int], str],
) -> tuple[int, DeterminationType]:
    if direct_ref_pos is None or direct_ref_pos == probe_ref_pos:
        return probe_ref_pos, DeterminationType.SNP_PROBE_ADJACENT

    alleles = {
        allele.upper()
        for allele in (allele1, allele2)
        if allele and len(allele) == 1
    }
    if not alleles:
        return probe_ref_pos, DeterminationType.SNP_AMBIGUOUS

    candidates: list[tuple[int, str, bool]] = []
    for pos in dict.fromkeys((probe_ref_pos, direct_ref_pos)):
        ref_base = fetch_ref_base(pos).upper()
        if ref_base:
            candidates.append((pos, ref_base, ref_base in alleles))

    compatible = [pos for pos, _base, is_compatible in candidates if is_compatible]
    if len(compatible) == 1:
        return compatible[0], DeterminationType.SNP_ALLELE_ASSISTED
    return probe_ref_pos, DeterminationType.SNP_AMBIGUOUS


def _refine_snp_with_local_alignment(
    record: ManifestRecord,
    query: str,
    n_pos: int,
    chrom: str,
    hit,
    strand: GenomicStrand,
    probe_side: Optional[ProbeSide],
    direct_ref_pos: Optional[int],
    allele1: Optional[str],
    allele2: Optional[str],
    fetch_ref: Callable[[str, int, int], str],
    placements: Optional[list[ProbePlacement]] = None,
) -> tuple[Optional[int], Optional[DeterminationType], Optional[str], Optional[str], int]:
    if probe_side is None:
        return None, None, None, None, hit.r_st

    oriented_query = query if strand == GenomicStrand.PLUS else reverse_complement(query)
    oriented_n_index = (
        n_pos if strand == GenomicStrand.PLUS else len(query) - 1 - n_pos
    )
    pad = 50
    window_start = max(0, hit.r_st - pad)
    window_end = hit.r_en + pad
    ref_window = fetch_ref(chrom, window_start, window_end)
    if not ref_window:
        return None, None, None, None, hit.r_st

    q_aln, r_aln, offset = _semiglobal_align_query_to_ref(oriented_query, ref_window)
    local_r_start = window_start + offset
    ref_probe_pos = _reference_probe_adjacent_ref_pos(
        record,
        q_aln,
        r_aln,
        local_r_start,
        probe_side,
        oriented_n_index,
        direct_ref_pos,
        placements=placements,
    )
    query_probe_pos = _probe_adjacent_ref_pos(
        record,
        q_aln,
        r_aln,
        local_r_start,
        probe_side,
        oriented_n_index,
        placements=placements,
    )

    probe_ref_pos = ref_probe_pos if ref_probe_pos is not None else query_probe_pos
    if probe_ref_pos is None:
        return None, None, q_aln, r_aln, local_r_start

    if ref_probe_pos is not None:
        return (
            ref_probe_pos,
            DeterminationType.SNP_PROBE_ADJACENT,
            q_aln,
            r_aln,
            local_r_start,
        )

    ref_pos, determination = _choose_probe_ref_pos(
        probe_ref_pos,
        direct_ref_pos,
        allele1,
        allele2,
        lambda pos: fetch_ref(chrom, pos, pos + 1),
    )
    return ref_pos, determination, q_aln, r_aln, local_r_start


def _cigar_query_to_ref(
    q_offset: int,
    r_start: int,
    cigar,
    probe_side: Optional[ProbeSide] = None,
) -> Optional[int]:
    """Walk a mappy cigar list to map a 0-based query offset to a 0-based reference position."""
    q = 0
    r = r_start
    prev_op = None
    prev_len = 0
    for length, op in cigar:
        if op in (4, 5):  # S, H
            prev_op = op
            prev_len = length
            continue
        if op in (0, 7, 8):  # M, =, X
            if q + length > q_offset:
                if probe_side == ProbeSide.LEFT and prev_op == 2:
                    gap_delta = q_offset - q
                    if 0 <= gap_delta <= prev_len:
                        return max(r_start, r - (prev_len - gap_delta))
                return r + (q_offset - q)
            q += length
            r += length
        elif op == 1:  # I – insertion in query
            if q + length > q_offset:
                if probe_side == ProbeSide.RIGHT:
                    return max(r_start, r - 1)
                return r
            q += length
        elif op == 2:  # D – deletion in query
            r += length
        prev_op = op
        prev_len = length
    return None


def _query_pos_to_ref_pos(
    n_pos: int,
    hit,
    probe_side: Optional[ProbeSide] = None,
) -> Optional[int]:
    """Map a 0-based query position to a 0-based reference position via CIGAR."""
    if hit.strand == 1:
        if not (hit.q_st <= n_pos < hit.q_en):
            return None
        return _cigar_query_to_ref(
            n_pos - hit.q_st,
            hit.r_st,
            hit.cigar,
            probe_side,
        )
    else:
        rc_offset = hit.q_en - 1 - n_pos
        if not (0 <= rc_offset < hit.q_en - hit.q_st):
            return None
        return _cigar_query_to_ref(
            rc_offset,
            hit.r_st,
            hit.cigar,
            probe_side,
        )


def _cigar_has_gap_near_query_offset(q_offset: int, cigar, radius: int = 6) -> bool:
    q = 0
    for length, op in cigar:
        if op in (4, 5):
            continue
        if op in (0, 7, 8):
            q += length
        elif op == 1:
            if q - radius <= q_offset <= q + length + radius:
                return True
            q += length
        elif op == 2:
            if q - radius <= q_offset <= q + radius:
                return True
    return False


def _hit_has_gap_near_variant(n_pos: int, hit, radius: int = 6) -> bool:
    q_offset = _aligned_query_variant_index(n_pos, hit)
    if q_offset is None:
        return False
    return _cigar_has_gap_near_query_offset(q_offset, hit.cigar, radius)


def _build_alignment_strings(query: str, ref_sub: str, hit) -> tuple[str, str]:
    """Reconstruct gapped alignment strings from mappy hit + reference sequence."""
    if hit.strand == 1:
        q_bases = query[hit.q_st : hit.q_en]
    else:
        q_bases = reverse_complement(query[hit.q_st : hit.q_en])

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


def _ungapped_alignment_map(aln: str) -> list[int]:
    return [col for col, base in enumerate(aln) if base != "-"]


def _format_probe_alignment_line(
    q_aln: str,
    r_aln: str,
    placements: list[ProbePlacement],
) -> Optional[str]:
    if not placements:
        return None

    candidates: list[tuple[int, int, str, list[int]]] = []
    for priority, aln in enumerate((r_aln, q_aln)):
        ungapped = aln.replace("-", "")
        col_map = _ungapped_alignment_map(aln)
        for placement in placements:
            seq = placement.sequence
            start = ungapped.find(seq)
            while start != -1:
                end = start + len(seq)
                if end <= len(col_map):
                    candidates.append((len(seq), -priority, seq, col_map[start:end]))
                start = ungapped.find(seq, start + 1)

    if not candidates:
        return None

    _length, _priority, seq, cols = max(candidates, key=lambda item: (item[0], item[1]))
    line = [" "] * len(q_aln)
    for base, col in zip(seq, cols):
        if 0 <= col < len(line):
            line[col] = base
    return "".join(line).rstrip()


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
    probe_placements: Optional[list[ProbePlacement]] = None,
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
    probe_line = _format_probe_alignment_line(
        q_aln,
        r_aln,
        probe_placements or [],
    )
    if probe_line:
        lines.append(f"{'PROBE':>{pad}} {probe_line}\n")
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


def _resolve_indel_vcf_from_fetch(
    record: ManifestRecord,
    strand: GenomicStrand,
    fetch_ref: Callable[[int, int], str],
    ref_pos_0: int,
) -> tuple[str, str, int, DeterminationType, Optional[bool]]:
    ref_base = fetch_ref(ref_pos_0, ref_pos_0 + 1) or "N"

    if ref_pos_0 == 0:
        return ref_base, ref_base, ref_pos_0, DeterminationType.INDEL_NO_ANCHOR, None

    anchor_pos_0 = ref_pos_0 - 1
    anchor = fetch_ref(anchor_pos_0, anchor_pos_0 + 1) or "N"

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
        vcf_ref = fetch_ref(true_anchor_pos_0, true_anchor_pos_0 + del_len + 1)
        if len(vcf_ref) != del_len + 1:
            return ref_base, ref_base, ref_pos_0, DeterminationType.INDEL_NO_ANCHOR, None
        return vcf_ref, vcf_ref[0], true_anchor_pos_0, DeterminationType.INDEL_DELETION, False

    return ref_base, ref_base, ref_pos_0, DeterminationType.INDEL_SITE, None


class VariantAligner:
    """Aligns variant flanking sequences to a reference genome using minimap2."""

    def __init__(
        self,
        reference_path: str,
        save_alignment: bool = False,
        progress: bool = False,
    ):
        self._ref_path = reference_path
        self._save_alignment = save_alignment
        if progress:
            sys.stderr.write(f"Loading reference index: {reference_path}\n")
            sys.stderr.flush()
        self._aligner = mappy.Aligner(reference_path, preset="sr", best_n=5)
        if not self._aligner:
            raise RuntimeError(f"Failed to load reference: {reference_path}")

    def _ref_slice(self, chrom: str, start: int, end: int) -> str:
        if start < 0 or end <= start:
            return ""
        seq = self._aligner.seq(chrom, start, end)
        return (seq or "").upper()

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

        probe_placements = _aligned_probe_placements(record, hit) if record.is_snp else []
        probe_side = None
        if probe_placements:
            probe_sides = {placement.side for placement in probe_placements}
            if len(probe_sides) == 1:
                probe_side = next(iter(probe_sides))
        ref_sub = self._ref_slice(chrom, hit.r_st, hit.r_en)
        if not ref_sub:
            return AlignmentResult(chromosome=chrom, strand=strand)

        direct_ref_pos_0 = _query_pos_to_ref_pos(
            n_pos,
            hit,
            probe_side=probe_side,
        )
        ref_pos_0 = direct_ref_pos_0
        snp_determination: Optional[DeterminationType] = None
        aln_r_start = hit.r_st
        q_aln = None
        r_aln = None
        if record.is_snp and probe_side is not None:
            q_aln, r_aln = _build_alignment_strings(query, ref_sub, hit)
            q_n_index = _aligned_query_variant_index(n_pos, hit)
            probe_ref_pos_0 = _probe_adjacent_ref_pos(
                record,
                q_aln,
                r_aln,
                hit.r_st,
                probe_side,
                q_n_index if q_n_index is not None else -1,
                placements=probe_placements,
            )
            if probe_ref_pos_0 is not None:
                ref_pos_0, snp_determination = _choose_probe_ref_pos(
                    probe_ref_pos_0,
                    direct_ref_pos_0,
                    allele1,
                    allele2,
                    lambda pos: self._ref_slice(chrom, pos, pos + 1),
                )
            needs_refinement = (
                probe_ref_pos_0 is not None
                and (
                    (
                        direct_ref_pos_0 is not None
                        and probe_ref_pos_0 != direct_ref_pos_0
                    )
                    or _hit_has_gap_near_variant(n_pos, hit)
                )
            )
            if needs_refinement:
                refined_pos, refined_det, refined_q, refined_r, refined_start = (
                    _refine_snp_with_local_alignment(
                        record,
                        query,
                        n_pos,
                        chrom,
                        hit,
                        strand,
                        probe_side,
                        direct_ref_pos_0,
                        allele1,
                        allele2,
                        lambda ref_chrom, start, end: self._ref_slice(ref_chrom, start, end),
                        placements=probe_placements,
                    )
                )
                if refined_pos is not None:
                    ref_pos_0 = refined_pos
                    snp_determination = refined_det
                    q_aln = refined_q
                    r_aln = refined_r
                    aln_r_start = refined_start
        if ref_pos_0 is None:
            return AlignmentResult(chromosome=chrom, strand=strand)

        ref_base = self._ref_slice(chrom, ref_pos_0, ref_pos_0 + 1)
        if not ref_base:
            return AlignmentResult(chromosome=chrom, strand=strand)

        if record.is_indel:
            vcf_ref, vcf_alt, pos_0, determination, indel_ref_is_del = _resolve_indel_vcf_from_fetch(
                record,
                strand,
                lambda start, end: self._ref_slice(chrom, start, end),
                ref_pos_0,
            )
            variant_type = "INDEL"
        else:
            if allele1 is not None and allele2 is not None:
                vcf_ref, vcf_alt, determination = _resolve_vcf(ref_base, allele1, allele2)
            else:
                vcf_ref, vcf_alt, determination = ref_base, None, DeterminationType.SNP_NO_ALLELES
            if snp_determination is not None:
                determination = snp_determination
            pos_0 = ref_pos_0
            indel_ref_is_del = None
            variant_type = "SNP"

        aln_text = None
        if self._save_alignment:
            if q_aln is None or r_aln is None:
                q_aln, r_aln = _build_alignment_strings(query, ref_sub, hit)
                aln_r_start = hit.r_st
            aln_text = _format_alignment(
                record=record,
                query=query,
                q_aln=q_aln,
                r_aln=r_aln,
                r_start=aln_r_start,
                ref_pos_0=ref_pos_0,
                ref_base=ref_base,
                allele1=allele1 or "",
                allele2=allele2 or "",
                vcf_ref=vcf_ref,
                vcf_alt=vcf_alt or "",
                variant_type=variant_type,
                determination_type=determination,
                probe_placements=probe_placements,
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


def _init_worker(
    reference_path: str,
    save_alignment: bool = False,
    progress: bool = False,
) -> None:
    global _worker_aligner
    _worker_aligner = VariantAligner(
        reference_path,
        save_alignment=save_alignment,
        progress=progress,
    )


def _align_record(record: ManifestRecord) -> AlignmentResult:
    assert _worker_aligner is not None
    return _worker_aligner.align(record)
