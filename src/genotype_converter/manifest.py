from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import Optional


class ManifestType(Enum):
    ILLUMINA = auto()
    AFFYMETRIX = auto()


class IllmnStrand(Enum):
    TOP = auto()
    BOT = auto()


@dataclass
class ManifestRecord:
    name: str
    alt_name: Optional[str]
    ab_a: str
    ab_b: str
    flanking: str
    first_allele: Optional[str]
    second_allele: Optional[str]
    manifest_type: ManifestType
    is_indel: bool
    ilmn_strand: Optional[IllmnStrand] = None
    source_strand: Optional[IllmnStrand] = None
    ilmn_id_tbpm_fru: Optional[str] = None
    allele_a_probe_seq: Optional[str] = None
    allele_b_probe_seq: Optional[str] = None

    @property
    def is_snp(self) -> bool:
        return not self.is_indel


def _clean(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    value = value.strip().strip('"\'')
    if not value:
        return None
    value = re.sub(r"\s+", "_", value)
    return value


def _find_header_and_rows(lines: list[str]) -> tuple[Optional[int], Optional[ManifestType]]:
    """Return (header_line_index, detected_format) or (None, None)."""
    for i, line in enumerate(lines):
        stripped = line.strip().strip('"')
        if stripped.startswith("IlmnID,Name,IlmnStrand"):
            return i, ManifestType.ILLUMINA
        if re.match(r'^"?Probe Set ID"?[,\t]"?Affy SNP ID"?', stripped):
            return i, ManifestType.AFFYMETRIX
        if stripped.startswith("Probe_Set_ID,Affy_SNP_ID"):
            return i, ManifestType.AFFYMETRIX
    return None, None


def _read_csv_rows(path: str) -> tuple[list[dict], Optional[ManifestType]]:
    text = Path(path).read_text(encoding="utf-8-sig")
    lines = text.splitlines()

    header_idx, fmt = _find_header_and_rows(lines)
    if header_idx is None:
        return [], None

    header_line = lines[header_idx]
    delim = "\t" if "\t" in header_line else ","

    data_lines: list[str] = []
    for line in lines[header_idx + 1 :]:
        if re.match(r"^\[Controls\]", line.strip()):
            break
        if line.strip():
            data_lines.append(line)

    reader = csv.DictReader([header_line] + data_lines, delimiter=delim)
    rows = []
    for row in reader:
        cleaned = {
            re.sub(r"\s+", "_", k.strip().strip('"')): _clean(v)
            for k, v in row.items()
            if k is not None
        }
        rows.append(cleaned)
    return rows, fmt


def _parse_alleles_from_snp(snp: str) -> tuple[str, str]:
    m = re.match(r"\[([^/\]]+)/([^/\]]+)\]", snp)
    if not m:
        raise ValueError(f"Cannot parse alleles from SNP field: {snp!r}")
    return m.group(1), m.group(2)


def _parse_alleles_from_flanking(flanking: str) -> tuple[Optional[str], Optional[str]]:
    m = re.search(r"\[([^/\]]*)/([^/\]]*)\]", flanking)
    if m:
        return m.group(1) or None, m.group(2) or None
    return None, None


def _parse_ilmn_strand(s: Optional[str]) -> Optional[IllmnStrand]:
    if s == "TOP":
        return IllmnStrand.TOP
    if s == "BOT":
        return IllmnStrand.BOT
    return None


def _parse_tbpm_fru(ilmn_id: str, name: str) -> Optional[str]:
    """Extract the T/B/M/P _ F/R/U probe strand code from IlmnID."""
    m = re.search(re.escape(name) + r"-\d_([TBMP]_[FRU])", ilmn_id)
    if m:
        return m.group(1)
    m = re.search(r"([TBMP]_[FRU])_\d+$", ilmn_id)
    if m:
        return m.group(1)
    m = re.search(r"([TBMP]_[FRU])", ilmn_id)
    if m:
        return m.group(1)
    return None


def _parse_illumina_rows(rows: list[dict]) -> list[ManifestRecord]:
    records = []
    for row in rows:
        source_seq = row.get("SourceSeq")
        if source_seq is None:
            continue

        name = row.get("Name")
        if not name:
            continue

        snp = row.get("SNP")
        if not snp:
            continue

        ilmn_id = row.get("IlmnID")
        ilmn_strand = _parse_ilmn_strand(row.get("IlmnStrand"))
        source_strand = _parse_ilmn_strand(row.get("SourceStrand"))

        ab_a, ab_b = _parse_alleles_from_snp(snp)
        first_allele, second_allele = _parse_alleles_from_flanking(source_seq)
        tbpm_fru = _parse_tbpm_fru(ilmn_id, name) if ilmn_id else None
        is_indel = ab_a in ("I", "D")

        records.append(
            ManifestRecord(
                name=name,
                alt_name=ilmn_id,
                ab_a=ab_a,
                ab_b=ab_b,
                flanking=source_seq,
                first_allele=first_allele,
                second_allele=second_allele,
                manifest_type=ManifestType.ILLUMINA,
                is_indel=is_indel,
                ilmn_strand=ilmn_strand,
                source_strand=source_strand,
                ilmn_id_tbpm_fru=tbpm_fru,
                allele_a_probe_seq=row.get("AlleleA_ProbeSeq"),
                allele_b_probe_seq=row.get("AlleleB_ProbeSeq"),
            )
        )
    return records


def _parse_affymetrix_rows(rows: list[dict]) -> list[ManifestRecord]:
    records = []
    for row in rows:
        affy_snp_id = row.get("Affy_SNP_ID")
        if affy_snp_id is None:
            continue

        probe_set_id = row.get("Probe_Set_ID") or affy_snp_id
        ab_a = row.get("Allele_A")
        ab_b = row.get("Allele_B")
        if ab_a is None or ab_b is None:
            continue

        flank = row.get("Flank") or ""
        first_allele, second_allele = _parse_alleles_from_flanking(flank)
        is_indel = ab_a in ("I", "D")

        records.append(
            ManifestRecord(
                name=probe_set_id,
                alt_name=affy_snp_id,
                ab_a=ab_a,
                ab_b=ab_b,
                flanking=flank,
                first_allele=first_allele,
                second_allele=second_allele,
                manifest_type=ManifestType.AFFYMETRIX,
                is_indel=is_indel,
            )
        )
    return records


def parse_manifest(path: str) -> list[ManifestRecord]:
    """Auto-detect Illumina vs Affymetrix and return parsed records."""
    rows, fmt = _read_csv_rows(path)
    if fmt == ManifestType.AFFYMETRIX:
        return _parse_affymetrix_rows(rows)
    if fmt == ManifestType.ILLUMINA:
        return _parse_illumina_rows(rows)
    raise ValueError(
        f"Manifest format not recognised for {path!r}. Expected an Illumina "
        "manifest with an IlmnID/Name/IlmnStrand header or an Affymetrix "
        "manifest with Probe Set ID/Affy SNP ID columns."
    )
