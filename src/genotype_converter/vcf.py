from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Any
from urllib.parse import quote

from .conversion import VariantResult


VCF_REQUIRED_COLUMNS = [
    "marker_name",
    "chromosome",
    "position",
    "ref_allele",
    "alt_allele",
]


@dataclass(frozen=True)
class SiteVcfStats:
    output_path: str
    records_total: int
    records_written: int
    records_skipped: int


def rows_from_results(results: Iterable[VariantResult]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for result in results:
        rows.append(
            {
                "marker_name": result.name,
                "alt_marker_name": result.alt_name or "",
                "chromosome": result.chromosome or "",
                "position": result.position or "",
                "ref_allele": result.vcf_ref or "",
                "alt_allele": result.vcf_alt or "",
                "determination_type": result.determination_type or "",
                "A_in_PLUS": result.plus_a or "",
                "B_in_PLUS": result.plus_b or "",
                "A_vcf": result.vcf_a or "",
                "B_vcf": result.vcf_b or "",
            }
        )
    return rows


def rows_from_lookup_csv(lookup_path: str) -> list[dict[str, str]]:
    with open(lookup_path, newline="") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    reader = csv.DictReader(lines)
    if reader.fieldnames is None:
        raise ValueError(f"Lookup file has no header row: {lookup_path}")
    missing = [column for column in VCF_REQUIRED_COLUMNS if column not in reader.fieldnames]
    if missing:
        raise ValueError(
            "Lookup file is missing required VCF export column(s): "
            f"{', '.join(missing)}"
        )
    return [dict(row) for row in reader if row.get("marker_name")]


def write_site_vcf(
    rows: Iterable[Mapping[str, Any]],
    path: str,
    *,
    source_name: str | None = None,
) -> SiteVcfStats:
    """Write a site-only VCF with nucleotide REF/ALT alleles and no samples."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    written = 0

    with output_path.open("w", newline="") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=genotype_converter\n")
        if source_name:
            handle.write(f"##genotype_converter_source={_meta_value(source_name)}\n")
        handle.write(
            '##INFO=<ID=ALT_ID,Number=1,Type=String,Description="Alternate marker identifier from the manifest">\n'
        )
        handle.write(
            '##INFO=<ID=DETERMINATION_TYPE,Number=1,Type=String,Description="How the marker position and alleles were resolved">\n'
        )
        handle.write(
            '##INFO=<ID=A_ALLELE,Number=1,Type=String,Description="Nucleotide allele represented by manifest allele A when resolvable">\n'
        )
        handle.write(
            '##INFO=<ID=B_ALLELE,Number=1,Type=String,Description="Nucleotide allele represented by manifest allele B when resolvable">\n'
        )
        handle.write(
            '##INFO=<ID=A_VCF,Number=1,Type=String,Description="Whether manifest allele A is REF or ALT">\n'
        )
        handle.write(
            '##INFO=<ID=B_VCF,Number=1,Type=String,Description="Whether manifest allele B is REF or ALT">\n'
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in rows:
            total += 1
            record = _site_record(row)
            if record is None:
                continue
            handle.write("\t".join(record) + "\n")
            written += 1

    return SiteVcfStats(
        output_path=str(output_path),
        records_total=total,
        records_written=written,
        records_skipped=total - written,
    )


def _site_record(row: Mapping[str, Any]) -> list[str] | None:
    chrom = _text(row.get("chromosome"))
    pos = _text(row.get("position"))
    marker_id = _text(row.get("marker_name")) or "."
    ref = _text(row.get("ref_allele")).upper()
    alt = _normalize_alt_alleles(ref, _text(row.get("alt_allele")))
    if not (chrom and pos and ref and alt):
        return None
    info = _info_field(
        {
            "ALT_ID": _text(row.get("alt_marker_name")),
            "DETERMINATION_TYPE": _text(row.get("determination_type")),
            "A_ALLELE": _manifest_allele(row, "A", ref, alt),
            "B_ALLELE": _manifest_allele(row, "B", ref, alt),
            "A_VCF": _text(row.get("A_vcf")),
            "B_VCF": _text(row.get("B_vcf")),
        }
    )
    return [chrom, pos, marker_id, ref, alt, ".", "PASS", info]


def _normalize_alt_alleles(ref: str, alt: str) -> str:
    if not alt:
        return ""
    raw_parts = []
    for comma_part in alt.replace("/", ",").split(","):
        part = comma_part.strip().upper()
        if part and part != ".":
            raw_parts.append(part)
    unique_parts: list[str] = []
    for part in raw_parts:
        if part == ref:
            continue
        if part not in unique_parts:
            unique_parts.append(part)
    return ",".join(unique_parts)


def _manifest_allele(row: Mapping[str, Any], allele: str, ref: str, alt: str) -> str:
    vcf_label = _text(row.get(f"{allele}_vcf")).upper()
    if vcf_label == "REF":
        return ref
    if vcf_label != "ALT":
        return ""
    alt_parts = [part for part in alt.split(",") if part]
    plus_value = _text(row.get(f"{allele}_in_PLUS")).upper()
    if plus_value and plus_value in alt_parts:
        return plus_value
    if len(alt_parts) == 1:
        return alt_parts[0]
    return ""


def _info_field(values: Mapping[str, str]) -> str:
    fields = [
        f"{key}={_info_value(value)}"
        for key, value in values.items()
        if value
    ]
    return ";".join(fields) if fields else "."


def _text(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _info_value(value: str) -> str:
    return quote(value, safe="._:-")


def _meta_value(value: str) -> str:
    return quote(value, safe="._:-/")
