from __future__ import annotations

import csv
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

from .conversion import VariantResult
from .manifest import ManifestType


@dataclass
class BuildStats:
    manifest_path: str
    manifest_sha256: str
    reference_path: str
    reference_sha256: str
    species: str
    workers: int
    total_markers: int
    n_snp: int
    n_indel: int
    by_manifest_type: dict = field(default_factory=dict)  # e.g. {"illumina": 50000}
    n_positioned: int = 0
    n_not_positioned: int = 0
    by_chromosome: dict = field(default_factory=dict)
    output_files: list = field(default_factory=list)


def _header_lines(info: list[str], file_type: str) -> list[str]:
    ts = datetime.now().strftime("%A %B %d %H:%M:%S %Y")
    lines = [f"#{i}" for i in info]
    lines += [
        "#",
        f"#{file_type} generated on {ts}.",
        "#Using genotype_converter, written by Paul Stothard, stothard@ualberta.ca.",
        "#",
    ]
    return lines


def _val(v: Optional[str]) -> str:
    return v if v is not None else ""


def write_position(
    results: list[VariantResult],
    path: str,
    info: list[str],
) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="") as f:
        for line in _header_lines(info, "Variant position file"):
            f.write(line + "\n")
        w = csv.writer(f)
        w.writerow(["marker_name", "alt_marker_name", "chromosome", "position", "VCF_REF", "VCF_ALT"])
        for r in results:
            w.writerow([
                _val(r.name), _val(r.alt_name),
                _val(r.chromosome), _val(str(r.position) if r.position else ""),
                _val(r.vcf_ref), _val(r.vcf_alt),
            ])


def write_conversion(
    results: list[VariantResult],
    path: str,
    info: list[str],
) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="") as f:
        for line in _header_lines(info, "Genotype conversion file"):
            f.write(line + "\n")
        w = csv.writer(f)
        w.writerow(["marker_name", "alt_marker_name", "AB", "TOP", "FORWARD", "DESIGN", "PLUS", "VCF"])
        for r in results:
            w.writerow([
                _val(r.name), _val(r.alt_name),
                "A", _val(r.top_a), _val(r.forward_a), _val(r.design_a), _val(r.plus_a), _val(r.vcf_a),
            ])
            w.writerow([
                _val(r.name), _val(r.alt_name),
                "B", _val(r.top_b), _val(r.forward_b), _val(r.design_b), _val(r.plus_b), _val(r.vcf_b),
            ])


def write_wide(
    results: list[VariantResult],
    path: str,
    info: list[str],
) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="") as f:
        for line in _header_lines(info, "Wide file"):
            f.write(line + "\n")
        w = csv.writer(f)
        w.writerow([
            "marker_name", "alt_marker_name", "chromosome", "position",
            "VCF_REF", "VCF_ALT",
            "AB_A", "AB_B", "TOP_A", "TOP_B",
            "FORWARD_A", "FORWARD_B", "DESIGN_A", "DESIGN_B",
            "PLUS_A", "PLUS_B", "VCF_A", "VCF_B",
        ])
        for r in results:
            w.writerow([
                _val(r.name), _val(r.alt_name),
                _val(r.chromosome), _val(str(r.position) if r.position else ""),
                _val(r.vcf_ref), _val(r.vcf_alt),
                _val(r.ab_a), _val(r.ab_b),
                _val(r.top_a), _val(r.top_b),
                _val(r.forward_a), _val(r.forward_b),
                _val(r.design_a), _val(r.design_b),
                _val(r.plus_a), _val(r.plus_b),
                _val(r.vcf_a), _val(r.vcf_b),
            ])


def write_alignment(
    results: list[VariantResult],
    path: str,
    info: list[str],
) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    sep = "=" * 88
    with open(p, "w") as f:
        for line in _header_lines(info, "Alignment file"):
            f.write(line + "\n")
        for r in results:
            f.write(sep + "\n")
            f.write(f"{_val(r.name)},{_val(r.alt_name)}\n")
            if r.alignment_text:
                f.write(r.alignment_text.rstrip() + "\n")
            else:
                f.write("No alignment obtained.\n")


def write_lookup(
    results: list[VariantResult],
    path: str,
    info: list[str],
) -> None:
    """
    Write a single self-describing lookup file (Option B combined format).

    One row per variant with clearly named columns. Column naming:
      A_in_TOP = what allele A looks like in TOP format, etc.
    This file can be used directly with the 'convert' subcommand.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="") as f:
        for line in _header_lines(info, "Genotype lookup file"):
            f.write(line + "\n")
        f.write(
            "# Columns: marker_name, alt_marker_name, chromosome, position (1-based),\n"
            "#   ref_allele (VCF REF), alt_allele (VCF ALT),\n"
            "#   A_in_<FORMAT> = nucleotide for allele A in that format encoding,\n"
            "#   B_in_<FORMAT> = nucleotide for allele B in that format encoding,\n"
            "#   A_vcf / B_vcf = whether allele A/B is REF or ALT in VCF notation.\n"
            "# Supported formats: AB, TOP, FORWARD, DESIGN, PLUS.\n"
            "#\n"
        )
        w = csv.writer(f)
        w.writerow([
            "marker_name", "alt_marker_name", "chromosome", "position",
            "ref_allele", "alt_allele",
            "A_in_AB", "B_in_AB",
            "A_in_TOP", "B_in_TOP",
            "A_in_FORWARD", "B_in_FORWARD",
            "A_in_DESIGN", "B_in_DESIGN",
            "A_in_PLUS", "B_in_PLUS",
            "A_vcf", "B_vcf",
        ])
        for r in results:
            w.writerow([
                _val(r.name), _val(r.alt_name),
                _val(r.chromosome),
                _val(str(r.position) if r.position else ""),
                _val(r.vcf_ref), _val(r.vcf_alt),
                _val(r.ab_a), _val(r.ab_b),
                _val(r.top_a), _val(r.top_b),
                _val(r.forward_a), _val(r.forward_b),
                _val(r.design_a), _val(r.design_b),
                _val(r.plus_a), _val(r.plus_b),
                _val(r.vcf_a), _val(r.vcf_b),
            ])


def write_lookup_parquet(
    results: list[VariantResult],
    path: str,
    info: list[str],
) -> None:
    """Write lookup data as a Parquet file (requires pyarrow)."""
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError:
        raise ImportError(
            "pyarrow is required for Parquet output: pip install pyarrow"
        )

    rows = []
    for r in results:
        rows.append({
            "marker_name": _val(r.name),
            "alt_marker_name": _val(r.alt_name),
            "chromosome": _val(r.chromosome),
            "position": r.position,
            "ref_allele": _val(r.vcf_ref),
            "alt_allele": _val(r.vcf_alt),
            "A_in_AB": _val(r.ab_a),
            "B_in_AB": _val(r.ab_b),
            "A_in_TOP": _val(r.top_a),
            "B_in_TOP": _val(r.top_b),
            "A_in_FORWARD": _val(r.forward_a),
            "B_in_FORWARD": _val(r.forward_b),
            "A_in_DESIGN": _val(r.design_a),
            "B_in_DESIGN": _val(r.design_b),
            "A_in_PLUS": _val(r.plus_a),
            "B_in_PLUS": _val(r.plus_b),
            "A_vcf": _val(r.vcf_a),
            "B_vcf": _val(r.vcf_b),
        })

    # Build metadata dict from info lines
    meta = {line.split("=", 1)[0].lstrip("#"): line.split("=", 1)[1]
            for line in info if "=" in line}

    table = pa.Table.from_pylist(rows)
    table = table.replace_schema_metadata(
        {k.encode(): v.encode() for k, v in meta.items()}
    )
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, path, compression="snappy")


def write_summary(stats: BuildStats, path: str) -> None:
    """Write a human-readable build summary with marker counts and file checksums."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%A %B %d %H:%M:%S %Y")

    pct = lambda n: f"{100 * n / stats.total_markers:.1f}%" if stats.total_markers else "n/a"

    lines: list[str] = [
        "genotype_converter build summary",
        "=" * 50,
        f"Generated:  {ts}",
        "",
        "Input files",
        "-----------",
        f"Manifest:   {stats.manifest_path}",
        f"  SHA256:   {stats.manifest_sha256}",
        f"Reference:  {stats.reference_path}",
        f"  SHA256:   {stats.reference_sha256}",
        "",
        "Run parameters",
        "--------------",
        f"Species:    {stats.species}",
        f"Workers:    {stats.workers}",
        "",
        "Marker counts",
        "-------------",
        f"Total:      {stats.total_markers}",
        f"  SNPs:     {stats.n_snp}  ({pct(stats.n_snp)})",
        f"  Indels:   {stats.n_indel}  ({pct(stats.n_indel)})",
    ]
    if stats.by_manifest_type:
        lines.append("  By manifest type:")
        for mtype, count in sorted(stats.by_manifest_type.items(), key=lambda x: x[0].name):
            lines.append(f"    {mtype.name.title()}: {count}  ({pct(count)})")
    lines += [
        "",
        "Alignment results",
        "-----------------",
        f"Positioned:     {stats.n_positioned}  ({pct(stats.n_positioned)})",
        f"Not positioned: {stats.n_not_positioned}  ({pct(stats.n_not_positioned)})",
    ]
    if stats.by_chromosome:
        lines += ["", "Per-chromosome marker counts", "---------------------------"]
        for chrom, count in sorted(
            stats.by_chromosome.items(),
            key=lambda x: (x[0].lstrip("chr").zfill(20)),
        ):
            lines.append(f"  {chrom:<20} {count}")
    lines += ["", "Output files", "------------"]
    for f in stats.output_files:
        lines.append(f"  {f}")

    with open(p, "w") as fh:
        fh.write("\n".join(lines) + "\n")
