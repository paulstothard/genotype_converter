from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Optional


MISSING = {"0", "00", "na", "n/a", "--", ".", "", "0/0", "00/00"}
FORMATS = {"AB", "TOP", "FORWARD", "DESIGN", "PLUS", "VCF"}


def _require_formats(from_fmt: str, to_fmt: str) -> None:
    unknown = [fmt for fmt in (from_fmt, to_fmt) if fmt not in FORMATS]
    if unknown:
        supported = ", ".join(sorted(FORMATS))
        raise ValueError(f"Unsupported genotype format(s): {', '.join(unknown)}. Supported: {supported}")


def _require_columns(fieldnames: list[str], required: list[str], input_path: str) -> None:
    missing = [col for col in required if col not in fieldnames]
    if missing:
        raise ValueError(
            f"{input_path} is missing required column(s): {', '.join(missing)}"
        )


def _require_markers(markers: list[str], table: dict, context: str) -> None:
    missing = sorted({marker for marker in markers if marker and marker not in table})
    if missing:
        preview = ", ".join(missing[:10])
        more = f" and {len(missing) - 10} more" if len(missing) > 10 else ""
        raise ValueError(
            f"{context} contains marker(s) not found in the lookup table: {preview}{more}"
        )


def load_lookup_table(lookup_path: str) -> dict:
    """
    Load a lookup.csv or conversion.csv and return:
      { marker_name: [ a_row_dict, b_row_dict ] }

    lookup.csv (one row per marker, columns A_in_TOP / B_in_TOP etc.) and
    conversion.csv (two rows per marker, AB column = 'A' or 'B') are both
    supported.
    """
    with open(lookup_path, newline="") as f:
        lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines)
    fieldnames = list(reader.fieldnames or [])
    rows = list(reader)
    if not fieldnames:
        raise ValueError(f"Lookup file has no header row: {lookup_path}")
    if "marker_name" not in fieldnames:
        raise ValueError(f"Lookup file is missing required column: marker_name")

    # Detect format by checking for A_in_<FORMAT> column pattern
    is_lookup_fmt = any(
        f"A_in_{fmt}" in fieldnames for fmt in ("TOP", "FORWARD", "PLUS", "DESIGN")
    )

    table: dict[str, list[dict]] = {}

    if is_lookup_fmt:
        # One row per marker: split into A-dict and B-dict
        for row in rows:
            name = row.get("marker_name", "")
            if not name:
                continue
            row_a: dict = {"AB": "A"}
            row_b: dict = {"AB": "B"}
            for fmt in ("TOP", "FORWARD", "DESIGN", "PLUS"):
                row_a[fmt] = row.get(f"A_in_{fmt}", "")
                row_b[fmt] = row.get(f"B_in_{fmt}", "")
            # AB-format values: A_in_AB / B_in_AB, fall back to "A"/"B"
            row_a["AB"] = row.get("A_in_AB", "A")
            row_b["AB"] = row.get("B_in_AB", "B")
            # VCF column uses A_vcf / B_vcf names
            row_a["VCF"] = row.get("A_vcf", row.get("A_in_VCF", ""))
            row_b["VCF"] = row.get("B_vcf", row.get("B_in_VCF", ""))
            table[name] = [row_a, row_b]
    else:
        # Two rows per marker (conversion.csv style)
        for row in rows:
            name = row.get("marker_name", "")
            if not name:
                continue
            entry: dict = {}
            for fmt in ("AB", "TOP", "FORWARD", "DESIGN", "PLUS", "VCF"):
                if fmt in row:
                    entry[fmt] = row[fmt]
            table.setdefault(name, []).append(entry)
        table = {k: v for k, v in table.items() if len(v) == 2}

    if not table:
        raise ValueError(f"No usable marker conversion rows found in lookup file: {lookup_path}")

    return table


def _split_genotype(gt: str, sep: Optional[str]) -> Optional[tuple[str, str]]:
    """Split a genotype string into two alleles. Returns None for missing data."""
    gt = gt.strip()
    if gt.lower() in MISSING:
        return None
    if sep:
        parts = gt.split(sep, 1)
    else:
        # Auto-detect
        for s in ("/", " ", "\t"):
            if s in gt:
                parts = gt.split(s, 1)
                break
        else:
            # Single character alleles concatenated: "AG" → ["A", "G"]
            if len(gt) == 2:
                parts = [gt[0], gt[1]]
            else:
                return None
    if len(parts) != 2:
        return None
    a, b = parts[0].strip(), parts[1].strip()
    if a.lower() in MISSING or b.lower() in MISSING:
        return None
    return a, b


def _convert_allele(
    allele: str, marker: str, from_fmt: str, to_fmt: str, table: dict
) -> str:
    """Return converted allele, or original if not found in table."""
    rows = table.get(marker)
    if not rows:
        return allele
    for row in rows:
        if row.get(from_fmt, "").upper() == allele.upper():
            return row.get(to_fmt, allele)
    return allele  # unknown allele, pass through


def _join_alleles(a: str, b: str, sep: str) -> str:
    return f"{a}{sep}{b}"


def convert_wide(
    input_path: str,
    output_path: str,
    table: dict,
    from_fmt: str,
    to_fmt: str,
    in_sep: Optional[str],
    out_sep: str,
    sample_col: str,
) -> None:
    """Convert a wide-format genotype file (samples × markers)."""
    _require_formats(from_fmt, to_fmt)
    with open(input_path, newline="") as f:
        lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines)
    if reader.fieldnames is None:
        raise ValueError("Input file has no header row")

    fieldnames = list(reader.fieldnames)
    _require_columns(fieldnames, [sample_col], input_path)
    marker_cols = [col for col in fieldnames if col != sample_col]
    _require_markers(marker_cols, table, input_path)
    rows = list(reader)

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            new_row = dict(row)
            for col in fieldnames:
                if col == sample_col:
                    continue
                gt = row.get(col, "")
                pair = _split_genotype(gt, in_sep)
                if pair is None:
                    continue  # missing – keep as-is
                a_out = _convert_allele(pair[0], col, from_fmt, to_fmt, table)
                b_out = _convert_allele(pair[1], col, from_fmt, to_fmt, table)
                new_row[col] = _join_alleles(a_out, b_out, out_sep)
            writer.writerow(new_row)


def convert_long(
    input_path: str,
    output_path: str,
    table: dict,
    from_fmt: str,
    to_fmt: str,
    in_sep: Optional[str],
    out_sep: str,
    sample_col: str,
    marker_col: str,
    genotype_col: str,
) -> None:
    """Convert a long-format genotype file (one row per sample×marker)."""
    _require_formats(from_fmt, to_fmt)
    with open(input_path, newline="") as f:
        lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines)
    if reader.fieldnames is None:
        raise ValueError("Input file has no header row")

    fieldnames = list(reader.fieldnames)
    _require_columns(fieldnames, [sample_col, marker_col, genotype_col], input_path)
    rows = list(reader)
    _require_markers([row.get(marker_col, "") for row in rows], table, input_path)

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            new_row = dict(row)
            marker = row.get(marker_col, "")
            gt = row.get(genotype_col, "")
            pair = _split_genotype(gt, in_sep)
            if pair is not None:
                a_out = _convert_allele(pair[0], marker, from_fmt, to_fmt, table)
                b_out = _convert_allele(pair[1], marker, from_fmt, to_fmt, table)
                new_row[genotype_col] = _join_alleles(a_out, b_out, out_sep)
            writer.writerow(new_row)
