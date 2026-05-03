from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


MISSING = {"0", "00", "na", "n/a", "--", ".", "", "0/0", "00/00"}
FORMATS = {"AB", "TOP", "FORWARD", "DESIGN", "PLUS", "VCF"}


@dataclass(frozen=True)
class ConvertStats:
    input_path: str
    output_path: str
    rows_total: int
    markers_total: int
    genotype_cells_total: int
    genotypes_parsed: int
    genotypes_changed: int
    missing_or_unparsed_genotypes: int
    alleles_changed: int
    unknown_alleles: int


@dataclass(frozen=True)
class BatchConvertStats:
    files_total: int
    files_converted: int
    stats: list[ConvertStats]
    summary_path: str

    @property
    def output_paths(self) -> list[str]:
        return [item.output_path for item in self.stats]


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
            for key in ("chromosome", "position"):
                row_a[key] = row.get(key, "")
                row_b[key] = row.get(key, "")
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


def _convert_allele_with_status(
    allele: str, marker: str, from_fmt: str, to_fmt: str, table: dict
) -> tuple[str, bool]:
    rows = table.get(marker)
    if not rows:
        return allele, False
    for row in rows:
        if row.get(from_fmt, "").upper() == allele.upper():
            return row.get(to_fmt, allele), True
    return allele, False


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
) -> ConvertStats:
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
    genotype_cells_total = len(rows) * len(marker_cols)
    genotypes_parsed = 0
    genotypes_changed = 0
    missing_or_unparsed_genotypes = 0
    alleles_changed = 0
    unknown_alleles = 0

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
                    missing_or_unparsed_genotypes += 1
                    continue  # missing – keep as-is
                genotypes_parsed += 1
                a_out, a_known = _convert_allele_with_status(pair[0], col, from_fmt, to_fmt, table)
                b_out, b_known = _convert_allele_with_status(pair[1], col, from_fmt, to_fmt, table)
                unknown_alleles += int(not a_known) + int(not b_known)
                alleles_changed += int(a_out != pair[0]) + int(b_out != pair[1])
                converted = _join_alleles(a_out, b_out, out_sep)
                if converted != gt:
                    genotypes_changed += 1
                new_row[col] = converted
            writer.writerow(new_row)
    return ConvertStats(
        input_path=str(input_path),
        output_path=str(output_path),
        rows_total=len(rows),
        markers_total=len(marker_cols),
        genotype_cells_total=genotype_cells_total,
        genotypes_parsed=genotypes_parsed,
        genotypes_changed=genotypes_changed,
        missing_or_unparsed_genotypes=missing_or_unparsed_genotypes,
        alleles_changed=alleles_changed,
        unknown_alleles=unknown_alleles,
    )


def convert_file(
    input_path: str,
    output_path: str,
    table: dict,
    from_fmt: str,
    to_fmt: str,
    layout: str,
    in_sep: Optional[str],
    out_sep: str,
    sample_col: str,
    marker_col: str,
    genotype_col: str,
) -> ConvertStats:
    """Convert one CSV genotype file in wide or long layout."""
    layout = layout.lower()
    if layout == "wide":
        return convert_wide(
            input_path=input_path,
            output_path=output_path,
            table=table,
            from_fmt=from_fmt,
            to_fmt=to_fmt,
            in_sep=in_sep,
            out_sep=out_sep,
            sample_col=sample_col,
        )
    if layout == "long":
        return convert_long(
            input_path=input_path,
            output_path=output_path,
            table=table,
            from_fmt=from_fmt,
            to_fmt=to_fmt,
            in_sep=in_sep,
            out_sep=out_sep,
            sample_col=sample_col,
            marker_col=marker_col,
            genotype_col=genotype_col,
        )
    raise ValueError(f"Unsupported genotype layout: {layout}")


def _write_csv_batch_summary(summary_path: Path, stats: list[ConvertStats]) -> None:
    fieldnames = [
        "input_path",
        "output_path",
        "rows_total",
        "markers_total",
        "genotype_cells_total",
        "genotypes_parsed",
        "genotypes_changed",
        "missing_or_unparsed_genotypes",
        "alleles_changed",
        "unknown_alleles",
    ]
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in stats:
            writer.writerow({field: getattr(item, field) for field in fieldnames})


def convert_batch(
    input_dir: str,
    output_dir: str,
    pattern: str,
    suffix: str,
    overwrite: bool,
    table: dict,
    from_fmt: str,
    to_fmt: str,
    layout: str,
    in_sep: Optional[str],
    out_sep: str,
    sample_col: str,
    marker_col: str,
    genotype_col: str,
) -> BatchConvertStats:
    """Convert every matching CSV genotype file in a directory."""
    source_dir = Path(input_dir)
    if not source_dir.is_dir():
        raise ValueError(f"Genotype input directory does not exist: {input_dir}")
    paths = sorted(path for path in source_dir.glob(pattern) if path.is_file())
    if not paths:
        raise ValueError(f"No genotype files matched pattern {pattern!r} in {input_dir}")

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "conversion_summary.csv"
    if summary_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output file already exists: {summary_path}. Use overwrite=True to replace it."
        )
    stats: list[ConvertStats] = []
    for path in paths:
        out_path = out_dir / f"{path.stem}{suffix}"
        if out_path.exists() and not overwrite:
            raise FileExistsError(
                f"Output file already exists: {out_path}. Use overwrite=True to replace it."
            )
        stats.append(
            convert_file(
                input_path=str(path),
                output_path=str(out_path),
                table=table,
                from_fmt=from_fmt,
                to_fmt=to_fmt,
                layout=layout,
                in_sep=in_sep,
                out_sep=out_sep,
                sample_col=sample_col,
                marker_col=marker_col,
                genotype_col=genotype_col,
            )
        )
    _write_csv_batch_summary(summary_path, stats)
    return BatchConvertStats(
        files_total=len(paths),
        files_converted=len(stats),
        stats=stats,
        summary_path=str(summary_path),
    )


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
) -> ConvertStats:
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
    genotypes_parsed = 0
    genotypes_changed = 0
    missing_or_unparsed_genotypes = 0
    alleles_changed = 0
    unknown_alleles = 0

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
                genotypes_parsed += 1
                a_out, a_known = _convert_allele_with_status(pair[0], marker, from_fmt, to_fmt, table)
                b_out, b_known = _convert_allele_with_status(pair[1], marker, from_fmt, to_fmt, table)
                unknown_alleles += int(not a_known) + int(not b_known)
                alleles_changed += int(a_out != pair[0]) + int(b_out != pair[1])
                converted = _join_alleles(a_out, b_out, out_sep)
                if converted != gt:
                    genotypes_changed += 1
                new_row[genotype_col] = converted
            else:
                missing_or_unparsed_genotypes += 1
            writer.writerow(new_row)
    return ConvertStats(
        input_path=str(input_path),
        output_path=str(output_path),
        rows_total=len(rows),
        markers_total=len({row.get(marker_col, "") for row in rows if row.get(marker_col, "")}),
        genotype_cells_total=len(rows),
        genotypes_parsed=genotypes_parsed,
        genotypes_changed=genotypes_changed,
        missing_or_unparsed_genotypes=missing_or_unparsed_genotypes,
        alleles_changed=alleles_changed,
        unknown_alleles=unknown_alleles,
    )
