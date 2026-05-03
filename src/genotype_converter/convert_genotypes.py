from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


MISSING = {"0", "00", "na", "n/a", "-", "--", ".", "", "0/0", "00/00"}
FORMATS = {"AB", "TOP", "FORWARD", "DESIGN", "PLUS", "VCF"}
LAYOUTS = {"wide", "long", "illumina-matrix", "illumina-long", "affymetrix-matrix"}
ILLUMINA_LONG_COLUMNS = {
    "TOP": ("Allele1 - Top", "Allele2 - Top"),
    "FORWARD": ("Allele1 - Forward", "Allele2 - Forward"),
    "AB": ("Allele1 - AB", "Allele2 - AB"),
    "DESIGN": ("Allele1 - Design", "Allele2 - Design"),
    "PLUS": ("Allele1 - Plus", "Allele2 - Plus"),
}


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


def _convert_pair(
    pair: tuple[str, str],
    marker: str,
    from_fmt: str,
    to_fmt: str,
    table: dict,
) -> tuple[tuple[str, str], int, int]:
    a_out, a_known = _convert_allele_with_status(pair[0], marker, from_fmt, to_fmt, table)
    b_out, b_known = _convert_allele_with_status(pair[1], marker, from_fmt, to_fmt, table)
    alleles_changed = int(a_out != pair[0]) + int(b_out != pair[1])
    unknown_alleles = int(not a_known) + int(not b_known)
    return (a_out, b_out), alleles_changed, unknown_alleles


def _read_gsgt_sections(input_path: str) -> tuple[list[str], list[str]]:
    lines = Path(input_path).read_text().splitlines()
    try:
        data_index = next(index for index, line in enumerate(lines) if line.strip() == "[Data]")
    except StopIteration as exc:
        raise ValueError(f"{input_path} is missing a [Data] section") from exc
    return lines[: data_index + 1], [line for line in lines[data_index + 1 :] if line.strip()]


def _write_lines(path: str, lines: list[str]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")


def _marker_count(markers: list[str]) -> int:
    return len({marker for marker in markers if marker})


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
                (a_out, b_out), pair_changed, pair_unknown = _convert_pair(
                    pair, col, from_fmt, to_fmt, table
                )
                unknown_alleles += pair_unknown
                alleles_changed += pair_changed
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


def convert_illumina_matrix(
    input_path: str,
    output_path: str,
    table: dict,
    from_fmt: str,
    to_fmt: str,
    in_sep: Optional[str],
    out_sep: str,
) -> ConvertStats:
    """Convert an Illumina GenomeStudio/GSGT matrix file."""
    _require_formats(from_fmt, to_fmt)
    header_lines, data_lines = _read_gsgt_sections(input_path)
    if not data_lines:
        raise ValueError(f"{input_path} has no data rows")
    sample_header = data_lines[0].split("\t")
    if len(sample_header) < 2:
        raise ValueError(f"{input_path} has no sample columns")
    marker_rows = [line.split("\t") for line in data_lines[1:]]
    markers = [row[0] for row in marker_rows if row and row[0]]
    _require_markers(markers, table, input_path)

    genotypes_parsed = 0
    genotypes_changed = 0
    missing_or_unparsed_genotypes = 0
    alleles_changed = 0
    unknown_alleles = 0
    out_lines = list(header_lines)
    out_lines.append(data_lines[0])

    for row in marker_rows:
        if not row:
            continue
        marker = row[0]
        new_row = [marker]
        for gt in row[1:]:
            pair = _split_genotype(gt, in_sep)
            if pair is None:
                missing_or_unparsed_genotypes += 1
                new_row.append(gt)
                continue
            genotypes_parsed += 1
            converted_pair, pair_changed, pair_unknown = _convert_pair(
                pair, marker, from_fmt, to_fmt, table
            )
            alleles_changed += pair_changed
            unknown_alleles += pair_unknown
            converted = _join_alleles(converted_pair[0], converted_pair[1], out_sep)
            if converted != gt:
                genotypes_changed += 1
            new_row.append(converted)
        out_lines.append("\t".join(new_row))

    _write_lines(output_path, out_lines)
    sample_count = max(len(sample_header) - 1, 0)
    return ConvertStats(
        input_path=str(input_path),
        output_path=str(output_path),
        rows_total=len(marker_rows),
        markers_total=_marker_count(markers),
        genotype_cells_total=len(marker_rows) * sample_count,
        genotypes_parsed=genotypes_parsed,
        genotypes_changed=genotypes_changed,
        missing_or_unparsed_genotypes=missing_or_unparsed_genotypes,
        alleles_changed=alleles_changed,
        unknown_alleles=unknown_alleles,
    )


def convert_illumina_long(
    input_path: str,
    output_path: str,
    table: dict,
    from_fmt: str,
    to_fmt: str,
) -> ConvertStats:
    """Convert an Illumina GenomeStudio/GSGT long report."""
    _require_formats(from_fmt, to_fmt)
    if from_fmt == "VCF" or to_fmt == "VCF":
        raise ValueError("Illumina long layout does not have VCF allele columns")
    from_cols = ILLUMINA_LONG_COLUMNS[from_fmt]
    to_cols = ILLUMINA_LONG_COLUMNS[to_fmt]
    header_lines, data_lines = _read_gsgt_sections(input_path)
    reader = csv.DictReader(data_lines, delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError(f"{input_path} has no data header row")
    fieldnames = list(reader.fieldnames)
    _require_columns(fieldnames, ["SNP Name", "Sample ID", *from_cols, *to_cols], input_path)
    rows = list(reader)
    _require_markers([row.get("SNP Name", "") for row in rows], table, input_path)

    genotypes_parsed = 0
    genotypes_changed = 0
    missing_or_unparsed_genotypes = 0
    alleles_changed = 0
    unknown_alleles = 0
    out_lines = list(header_lines)
    out_lines.append("\t".join(fieldnames))

    for row in rows:
        marker = row.get("SNP Name", "")
        pair = (row.get(from_cols[0], ""), row.get(from_cols[1], ""))
        if pair[0].lower() in MISSING or pair[1].lower() in MISSING:
            missing_or_unparsed_genotypes += 1
        else:
            genotypes_parsed += 1
            converted_pair, pair_changed, pair_unknown = _convert_pair(
                pair, marker, from_fmt, to_fmt, table
            )
            alleles_changed += pair_changed
            unknown_alleles += pair_unknown
            if (row.get(to_cols[0], ""), row.get(to_cols[1], "")) != converted_pair:
                genotypes_changed += 1
            row[to_cols[0]], row[to_cols[1]] = converted_pair
        out_lines.append("\t".join(row.get(field, "") for field in fieldnames))

    _write_lines(output_path, out_lines)
    return ConvertStats(
        input_path=str(input_path),
        output_path=str(output_path),
        rows_total=len(rows),
        markers_total=_marker_count([row.get("SNP Name", "") for row in rows]),
        genotype_cells_total=len(rows),
        genotypes_parsed=genotypes_parsed,
        genotypes_changed=genotypes_changed,
        missing_or_unparsed_genotypes=missing_or_unparsed_genotypes,
        alleles_changed=alleles_changed,
        unknown_alleles=unknown_alleles,
    )


def convert_affymetrix_matrix(
    input_path: str,
    output_path: str,
    table: dict,
    from_fmt: str,
    to_fmt: str,
) -> ConvertStats:
    """Convert an Affymetrix/Axiom matrix with paired AB and native-call columns."""
    _require_formats(from_fmt, to_fmt)
    if from_fmt == "VCF" or to_fmt == "VCF":
        raise ValueError("Affymetrix/Axiom matrix layout does not have VCF allele columns")
    lines = Path(input_path).read_text().splitlines()
    if not lines:
        raise ValueError(f"Input file has no header row: {input_path}")
    header = lines[0].split("\t")
    if len(header) < 3 or (len(header) - 1) % 2 != 0:
        raise ValueError(
            f"{input_path} should have probeset_id followed by paired sample columns"
        )
    if header[0] != "probeset_id":
        raise ValueError(f"{input_path} first column should be probeset_id")
    rows = [line.split("\t") for line in lines[1:] if line.strip()]
    markers = [row[0] for row in rows if row]
    _require_markers(markers, table, input_path)

    source_offset = 1 if from_fmt == "AB" else 2
    target_offset = 1 if to_fmt == "AB" else 2
    genotypes_parsed = 0
    genotypes_changed = 0
    missing_or_unparsed_genotypes = 0
    alleles_changed = 0
    unknown_alleles = 0
    out_lines = ["\t".join(header)]

    for row in rows:
        if len(row) != len(header):
            raise ValueError(
                f"{input_path} row for marker {row[0] if row else '<blank>'!r} "
                f"has {len(row)} columns; expected {len(header)}"
            )
        marker = row[0]
        new_row = list(row)
        for index in range(1, len(row), 2):
            gt = row[index if source_offset == 1 else index + 1]
            pair = _split_genotype(gt, None)
            if pair is None:
                missing_or_unparsed_genotypes += 1
                continue
            genotypes_parsed += 1
            converted_pair, pair_changed, pair_unknown = _convert_pair(
                pair, marker, from_fmt, to_fmt, table
            )
            alleles_changed += pair_changed
            unknown_alleles += pair_unknown
            converted = "".join(converted_pair)
            target_index = index if target_offset == 1 else index + 1
            if converted != row[target_index]:
                genotypes_changed += 1
            new_row[target_index] = converted
        out_lines.append("\t".join(new_row))

    _write_lines(output_path, out_lines)
    sample_count = (len(header) - 1) // 2
    return ConvertStats(
        input_path=str(input_path),
        output_path=str(output_path),
        rows_total=len(rows),
        markers_total=_marker_count(markers),
        genotype_cells_total=len(rows) * sample_count,
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
    """Convert one genotype file in a supported text layout."""
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
    if layout == "illumina-matrix":
        return convert_illumina_matrix(
            input_path=input_path,
            output_path=output_path,
            table=table,
            from_fmt=from_fmt,
            to_fmt=to_fmt,
            in_sep=in_sep,
            out_sep=out_sep,
        )
    if layout == "illumina-long":
        return convert_illumina_long(
            input_path=input_path,
            output_path=output_path,
            table=table,
            from_fmt=from_fmt,
            to_fmt=to_fmt,
        )
    if layout == "affymetrix-matrix":
        return convert_affymetrix_matrix(
            input_path=input_path,
            output_path=output_path,
            table=table,
            from_fmt=from_fmt,
            to_fmt=to_fmt,
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
