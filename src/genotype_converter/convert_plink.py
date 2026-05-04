from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .convert_genotypes import FORMATS, _convert_allele_with_status


PLINK_FORMATS = FORMATS - {"VCF"}
UNCONVERTIBLE_POLICIES = {"exclude", "fail", "keep"}


@dataclass
class PlinkConvertStats:
    variants_total: int
    variants_converted: int
    variants_missing_lookup: int
    variants_incomplete_mapping: int
    variants_excluded: int
    alleles_changed: int
    genotype_path: str
    variant_path: str
    sample_path: str
    marker_report_path: str
    exclude_marker_path: str

    @property
    def bed_path(self) -> str:
        return self.genotype_path

    @property
    def bim_path(self) -> str:
        return self.variant_path

    @property
    def fam_path(self) -> str:
        return self.sample_path


@dataclass
class PlinkBatchConvertStats:
    filesets_total: int
    filesets_converted: int
    stats: list[PlinkConvertStats]
    summary_path: str

    @property
    def output_prefixes(self) -> list[str]:
        return [
            str(Path(item.genotype_path).with_suffix(""))
            for item in self.stats
        ]


def _prefix_path(prefix: str, suffix: str) -> Path:
    return Path(f"{prefix}.{suffix}")


def _require_plink_formats(from_fmt: str, to_fmt: str) -> None:
    unknown = [fmt for fmt in (from_fmt, to_fmt) if fmt not in PLINK_FORMATS]
    if unknown:
        supported = ", ".join(sorted(PLINK_FORMATS))
        raise ValueError(
            "PLINK conversion supports allele-label formats only. "
            f"Unsupported format(s): {', '.join(unknown)}. Supported: {supported}"
        )


def _require_bfile(prefix: str) -> tuple[Path, Path, Path]:
    bed = _prefix_path(prefix, "bed")
    bim = _prefix_path(prefix, "bim")
    fam = _prefix_path(prefix, "fam")
    missing = [str(path) for path in (bed, bim, fam) if not path.exists()]
    if missing:
        raise ValueError(
            "PLINK binary conversion requires all three files: .bed, .bim, and .fam. "
            f"Missing: {', '.join(missing)}"
        )
    return bed, bim, fam


def _require_pfile(prefix: str) -> tuple[Path, Path, Path]:
    pgen = _prefix_path(prefix, "pgen")
    pvar = _prefix_path(prefix, "pvar")
    psam = _prefix_path(prefix, "psam")
    missing = [str(path) for path in (pgen, pvar, psam) if not path.exists()]
    if missing:
        raise ValueError(
            "PLINK 2 conversion requires all three files: .pgen, .pvar, and .psam. "
            f"Missing: {', '.join(missing)}"
        )
    return pgen, pvar, psam


def _discover_prefixes(input_dir: str, pattern: str, suffix: str, label: str) -> list[Path]:
    source_dir = Path(input_dir)
    if not source_dir.is_dir():
        raise ValueError(f"{label} input directory does not exist: {input_dir}")
    paths = sorted(path for path in source_dir.glob(pattern) if path.is_file())
    prefixes = [path.with_suffix("") for path in paths if path.suffix == f".{suffix}"]
    if not prefixes:
        raise ValueError(
            f"No {label} filesets matched pattern {pattern!r} in {input_dir}. "
            f"The pattern should match .{suffix} files."
        )
    return prefixes


def _require_batch_outputs_available(
    output_prefixes: list[Path],
    suffixes: tuple[str, str, str],
    overwrite: bool,
) -> None:
    if overwrite:
        return
    existing: list[str] = []
    for prefix in output_prefixes:
        existing.extend(str(_prefix_path(str(prefix), suffix)) for suffix in suffixes
                        if _prefix_path(str(prefix), suffix).exists())
    if existing:
        preview = ", ".join(existing[:10])
        more = f" and {len(existing) - 10} more" if len(existing) > 10 else ""
        raise FileExistsError(
            f"Output file(s) already exist: {preview}{more}. "
            "Use overwrite=True to replace them."
        )


def _write_plink_batch_summary(
    summary_path: Path,
    input_prefixes: list[Path],
    stats: list[PlinkConvertStats],
) -> None:
    fieldnames = [
        "input_prefix",
        "output_prefix",
        "variants_total",
        "variants_converted",
        "variants_missing_lookup",
        "variants_incomplete_mapping",
        "variants_excluded",
        "alleles_changed",
        "genotype_path",
        "variant_path",
        "sample_path",
        "marker_report_path",
        "exclude_marker_path",
    ]
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for input_prefix, item in zip(input_prefixes, stats):
            writer.writerow(
                {
                    "input_prefix": str(input_prefix),
                    "output_prefix": str(Path(item.genotype_path).with_suffix("")),
                    "variants_total": item.variants_total,
                    "variants_converted": item.variants_converted,
                    "variants_missing_lookup": item.variants_missing_lookup,
                    "variants_incomplete_mapping": item.variants_incomplete_mapping,
                    "variants_excluded": item.variants_excluded,
                    "alleles_changed": item.alleles_changed,
                    "genotype_path": item.genotype_path,
                    "variant_path": item.variant_path,
                    "sample_path": item.sample_path,
                    "marker_report_path": item.marker_report_path,
                    "exclude_marker_path": item.exclude_marker_path,
                }
            )


def _preview_markers(markers: list[str]) -> str:
    unique_markers = sorted(set(markers))
    preview = ", ".join(unique_markers[:10])
    more = f" and {len(unique_markers) - 10} more" if len(unique_markers) > 10 else ""
    return f"{preview}{more}"


def _convert_plink_allele(
    allele: str,
    marker: str,
    from_fmt: str,
    to_fmt: str,
    table: dict,
) -> tuple[str, bool]:
    if allele == "0":
        return allele, True
    return _convert_allele_with_status(allele, marker, from_fmt, to_fmt, table)


def _require_unconvertible_policy(policy: str) -> None:
    if policy not in UNCONVERTIBLE_POLICIES:
        supported = ", ".join(sorted(UNCONVERTIBLE_POLICIES))
        raise ValueError(
            f"Unsupported unconvertible marker policy: {policy}. Supported: {supported}"
        )


def _marker_report_path(output_prefix: str) -> Path:
    return _prefix_path(output_prefix, "marker_conversion_report.csv")


def _exclude_marker_path(output_prefix: str) -> Path:
    return _prefix_path(output_prefix, "exclude_markers.txt")


def _marker_issue(
    marker: str,
    alleles: tuple[str, str],
    table: dict,
    from_fmt: str,
    to_fmt: str,
) -> str:
    if marker not in table:
        return "missing_lookup"
    rows = table.get(marker) or []
    if not rows:
        return "empty_lookup_rule"

    reasons: list[str] = []
    for allele in alleles:
        if allele == "0":
            continue
        matched = None
        for row in rows:
            if row.get(from_fmt, "").upper() == allele.upper():
                matched = row
                break
        if matched is None:
            reasons.append(f"input_allele_not_in_{from_fmt}")
        elif not matched.get(to_fmt):
            reasons.append(f"missing_{to_fmt}_allele")
    return ";".join(sorted(set(reasons)))


def _write_marker_report(report_path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "marker_name",
        "action",
        "reason",
        "from_format",
        "to_format",
        "allele1",
        "allele2",
        "chromosome",
        "position",
        "message",
    ]
    with report_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_exclude_markers(path: Path, rows: list[dict[str, str]]) -> None:
    path.write_text(
        "".join(f"{row['marker_name']}\n" for row in rows if row["action"] == "exclude")
    )


def _scan_bim_marker_issues(
    input_bim: Path,
    table: dict,
    from_fmt: str,
    to_fmt: str,
    action: str,
) -> tuple[int, list[dict[str, str]]]:
    total = 0
    rows: list[dict[str, str]] = []
    with input_bim.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            fields = stripped.split()
            if len(fields) != 6:
                raise ValueError(
                    f"{input_bim} line {line_number} has {len(fields)} fields; expected 6"
                )
            chrom, marker, _cm, pos, allele1, allele2 = fields
            total += 1
            reason = _marker_issue(marker, (allele1, allele2), table, from_fmt, to_fmt)
            if reason:
                rows.append(
                    {
                        "marker_name": marker,
                        "action": action,
                        "reason": reason,
                        "from_format": from_fmt,
                        "to_format": to_fmt,
                        "allele1": allele1,
                        "allele2": allele2,
                        "chromosome": chrom,
                        "position": pos,
                        "message": "Variant cannot be converted completely.",
                    }
                )
    return total, rows


def _scan_pvar_marker_issues(
    input_pvar: Path,
    table: dict,
    from_fmt: str,
    to_fmt: str,
    action: str,
) -> tuple[int, list[dict[str, str]]]:
    total = 0
    rows: list[dict[str, str]] = []
    indexes: dict[str, int] | None = None
    with input_pvar.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.rstrip("\n")
            if not stripped or stripped.startswith("##"):
                continue
            if stripped.startswith("#CHROM"):
                indexes = _pvar_header_indexes(stripped.split())
                continue
            if indexes is None:
                raise ValueError(
                    f"{input_pvar} line {line_number} appears before a #CHROM header. "
                    "Only headered .pvar files are supported."
                )
            fields = stripped.split()
            max_index = max(indexes.values())
            if len(fields) <= max_index:
                raise ValueError(
                    f"{input_pvar} line {line_number} has {len(fields)} fields; "
                    f"expected at least {max_index + 1}"
                )
            marker = fields[indexes["ID"]]
            alt = fields[indexes["ALT"]]
            if "," in alt:
                raise ValueError(
                    f"{input_pvar} line {line_number} marker {marker!r} is multiallelic; "
                    "PLINK 2 conversion currently supports biallelic records only"
                )
            ref = fields[indexes["REF"]]
            total += 1
            reason = _marker_issue(marker, (ref, alt), table, from_fmt, to_fmt)
            if reason:
                rows.append(
                    {
                        "marker_name": marker,
                        "action": action,
                        "reason": reason,
                        "from_format": from_fmt,
                        "to_format": to_fmt,
                        "allele1": ref,
                        "allele2": alt,
                        "chromosome": fields[indexes["CHROM"]],
                        "position": fields[indexes["POS"]],
                        "message": "Variant cannot be converted completely.",
                    }
                )
    if indexes is None:
        raise ValueError(f"{input_pvar} is missing a #CHROM header")
    return total, rows


def _issue_counts(rows: list[dict[str, str]]) -> tuple[int, int]:
    missing = sum(1 for row in rows if "missing_lookup" in row["reason"].split(";"))
    incomplete = len(rows) - missing
    return missing, incomplete


def _run_plink_exclude(
    *,
    command: str,
    mode: str,
    input_prefix: str,
    exclude_path: Path,
    output_prefix: Path,
) -> None:
    if mode == "bfile":
        args = [
            command,
            "--bfile",
            input_prefix,
            "--exclude",
            str(exclude_path),
            "--make-bed",
            "--out",
            str(output_prefix),
        ]
    else:
        args = [
            command,
            "--pfile",
            input_prefix,
            "--exclude",
            str(exclude_path),
            "--make-pgen",
            "--out",
            str(output_prefix),
        ]
    try:
        result = subprocess.run(
            args,
            check=False,
            text=True,
            capture_output=True,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"{command!r} was not found. Install PLINK/PLINK2 or rerun with "
            "--on-unconvertible-marker fail or --on-unconvertible-marker keep."
        ) from exc
    if result.returncode != 0:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(
            f"{command} failed while excluding unconvertible markers"
            + (f": {detail}" if detail else ".")
        )


def convert_plink_bfile(
    input_prefix: str,
    output_prefix: str,
    table: dict,
    from_fmt: str,
    to_fmt: str,
    *,
    update_position: bool = False,
    on_unconvertible_marker: str = "exclude",
    plink_command: str = "plink",
) -> PlinkConvertStats:
    """
    Convert allele labels in a PLINK binary fileset.

    The .bim allele columns are rewritten according to the lookup table. When
    on_unconvertible_marker="exclude", PLINK removes unconvertible variants
    before the final .bed/.bim/.fam files are written.
    """
    _require_plink_formats(from_fmt, to_fmt)
    _require_unconvertible_policy(on_unconvertible_marker)
    input_bed, input_bim, input_fam = _require_bfile(input_prefix)

    output_bed = _prefix_path(output_prefix, "bed")
    output_bim = _prefix_path(output_prefix, "bim")
    output_fam = _prefix_path(output_prefix, "fam")
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    output_bim.parent.mkdir(parents=True, exist_ok=True)
    output_fam.parent.mkdir(parents=True, exist_ok=True)
    report_path = _marker_report_path(output_prefix)
    exclude_path = _exclude_marker_path(output_prefix)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    exclude_path.parent.mkdir(parents=True, exist_ok=True)

    action = "exclude" if on_unconvertible_marker == "exclude" else on_unconvertible_marker
    variants_total, issue_rows = _scan_bim_marker_issues(
        input_bim, table, from_fmt, to_fmt, action
    )
    variants_missing_lookup, variants_incomplete_mapping = _issue_counts(issue_rows)
    _write_marker_report(report_path, issue_rows)
    if on_unconvertible_marker == "exclude" and issue_rows:
        _write_exclude_markers(exclude_path, issue_rows)
    else:
        exclude_path.unlink(missing_ok=True)

    if issue_rows and on_unconvertible_marker == "fail":
        raise ValueError(
            "PLINK .bim contains marker(s) that cannot be fully converted from "
            f"{from_fmt} to {to_fmt}: "
            f"{_preview_markers([row['marker_name'] for row in issue_rows])}. "
            f"See {report_path}."
        )

    conversion_prefix = input_prefix
    conversion_bed = input_bed
    conversion_bim = input_bim
    conversion_fam = input_fam
    variants_excluded = 0

    temp_dir: tempfile.TemporaryDirectory[str] | None = None
    try:
        if issue_rows and on_unconvertible_marker == "exclude":
            variants_excluded = len(issue_rows)
            temp_dir = tempfile.TemporaryDirectory(prefix="plink_exclude_", dir=output_bed.parent)
            filtered_prefix = Path(temp_dir.name) / "filtered"
            _run_plink_exclude(
                command=plink_command,
                mode="bfile",
                input_prefix=input_prefix,
                exclude_path=exclude_path,
                output_prefix=filtered_prefix,
            )
            conversion_prefix = str(filtered_prefix)
            conversion_bed, conversion_bim, conversion_fam = _require_bfile(conversion_prefix)

        variants_converted = 0
        alleles_changed = 0

        with conversion_bim.open() as in_handle, output_bim.open("w") as out_handle:
            for line_number, line in enumerate(in_handle, start=1):
                stripped = line.strip()
                if not stripped:
                    continue
                fields = stripped.split()
                if len(fields) != 6:
                    raise ValueError(
                        f"{conversion_bim} line {line_number} has {len(fields)} fields; expected 6"
                    )

                chrom, marker, cm, pos, allele1, allele2 = fields
                if marker not in table:
                    if on_unconvertible_marker != "keep":
                        raise ValueError(
                            f"Marker {marker!r} remained after filtering but is absent from lookup"
                        )
                    out_handle.write("\t".join(fields) + "\n")
                    continue
                rows = table.get(marker) or []
                if not rows:
                    if on_unconvertible_marker != "keep":
                        raise ValueError(
                            f"Marker {marker!r} remained after filtering but has no lookup rule"
                        )
                    out_handle.write("\t".join(fields) + "\n")
                    continue

                allele1_out, allele1_complete = _convert_plink_allele(
                    allele1, marker, from_fmt, to_fmt, table
                )
                allele2_out, allele2_complete = _convert_plink_allele(
                    allele2, marker, from_fmt, to_fmt, table
                )
                if not allele1_complete or not allele2_complete:
                    if on_unconvertible_marker != "keep":
                        raise ValueError(
                            f"Marker {marker!r} remained after filtering but cannot be fully converted"
                        )
                if allele1_out != allele1:
                    alleles_changed += 1
                if allele2_out != allele2:
                    alleles_changed += 1

                lookup_row = rows[0]
                if update_position:
                    chrom = lookup_row.get("chromosome", chrom) or chrom
                    pos = lookup_row.get("position", pos) or pos

                out_handle.write(
                    "\t".join([chrom, marker, cm, pos, allele1_out, allele2_out]) + "\n"
                )
                variants_converted += 1

        shutil.copyfile(conversion_bed, output_bed)
        shutil.copyfile(conversion_fam, output_fam)
    finally:
        if temp_dir is not None:
            temp_dir.cleanup()

    return PlinkConvertStats(
        variants_total=variants_total,
        variants_converted=variants_converted,
        variants_missing_lookup=variants_missing_lookup,
        variants_incomplete_mapping=variants_incomplete_mapping,
        variants_excluded=variants_excluded,
        alleles_changed=alleles_changed,
        genotype_path=str(output_bed),
        variant_path=str(output_bim),
        sample_path=str(output_fam),
        marker_report_path=str(report_path),
        exclude_marker_path=str(exclude_path) if exclude_path.exists() else "",
    )


def convert_plink_bfile_batch(
    input_dir: str,
    output_dir: str,
    pattern: str,
    suffix: str,
    overwrite: bool,
    table: dict,
    from_fmt: str,
    to_fmt: str,
    *,
    update_position: bool = False,
    on_unconvertible_marker: str = "exclude",
    plink_command: str = "plink",
) -> PlinkBatchConvertStats:
    """Convert every matching PLINK 1 binary fileset in a directory."""
    input_prefixes = _discover_prefixes(input_dir, pattern, "bed", "PLINK binary")
    out_dir = Path(output_dir)
    output_prefixes = [out_dir / f"{prefix.name}{suffix}" for prefix in input_prefixes]
    _require_batch_outputs_available(output_prefixes, ("bed", "bim", "fam"), overwrite)
    summary_path = out_dir / "conversion_summary.csv"
    if summary_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output file already exists: {summary_path}. Use overwrite=True to replace it."
        )

    stats: list[PlinkConvertStats] = []
    for input_prefix, output_prefix in zip(input_prefixes, output_prefixes):
        stats.append(
            convert_plink_bfile(
                input_prefix=str(input_prefix),
                output_prefix=str(output_prefix),
                table=table,
                from_fmt=from_fmt,
                to_fmt=to_fmt,
                update_position=update_position,
                on_unconvertible_marker=on_unconvertible_marker,
                plink_command=plink_command,
            )
        )
    _write_plink_batch_summary(summary_path, input_prefixes, stats)
    return PlinkBatchConvertStats(
        filesets_total=len(input_prefixes),
        filesets_converted=len(stats),
        stats=stats,
        summary_path=str(summary_path),
    )


def _pvar_header_indexes(header: list[str]) -> dict[str, int]:
    normalized = ["CHROM" if col == "#CHROM" else col for col in header]
    required = ["CHROM", "POS", "ID", "REF", "ALT"]
    missing = [col for col in required if col not in normalized]
    if missing:
        raise ValueError(f"PLINK 2 .pvar header is missing column(s): {', '.join(missing)}")
    return {col: normalized.index(col) for col in required}


def convert_plink_pfile(
    input_prefix: str,
    output_prefix: str,
    table: dict,
    from_fmt: str,
    to_fmt: str,
    *,
    update_position: bool = False,
    on_unconvertible_marker: str = "exclude",
    plink2_command: str = "plink2",
) -> PlinkConvertStats:
    """
    Convert allele labels in a PLINK 2 pfile.

    The .pvar REF/ALT columns are rewritten according to the lookup table. When
    on_unconvertible_marker="exclude", PLINK2 removes unconvertible variants
    before the final .pgen/.pvar/.psam files are written.
    """
    _require_plink_formats(from_fmt, to_fmt)
    _require_unconvertible_policy(on_unconvertible_marker)
    input_pgen, input_pvar, input_psam = _require_pfile(input_prefix)

    output_pgen = _prefix_path(output_prefix, "pgen")
    output_pvar = _prefix_path(output_prefix, "pvar")
    output_psam = _prefix_path(output_prefix, "psam")
    output_pgen.parent.mkdir(parents=True, exist_ok=True)
    output_pvar.parent.mkdir(parents=True, exist_ok=True)
    output_psam.parent.mkdir(parents=True, exist_ok=True)
    report_path = _marker_report_path(output_prefix)
    exclude_path = _exclude_marker_path(output_prefix)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    exclude_path.parent.mkdir(parents=True, exist_ok=True)

    action = "exclude" if on_unconvertible_marker == "exclude" else on_unconvertible_marker
    variants_total, issue_rows = _scan_pvar_marker_issues(
        input_pvar, table, from_fmt, to_fmt, action
    )
    variants_missing_lookup, variants_incomplete_mapping = _issue_counts(issue_rows)
    _write_marker_report(report_path, issue_rows)
    if on_unconvertible_marker == "exclude" and issue_rows:
        _write_exclude_markers(exclude_path, issue_rows)
    else:
        exclude_path.unlink(missing_ok=True)

    if issue_rows and on_unconvertible_marker == "fail":
        raise ValueError(
            "PLINK 2 .pvar contains marker(s) that cannot be fully converted from "
            f"{from_fmt} to {to_fmt}: "
            f"{_preview_markers([row['marker_name'] for row in issue_rows])}. "
            f"See {report_path}."
        )

    conversion_prefix = input_prefix
    conversion_pgen = input_pgen
    conversion_pvar = input_pvar
    conversion_psam = input_psam
    variants_excluded = 0

    temp_dir: tempfile.TemporaryDirectory[str] | None = None
    try:
        if issue_rows and on_unconvertible_marker == "exclude":
            variants_excluded = len(issue_rows)
            temp_dir = tempfile.TemporaryDirectory(prefix="plink2_exclude_", dir=output_pgen.parent)
            filtered_prefix = Path(temp_dir.name) / "filtered"
            _run_plink_exclude(
                command=plink2_command,
                mode="pfile",
                input_prefix=input_prefix,
                exclude_path=exclude_path,
                output_prefix=filtered_prefix,
            )
            conversion_prefix = str(filtered_prefix)
            conversion_pgen, conversion_pvar, conversion_psam = _require_pfile(conversion_prefix)

        variants_converted = 0
        alleles_changed = 0
        indexes: dict[str, int] | None = None

        with conversion_pvar.open() as in_handle, output_pvar.open("w") as out_handle:
            for line_number, line in enumerate(in_handle, start=1):
                stripped = line.rstrip("\n")
                if not stripped:
                    out_handle.write(line)
                    continue
                if stripped.startswith("##"):
                    out_handle.write(line)
                    continue
                if stripped.startswith("#CHROM"):
                    header = stripped.split()
                    indexes = _pvar_header_indexes(header)
                    out_handle.write("\t".join(header) + "\n")
                    continue
                if indexes is None:
                    raise ValueError(
                        f"{conversion_pvar} line {line_number} appears before a #CHROM header. "
                        "Only headered .pvar files are supported."
                    )

                fields = stripped.split()
                max_index = max(indexes.values())
                if len(fields) <= max_index:
                    raise ValueError(
                        f"{conversion_pvar} line {line_number} has {len(fields)} fields; "
                        f"expected at least {max_index + 1}"
                    )

                marker = fields[indexes["ID"]]
                alt = fields[indexes["ALT"]]
                if "," in alt:
                    raise ValueError(
                        f"{conversion_pvar} line {line_number} marker {marker!r} is multiallelic; "
                        "PLINK 2 conversion currently supports biallelic records only"
                    )

                if marker not in table:
                    if on_unconvertible_marker != "keep":
                        raise ValueError(
                            f"Marker {marker!r} remained after filtering but is absent from lookup"
                        )
                    out_handle.write("\t".join(fields) + "\n")
                    continue
                rows = table.get(marker) or []
                if not rows:
                    if on_unconvertible_marker != "keep":
                        raise ValueError(
                            f"Marker {marker!r} remained after filtering but has no lookup rule"
                        )
                    out_handle.write("\t".join(fields) + "\n")
                    continue

                ref = fields[indexes["REF"]]
                ref_out, ref_complete = _convert_plink_allele(
                    ref, marker, from_fmt, to_fmt, table
                )
                alt_out, alt_complete = _convert_plink_allele(
                    alt, marker, from_fmt, to_fmt, table
                )
                if not ref_complete or not alt_complete:
                    if on_unconvertible_marker != "keep":
                        raise ValueError(
                            f"Marker {marker!r} remained after filtering but cannot be fully converted"
                        )
                if ref_out != ref:
                    alleles_changed += 1
                if alt_out != alt:
                    alleles_changed += 1
                fields[indexes["REF"]] = ref_out
                fields[indexes["ALT"]] = alt_out

                lookup_row = rows[0]
                if update_position:
                    fields[indexes["CHROM"]] = lookup_row.get("chromosome", fields[indexes["CHROM"]]) or fields[indexes["CHROM"]]
                    fields[indexes["POS"]] = lookup_row.get("position", fields[indexes["POS"]]) or fields[indexes["POS"]]

                out_handle.write("\t".join(fields) + "\n")
                variants_converted += 1

        if indexes is None:
            output_pvar.unlink(missing_ok=True)
            raise ValueError(f"{conversion_pvar} is missing a #CHROM header")

        shutil.copyfile(conversion_pgen, output_pgen)
        shutil.copyfile(conversion_psam, output_psam)
    finally:
        if temp_dir is not None:
            temp_dir.cleanup()

    return PlinkConvertStats(
        variants_total=variants_total,
        variants_converted=variants_converted,
        variants_missing_lookup=variants_missing_lookup,
        variants_incomplete_mapping=variants_incomplete_mapping,
        variants_excluded=variants_excluded,
        alleles_changed=alleles_changed,
        genotype_path=str(output_pgen),
        variant_path=str(output_pvar),
        sample_path=str(output_psam),
        marker_report_path=str(report_path),
        exclude_marker_path=str(exclude_path) if exclude_path.exists() else "",
    )


def convert_plink_pfile_batch(
    input_dir: str,
    output_dir: str,
    pattern: str,
    suffix: str,
    overwrite: bool,
    table: dict,
    from_fmt: str,
    to_fmt: str,
    *,
    update_position: bool = False,
    on_unconvertible_marker: str = "exclude",
    plink2_command: str = "plink2",
) -> PlinkBatchConvertStats:
    """Convert every matching PLINK 2 fileset in a directory."""
    input_prefixes = _discover_prefixes(input_dir, pattern, "pgen", "PLINK 2")
    out_dir = Path(output_dir)
    output_prefixes = [out_dir / f"{prefix.name}{suffix}" for prefix in input_prefixes]
    _require_batch_outputs_available(output_prefixes, ("pgen", "pvar", "psam"), overwrite)
    summary_path = out_dir / "conversion_summary.csv"
    if summary_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output file already exists: {summary_path}. Use overwrite=True to replace it."
        )

    stats: list[PlinkConvertStats] = []
    for input_prefix, output_prefix in zip(input_prefixes, output_prefixes):
        stats.append(
            convert_plink_pfile(
                input_prefix=str(input_prefix),
                output_prefix=str(output_prefix),
                table=table,
                from_fmt=from_fmt,
                to_fmt=to_fmt,
                update_position=update_position,
                on_unconvertible_marker=on_unconvertible_marker,
                plink2_command=plink2_command,
            )
        )
    _write_plink_batch_summary(summary_path, input_prefixes, stats)
    return PlinkBatchConvertStats(
        filesets_total=len(input_prefixes),
        filesets_converted=len(stats),
        stats=stats,
        summary_path=str(summary_path),
    )
