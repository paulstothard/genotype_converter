#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path


WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = WORKSPACE_ROOT / "sources"
OUT_ROOT = WORKSPACE_ROOT / "converted" / "example_conversions"
REPORT_ROOT = WORKSPACE_ROOT / "reports" / "example_conversions"
SUMMARY_CSV = WORKSPACE_ROOT / "reports" / "example_conversion_check.csv"
SUMMARY_MD = WORKSPACE_ROOT / "reports" / "example_conversion_check.md"


@dataclass(frozen=True)
class CheckResult:
    species: str
    assembly: str
    format_name: str
    status: str
    input_records: int
    output_records: int
    resolution_rows: int
    skipped_or_ambiguous: int
    message: str


def data_lines(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines() if line.strip()]


def count_csv_records(path: Path) -> int:
    with path.open(newline="") as handle:
        return sum(1 for _row in csv.DictReader(handle))


def count_csv_wide_markers(path: Path) -> int:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return 0
        return max(len(reader.fieldnames) - 1, 0)


def count_csv_long_markers(path: Path) -> int:
    with path.open(newline="") as handle:
        return len(
            {
                row.get("marker_name", "")
                for row in csv.DictReader(handle)
                if row.get("marker_name")
            }
        )


def count_gsgt_matrix_rows(path: Path) -> int:
    lines = data_lines(path)
    data_index = lines.index("[Data]")
    return max(len(lines[data_index + 2 :]), 0)


def count_gsgt_long_rows(path: Path) -> int:
    lines = data_lines(path)
    data_index = lines.index("[Data]")
    return max(len(lines[data_index + 2 :]), 0)


def count_gsgt_long_markers(path: Path) -> int:
    lines = data_lines(path)
    data_index = lines.index("[Data]")
    reader = csv.DictReader(lines[data_index + 1 :], delimiter="\t")
    return len(
        {
            row.get("SNP Name", "")
            for row in reader
            if row.get("SNP Name")
        }
    )


def count_tabular_body_rows(path: Path) -> int:
    return max(len(data_lines(path)) - 1, 0)


def count_plink_bim(path: Path) -> int:
    return len(data_lines(path))


def count_pvar(path: Path) -> int:
    return sum(
        1
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    )


def resolution_counts(path: Path) -> tuple[int, int]:
    if not path.exists():
        return 0, 0
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    flagged = sum(
        1
        for row in rows
        if row.get("status") in {"ambiguous", "missing", "skipped"}
        or row.get("selected_manifest_name", "") == ""
    )
    return len(rows), flagged


def check_one(
    species: str,
    assembly: str,
    format_name: str,
    input_path: Path,
    output_path: Path,
    report_path: Path,
    record_counter,
    marker_counter,
) -> CheckResult:
    if not input_path.exists():
        return CheckResult(species, assembly, format_name, "missing-input", 0, 0, 0, 0, str(input_path))
    if not output_path.exists():
        return CheckResult(species, assembly, format_name, "missing-output", record_counter(input_path), 0, 0, 0, str(output_path))
    input_records = record_counter(input_path)
    output_records = record_counter(output_path)
    expected_resolution_rows = marker_counter(input_path)
    resolution_rows, flagged = resolution_counts(report_path)
    status = (
        "pass"
        if input_records == output_records
        and resolution_rows >= expected_resolution_rows
        and flagged == 0
        else "warn"
    )
    message = ""
    if input_records != output_records:
        message = "input/output record counts differ"
    elif resolution_rows < expected_resolution_rows:
        message = "resolution report has fewer rows than distinct markers"
    elif flagged:
        message = "resolution report contains missing or ambiguous markers"
    return CheckResult(
        species,
        assembly,
        format_name,
        status,
        input_records,
        output_records,
        resolution_rows,
        flagged,
        message,
    )


def species_assemblies() -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for species_dir in sorted(path for path in SOURCE_ROOT.iterdir() if path.is_dir()):
        references_dir = species_dir / "references"
        if not references_dir.is_dir():
            continue
        for assembly_dir in sorted(path for path in references_dir.iterdir() if path.is_dir()):
            pairs.append((species_dir.name, assembly_dir.name))
    return pairs


def collect_results() -> list[CheckResult]:
    results: list[CheckResult] = []
    for species, assembly in species_assemblies():
        genotype_dir = SOURCE_ROOT / species / "genotypes" / "synthetic_mixed_manifest"
        out_dir = OUT_ROOT / species / assembly
        report_dir = REPORT_ROOT / species / assembly
        checks = [
            (
                "csv-wide",
                genotype_dir / "mixed_manifest_wide_ab.csv",
                out_dir / "mixed_manifest_wide_plus.csv",
                report_dir / "mixed_manifest_wide_resolution.csv",
                count_csv_records,
                count_csv_wide_markers,
            ),
            (
                "csv-long",
                genotype_dir / "mixed_manifest_long_ab.csv",
                out_dir / "mixed_manifest_long_plus.csv",
                report_dir / "mixed_manifest_long_resolution.csv",
                count_csv_records,
                count_csv_long_markers,
            ),
            (
                "illumina-matrix",
                genotype_dir / "illumina_gsgt_matrix_ab.txt",
                out_dir / "illumina_gsgt_matrix_plus.txt",
                report_dir / "illumina_gsgt_matrix_resolution.csv",
                count_gsgt_matrix_rows,
                count_gsgt_matrix_rows,
            ),
            (
                "illumina-long",
                genotype_dir / "illumina_gsgt_long_multiformat.txt",
                out_dir / "illumina_gsgt_long_plus.txt",
                report_dir / "illumina_gsgt_long_resolution.csv",
                count_gsgt_long_rows,
                count_gsgt_long_markers,
            ),
            (
                "affymetrix-matrix",
                genotype_dir / "affymetrix_axiom_dual_call_matrix.txt",
                out_dir / "affymetrix_axiom_dual_call_matrix_plus.txt",
                report_dir / "affymetrix_axiom_dual_call_matrix_resolution.csv",
                count_tabular_body_rows,
                count_tabular_body_rows,
            ),
            (
                "plink1",
                genotype_dir / "mixed_manifest_plink1_ab.bim",
                out_dir / "mixed_manifest_plink1_plus.bim",
                report_dir / "mixed_manifest_plink1_resolution.csv",
                count_plink_bim,
                count_plink_bim,
            ),
            (
                "plink2",
                genotype_dir / "mixed_manifest_plink2_ab.pvar",
                out_dir / "mixed_manifest_plink2_plus.pvar",
                report_dir / "mixed_manifest_plink2_resolution.csv",
                count_pvar,
                count_pvar,
            ),
        ]
        for item in checks:
            results.append(check_one(species, assembly, *item))
    return results


def write_reports(results: list[CheckResult]) -> None:
    SUMMARY_CSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "species",
        "assembly",
        "format_name",
        "status",
        "input_records",
        "output_records",
        "resolution_rows",
        "skipped_or_ambiguous",
        "message",
    ]
    with SUMMARY_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow({field: getattr(result, field) for field in fieldnames})

    failures = [result for result in results if result.status != "pass"]
    lines = [
        "# Example Conversion Check",
        "",
        f"Checks: {len(results)}",
        f"Passed: {len(results) - len(failures)}",
        f"Warnings or failures: {len(failures)}",
        "",
        "| Species | Assembly | Format | Status | Records | Resolution rows | Skipped/Ambiguous |",
        "| --- | --- | --- | --- | ---: | ---: | ---: |",
    ]
    for result in results:
        lines.append(
            "| "
            + " | ".join(
                [
                    result.species,
                    result.assembly,
                    result.format_name,
                    result.status,
                    str(result.output_records),
                    str(result.resolution_rows),
                    str(result.skipped_or_ambiguous),
                ]
            )
            + " |"
        )
    SUMMARY_MD.write_text("\n".join(lines) + "\n")


def main() -> int:
    results = collect_results()
    write_reports(results)
    failures = [result for result in results if result.status != "pass"]
    print(f"Wrote {SUMMARY_CSV}")
    print(f"Wrote {SUMMARY_MD}")
    if failures:
        print(f"{len(failures)} check(s) require review")
        return 1
    print(f"All {len(results)} example conversion checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
