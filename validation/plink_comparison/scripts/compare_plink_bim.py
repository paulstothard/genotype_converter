#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
STATUS_DEFINITIONS = {
    "exact": "The allele pair is identical in both files, in the same allele-column order.",
    "swapped": "The same two allele labels are present, but allele1 and allele2 are reversed.",
    "complement": "Both allele labels in the second file are DNA complements of the first file, in the same allele-column order.",
    "swapped_complement": "Both allele labels are DNA complements and allele1/allele2 are also reversed.",
    "mismatch": "The allele labels are not explained by exact, swapped, complement, or swapped-complement relationships.",
    "missing_from_genotype_converter_plus": "The marker is present in Mike's PLUS file but absent from genotype_converter PLUS.",
    "extra_in_genotype_converter_plus": "The marker is present in genotype_converter PLUS but absent from Mike's PLUS file.",
}
CSV_FIELD_MAP = [
    ("marker_name", "marker_name"),
    ("status", "status"),
    ("mike_vs_original", "mike_plus_vs_mike_top"),
    ("genotype_converter_vs_original", "genotype_converter_plus_vs_mike_top"),
    ("original_chrom", "mike_top_chrom"),
    ("original_position", "mike_top_position"),
    ("original_allele1", "mike_top_allele1"),
    ("original_allele2", "mike_top_allele2"),
    ("expected_chrom", "mike_plus_chrom"),
    ("expected_position", "mike_plus_position"),
    ("expected_allele1", "mike_plus_allele1"),
    ("expected_allele2", "mike_plus_allele2"),
    ("actual_chrom", "genotype_converter_plus_chrom"),
    ("actual_position", "genotype_converter_plus_position"),
    ("actual_allele1", "genotype_converter_plus_allele1"),
    ("actual_allele2", "genotype_converter_plus_allele2"),
    ("position_status", "position_status"),
]


@dataclass(frozen=True)
class BimRow:
    chrom: str
    marker: str
    cm: str
    pos: str
    allele1: str
    allele2: str


def read_bim(path: Path) -> dict[str, BimRow]:
    rows: dict[str, BimRow] = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split()
            if len(fields) != 6:
                raise ValueError(f"{path}:{line_number} does not have 6 BIM fields")
            row = BimRow(*fields)
            if row.marker in rows:
                raise ValueError(f"{path}:{line_number} duplicate marker {row.marker}")
            rows[row.marker] = row
    return rows


def complement(value: str) -> str:
    if len(value) != 1 or value.upper() not in {"A", "C", "G", "T"}:
        return value
    return value.translate(COMPLEMENT).upper()


def classify(source: BimRow, target: BimRow) -> str:
    source_pair = (source.allele1, source.allele2)
    target_pair = (target.allele1, target.allele2)
    if target_pair == source_pair:
        return "exact"
    if target_pair == (source.allele2, source.allele1):
        return "swapped"
    complement_pair = (complement(source.allele1), complement(source.allele2))
    if target_pair == complement_pair:
        return "complement"
    if target_pair == (complement_pair[1], complement_pair[0]):
        return "swapped_complement"
    return "mismatch"


def compare(
    mike_plus_path: Path,
    genotype_converter_plus_path: Path,
    *,
    mike_top_path: Path | None = None,
) -> list[dict[str, str]]:
    expected = read_bim(mike_plus_path)
    actual = read_bim(genotype_converter_plus_path)
    original = read_bim(mike_top_path) if mike_top_path else {}
    rows: list[dict[str, str]] = []
    for marker in expected:
        expected_row = expected[marker]
        actual_row = actual.get(marker)
        original_row = original.get(marker)
        if actual_row is None:
            rows.append({
                "marker_name": marker,
                "status": "missing_from_genotype_converter_plus",
                "mike_vs_original": classify(original_row, expected_row) if original_row else "",
                "genotype_converter_vs_original": "",
                "expected_chrom": expected_row.chrom,
                "expected_position": expected_row.pos,
                "expected_allele1": expected_row.allele1,
                "expected_allele2": expected_row.allele2,
                "original_chrom": original_row.chrom if original_row else "",
                "original_position": original_row.pos if original_row else "",
                "original_allele1": original_row.allele1 if original_row else "",
                "original_allele2": original_row.allele2 if original_row else "",
                "actual_chrom": "",
                "actual_position": "",
                "actual_allele1": "",
                "actual_allele2": "",
                "position_status": "",
            })
            continue
        position_status = (
            "same"
            if (expected_row.chrom, expected_row.pos) == (actual_row.chrom, actual_row.pos)
            else "different"
        )
        rows.append({
            "marker_name": marker,
            "status": classify(expected_row, actual_row),
            "mike_vs_original": classify(original_row, expected_row) if original_row else "",
            "genotype_converter_vs_original": classify(original_row, actual_row) if original_row else "",
            "original_chrom": original_row.chrom if original_row else "",
            "original_position": original_row.pos if original_row else "",
            "original_allele1": original_row.allele1 if original_row else "",
            "original_allele2": original_row.allele2 if original_row else "",
            "expected_chrom": expected_row.chrom,
            "expected_position": expected_row.pos,
            "expected_allele1": expected_row.allele1,
            "expected_allele2": expected_row.allele2,
            "actual_chrom": actual_row.chrom,
            "actual_position": actual_row.pos,
            "actual_allele1": actual_row.allele1,
            "actual_allele2": actual_row.allele2,
            "position_status": position_status,
        })
    for marker in actual.keys() - expected.keys():
        actual_row = actual[marker]
        rows.append({
            "marker_name": marker,
            "status": "extra_in_genotype_converter_plus",
            "mike_vs_original": "",
            "genotype_converter_vs_original": classify(original[marker], actual_row) if marker in original else "",
            "original_chrom": original[marker].chrom if marker in original else "",
            "original_position": original[marker].pos if marker in original else "",
            "original_allele1": original[marker].allele1 if marker in original else "",
            "original_allele2": original[marker].allele2 if marker in original else "",
            "expected_chrom": "",
            "expected_position": "",
            "expected_allele1": "",
            "expected_allele2": "",
            "actual_chrom": actual_row.chrom,
            "actual_position": actual_row.pos,
            "actual_allele1": actual_row.allele1,
            "actual_allele2": actual_row.allele2,
            "position_status": "",
        })
    return rows


def write_reports(rows: list[dict[str, str]], report_prefix: Path) -> None:
    report_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = report_prefix.with_suffix(".csv")
    md_path = report_prefix.with_suffix(".md")
    fieldnames = [output_name for _input_name, output_name in CSV_FIELD_MAP]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {output_name: row.get(input_name, "") for input_name, output_name in CSV_FIELD_MAP}
            )

    counts = Counter(row["status"] for row in rows)
    transition_counts = Counter(
        (row["mike_vs_original"], row["genotype_converter_vs_original"], row["status"])
        for row in rows
    )
    position_differences = sum(1 for row in rows if row["position_status"] == "different")
    lines = [
        "# PLINK BIM Comparison",
        "",
        f"Markers compared: {len(rows)}",
        f"Position differences: {position_differences}",
        "",
        "| Status | Count |",
        "| --- | ---: |",
    ]
    for status, count in sorted(counts.items()):
        lines.append(f"| {status} | {count} |")
    lines.extend([
        "",
        "## Status Definitions",
        "",
        "| Status | Meaning |",
        "| --- | --- |",
    ])
    for status in sorted(counts):
        lines.append(f"| {status} | {STATUS_DEFINITIONS.get(status, '')} |")
    if any(row["mike_vs_original"] or row["genotype_converter_vs_original"] for row in rows):
        lines.extend([
            "",
            "## Conversion Relationships",
            "",
            "These rows compare Mike's PLUS file and genotype_converter PLUS to Mike's TOP file, then compare genotype_converter PLUS to Mike's PLUS file.",
            "",
            "| Mike PLUS vs Mike TOP | genotype_converter PLUS vs Mike TOP | genotype_converter PLUS vs Mike PLUS | Count |",
            "| --- | --- | --- | ---: |",
        ])
        for (mike_status, genotype_converter_status, comparison_status), count in sorted(
            transition_counts.items(),
            key=lambda item: (-item[1], item[0]),
        ):
            lines.append(
                f"| {mike_status} | {genotype_converter_status} | {comparison_status} | {count} |"
            )
    mismatch_examples = [
        row for row in rows
        if row["status"] not in {"exact", "swapped"}
        or row["position_status"] == "different"
    ][:20]
    if mismatch_examples:
        lines.extend([
            "",
            "## Review Examples",
            "",
            "| Marker | Status | Mike TOP | Mike PLUS | genotype_converter PLUS | Position |",
            "| --- | --- | --- | --- | --- | --- |",
        ])
        for row in mismatch_examples:
            original = (
                f"{row['original_chrom']}:{row['original_position']} "
                f"{row['original_allele1']}/{row['original_allele2']}"
            ).strip()
            expected = (
                f"{row['expected_chrom']}:{row['expected_position']} "
                f"{row['expected_allele1']}/{row['expected_allele2']}"
            ).strip()
            actual = (
                f"{row['actual_chrom']}:{row['actual_position']} "
                f"{row['actual_allele1']}/{row['actual_allele2']}"
            ).strip()
            lines.append(
                f"| {row['marker_name']} | {row['status']} | "
                f"{original} | {expected} | {actual} | {row['position_status']} |"
            )
    md_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compare Mike PLUS and genotype_converter PLUS PLINK BIM files and write reports. "
            "Differences are reported, not treated as command failure by default."
        )
    )
    parser.add_argument("--mike-plus", dest="mike_plus", type=Path, help="Mike PLUS .bim")
    parser.add_argument(
        "--genotype-converter-plus",
        dest="genotype_converter_plus",
        type=Path,
        help="genotype_converter PLUS .bim",
    )
    parser.add_argument("--mike-top", dest="mike_top", type=Path, help="Mike TOP .bim")
    parser.add_argument("--expected", dest="mike_plus", type=Path, help=argparse.SUPPRESS)
    parser.add_argument(
        "--actual",
        dest="genotype_converter_plus",
        type=Path,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--original", dest="mike_top", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--report-prefix", required=True, type=Path, help="Output report prefix")
    parser.add_argument(
        "--fail-on-differences",
        action="store_true",
        help="Return exit code 1 when non-exact differences are found.",
    )
    args = parser.parse_args()
    if args.mike_plus is None or args.genotype_converter_plus is None:
        parser.error("--mike-plus and --genotype-converter-plus are required")

    rows = compare(
        args.mike_plus,
        args.genotype_converter_plus,
        mike_top_path=args.mike_top,
    )
    write_reports(rows, args.report_prefix)
    failing = [
        row for row in rows
        if row["status"] not in {"exact", "swapped"}
        or row["position_status"] == "different"
    ]
    if args.fail_on_differences and failing:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
