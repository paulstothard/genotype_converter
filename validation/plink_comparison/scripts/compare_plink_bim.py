#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


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
    expected_path: Path,
    actual_path: Path,
    *,
    original_path: Path | None = None,
) -> list[dict[str, str]]:
    expected = read_bim(expected_path)
    actual = read_bim(actual_path)
    original = read_bim(original_path) if original_path else {}
    rows: list[dict[str, str]] = []
    for marker in expected:
        expected_row = expected[marker]
        actual_row = actual.get(marker)
        original_row = original.get(marker)
        if actual_row is None:
            rows.append({
                "marker_name": marker,
                "status": "missing_in_actual",
                "mike_vs_original": classify(original_row, expected_row) if original_row else "",
                "ours_vs_original": "",
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
            "ours_vs_original": classify(original_row, actual_row) if original_row else "",
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
            "status": "extra_in_actual",
            "mike_vs_original": "",
            "ours_vs_original": classify(original[marker], actual_row) if marker in original else "",
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
    fieldnames = [
        "marker_name", "status", "mike_vs_original", "ours_vs_original",
        "original_chrom", "original_position", "original_allele1",
        "original_allele2", "expected_chrom", "expected_position",
        "expected_allele1", "expected_allele2", "actual_chrom",
        "actual_position", "actual_allele1", "actual_allele2",
        "position_status",
    ]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    counts = Counter(row["status"] for row in rows)
    transition_counts = Counter(
        (row["mike_vs_original"], row["ours_vs_original"], row["status"])
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
    if any(row["mike_vs_original"] or row["ours_vs_original"] for row in rows):
        lines.extend([
            "",
            "## Original-To-Converted Relationships",
            "",
            "| Mike vs original | Ours vs original | Ours vs Mike | Count |",
            "| --- | --- | --- | ---: |",
        ])
        for (mike_status, ours_status, comparison_status), count in sorted(
            transition_counts.items(),
            key=lambda item: (-item[1], item[0]),
        ):
            lines.append(
                f"| {mike_status} | {ours_status} | {comparison_status} | {count} |"
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
            "| Marker | Status | Original | Mike expected | Ours | Position |",
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
    parser = argparse.ArgumentParser(description="Compare expected and actual PLINK BIM files.")
    parser.add_argument("--expected", required=True, type=Path, help="Expected converted .bim")
    parser.add_argument("--actual", required=True, type=Path, help="Actual converted .bim")
    parser.add_argument("--original", required=False, type=Path, help="Original pre-conversion .bim")
    parser.add_argument("--report-prefix", required=True, type=Path, help="Output report prefix")
    args = parser.parse_args()

    rows = compare(args.expected, args.actual, original_path=args.original)
    write_reports(rows, args.report_prefix)
    failing = [
        row for row in rows
        if row["status"] not in {"exact", "swapped"}
        or row["position_status"] == "different"
    ]
    return 1 if failing else 0


if __name__ == "__main__":
    raise SystemExit(main())
