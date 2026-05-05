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
    "missing_from_genotype_converter_plus": "The marker is present in the collaborator PLUS file but absent from genotype_converter PLUS.",
    "extra_in_genotype_converter_plus": "The marker is present in genotype_converter PLUS but absent from the collaborator PLUS file.",
}
CSV_FIELD_MAP = [
    ("marker_name", "marker_name"),
    ("status", "status"),
    ("collaborator_plus_vs_top", "collaborator_plus_vs_top"),
    ("genotype_converter_plus_vs_top", "genotype_converter_plus_vs_collaborator_top"),
    ("collaborator_top_chrom", "collaborator_top_chrom"),
    ("collaborator_top_position", "collaborator_top_position"),
    ("collaborator_top_allele1", "collaborator_top_allele1"),
    ("collaborator_top_allele2", "collaborator_top_allele2"),
    ("collaborator_plus_chrom", "collaborator_plus_chrom"),
    ("collaborator_plus_position", "collaborator_plus_position"),
    ("collaborator_plus_allele1", "collaborator_plus_allele1"),
    ("collaborator_plus_allele2", "collaborator_plus_allele2"),
    ("genotype_converter_plus_chrom", "genotype_converter_plus_chrom"),
    ("genotype_converter_plus_position", "genotype_converter_plus_position"),
    ("genotype_converter_plus_allele1", "genotype_converter_plus_allele1"),
    ("genotype_converter_plus_allele2", "genotype_converter_plus_allele2"),
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
    collaborator_plus_path: Path,
    genotype_converter_plus_path: Path,
    *,
    collaborator_top_path: Path | None = None,
) -> list[dict[str, str]]:
    collaborator_plus = read_bim(collaborator_plus_path)
    genotype_converter_plus = read_bim(genotype_converter_plus_path)
    collaborator_top = read_bim(collaborator_top_path) if collaborator_top_path else {}
    rows: list[dict[str, str]] = []
    for marker in collaborator_plus:
        collaborator_plus_row = collaborator_plus[marker]
        genotype_converter_plus_row = genotype_converter_plus.get(marker)
        collaborator_top_row = collaborator_top.get(marker)
        if genotype_converter_plus_row is None:
            rows.append({
                "marker_name": marker,
                "status": "missing_from_genotype_converter_plus",
                "collaborator_plus_vs_top": classify(collaborator_top_row, collaborator_plus_row) if collaborator_top_row else "",
                "genotype_converter_plus_vs_top": "",
                "collaborator_plus_chrom": collaborator_plus_row.chrom,
                "collaborator_plus_position": collaborator_plus_row.pos,
                "collaborator_plus_allele1": collaborator_plus_row.allele1,
                "collaborator_plus_allele2": collaborator_plus_row.allele2,
                "collaborator_top_chrom": collaborator_top_row.chrom if collaborator_top_row else "",
                "collaborator_top_position": collaborator_top_row.pos if collaborator_top_row else "",
                "collaborator_top_allele1": collaborator_top_row.allele1 if collaborator_top_row else "",
                "collaborator_top_allele2": collaborator_top_row.allele2 if collaborator_top_row else "",
                "genotype_converter_plus_chrom": "",
                "genotype_converter_plus_position": "",
                "genotype_converter_plus_allele1": "",
                "genotype_converter_plus_allele2": "",
                "position_status": "",
            })
            continue
        position_status = (
            "same"
            if (collaborator_plus_row.chrom, collaborator_plus_row.pos) == (genotype_converter_plus_row.chrom, genotype_converter_plus_row.pos)
            else "different"
        )
        rows.append({
            "marker_name": marker,
            "status": classify(collaborator_plus_row, genotype_converter_plus_row),
            "collaborator_plus_vs_top": classify(collaborator_top_row, collaborator_plus_row) if collaborator_top_row else "",
            "genotype_converter_plus_vs_top": classify(collaborator_top_row, genotype_converter_plus_row) if collaborator_top_row else "",
            "collaborator_top_chrom": collaborator_top_row.chrom if collaborator_top_row else "",
            "collaborator_top_position": collaborator_top_row.pos if collaborator_top_row else "",
            "collaborator_top_allele1": collaborator_top_row.allele1 if collaborator_top_row else "",
            "collaborator_top_allele2": collaborator_top_row.allele2 if collaborator_top_row else "",
            "collaborator_plus_chrom": collaborator_plus_row.chrom,
            "collaborator_plus_position": collaborator_plus_row.pos,
            "collaborator_plus_allele1": collaborator_plus_row.allele1,
            "collaborator_plus_allele2": collaborator_plus_row.allele2,
            "genotype_converter_plus_chrom": genotype_converter_plus_row.chrom,
            "genotype_converter_plus_position": genotype_converter_plus_row.pos,
            "genotype_converter_plus_allele1": genotype_converter_plus_row.allele1,
            "genotype_converter_plus_allele2": genotype_converter_plus_row.allele2,
            "position_status": position_status,
        })
    for marker in genotype_converter_plus.keys() - collaborator_plus.keys():
        genotype_converter_plus_row = genotype_converter_plus[marker]
        rows.append({
            "marker_name": marker,
            "status": "extra_in_genotype_converter_plus",
            "collaborator_plus_vs_top": "",
            "genotype_converter_plus_vs_top": classify(collaborator_top[marker], genotype_converter_plus_row) if marker in collaborator_top else "",
            "collaborator_top_chrom": collaborator_top[marker].chrom if marker in collaborator_top else "",
            "collaborator_top_position": collaborator_top[marker].pos if marker in collaborator_top else "",
            "collaborator_top_allele1": collaborator_top[marker].allele1 if marker in collaborator_top else "",
            "collaborator_top_allele2": collaborator_top[marker].allele2 if marker in collaborator_top else "",
            "collaborator_plus_chrom": "",
            "collaborator_plus_position": "",
            "collaborator_plus_allele1": "",
            "collaborator_plus_allele2": "",
            "genotype_converter_plus_chrom": genotype_converter_plus_row.chrom,
            "genotype_converter_plus_position": genotype_converter_plus_row.pos,
            "genotype_converter_plus_allele1": genotype_converter_plus_row.allele1,
            "genotype_converter_plus_allele2": genotype_converter_plus_row.allele2,
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
        (row["collaborator_plus_vs_top"], row["genotype_converter_plus_vs_top"], row["status"])
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
    if any(row["collaborator_plus_vs_top"] or row["genotype_converter_plus_vs_top"] for row in rows):
        lines.extend([
            "",
            "## Conversion Relationships",
            "",
            "These rows compare collaborator PLUS and genotype_converter PLUS to collaborator TOP, then compare genotype_converter PLUS to collaborator PLUS.",
            "",
            "| collaborator PLUS vs TOP | genotype_converter PLUS vs TOP | genotype_converter PLUS vs collaborator PLUS | Count |",
            "| --- | --- | --- | ---: |",
        ])
        for (collaborator_status, genotype_converter_status, comparison_status), count in sorted(
            transition_counts.items(),
            key=lambda item: (-item[1], item[0]),
        ):
            lines.append(
                f"| {collaborator_status} | {genotype_converter_status} | {comparison_status} | {count} |"
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
            "| Marker | Status | collaborator TOP | collaborator PLUS | genotype_converter PLUS | Position |",
            "| --- | --- | --- | --- | --- | --- |",
        ])
        for row in mismatch_examples:
            collaborator_top = (
                f"{row['collaborator_top_chrom']}:{row['collaborator_top_position']} "
                f"{row['collaborator_top_allele1']}/{row['collaborator_top_allele2']}"
            ).strip()
            collaborator_plus = (
                f"{row['collaborator_plus_chrom']}:{row['collaborator_plus_position']} "
                f"{row['collaborator_plus_allele1']}/{row['collaborator_plus_allele2']}"
            ).strip()
            genotype_converter_plus = (
                f"{row['genotype_converter_plus_chrom']}:{row['genotype_converter_plus_position']} "
                f"{row['genotype_converter_plus_allele1']}/{row['genotype_converter_plus_allele2']}"
            ).strip()
            lines.append(
                f"| {row['marker_name']} | {row['status']} | "
                f"{collaborator_top} | {collaborator_plus} | {genotype_converter_plus} | {row['position_status']} |"
            )
    md_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compare collaborator PLUS and genotype_converter PLUS PLINK BIM files and write reports. "
            "Differences are reported, not treated as command failure by default."
        )
    )
    parser.add_argument("--collaborator-plus", dest="collaborator_plus", type=Path, help="Collaborator PLUS .bim")
    parser.add_argument(
        "--genotype-converter-plus",
        dest="genotype_converter_plus",
        type=Path,
        help="genotype_converter PLUS .bim",
    )
    parser.add_argument("--collaborator-top", dest="collaborator_top", type=Path, help="Collaborator TOP .bim")
    parser.add_argument("--report-prefix", required=True, type=Path, help="Output report prefix")
    parser.add_argument(
        "--fail-on-differences",
        action="store_true",
        help="Return exit code 1 when non-exact differences are found.",
    )
    args = parser.parse_args()
    if args.collaborator_plus is None or args.genotype_converter_plus is None:
        parser.error("--collaborator-plus and --genotype-converter-plus are required")

    rows = compare(
        args.collaborator_plus,
        args.genotype_converter_plus,
        collaborator_top_path=args.collaborator_top,
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
