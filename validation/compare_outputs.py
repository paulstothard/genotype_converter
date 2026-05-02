from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
from typing import Iterable


DEFAULT_OLD_POSITION = Path("validation/old_pipeline/bovinehd-manifest-b.ARS_UCD_v2_0.position.csv")
DEFAULT_OLD_CONVERSION = Path("validation/old_pipeline/bovinehd-manifest-b.ARS_UCD_v2_0.conversion.csv")
DEFAULT_NEW_DIR = Path("validation/new_pipeline/bos_taurus/ARS_UCD_v2_0")
DEFAULT_NEW_POSITION = DEFAULT_NEW_DIR / "bovinehd-manifest-b.ARS_UCD_v2_0.position.csv"
DEFAULT_NEW_CONVERSION = DEFAULT_NEW_DIR / "bovinehd-manifest-b.ARS_UCD_v2_0.conversion.csv"
DEFAULT_REPORT_DIR = Path("validation/reports")


def _data_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    return list(csv.DictReader(lines))


def _clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def _by_marker(rows: Iterable[dict[str, str]]) -> dict[str, dict[str, str]]:
    result = {}
    for row in rows:
        marker = _clean(row.get("marker_name"))
        if marker:
            result[marker] = row
    return result


def _by_marker_ab(rows: Iterable[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    result = {}
    for row in rows:
        marker = _clean(row.get("marker_name"))
        ab = _clean(row.get("AB"))
        if marker and ab:
            result[(marker, ab)] = row
    return result


def _compare_tables(
    old: dict,
    new: dict,
    columns: list[str],
) -> tuple[set, set, list[dict[str, str]], Counter]:
    old_keys = set(old)
    new_keys = set(new)
    missing_in_new = old_keys - new_keys
    extra_in_new = new_keys - old_keys
    mismatches: list[dict[str, str]] = []
    mismatch_columns: Counter = Counter()

    for key in sorted(old_keys & new_keys):
        old_row = old[key]
        new_row = new[key]
        for column in columns:
            old_value = _clean(old_row.get(column))
            new_value = _clean(new_row.get(column))
            if old_value != new_value:
                mismatch_columns[column] += 1
                if isinstance(key, tuple):
                    marker, ab = key
                else:
                    marker, ab = key, ""
                mismatches.append({
                    "marker_name": marker,
                    "AB": ab,
                    "column": column,
                    "old_value": old_value,
                    "new_value": new_value,
                })

    return missing_in_new, extra_in_new, mismatches, mismatch_columns


def _write_mismatches(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["marker_name", "AB", "column", "old_value", "new_value"],
        )
        writer.writeheader()
        writer.writerows(rows)


def _value_state(old_value: str, new_value: str) -> str:
    if old_value and new_value:
        return "old_set/new_set"
    if old_value:
        return "old_set/new_blank"
    if new_value:
        return "old_blank/new_set"
    return "old_blank/new_blank"


def _position_breakdown(
    old_position: dict[str, dict[str, str]],
    new_position: dict[str, dict[str, str]],
    mismatches: list[dict[str, str]],
) -> list[str]:
    mismatch_markers = {row["marker_name"] for row in mismatches}
    state_counts: Counter = Counter()
    position_deltas: Counter = Counter()
    ref_alt_swaps = 0

    for row in mismatches:
        state_counts[_value_state(row["old_value"], row["new_value"])] += 1
        if row["column"] != "position":
            continue
        try:
            delta = int(row["new_value"]) - int(row["old_value"])
        except ValueError:
            continue
        position_deltas[delta] += 1

    for marker in mismatch_markers:
        old_row = old_position.get(marker, {})
        new_row = new_position.get(marker, {})
        if (
            _clean(old_row.get("VCF_REF"))
            and _clean(old_row.get("VCF_ALT"))
            and _clean(old_row.get("VCF_REF")) == _clean(new_row.get("VCF_ALT"))
            and _clean(old_row.get("VCF_ALT")) == _clean(new_row.get("VCF_REF"))
        ):
            ref_alt_swaps += 1

    return [
        f"- Markers with any position-output mismatch: {len(mismatch_markers)}",
        f"- Value-state counts: {dict(state_counts)}",
        f"- Position deltas, new minus old: {dict(position_deltas.most_common(10))}",
        f"- Markers with old/new VCF REF/ALT swapped: {ref_alt_swaps}",
    ]


def _conversion_breakdown(mismatches: list[dict[str, str]]) -> list[str]:
    mismatch_markers = {row["marker_name"] for row in mismatches}
    state_counts: Counter = Counter(
        _value_state(row["old_value"], row["new_value"]) for row in mismatches
    )
    marker_column_counts: Counter = Counter(
        row["column"] for row in mismatches if row["column"] in {"PLUS", "VCF"}
    )
    return [
        f"- Markers with any conversion-output mismatch: {len(mismatch_markers)}",
        f"- Value-state counts: {dict(state_counts)}",
        f"- PLUS/VCF mismatch counts: {dict(marker_column_counts)}",
    ]


def _sample_keys(keys: set, limit: int = 20) -> str:
    if not keys:
        return "None"
    return ", ".join(str(key) for key in sorted(keys)[:limit])


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare old and new validation outputs.")
    parser.add_argument("--old-position", type=Path, default=DEFAULT_OLD_POSITION)
    parser.add_argument("--old-conversion", type=Path, default=DEFAULT_OLD_CONVERSION)
    parser.add_argument("--new-position", type=Path, default=DEFAULT_NEW_POSITION)
    parser.add_argument("--new-conversion", type=Path, default=DEFAULT_NEW_CONVERSION)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    args = parser.parse_args()

    old_position = _by_marker(_data_rows(args.old_position))
    new_position = _by_marker(_data_rows(args.new_position))
    old_conversion = _by_marker_ab(_data_rows(args.old_conversion))
    new_conversion = _by_marker_ab(_data_rows(args.new_conversion))

    pos_missing, pos_extra, pos_mismatches, pos_columns = _compare_tables(
        old_position,
        new_position,
        ["chromosome", "position", "VCF_REF", "VCF_ALT"],
    )
    conv_missing, conv_extra, conv_mismatches, conv_columns = _compare_tables(
        old_conversion,
        new_conversion,
        ["TOP", "FORWARD", "DESIGN", "PLUS", "VCF"],
    )

    args.report_dir.mkdir(parents=True, exist_ok=True)
    pos_mismatch_path = args.report_dir / "position_mismatches.csv"
    conv_mismatch_path = args.report_dir / "conversion_mismatches.csv"
    _write_mismatches(pos_mismatch_path, pos_mismatches)
    _write_mismatches(conv_mismatch_path, conv_mismatches)
    pos_breakdown = _position_breakdown(old_position, new_position, pos_mismatches)
    conv_breakdown = _conversion_breakdown(conv_mismatches)

    report = args.report_dir / "comparison_summary.md"
    report.write_text(
        "\n".join([
            "# Validation Comparison Summary",
            "",
            "## Inputs",
            "",
            f"- Old position: `{args.old_position}`",
            f"- New position: `{args.new_position}`",
            f"- Old conversion: `{args.old_conversion}`",
            f"- New conversion: `{args.new_conversion}`",
            "",
            "## Position Output",
            "",
            f"- Old markers: {len(old_position)}",
            f"- New markers: {len(new_position)}",
            f"- Missing in new: {len(pos_missing)}",
            f"- Extra in new: {len(pos_extra)}",
            f"- Cell mismatches: {len(pos_mismatches)}",
            f"- Mismatches by column: {dict(pos_columns)}",
            f"- Missing sample: {_sample_keys(pos_missing)}",
            f"- Extra sample: {_sample_keys(pos_extra)}",
            f"- Mismatch details: `{pos_mismatch_path}`",
            "",
            "### Position Breakdown",
            "",
            *pos_breakdown,
            "",
            "## Conversion Output",
            "",
            f"- Old marker/AB rows: {len(old_conversion)}",
            f"- New marker/AB rows: {len(new_conversion)}",
            f"- Missing in new: {len(conv_missing)}",
            f"- Extra in new: {len(conv_extra)}",
            f"- Cell mismatches: {len(conv_mismatches)}",
            f"- Mismatches by column: {dict(conv_columns)}",
            f"- Missing sample: {_sample_keys(conv_missing)}",
            f"- Extra sample: {_sample_keys(conv_extra)}",
            f"- Mismatch details: `{conv_mismatch_path}`",
            "",
            "### Conversion Breakdown",
            "",
            *conv_breakdown,
            "",
        ]) + "\n",
        encoding="utf-8",
    )

    print(report)


if __name__ == "__main__":
    main()
