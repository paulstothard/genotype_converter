#!/usr/bin/env python3
from __future__ import annotations

import csv
import shutil
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = WORKSPACE_ROOT / "sources"
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from genotype_converter.manifest import ManifestRecord, parse_manifest  # noqa: E402


WINDOW_SIZE = 16
MAX_MARKERS = 28


@dataclass(frozen=True)
class Occurrence:
    species: str
    manifest_name: str
    order_index: int
    record: ManifestRecord


def manifest_paths(species_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in (species_dir / "manifests").iterdir()
        if path.is_file() and path.suffix.lower() in {".csv", ".txt"}
    )


def load_occurrences(species: str, species_dir: Path) -> list[Occurrence]:
    occurrences: list[Occurrence] = []
    for manifest_path in manifest_paths(species_dir):
        records = parse_manifest(str(manifest_path))
        for order_index, record in enumerate(records):
            occurrences.append(
                Occurrence(
                    species=species,
                    manifest_name=manifest_path.stem,
                    order_index=order_index,
                    record=record,
                )
            )
    return occurrences


def marker_index(occurrences: list[Occurrence]) -> dict[str, list[Occurrence]]:
    by_marker: dict[str, list[Occurrence]] = defaultdict(list)
    for occurrence in occurrences:
        by_marker[occurrence.record.name].append(occurrence)
    return dict(by_marker)


def best_marker_window(
    occurrences: list[Occurrence],
    by_marker: dict[str, list[Occurrence]],
) -> list[Occurrence]:
    by_manifest: dict[str, list[Occurrence]] = defaultdict(list)
    for occurrence in occurrences:
        by_manifest[occurrence.manifest_name].append(occurrence)

    best: tuple[int, int, str, int, list[Occurrence]] | None = None
    for manifest_name, manifest_occurrences in by_manifest.items():
        manifest_occurrences.sort(key=lambda occurrence: occurrence.order_index)
        if len(manifest_occurrences) < WINDOW_SIZE:
            windows = [manifest_occurrences]
        else:
            windows = [
                manifest_occurrences[start : start + WINDOW_SIZE]
                for start in range(0, len(manifest_occurrences) - WINDOW_SIZE + 1)
            ]
        for window in windows:
            shared = sum(1 for item in window if len(by_marker[item.record.name]) > 1)
            unique = sum(1 for item in window if len(by_marker[item.record.name]) == 1)
            score = (shared * 10) + unique
            candidate = (score, shared, manifest_name, window[0].order_index, window)
            if best is None or candidate[:4] > best[:4]:
                best = candidate
    if best is None:
        return []
    return best[4]


def choose_markers(occurrences: list[Occurrence]) -> list[Occurrence]:
    by_marker = marker_index(occurrences)
    chosen: dict[str, Occurrence] = {}

    for occurrence in best_marker_window(occurrences, by_marker):
        chosen.setdefault(occurrence.record.name, occurrence)

    shared = sorted(
        (
            occurrences_for_marker
            for occurrences_for_marker in by_marker.values()
            if len({item.manifest_name for item in occurrences_for_marker}) > 1
        ),
        key=lambda items: (-len({item.manifest_name for item in items}), items[0].record.name),
    )
    for occurrences_for_marker in shared:
        chosen.setdefault(occurrences_for_marker[0].record.name, occurrences_for_marker[0])
        if sum(1 for marker in chosen if len(by_marker[marker]) > 1) >= 10:
            break

    unique = sorted(
        (
            occurrences_for_marker[0]
            for occurrences_for_marker in by_marker.values()
            if len({item.manifest_name for item in occurrences_for_marker}) == 1
        ),
        key=lambda item: (item.manifest_name, item.order_index, item.record.name),
    )
    for occurrence in unique:
        chosen.setdefault(occurrence.record.name, occurrence)
        if len(chosen) >= MAX_MARKERS:
            break

    return sorted(
        chosen.values(),
        key=lambda occurrence: (
            occurrence.manifest_name,
            occurrence.order_index,
            occurrence.record.name,
        ),
    )[:MAX_MARKERS]


def genotype_for(sample_index: int, marker_index: int) -> str:
    patterns = ["A/A", "A/B", "B/B", "0", "A/B"]
    return patterns[(sample_index + marker_index) % len(patterns)]


def compact_ab_genotype(sample_index: int, marker_index: int) -> str:
    return genotype_for(sample_index, marker_index).replace("/", "")


def native_genotype(record: ManifestRecord, compact_ab: str) -> str:
    if compact_ab == "0":
        return "---"
    alleles = []
    for allele in compact_ab:
        if allele == "A":
            alleles.append(record.ab_a)
        elif allele == "B":
            alleles.append(record.ab_b)
        else:
            alleles.append("-")
    return "".join(alleles)


def write_gsgt_header(
    handle,
    *,
    content: str,
    marker_count: int,
    sample_count: int,
) -> None:
    handle.write("[Header]\n")
    handle.write("GSGT Version\t2.0.4\n")
    handle.write("Processing Date\tvalidation synthetic\n")
    handle.write(f"Content\t\t{content}\n")
    handle.write(f"Num SNPs\t{marker_count}\n")
    handle.write(f"Total SNPs\t{marker_count}\n")
    handle.write(f"Num Samples\t{sample_count}\n")
    handle.write(f"Total Samples\t{sample_count}\n")
    handle.write("[Data]\n")


def write_wide(path: Path, markers: list[Occurrence]) -> None:
    samples = ["SAMPLE001", "SAMPLE002", "SAMPLE003"]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", *[marker.record.name for marker in markers]])
        for sample_index, sample in enumerate(samples):
            writer.writerow(
                [sample]
                + [
                    genotype_for(sample_index, marker_index)
                    for marker_index, _marker in enumerate(markers)
                ]
            )


def write_long(path: Path, markers: list[Occurrence]) -> None:
    samples = ["SAMPLE001", "SAMPLE002", "SAMPLE003"]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", "marker_name", "genotype"])
        for sample_index, sample in enumerate(samples):
            for marker_index, marker in enumerate(markers):
                writer.writerow(
                    [
                        sample,
                        marker.record.name,
                        genotype_for(sample_index, marker_index),
                    ]
                )


def write_illumina_gsgt_matrix(path: Path, markers: list[Occurrence]) -> None:
    samples = ["SAMPLE001", "SAMPLE002", "SAMPLE003"]
    with path.open("w") as handle:
        write_gsgt_header(
            handle,
            content="synthetic_mixed_manifest.bpm",
            marker_count=len(markers),
            sample_count=len(samples),
        )
        handle.write("\t" + "\t".join(samples) + "\n")
        for marker_index, marker in enumerate(markers):
            values = [
                compact_ab_genotype(sample_index, marker_index)
                for sample_index, _sample in enumerate(samples)
            ]
            handle.write(marker.record.name + "\t" + "\t".join(values) + "\n")


def write_illumina_gsgt_long(path: Path, markers: list[Occurrence]) -> None:
    samples = ["SAMPLE001", "SAMPLE002", "SAMPLE003"]
    columns = [
        "SNP Name",
        "Sample ID",
        "Allele1 - Top",
        "Allele2 - Top",
        "Allele1 - Forward",
        "Allele2 - Forward",
        "Allele1 - AB",
        "Allele2 - AB",
        "Allele1 - Design",
        "Allele2 - Design",
        "Allele1 - Plus",
        "Allele2 - Plus",
        "GC Score",
    ]
    with path.open("w") as handle:
        write_gsgt_header(
            handle,
            content="synthetic_mixed_manifest.bpm",
            marker_count=len(markers),
            sample_count=len(samples),
        )
        handle.write("\t".join(columns) + "\n")
        for sample_index, sample in enumerate(samples):
            for marker_index, marker in enumerate(markers):
                ab = compact_ab_genotype(sample_index, marker_index)
                if ab == "0":
                    allele1 = allele2 = "-"
                    native1 = native2 = "-"
                else:
                    allele1, allele2 = list(ab)
                    native = native_genotype(marker.record, ab)
                    native1, native2 = list(native)
                handle.write(
                    "\t".join(
                        [
                            marker.record.name,
                            sample,
                            native1,
                            native2,
                            native1,
                            native2,
                            allele1,
                            allele2,
                            native1,
                            native2,
                            "-",
                            "-",
                            "0.99",
                        ]
                    )
                    + "\n"
                )


def write_affymetrix_axiom_dual_call_matrix(path: Path, markers: list[Occurrence]) -> None:
    samples = ["SAMPLE001", "SAMPLE002", "SAMPLE003"]
    with path.open("w") as handle:
        header = ["probeset_id"]
        for sample in samples:
            header.extend([sample, sample])
        handle.write("\t".join(header) + "\n")
        for marker_index, marker in enumerate(markers):
            values = [marker.record.name]
            for sample_index, _sample in enumerate(samples):
                ab = compact_ab_genotype(sample_index, marker_index)
                values.extend([("NoCall" if ab == "0" else ab), native_genotype(marker.record, ab)])
            handle.write("\t".join(values) + "\n")


def write_plink1(prefix: Path, markers: list[Occurrence]) -> None:
    prefix.with_suffix(".bed").write_bytes(b"synthetic validation placeholder\n")
    prefix.with_suffix(".fam").write_text(
        "\n".join(
            [
                "FAM001 SAMPLE001 0 0 0 -9",
                "FAM001 SAMPLE002 0 0 0 -9",
                "FAM001 SAMPLE003 0 0 0 -9",
            ]
        )
        + "\n"
    )
    with prefix.with_suffix(".bim").open("w") as handle:
        for index, marker in enumerate(markers, start=1):
            handle.write(f"0\t{marker.record.name}\t0\t{index}\tA\tB\n")


def write_plink2(prefix: Path, markers: list[Occurrence]) -> None:
    prefix.with_suffix(".pgen").write_bytes(b"synthetic validation placeholder\n")
    prefix.with_suffix(".psam").write_text(
        "#IID\nSAMPLE001\nSAMPLE002\nSAMPLE003\n"
    )
    with prefix.with_suffix(".pvar").open("w") as handle:
        handle.write("##fileformat=PLINK2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\n")
        for index, marker in enumerate(markers, start=1):
            handle.write(f"0\t{index}\t{marker.record.name}\tA\tB\n")


def write_selection_report(
    path: Path,
    selected: list[Occurrence],
    by_marker: dict[str, list[Occurrence]],
) -> None:
    with path.open("w", newline="") as handle:
        fieldnames = [
            "species",
            "marker_name",
            "category",
            "selected_manifest",
            "selected_manifest_order_index",
            "manifest_count",
            "manifests",
            "ab_a",
            "ab_b",
            "is_indel",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for occurrence in selected:
            manifests = sorted(
                {item.manifest_name for item in by_marker[occurrence.record.name]}
            )
            writer.writerow(
                {
                    "species": occurrence.species,
                    "marker_name": occurrence.record.name,
                    "category": "shared" if len(manifests) > 1 else "unique",
                    "selected_manifest": occurrence.manifest_name,
                    "selected_manifest_order_index": occurrence.order_index,
                    "manifest_count": len(manifests),
                    "manifests": ";".join(manifests),
                    "ab_a": occurrence.record.ab_a,
                    "ab_b": occurrence.record.ab_b,
                    "is_indel": str(occurrence.record.is_indel).lower(),
                }
            )


def write_examples(species_dir: Path) -> None:
    species = species_dir.name
    occurrences = load_occurrences(species, species_dir)
    if not occurrences:
        print(f"No manifest markers found for {species}", flush=True)
        return
    by_marker = marker_index(occurrences)
    selected = choose_markers(occurrences)
    if not selected:
        print(f"No markers selected for {species}", flush=True)
        return

    genotype_dir = species_dir / "genotypes" / "synthetic_mixed_manifest"
    if genotype_dir.exists():
        shutil.rmtree(genotype_dir)
    genotype_dir.mkdir(parents=True)

    write_selection_report(genotype_dir / "marker_selection.csv", selected, by_marker)
    write_wide(genotype_dir / "mixed_manifest_wide_ab.csv", selected)
    write_long(genotype_dir / "mixed_manifest_long_ab.csv", selected)
    write_illumina_gsgt_matrix(genotype_dir / "illumina_gsgt_matrix_ab.txt", selected)
    write_illumina_gsgt_long(genotype_dir / "illumina_gsgt_long_multiformat.txt", selected)
    write_affymetrix_axiom_dual_call_matrix(
        genotype_dir / "affymetrix_axiom_dual_call_matrix.txt", selected
    )
    write_plink1(genotype_dir / "mixed_manifest_plink1_ab", selected)
    write_plink2(genotype_dir / "mixed_manifest_plink2_ab", selected)
    shared_count = sum(1 for marker in selected if len(by_marker[marker.record.name]) > 1)
    print(
        f"Wrote {len(selected)} markers for {species} "
        f"({shared_count} shared, {len(selected) - shared_count} unique) to {genotype_dir}",
        flush=True,
    )


def main() -> int:
    for species_dir in sorted(path for path in SOURCE_ROOT.iterdir() if path.is_dir()):
        write_examples(species_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
