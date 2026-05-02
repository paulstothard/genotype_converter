from __future__ import annotations

import hashlib
import multiprocessing
import sys
import time
from collections import Counter
from pathlib import Path

from .aligner import VariantAligner, _init_worker, _align_record
from .conversion import compute_conversion
from .manifest import parse_manifest
from .output import (
    BuildStats,
    write_alignment,
    write_conversion,
    write_lookup,
    write_lookup_parquet,
    write_position,
    write_summary,
    write_wide,
)


class _Progress:
    def __init__(self, total: int, label: str, enabled: bool):
        self.total = total
        self.label = label
        self.enabled = enabled and total > 0
        self.count = 0
        self.started = time.monotonic()
        self.last_draw = 0.0

    def update(self, step: int = 1) -> None:
        self.count += step
        if not self.enabled:
            return
        now = time.monotonic()
        if self.count < self.total and now - self.last_draw < 0.2:
            return
        self.last_draw = now
        width = 30
        done = int(width * self.count / self.total)
        bar = "#" * done + "-" * (width - done)
        pct = 100 * self.count / self.total
        elapsed = max(now - self.started, 0.001)
        rate = self.count / elapsed
        remaining = (self.total - self.count) / rate if rate else 0
        sys.stderr.write(
            f"\r{self.label}: [{bar}] {self.count}/{self.total} "
            f"({pct:5.1f}%) ETA {remaining:,.0f}s"
        )
        sys.stderr.flush()

    def finish(self) -> None:
        if self.enabled:
            if self.count < self.total:
                self.update(0)
            sys.stderr.write("\n")
            sys.stderr.flush()


def _progress_message(enabled: bool, message: str) -> None:
    if enabled:
        sys.stderr.write(message + "\n")
        sys.stderr.flush()


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _align_records(
    records,
    reference_path: str,
    workers: int,
    save_alignment: bool,
    progress: bool,
):
    tracker = _Progress(len(records), "Aligning variants", progress)
    try:
        if workers == 1:
            aligner = VariantAligner(
                reference_path,
                save_alignment=save_alignment,
                progress=progress,
            )
            alignments = []
            for record in records:
                alignments.append(aligner.align(record))
                tracker.update()
            return alignments

        chunksize = min(1000, max(1, len(records) // (workers * 20))) if records else 1
        with multiprocessing.Pool(
            workers,
            initializer=_init_worker,
            initargs=(reference_path, save_alignment, progress),
        ) as pool:
            alignments = []
            for alignment in pool.imap(_align_record, records, chunksize=chunksize):
                alignments.append(alignment)
                tracker.update()
            return alignments
    finally:
        tracker.finish()


def run(
    manifest_path: str,
    reference_path: str,
    outdir: str = "output",
    species: str = "all",
    workers: int = 1,
    save_alignment: bool = False,
    save_parquet: bool = False,
    progress: bool = False,
) -> BuildStats:
    if workers < 1:
        workers = 1

    _progress_message(progress, f"Hashing manifest: {manifest_path}")
    manifest_sha = _sha256(manifest_path)
    _progress_message(progress, f"Hashing reference: {reference_path}")
    reference_sha = _sha256(reference_path)

    ref_name = Path(reference_path).stem.replace(".", "_")
    panel_name = Path(manifest_path).stem.replace(".", "_")
    output_name = f"{panel_name}.{ref_name}"
    out_dir = Path(outdir) / species / ref_name
    _progress_message(progress, f"Preparing output directory: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    _progress_message(progress, f"Parsing manifest: {manifest_path}")
    records = parse_manifest(manifest_path)
    _progress_message(progress, f"Loaded {len(records)} markers")
    alignments = _align_records(records, reference_path, workers, save_alignment, progress)

    _progress_message(progress, "Computing conversions")
    results = [compute_conversion(r, a) for r, a in zip(records, alignments)]

    info = [
        f"SPECIES={species}",
        f"REF={ref_name}",
        f"PANEL={panel_name}",
    ]

    output_files: list[str] = []

    path_position = str(out_dir / f"{output_name}.position.csv")
    path_conversion = str(out_dir / f"{output_name}.conversion.csv")
    path_wide = str(out_dir / f"{output_name}.wide.csv")
    path_lookup = str(out_dir / f"{output_name}.lookup.csv")
    path_summary = str(out_dir / f"{output_name}.summary.txt")

    _progress_message(progress, f"Writing position file: {path_position}")
    write_position(results, path_position, info)
    _progress_message(progress, f"Writing conversion file: {path_conversion}")
    write_conversion(results, path_conversion, info)
    _progress_message(progress, f"Writing wide file: {path_wide}")
    write_wide(results, path_wide, info)
    _progress_message(progress, f"Writing lookup file: {path_lookup}")
    write_lookup(results, path_lookup, info)
    output_files += [path_position, path_conversion, path_wide, path_lookup]

    if save_alignment:
        path_aln = str(out_dir / f"{output_name}.alignment.txt")
        _progress_message(progress, f"Writing alignment file: {path_aln}")
        write_alignment(results, path_aln, info)
        output_files.append(path_aln)
    if save_parquet:
        path_parquet = str(out_dir / f"{output_name}.lookup.parquet")
        _progress_message(progress, f"Writing Parquet lookup file: {path_parquet}")
        write_lookup_parquet(results, path_parquet, info)
        output_files.append(path_parquet)

    n_positioned = sum(1 for a in alignments if a.position is not None)
    by_chromosome: Counter = Counter(
        r.chromosome for r in results if r.chromosome is not None
    )

    by_type: Counter = Counter(r.manifest_type for r in records)
    by_determination_type: Counter = Counter(
        a.determination_type.name for a in alignments if a.determination_type is not None
    )

    stats = BuildStats(
        manifest_path=str(Path(manifest_path).resolve()),
        manifest_sha256=manifest_sha,
        reference_path=str(Path(reference_path).resolve()),
        reference_sha256=reference_sha,
        species=species,
        workers=workers,
        total_markers=len(records),
        n_snp=sum(1 for r in records if r.is_snp),
        n_indel=sum(1 for r in records if r.is_indel),
        by_manifest_type=dict(by_type),
        n_positioned=n_positioned,
        n_not_positioned=len(alignments) - n_positioned,
        by_chromosome=dict(by_chromosome),
        by_determination_type=dict(by_determination_type),
        output_files=output_files + [path_summary],
    )

    _progress_message(progress, f"Writing summary: {path_summary}")
    write_summary(stats, path_summary)

    return stats
