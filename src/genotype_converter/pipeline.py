from __future__ import annotations

import hashlib
import multiprocessing
import os
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


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def run(
    manifest_path: str,
    reference_path: str,
    outdir: str = "output",
    species: str = "all",
    workers: int = 0,
    save_alignment: bool = False,
    save_parquet: bool = False,
) -> BuildStats:
    if workers <= 0:
        workers = max(1, os.cpu_count() or 1)

    manifest_sha = _sha256(manifest_path)
    reference_sha = _sha256(reference_path)

    ref_name = Path(reference_path).stem.replace(".", "_")
    panel_name = Path(manifest_path).stem.replace(".", "_")
    output_name = f"{panel_name}.{ref_name}"
    out_dir = Path(outdir) / species / ref_name
    out_dir.mkdir(parents=True, exist_ok=True)

    records = parse_manifest(manifest_path)

    if workers == 1:
        aligner = VariantAligner(reference_path)
        alignments = [aligner.align(r) for r in records]
    else:
        with multiprocessing.Pool(
            workers,
            initializer=_init_worker,
            initargs=(reference_path,),
        ) as pool:
            alignments = pool.map(_align_record, records)

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

    write_position(results, path_position, info)
    write_conversion(results, path_conversion, info)
    write_wide(results, path_wide, info)
    write_lookup(results, path_lookup, info)
    output_files += [path_position, path_conversion, path_wide, path_lookup]

    if save_alignment:
        path_aln = str(out_dir / f"{output_name}.alignment.txt")
        write_alignment(results, path_aln, info)
        output_files.append(path_aln)
    if save_parquet:
        path_parquet = str(out_dir / f"{output_name}.lookup.parquet")
        write_lookup_parquet(results, path_parquet, info)
        output_files.append(path_parquet)

    n_positioned = sum(1 for a in alignments if a.position is not None)
    by_chromosome: Counter = Counter(
        r.chromosome for r in results if r.chromosome is not None
    )

    by_type: Counter = Counter(r.manifest_type for r in records)

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
        output_files=output_files + [path_summary],
    )

    write_summary(stats, path_summary)

    return stats
