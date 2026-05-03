from __future__ import annotations

import csv
import hashlib
import json
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = 1

LOOKUP_COLUMNS = [
    "marker_name",
    "alt_marker_name",
    "chromosome",
    "position",
    "ref_allele",
    "alt_allele",
    "determination_type",
    "A_in_AB",
    "B_in_AB",
    "A_in_TOP",
    "B_in_TOP",
    "A_in_FORWARD",
    "B_in_FORWARD",
    "A_in_DESIGN",
    "B_in_DESIGN",
    "A_in_PLUS",
    "B_in_PLUS",
    "A_vcf",
    "B_vcf",
]


@dataclass(frozen=True)
class LookupImportStats:
    database_path: str
    lookup_path: str
    source_id: int
    rows_imported: int
    rows_replaced: int
    rows_skipped_duplicate_marker: int = 0
    duplicate_marker_names: tuple[str, ...] = ()


@dataclass(frozen=True)
class SourceFolder:
    species: str
    assembly: str
    root_path: str
    manifest_root_path: str
    reference_root_path: str
    manifest_paths: list[str]
    reference_paths: list[str]


@dataclass(frozen=True)
class InferredLookupSource:
    source_id: int
    species: str
    assembly: str
    manifest_name: str
    markers_matched: int
    markers_requested: int
    source_markers_total: int


@dataclass(frozen=True)
class MarkerResolution:
    marker_name: str
    status: str
    selected_manifest_name: str
    selected_source_id: int | None
    candidate_manifest_names: str
    candidate_count: int
    reason: str


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def connect(database_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(database_path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def init_database(database_path: str) -> None:
    path = Path(database_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with connect(str(path)) as conn:
        conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS lookup_sources (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                species TEXT NOT NULL,
                assembly TEXT NOT NULL,
                manifest_name TEXT NOT NULL,
                manifest_path TEXT,
                manifest_sha256 TEXT,
                reference_name TEXT,
                reference_path TEXT,
                reference_sha256 TEXT,
                lookup_path TEXT NOT NULL,
                lookup_sha256 TEXT NOT NULL,
                imported_at TEXT NOT NULL,
                tool_version TEXT,
                notes TEXT,
                UNIQUE(species, assembly, manifest_name, lookup_sha256)
            );

            CREATE TABLE IF NOT EXISTS marker_rules (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                source_id INTEGER NOT NULL REFERENCES lookup_sources(id) ON DELETE CASCADE,
                marker_name TEXT NOT NULL,
                alt_marker_name TEXT,
                chromosome TEXT,
                position INTEGER,
                ref_allele TEXT,
                alt_allele TEXT,
                determination_type TEXT,
                A_in_AB TEXT,
                B_in_AB TEXT,
                A_in_TOP TEXT,
                B_in_TOP TEXT,
                A_in_FORWARD TEXT,
                B_in_FORWARD TEXT,
                A_in_DESIGN TEXT,
                B_in_DESIGN TEXT,
                A_in_PLUS TEXT,
                B_in_PLUS TEXT,
                A_vcf TEXT,
                B_vcf TEXT,
                UNIQUE(source_id, marker_name)
            );

            CREATE INDEX IF NOT EXISTS idx_marker_rules_marker
                ON marker_rules(marker_name);
            CREATE INDEX IF NOT EXISTS idx_marker_rules_source_marker
                ON marker_rules(source_id, marker_name);
            CREATE INDEX IF NOT EXISTS idx_lookup_sources_context
                ON lookup_sources(species, assembly, manifest_name);
            """
        )
        conn.execute(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
            ("schema_version", str(SCHEMA_VERSION)),
        )


def _lookup_rows(lookup_path: str) -> list[dict[str, str]]:
    with open(lookup_path, newline="") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    reader = csv.DictReader(lines)
    if reader.fieldnames is None:
        raise ValueError(f"Lookup file has no header row: {lookup_path}")
    missing = [column for column in LOOKUP_COLUMNS if column not in reader.fieldnames]
    if missing:
        raise ValueError(
            f"Lookup file is missing required column(s): {', '.join(missing)}"
        )
    return [dict(row) for row in reader if row.get("marker_name")]


def _deduplicate_lookup_rows(
    rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], int, tuple[str, ...]]:
    by_marker: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_marker.setdefault(row["marker_name"], []).append(row)

    deduplicated: list[dict[str, str]] = []
    skipped = 0
    conflicting_markers: list[str] = []
    for marker_name, marker_rows in by_marker.items():
        if len(marker_rows) == 1:
            deduplicated.append(marker_rows[0])
            continue

        signatures = {
            tuple(str(row.get(column) or "") for column in LOOKUP_COLUMNS[1:])
            for row in marker_rows
        }
        if len(signatures) == 1:
            deduplicated.append(marker_rows[0])
            continue

        skipped += len(marker_rows)
        conflicting_markers.append(marker_name)

    return deduplicated, skipped, tuple(sorted(conflicting_markers))


def import_lookup(
    database_path: str,
    lookup_path: str,
    species: str,
    assembly: str,
    manifest_name: str,
    *,
    manifest_path: str | None = None,
    reference_name: str | None = None,
    reference_path: str | None = None,
    tool_version: str | None = None,
    notes: str | None = None,
    replace: bool = False,
) -> LookupImportStats:
    init_database(database_path)
    rows = _lookup_rows(lookup_path)
    rows, rows_skipped_duplicate_marker, duplicate_marker_names = (
        _deduplicate_lookup_rows(rows)
    )
    lookup_sha = sha256_file(lookup_path)
    manifest_sha = sha256_file(manifest_path) if manifest_path else None
    reference_sha = sha256_file(reference_path) if reference_path else None
    if reference_name is None and reference_path:
        reference_name = Path(reference_path).name

    with connect(database_path) as conn:
        existing = conn.execute(
            """
            SELECT id FROM lookup_sources
            WHERE species = ? AND assembly = ? AND manifest_name = ? AND lookup_sha256 = ?
            """,
            (species, assembly, manifest_name, lookup_sha),
        ).fetchone()
        rows_replaced = 0
        if existing and not replace:
            raise ValueError(
                "Lookup source already exists in database. Use replace=True to refresh it."
            )
        if existing and replace:
            conn.execute("DELETE FROM lookup_sources WHERE id = ?", (existing["id"],))
            rows_replaced = 1

        cursor = conn.execute(
            """
            INSERT INTO lookup_sources (
                species, assembly, manifest_name, manifest_path, manifest_sha256,
                reference_name, reference_path, reference_sha256,
                lookup_path, lookup_sha256, imported_at, tool_version, notes
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                species,
                assembly,
                manifest_name,
                manifest_path,
                manifest_sha,
                reference_name,
                reference_path,
                reference_sha,
                lookup_path,
                lookup_sha,
                _utc_now(),
                tool_version,
                notes,
            ),
        )
        source_id = int(cursor.lastrowid)
        for row in rows:
            conn.execute(
                """
                INSERT INTO marker_rules (
                    source_id, marker_name, alt_marker_name, chromosome, position,
                    ref_allele, alt_allele, determination_type,
                    A_in_AB, B_in_AB, A_in_TOP, B_in_TOP,
                    A_in_FORWARD, B_in_FORWARD, A_in_DESIGN, B_in_DESIGN,
                    A_in_PLUS, B_in_PLUS, A_vcf, B_vcf
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    source_id,
                    row.get("marker_name", ""),
                    row.get("alt_marker_name", ""),
                    row.get("chromosome", ""),
                    _int_or_none(row.get("position", "")),
                    row.get("ref_allele", ""),
                    row.get("alt_allele", ""),
                    row.get("determination_type", ""),
                    row.get("A_in_AB", ""),
                    row.get("B_in_AB", ""),
                    row.get("A_in_TOP", ""),
                    row.get("B_in_TOP", ""),
                    row.get("A_in_FORWARD", ""),
                    row.get("B_in_FORWARD", ""),
                    row.get("A_in_DESIGN", ""),
                    row.get("B_in_DESIGN", ""),
                    row.get("A_in_PLUS", ""),
                    row.get("B_in_PLUS", ""),
                    row.get("A_vcf", ""),
                    row.get("B_vcf", ""),
                ),
            )

    return LookupImportStats(
        database_path=database_path,
        lookup_path=lookup_path,
        source_id=source_id,
        rows_imported=len(rows),
        rows_replaced=rows_replaced,
        rows_skipped_duplicate_marker=rows_skipped_duplicate_marker,
        duplicate_marker_names=duplicate_marker_names,
    )


def _int_or_none(value: str | None) -> int | None:
    if value is None or value == "":
        return None
    return int(value)


def list_species(database_path: str) -> list[str]:
    with connect(database_path) as conn:
        rows = conn.execute(
            "SELECT DISTINCT species FROM lookup_sources ORDER BY species"
        ).fetchall()
    return [row["species"] for row in rows]


def list_assemblies(database_path: str, species: str) -> list[str]:
    with connect(database_path) as conn:
        rows = conn.execute(
            """
            SELECT DISTINCT assembly FROM lookup_sources
            WHERE species = ?
            ORDER BY assembly
            """,
            (species,),
        ).fetchall()
    return [row["assembly"] for row in rows]


def list_manifests(
    database_path: str,
    *,
    species: str | None = None,
    assembly: str | None = None,
) -> list[dict[str, Any]]:
    query = "SELECT * FROM lookup_sources"
    clauses: list[str] = []
    params: list[str] = []
    if species:
        clauses.append("species = ?")
        params.append(species)
    if assembly:
        clauses.append("assembly = ?")
        params.append(assembly)
    if clauses:
        query += " WHERE " + " AND ".join(clauses)
    query += " ORDER BY species, assembly, manifest_name, id"
    with connect(database_path) as conn:
        rows = conn.execute(query, params).fetchall()
    return [dict(row) for row in rows]


def query_marker(
    database_path: str,
    species: str,
    assembly: str,
    marker_name: str,
    *,
    manifest_name: str | None = None,
) -> list[dict[str, Any]]:
    clauses = [
        "s.species = ?",
        "s.assembly = ?",
        "r.marker_name = ?",
    ]
    params = [species, assembly, marker_name]
    if manifest_name:
        clauses.append("s.manifest_name = ?")
        params.append(manifest_name)
    query = f"""
        SELECT
            s.species, s.assembly, s.manifest_name, s.lookup_path,
            r.*
        FROM marker_rules r
        JOIN lookup_sources s ON s.id = r.source_id
        WHERE {" AND ".join(clauses)}
        ORDER BY s.manifest_name, r.id
    """
    with connect(database_path) as conn:
        rows = conn.execute(query, params).fetchall()
    return [dict(row) for row in rows]


def load_lookup_table_from_database(
    database_path: str,
    species: str,
    assembly: str,
    manifest_name: str,
) -> dict[str, list[dict[str, str]]]:
    source = _single_source(database_path, species, assembly, manifest_name)
    return load_lookup_table_from_database_source(database_path, int(source["id"]))


def load_lookup_table_from_database_source(
    database_path: str,
    source_id: int,
) -> dict[str, list[dict[str, str]]]:
    with connect(database_path) as conn:
        rows = conn.execute(
            """
            SELECT * FROM marker_rules
            WHERE source_id = ?
            ORDER BY marker_name
            """,
            (source_id,),
        ).fetchall()
    if not rows:
        raise ValueError(f"No marker rules found for lookup source id {source_id}")

    table: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        row_dict = dict(row)
        marker = row_dict["marker_name"]
        row_a = {
            "AB": row_dict.get("A_in_AB") or "A",
            "TOP": row_dict.get("A_in_TOP") or "",
            "FORWARD": row_dict.get("A_in_FORWARD") or "",
            "DESIGN": row_dict.get("A_in_DESIGN") or "",
            "PLUS": row_dict.get("A_in_PLUS") or "",
            "VCF": row_dict.get("A_vcf") or "",
            "chromosome": row_dict.get("chromosome") or "",
            "position": str(row_dict.get("position") or ""),
        }
        row_b = {
            "AB": row_dict.get("B_in_AB") or "B",
            "TOP": row_dict.get("B_in_TOP") or "",
            "FORWARD": row_dict.get("B_in_FORWARD") or "",
            "DESIGN": row_dict.get("B_in_DESIGN") or "",
            "PLUS": row_dict.get("B_in_PLUS") or "",
            "VCF": row_dict.get("B_vcf") or "",
            "chromosome": row_dict.get("chromosome") or "",
            "position": str(row_dict.get("position") or ""),
        }
        table[marker] = [row_a, row_b]
    return table


def _rule_pair_from_row(row_dict: dict[str, Any]) -> list[dict[str, str]]:
    row_a = {
        "AB": row_dict.get("A_in_AB") or "A",
        "TOP": row_dict.get("A_in_TOP") or "",
        "FORWARD": row_dict.get("A_in_FORWARD") or "",
        "DESIGN": row_dict.get("A_in_DESIGN") or "",
        "PLUS": row_dict.get("A_in_PLUS") or "",
        "VCF": row_dict.get("A_vcf") or "",
        "chromosome": row_dict.get("chromosome") or "",
        "position": str(row_dict.get("position") or ""),
    }
    row_b = {
        "AB": row_dict.get("B_in_AB") or "B",
        "TOP": row_dict.get("B_in_TOP") or "",
        "FORWARD": row_dict.get("B_in_FORWARD") or "",
        "DESIGN": row_dict.get("B_in_DESIGN") or "",
        "PLUS": row_dict.get("B_in_PLUS") or "",
        "VCF": row_dict.get("B_vcf") or "",
        "chromosome": row_dict.get("chromosome") or "",
        "position": str(row_dict.get("position") or ""),
    }
    return [row_a, row_b]


def _rule_signature(row_dict: dict[str, Any]) -> tuple[str, ...]:
    return tuple(str(row_dict.get(column) or "") for column in LOOKUP_COLUMNS[3:])


def load_mixed_lookup_table_from_database(
    database_path: str,
    species: str,
    assembly: str,
    marker_names: list[str],
    *,
    window_size: int = 10,
    on_ambiguous: str = "fail",
) -> tuple[dict[str, list[dict[str, str]]], list[MarkerResolution]]:
    ordered_markers = [marker for marker in marker_names if marker]
    unique_markers = list(dict.fromkeys(ordered_markers))
    if not unique_markers:
        raise ValueError("Cannot resolve mixed manifests because no marker names were found.")
    if on_ambiguous not in {"fail", "skip"}:
        raise ValueError("on_ambiguous must be 'fail' or 'skip'")

    placeholders = ",".join("?" for _ in unique_markers)
    with connect(database_path) as conn:
        rows = conn.execute(
            f"""
            SELECT
                s.id AS source_id,
                s.manifest_name,
                r.*
            FROM marker_rules r
            JOIN lookup_sources s ON s.id = r.source_id
            WHERE s.species = ?
              AND s.assembly = ?
              AND r.marker_name IN ({placeholders})
            ORDER BY r.marker_name, s.manifest_name, s.id
            """,
            (species, assembly, *unique_markers),
        ).fetchall()

    by_marker: dict[str, list[dict[str, Any]]] = {}
    markers_by_source: dict[int, set[str]] = {}
    for row in rows:
        row_dict = dict(row)
        marker = row_dict["marker_name"]
        by_marker.setdefault(marker, []).append(row_dict)
        markers_by_source.setdefault(int(row_dict["source_id"]), set()).add(marker)

    first_index: dict[str, int] = {}
    for index, marker in enumerate(ordered_markers):
        first_index.setdefault(marker, index)

    table: dict[str, list[dict[str, str]]] = {}
    resolutions: list[MarkerResolution] = []
    ambiguous: list[str] = []
    for marker in unique_markers:
        candidates = by_marker.get(marker, [])
        candidate_names = ";".join(row["manifest_name"] for row in candidates)
        if not candidates:
            resolutions.append(
                MarkerResolution(
                    marker_name=marker,
                    status="missing",
                    selected_manifest_name="",
                    selected_source_id=None,
                    candidate_manifest_names="",
                    candidate_count=0,
                    reason="no_source_contains_marker",
                )
            )
            table[marker] = []
            continue

        selected: dict[str, Any] | None = None
        reason = ""
        status = "resolved"
        signatures = {_rule_signature(row) for row in candidates}
        if len(candidates) == 1:
            selected = candidates[0]
            reason = "single_candidate"
        elif len(signatures) == 1:
            selected = candidates[0]
            reason = "duplicate_identical_rules"
        else:
            center = first_index[marker]
            start = max(0, center - window_size)
            end = min(len(ordered_markers), center + window_size + 1)
            window_markers = set(ordered_markers[start:end])
            scores = {
                int(row["source_id"]): len(markers_by_source[int(row["source_id"])] & window_markers)
                for row in candidates
            }
            best_score = max(scores.values())
            best_ids = [source_id for source_id, score in scores.items() if score == best_score]
            if len(best_ids) == 1:
                selected = next(row for row in candidates if int(row["source_id"]) == best_ids[0])
                reason = f"best_local_window_match:{best_score}/{len(window_markers)}"
            else:
                status = "ambiguous"
                reason = "tied_local_window_match"

        if selected is None:
            ambiguous.append(marker)
            resolutions.append(
                MarkerResolution(
                    marker_name=marker,
                    status=status,
                    selected_manifest_name="",
                    selected_source_id=None,
                    candidate_manifest_names=candidate_names,
                    candidate_count=len(candidates),
                    reason=reason,
                )
            )
            table[marker] = []
            continue

        table[marker] = _rule_pair_from_row(selected)
        resolutions.append(
            MarkerResolution(
                marker_name=marker,
                status=status,
                selected_manifest_name=str(selected["manifest_name"]),
                selected_source_id=int(selected["source_id"]),
                candidate_manifest_names=candidate_names,
                candidate_count=len(candidates),
                reason=reason,
            )
        )

    if ambiguous and on_ambiguous == "fail":
        preview = ", ".join(ambiguous[:10])
        more = f" and {len(ambiguous) - 10} more" if len(ambiguous) > 10 else ""
        raise ValueError(
            "Ambiguous marker rule(s) could not be resolved: "
            f"{preview}{more}. Use --on-ambiguous-marker skip to leave them unchanged."
        )
    return table, resolutions


def infer_lookup_source_for_markers(
    database_path: str,
    species: str,
    assembly: str,
    marker_names: list[str],
) -> InferredLookupSource:
    ordered_markers = [marker for marker in dict.fromkeys(marker_names) if marker]
    if not ordered_markers:
        raise ValueError("Cannot infer a manifest because no marker names were found.")

    placeholders = ",".join("?" for _ in ordered_markers)
    with connect(database_path) as conn:
        rows = conn.execute(
            f"""
            SELECT
                s.id AS source_id,
                s.species,
                s.assembly,
                s.manifest_name,
                COUNT(r.marker_name) AS markers_matched,
                (
                    SELECT COUNT(*)
                    FROM marker_rules all_rules
                    WHERE all_rules.source_id = s.id
                ) AS source_markers_total
            FROM lookup_sources s
            LEFT JOIN marker_rules r
                ON r.source_id = s.id
                AND r.marker_name IN ({placeholders})
            WHERE s.species = ? AND s.assembly = ?
            GROUP BY s.id
            HAVING markers_matched > 0
            ORDER BY markers_matched DESC, s.manifest_name, s.id
            """,
            (*ordered_markers, species, assembly),
        ).fetchall()

    if not rows:
        preview = ", ".join(ordered_markers[:10])
        more = f" and {len(ordered_markers) - 10} more" if len(ordered_markers) > 10 else ""
        raise ValueError(
            "Could not infer a manifest for "
            f"species={species!r}, assembly={assembly!r}; no source matched input marker(s): "
            f"{preview}{more}"
        )

    best = dict(rows[0])
    tied = [
        dict(row) for row in rows
        if row["markers_matched"] == best["markers_matched"]
    ]
    if len(tied) > 1:
        names = ", ".join(row["manifest_name"] for row in tied)
        raise ValueError(
            "Could not infer a unique manifest for "
            f"species={species!r}, assembly={assembly!r}; tied sources matched "
            f"{best['markers_matched']}/{len(ordered_markers)} input marker(s): {names}. "
            "Provide --manifest-name."
        )

    return InferredLookupSource(
        source_id=int(best["source_id"]),
        species=str(best["species"]),
        assembly=str(best["assembly"]),
        manifest_name=str(best["manifest_name"]),
        markers_matched=int(best["markers_matched"]),
        markers_requested=len(ordered_markers),
        source_markers_total=int(best["source_markers_total"]),
    )


def _single_source(
    database_path: str,
    species: str,
    assembly: str,
    manifest_name: str,
) -> dict[str, Any]:
    with connect(database_path) as conn:
        rows = conn.execute(
            """
            SELECT * FROM lookup_sources
            WHERE species = ? AND assembly = ? AND manifest_name = ?
            ORDER BY id
            """,
            (species, assembly, manifest_name),
        ).fetchall()
    if not rows:
        raise ValueError(
            "No lookup source found for "
            f"species={species!r}, assembly={assembly!r}, manifest_name={manifest_name!r}"
        )
    if len(rows) > 1:
        raise ValueError(
            "Multiple lookup sources found for "
            f"species={species!r}, assembly={assembly!r}, manifest_name={manifest_name!r}. "
            "Use a distinct manifest name or replace stale imports before converting."
        )
    return dict(rows[0])


def discover_source_folders(source_root: str) -> list[SourceFolder]:
    root = Path(source_root)
    if not root.is_dir():
        raise ValueError(f"Database source root does not exist: {source_root}")

    discovered: list[SourceFolder] = []
    for species_dir in sorted(path for path in root.iterdir() if path.is_dir()):
        shared_manifests_dir = species_dir / "manifests"
        shared_references_dir = species_dir / "references"
        shared_manifest_paths = _discover_files(shared_manifests_dir, {".csv", ".txt"})
        uses_shared_layout = shared_manifests_dir.is_dir() or shared_references_dir.is_dir()
        if shared_references_dir.is_dir():
            for reference_dir in sorted(
                path for path in shared_references_dir.iterdir() if path.is_dir()
            ):
                reference_paths = _discover_files(
                    reference_dir,
                    {".fa", ".fasta", ".fna", ".fas"},
                )
                discovered.append(
                    SourceFolder(
                        species=species_dir.name,
                        assembly=reference_dir.name,
                        root_path=str(species_dir),
                        manifest_root_path=str(shared_manifests_dir),
                        reference_root_path=str(reference_dir),
                        manifest_paths=[str(path) for path in shared_manifest_paths],
                        reference_paths=[str(path) for path in reference_paths],
                    )
                )
        if uses_shared_layout:
            continue

        legacy_skip = {"manifests", "references", "genotypes", "expected"}
        for assembly_dir in sorted(
            path
            for path in species_dir.iterdir()
            if path.is_dir() and path.name not in legacy_skip
        ):
            manifests_dir = assembly_dir / "manifests"
            references_dir = assembly_dir / "references"
            if not manifests_dir.is_dir() and not references_dir.is_dir():
                continue
            manifest_paths = _discover_files(manifests_dir, {".csv", ".txt"})
            reference_paths = _discover_files(
                references_dir,
                {".fa", ".fasta", ".fna", ".fas"},
            )
            discovered.append(
                SourceFolder(
                    species=species_dir.name,
                    assembly=assembly_dir.name,
                    root_path=str(assembly_dir),
                    manifest_root_path=str(manifests_dir),
                    reference_root_path=str(references_dir),
                    manifest_paths=[str(path) for path in manifest_paths],
                    reference_paths=[str(path) for path in reference_paths],
                )
            )
    return discovered


def _discover_files(directory: Path, suffixes: set[str]) -> list[Path]:
    if not directory.is_dir():
        return []
    return sorted(
        path
        for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() in suffixes
    )


def rows_to_csv(rows: list[dict[str, Any]], fieldnames: list[str]) -> str:
    from io import StringIO

    handle = StringIO()
    writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    return handle.getvalue().rstrip("\r\n")


def rows_to_json(rows: list[dict[str, Any]]) -> str:
    return json.dumps(rows, indent=2, sort_keys=True)
