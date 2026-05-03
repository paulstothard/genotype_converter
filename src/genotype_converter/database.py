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


@dataclass(frozen=True)
class SourceFolder:
    species: str
    assembly: str
    root_path: str
    manifest_paths: list[str]
    reference_paths: list[str]


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
    with connect(database_path) as conn:
        rows = conn.execute(
            """
            SELECT * FROM marker_rules
            WHERE source_id = ?
            ORDER BY marker_name
            """,
            (source["id"],),
        ).fetchall()
    if not rows:
        raise ValueError(
            "No marker rules found for "
            f"species={species!r}, assembly={assembly!r}, manifest_name={manifest_name!r}"
        )

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
        for assembly_dir in sorted(path for path in species_dir.iterdir() if path.is_dir()):
            manifests_dir = assembly_dir / "manifests"
            references_dir = assembly_dir / "references"
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
