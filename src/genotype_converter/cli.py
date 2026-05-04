from __future__ import annotations

import csv
from pathlib import Path

import click

from .pipeline import run
from .convert_genotypes import load_lookup_table, convert_file, convert_batch
from .convert_plink import (
    convert_plink_bfile,
    convert_plink_bfile_batch,
    convert_plink_pfile,
    convert_plink_pfile_batch,
)
from .database import (
    database_stats,
    discover_source_folders,
    duplicate_marker_report,
    infer_lookup_source_for_markers,
    import_lookup,
    init_database,
    list_import_warnings,
    list_assemblies,
    list_manifests,
    list_species,
    load_mixed_lookup_table_from_database,
    load_lookup_table_from_database,
    load_lookup_table_from_database_source,
    query_marker,
    remove_source,
    rows_to_csv,
    rows_to_json,
    source_details,
    unresolved_marker_report,
    vacuum_database,
    validate_database,
)


@click.group()
def main():
    """Build genotype conversion files and convert genotype data between formats."""


def _require_one_input(single_value, batch_value, single_flag: str, batch_flag: str) -> None:
    if bool(single_value) == bool(batch_value):
        raise click.UsageError(f"Provide exactly one of {single_flag} or {batch_flag}.")


def _require_output_for_mode(output_value, output_flag: str, input_flag: str) -> None:
    if not output_value:
        raise click.UsageError(f"{output_flag} is required with {input_flag}.")


def _echo_csv_batch_summary(stats, outdir) -> None:
    click.echo(
        f"Converted {stats.files_converted}/{stats.files_total} genotype files "
        f"to {Path(outdir)}"
    )
    for output_path in stats.output_paths:
        click.echo(f"Wrote {output_path}")
    click.echo(f"Batch summary: {stats.summary_path}")


def _echo_csv_single_summary(stats, output_path: str) -> None:
    click.echo(
        f"Converted genotypes written to {output_path}: "
        f"{stats.genotypes_changed}/{stats.genotype_cells_total} genotype cells changed, "
        f"{stats.alleles_changed} allele labels changed."
    )
    if stats.missing_or_unparsed_genotypes:
        click.echo(
            "Missing or unparsed genotype cells left unchanged: "
            f"{stats.missing_or_unparsed_genotypes}"
        )
    if stats.unknown_alleles:
        click.echo(f"Unknown allele labels left unchanged: {stats.unknown_alleles}")
    if stats.markers_missing_lookup:
        click.echo(f"Markers missing from lookup: {stats.markers_missing_lookup}")
    if stats.markers_incomplete_mapping:
        click.echo(f"Markers with incomplete allele mapping: {stats.markers_incomplete_mapping}")
    if stats.markers_excluded:
        click.echo(
            f"Markers excluded from output: {stats.markers_excluded} "
            f"({stats.genotypes_excluded} genotype cell/row value(s))"
        )
    if stats.marker_report_path:
        click.echo(f"Marker conversion report: {stats.marker_report_path}")
    if stats.exclude_marker_path:
        click.echo(f"Exclude marker list: {stats.exclude_marker_path}")


def _echo_plink_single_summary(stats, label: str) -> None:
    click.echo(
        f"Converted {label}: {stats.variants_converted}/{stats.variants_total} "
        f"variants, {stats.alleles_changed} allele labels changed."
    )
    if stats.variants_missing_lookup:
        click.echo(f"Markers missing from lookup: {stats.variants_missing_lookup}")
    if stats.variants_incomplete_mapping:
        click.echo(
            "Markers with incomplete allele mapping: "
            f"{stats.variants_incomplete_mapping}"
        )
    if stats.variants_excluded:
        click.echo(f"Markers excluded from output: {stats.variants_excluded}")
    if stats.marker_report_path:
        click.echo(f"Marker conversion report: {stats.marker_report_path}")
    if stats.exclude_marker_path:
        click.echo(f"PLINK exclude list: {stats.exclude_marker_path}")
    click.echo(f"Wrote {stats.genotype_path}")
    click.echo(f"Wrote {stats.variant_path}")
    click.echo(f"Wrote {stats.sample_path}")


def _echo_plink_batch_summary(stats, outdir, label: str) -> None:
    click.echo(
        f"Converted {stats.filesets_converted}/{stats.filesets_total} "
        f"{label} filesets to {Path(outdir)}"
    )
    for item in stats.stats:
        if item.variants_excluded:
            click.echo(
                f"{Path(item.genotype_path).with_suffix('')}: "
                f"excluded {item.variants_excluded} unconvertible marker(s)"
            )
        click.echo(f"Wrote {item.genotype_path}")
        click.echo(f"Wrote {item.variant_path}")
        click.echo(f"Wrote {item.sample_path}")
    click.echo(f"Batch summary: {stats.summary_path}")


def _echo_rows(rows, fieldnames: list[str], output_format: str) -> None:
    if output_format == "json":
        click.echo(rows_to_json(rows))
        return
    if output_format == "csv":
        click.echo(rows_to_csv(rows, fieldnames))
        return
    click.echo("\t".join(fieldnames))
    for row in rows:
        click.echo(
            "\t".join(
                "" if row.get(field) is None else str(row.get(field, ""))
                for field in fieldnames
            )
        )


def _manifest_name_from_path(path: str) -> str:
    return Path(path).stem.replace(".", "_")


def _write_marker_resolution_report(path: str, resolutions, marker_input_paths=None) -> None:
    fieldnames = [
        "input_path",
        "marker_name",
        "status",
        "selected_manifest_name",
        "selected_source_id",
        "candidate_manifest_names",
        "candidate_count",
        "reason",
    ]
    report_path = Path(path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for item in resolutions:
            row = {field: getattr(item, field) for field in fieldnames if field != "input_path"}
            row["input_path"] = ";".join(
                sorted((marker_input_paths or {}).get(item.marker_name, []))
            )
            writer.writerow(row)


def _default_resolution_report_path(
    output,
    outdir,
    output_prefix,
    genotypes,
    bfile,
    pfile,
    command: str,
) -> str:
    if output:
        return str(Path(output).with_suffix(".manifest_resolution.csv"))
    if output_prefix:
        return f"{output_prefix}.manifest_resolution.csv"
    if outdir:
        return str(Path(outdir) / "manifest_resolution.csv")
    if genotypes:
        return str(Path(genotypes).with_suffix(".manifest_resolution.csv"))
    if bfile:
        return f"{bfile}.manifest_resolution.csv"
    if pfile:
        return f"{pfile}.manifest_resolution.csv"
    return f"{command}.manifest_resolution.csv"


def _load_conversion_table(
    lookup,
    database,
    species,
    assembly,
    manifest_name,
    context: str,
    marker_names: list[str] | None = None,
    marker_input_paths: dict[str, set[str]] | None = None,
    resolve_mixed_manifests: bool = False,
    on_ambiguous_marker: str = "fail",
    resolution_report: str | None = None,
):
    _require_one_input(lookup, database, "--lookup", "--database")
    if lookup:
        return load_lookup_table(lookup)
    missing = [
        flag for flag, value in [
            ("--species", species),
            ("--assembly", assembly),
        ] if not value
    ]
    if missing:
        raise click.UsageError(
            f"{', '.join(missing)} required with --database for {context} conversion."
        )
    if resolve_mixed_manifests and manifest_name:
        raise click.UsageError(
            "Use either --manifest-name or --resolve-mixed-manifests, not both."
        )
    try:
        if resolve_mixed_manifests:
            if marker_names is None:
                raise click.UsageError(
                    f"Input marker names are required for mixed-manifest {context} conversion."
                )
            table, resolutions = load_mixed_lookup_table_from_database(
                database_path=database,
                species=species,
                assembly=assembly,
                marker_names=marker_names,
                on_ambiguous=on_ambiguous_marker,
            )
            if resolution_report:
                _write_marker_resolution_report(
                    resolution_report,
                    resolutions,
                    marker_input_paths=marker_input_paths,
                )
                click.echo(f"Marker resolution report: {resolution_report}")
            resolved = sum(1 for item in resolutions if item.status == "resolved")
            ambiguous = sum(1 for item in resolutions if item.status == "ambiguous")
            missing_count = sum(1 for item in resolutions if item.status == "missing")
            click.echo(
                "Resolved mixed manifests for "
                f"{resolved}/{len(resolutions)} marker(s); "
                f"{ambiguous} ambiguous, {missing_count} missing."
            )
            return table
        if not manifest_name:
            if marker_names is None:
                raise click.UsageError(
                    f"--manifest-name is required with --database for {context} conversion "
                    "unless input marker names are available for manifest inference."
                )
            inferred = infer_lookup_source_for_markers(
                database_path=database,
                species=species,
                assembly=assembly,
                marker_names=marker_names,
            )
            click.echo(
                "Inferred manifest "
                f"{inferred.manifest_name!r} from "
                f"{inferred.markers_matched}/{inferred.markers_requested} input marker(s)."
            )
            return load_lookup_table_from_database_source(
                database_path=database,
                source_id=inferred.source_id,
            )
        return load_lookup_table_from_database(
            database_path=database,
            species=species,
            assembly=assembly,
            manifest_name=manifest_name,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc


def _csv_input_paths(genotypes, genotypes_dir, pattern) -> list[Path]:
    if genotypes:
        return [Path(genotypes)]
    source_dir = Path(genotypes_dir)
    paths = sorted(path for path in source_dir.glob(pattern) if path.is_file())
    if not paths:
        raise click.ClickException(
            f"No genotype files matched pattern {pattern!r} in {genotypes_dir}"
        )
    return paths


def _csv_marker_names(paths: list[Path], layout: str, sample_col: str, marker_col: str) -> list[str]:
    return list(_csv_marker_source_map(paths, layout, sample_col, marker_col))


def _csv_marker_source_map(
    paths: list[Path],
    layout: str,
    sample_col: str,
    marker_col: str,
) -> dict[str, set[str]]:
    sources: dict[str, set[str]] = {}
    for path in paths:
        markers = _genotype_text_markers(path, layout, sample_col, marker_col)
        for marker in markers:
            if marker:
                sources.setdefault(marker, set()).add(str(path))
    return sources


def _gsgt_data_lines(path: Path) -> list[str]:
    lines = path.read_text().splitlines()
    try:
        data_index = next(index for index, line in enumerate(lines) if line.strip() == "[Data]")
    except StopIteration as exc:
        raise click.ClickException(f"{path} is missing a [Data] section") from exc
    return [line for line in lines[data_index + 1 :] if line.strip()]


def _genotype_text_markers(
    path: Path,
    layout: str,
    sample_col: str,
    marker_col: str,
) -> list[str]:
    layout = layout.lower()
    if layout in {"wide", "long"}:
        with path.open(newline="") as handle:
            lines = [line for line in handle if not line.startswith("#")]
        reader = csv.DictReader(lines)
        if reader.fieldnames is None:
            raise click.ClickException(f"Input file has no header row: {path}")
        if layout == "wide":
            return [col for col in reader.fieldnames if col != sample_col]
        if marker_col not in reader.fieldnames:
            raise click.ClickException(f"{path} is missing marker column {marker_col!r}")
        return [row.get(marker_col, "") for row in reader]
    if layout == "illumina-matrix":
        data_lines = _gsgt_data_lines(path)
        return [line.split("\t", 1)[0] for line in data_lines[1:]]
    if layout == "illumina-long":
        data_lines = _gsgt_data_lines(path)
        reader = csv.DictReader(data_lines, delimiter="\t")
        if reader.fieldnames is None:
            raise click.ClickException(f"{path} has no data header row")
        if "SNP Name" not in reader.fieldnames:
            raise click.ClickException(f"{path} is missing SNP Name column")
        return [row.get("SNP Name", "") for row in reader]
    if layout == "affymetrix-matrix":
        lines = [line for line in path.read_text().splitlines() if line.strip()]
        if not lines:
            raise click.ClickException(f"Input file has no header row: {path}")
        return [line.split("\t", 1)[0] for line in lines[1:]]
    raise click.ClickException(f"Unsupported genotype layout: {layout}")


def _plink_bfile_input_prefixes(bfile, bfile_dir, pattern) -> list[Path]:
    if bfile:
        return [Path(bfile)]
    paths = sorted(path for path in Path(bfile_dir).glob(pattern) if path.is_file())
    prefixes = [path.with_suffix("") for path in paths if path.suffix == ".bed"]
    if not prefixes:
        raise click.ClickException(
            f"No PLINK binary filesets matched pattern {pattern!r} in {bfile_dir}. "
            "The pattern should match .bed files."
        )
    return prefixes


def _plink_bfile_marker_names(prefixes: list[Path]) -> list[str]:
    return list(_plink_bfile_marker_source_map(prefixes))


def _plink_bfile_marker_source_map(prefixes: list[Path]) -> dict[str, set[str]]:
    markers: list[str] = []
    sources: dict[str, set[str]] = {}
    for prefix in prefixes:
        bim = prefix.with_suffix(".bim")
        if not bim.exists():
            raise click.ClickException(
                "PLINK binary conversion requires all three files: .bed, .bim, and .fam. "
                f"Missing: {bim}"
            )
        with bim.open() as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not stripped:
                    continue
                fields = stripped.split()
                if len(fields) != 6:
                    raise click.ClickException(
                        f"{bim} line {line_number} has {len(fields)} fields; expected 6"
                    )
                markers.append(fields[1])
                sources.setdefault(fields[1], set()).add(str(prefix))
    return sources


def _plink_pfile_input_prefixes(pfile, pfile_dir, pattern) -> list[Path]:
    if pfile:
        return [Path(pfile)]
    paths = sorted(path for path in Path(pfile_dir).glob(pattern) if path.is_file())
    prefixes = [path.with_suffix("") for path in paths if path.suffix == ".pgen"]
    if not prefixes:
        raise click.ClickException(
            f"No PLINK 2 filesets matched pattern {pattern!r} in {pfile_dir}. "
            "The pattern should match .pgen files."
        )
    return prefixes


def _plink_pfile_marker_names(prefixes: list[Path]) -> list[str]:
    return list(_plink_pfile_marker_source_map(prefixes))


def _plink_pfile_marker_source_map(prefixes: list[Path]) -> dict[str, set[str]]:
    markers: list[str] = []
    sources: dict[str, set[str]] = {}
    for prefix in prefixes:
        pvar = prefix.with_suffix(".pvar")
        if not pvar.exists():
            raise click.ClickException(
                "PLINK 2 conversion requires all three files: .pgen, .pvar, and .psam. "
                f"Missing: {pvar}"
            )
        indexes = None
        with pvar.open() as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not stripped or stripped.startswith("##"):
                    continue
                if stripped.startswith("#CHROM"):
                    header = stripped.split()
                    normalized = ["CHROM" if col == "#CHROM" else col for col in header]
                    if "ID" not in normalized:
                        raise click.ClickException(
                            f"{pvar} #CHROM header is missing ID column"
                        )
                    indexes = {"ID": normalized.index("ID")}
                    continue
                if indexes is None:
                    raise click.ClickException(
                        f"{pvar} line {line_number} appears before a #CHROM header."
                    )
                fields = stripped.split()
                if len(fields) <= indexes["ID"]:
                    raise click.ClickException(
                        f"{pvar} line {line_number} has too few fields for ID column"
                    )
                markers.append(fields[indexes["ID"]])
                sources.setdefault(fields[indexes["ID"]], set()).add(str(prefix))
    return sources


@main.command("build")
@click.option("--manifest", required=True, type=click.Path(exists=True),
              help="Illumina or Affymetrix manifest CSV")
@click.option("--reference", required=True, type=click.Path(exists=True),
              help="Reference genome FASTA")
@click.option("--outdir", default="output", show_default=True,
              help="Output directory")
@click.option("--species", default="all", show_default=True,
              help="Species name for output subdirectory")
@click.option("--workers", default=1, show_default=True,
              help="Worker processes. Each worker loads the reference index, so increase carefully for large genomes.")
@click.option("--align/--no-align", default=False, show_default=True,
              help="Write alignment.txt debug file")
@click.option("--parquet/--no-parquet", default=False, show_default=True,
              help="Also write lookup.parquet (requires pyarrow)")
@click.option("--progress/--no-progress", default=True, show_default=True,
              help="Show build progress while aligning variants")
def build_cmd(manifest, reference, outdir, species, workers, align, parquet, progress):
    """Align manifest variants and build conversion/position files."""
    stats = run(
        manifest_path=manifest,
        reference_path=reference,
        outdir=outdir,
        species=species,
        workers=workers,
        save_alignment=align,
        save_parquet=parquet,
        progress=progress,
    )
    click.echo(
        f"Build complete: {stats.total_markers} markers "
        f"({stats.n_snp} SNPs, {stats.n_indel} indels), "
        f"{stats.n_positioned} positioned, {stats.n_not_positioned} not positioned."
    )
    summary_path = next((f for f in stats.output_files if f.endswith(".summary.txt")), None)
    if summary_path:
        click.echo(f"Full summary: {summary_path}")


@main.group("db")
def db_cmd():
    """Maintain and inspect an optional SQLite conversion database."""


@db_cmd.command("init")
@click.option("--database", required=True, help="SQLite database path")
def db_init_cmd(database):
    """Create or update database tables."""
    init_database(database)
    click.echo(f"Initialized database: {database}")


@db_cmd.command("import-lookup")
@click.option("--database", required=True, help="SQLite database path")
@click.option("--lookup", required=True, type=click.Path(exists=True),
              help="Lookup CSV produced by the build command")
@click.option("--species", required=True, help="Species name")
@click.option("--assembly", required=True, help="Reference assembly name")
@click.option("--manifest-name", required=True,
              help="Manifest or panel name used to qualify marker rules")
@click.option("--manifest-path", required=False, type=click.Path(exists=True),
              help="Original manifest path, used for provenance")
@click.option("--reference-name", required=False,
              help="Reference name, if different from the reference filename")
@click.option("--reference-path", required=False, type=click.Path(exists=True),
              help="Reference FASTA path, used for provenance")
@click.option("--tool-version", required=False, help="Tool version or build label")
@click.option("--notes", required=False, help="Free-text import note")
@click.option("--replace/--no-replace", default=False, show_default=True,
              help="Replace an existing import with the same species, assembly, manifest, and lookup checksum")
@click.option("--replace-context/--no-replace-context", default=False, show_default=True,
              help="Replace all existing imports with the same species, assembly, and manifest")
def db_import_lookup_cmd(database, lookup, species, assembly, manifest_name,
                         manifest_path, reference_name, reference_path,
                         tool_version, notes, replace, replace_context):
    """Import a built lookup CSV into the database."""
    stats = import_lookup(
        database_path=database,
        lookup_path=lookup,
        species=species,
        assembly=assembly,
        manifest_name=manifest_name,
        manifest_path=manifest_path,
        reference_name=reference_name,
        reference_path=reference_path,
        tool_version=tool_version,
        notes=notes,
        replace=replace,
        replace_context=replace_context,
    )
    click.echo(
        f"Imported {stats.rows_imported} marker rules from {lookup} "
        f"as source {stats.source_id}."
    )
    if stats.rows_skipped_duplicate_marker:
        preview = ", ".join(stats.duplicate_marker_names[:10])
        more = (
            f" and {len(stats.duplicate_marker_names) - 10} more"
            if len(stats.duplicate_marker_names) > 10
            else ""
        )
        click.echo(
            "Skipped "
            f"{stats.rows_skipped_duplicate_marker} row(s) with conflicting "
            f"duplicate marker name(s): {preview}{more}."
        )
    if stats.rows_replaced:
        click.echo(f"Replaced {stats.rows_replaced} existing source(s).")


@db_cmd.command("build")
@click.option("--source-root", required=True, type=click.Path(exists=True, file_okay=False),
              help="Root folder organized as species/manifests and species/references/<assembly>")
@click.option("--database", required=True, help="SQLite database path")
@click.option("--build-outdir", default="database_build", show_default=True,
              help="Directory for generated build outputs")
@click.option("--workers", default=1, show_default=True,
              help="Worker processes for each build. Use 1 for large references unless memory is available.")
@click.option("--replace/--no-replace", default=False, show_default=True,
              help="Replace existing imports with the same species, assembly, manifest, and lookup checksum")
@click.option("--replace-context/--no-replace-context", default=False, show_default=True,
              help="Replace existing imports with the same species, assembly, and manifest")
@click.option("--progress/--no-progress", default=True, show_default=True,
              help="Show build progress")
def db_build_cmd(source_root, database, build_outdir, workers, replace, replace_context, progress):
    """Build lookup files from source folders and import them into SQLite."""
    folders = discover_source_folders(source_root)
    if not folders:
        raise click.ClickException(f"No source folders found under {source_root}")

    total_imported = 0
    total_rules = 0
    for folder in folders:
        if not folder.manifest_paths:
            raise click.ClickException(
                f"No manifest files found in {folder.manifest_root_path}"
            )
        if len(folder.reference_paths) != 1:
            raise click.ClickException(
                f"Expected exactly one reference file in "
                f"{folder.reference_root_path}, found {len(folder.reference_paths)}"
            )
        reference_path = folder.reference_paths[0]
        for manifest_path in folder.manifest_paths:
            manifest_name = _manifest_name_from_path(manifest_path)
            source_build_outdir = str(Path(build_outdir) / folder.species / folder.assembly)
            click.echo(
                f"Building {folder.species}/{folder.assembly}/{manifest_name}"
            )
            build_stats = run(
                manifest_path=manifest_path,
                reference_path=reference_path,
                outdir=source_build_outdir,
                species=folder.species,
                workers=workers,
                save_alignment=False,
                save_parquet=False,
                progress=progress,
            )
            lookup_path = next(
                (path for path in build_stats.output_files if path.endswith(".lookup.csv")),
                None,
            )
            if lookup_path is None:
                raise click.ClickException(
                    f"Build did not produce a lookup CSV for {manifest_path}"
                )
            import_stats = import_lookup(
                database_path=database,
                lookup_path=lookup_path,
                species=folder.species,
                assembly=folder.assembly,
                manifest_name=manifest_name,
                manifest_path=manifest_path,
                reference_name=Path(reference_path).name,
                reference_path=reference_path,
                replace=replace,
                replace_context=replace_context,
            )
            total_imported += 1
            total_rules += import_stats.rows_imported
            click.echo(
                f"Imported {import_stats.rows_imported} marker rules from {lookup_path}"
            )
            if import_stats.rows_skipped_duplicate_marker:
                preview = ", ".join(import_stats.duplicate_marker_names[:10])
                more = (
                    f" and {len(import_stats.duplicate_marker_names) - 10} more"
                    if len(import_stats.duplicate_marker_names) > 10
                    else ""
                )
                click.echo(
                    "Skipped "
                    f"{import_stats.rows_skipped_duplicate_marker} row(s) with "
                    f"conflicting duplicate marker name(s): {preview}{more}"
                )

    click.echo(
        f"Database build complete: {total_imported} lookup source(s), "
        f"{total_rules} marker rules imported."
    )


@db_cmd.command("list-species")
@click.option("--database", required=True, type=click.Path(exists=True),
              help="SQLite database path")
def db_list_species_cmd(database):
    """List species in the database."""
    for species in list_species(database):
        click.echo(species)


@db_cmd.command("list-assemblies")
@click.option("--database", required=True, type=click.Path(exists=True),
              help="SQLite database path")
@click.option("--species", required=True, help="Species name")
def db_list_assemblies_cmd(database, species):
    """List assemblies for a species."""
    for assembly in list_assemblies(database, species):
        click.echo(assembly)


@db_cmd.command("list-manifests")
@click.option("--database", required=True, type=click.Path(exists=True),
              help="SQLite database path")
@click.option("--species", required=False, help="Species name")
@click.option("--assembly", required=False, help="Reference assembly name")
def db_list_manifests_cmd(database, species, assembly):
    """List imported manifests."""
    for row in list_manifests(database, species=species, assembly=assembly):
        click.echo(
            f"{row['id']}\t{row['species']}\t{row['assembly']}\t"
            f"{row['manifest_name']}\t{row['lookup_path']}"
        )


@db_cmd.command("marker")
@click.option("--database", required=True, type=click.Path(exists=True),
              help="SQLite database path")
@click.option("--species", required=True, help="Species name")
@click.option("--assembly", required=True, help="Reference assembly name")
@click.option("--marker", "marker_name", required=True, help="Marker name")
@click.option("--manifest-name", required=False,
              help="Optional manifest or panel name to disambiguate duplicate marker names")
@click.option("--format", "output_format", default="table", show_default=True,
              type=click.Choice(["table", "csv", "json"]),
              help="Output format")
def db_marker_cmd(database, species, assembly, marker_name, manifest_name, output_format):
    """Show database rule rows for one marker."""
    rows = query_marker(
        database_path=database,
        species=species,
        assembly=assembly,
        marker_name=marker_name,
        manifest_name=manifest_name,
    )
    if not rows:
        raise click.ClickException("No marker rule found.")
    fieldnames = [
        "species", "assembly", "manifest_name", "marker_name", "chromosome",
        "position", "ref_allele", "alt_allele", "determination_type",
        "A_in_TOP", "B_in_TOP", "A_in_PLUS", "B_in_PLUS", "A_vcf", "B_vcf",
    ]
    _echo_rows(rows, fieldnames, output_format)


@db_cmd.command("discover-sources")
@click.option("--source-root", required=True, type=click.Path(exists=True, file_okay=False),
              help="Root folder organized as species/manifests and species/references/<assembly>")
@click.option("--format", "output_format", default="table", show_default=True,
              type=click.Choice(["table", "csv", "json"]),
              help="Output format")
def db_discover_sources_cmd(source_root, output_format):
    """List database source folders without running build."""
    folders = discover_source_folders(source_root)
    rows = [
        {
            "species": folder.species,
            "assembly": folder.assembly,
            "root_path": folder.root_path,
            "manifest_root_path": folder.manifest_root_path,
            "reference_root_path": folder.reference_root_path,
            "manifest_count": len(folder.manifest_paths),
            "reference_count": len(folder.reference_paths),
            "manifest_paths": ";".join(folder.manifest_paths),
            "reference_paths": ";".join(folder.reference_paths),
        }
        for folder in folders
    ]
    fieldnames = [
        "species", "assembly", "root_path", "manifest_root_path",
        "reference_root_path", "manifest_count", "reference_count",
        "manifest_paths", "reference_paths",
    ]
    _echo_rows(rows, fieldnames, output_format)


@db_cmd.command("stats")
@click.option("--database", required=True, type=click.Path(exists=True),
              help="SQLite database path")
@click.option("--format", "output_format", default="table", show_default=True,
              type=click.Choice(["table", "csv", "json"]),
              help="Output format")
def db_stats_cmd(database, output_format):
    """Summarize database contents."""
    row = database_stats(database)
    fieldnames = [
        "database_path", "schema_version", "species", "assemblies",
        "lookup_sources", "marker_rules", "import_warnings",
        "unpositioned_marker_rules",
    ]
    _echo_rows([row], fieldnames, output_format)


@db_cmd.command("source")
@click.option("--database", required=True, type=click.Path(exists=True),
              help="SQLite database path")
@click.option("--source-id", required=True, type=int, help="Lookup source id")
@click.option("--format", "output_format", default="table", show_default=True,
              type=click.Choice(["table", "csv", "json"]),
              help="Output format")
def db_source_cmd(database, source_id, output_format):
    """Show details for one lookup source."""
    try:
        row = source_details(database, source_id)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc
    fieldnames = [
        "id", "species", "assembly", "manifest_name", "manifest_path",
        "manifest_sha256", "reference_name", "reference_path",
        "reference_sha256", "lookup_path", "lookup_sha256", "imported_at",
        "tool_version", "notes", "marker_rules", "import_warnings",
    ]
    _echo_rows([row], fieldnames, output_format)


@db_cmd.command("remove-source")
@click.option("--database", required=True, type=click.Path(exists=True),
              help="SQLite database path")
@click.option("--source-id", required=False, type=int, help="Lookup source id")
@click.option("--species", required=False, help="Species name")
@click.option("--assembly", required=False, help="Reference assembly name")
@click.option("--manifest-name", required=False, help="Manifest or panel name")
@click.option("--yes", is_flag=True,
              help="Confirm removal")
def db_remove_source_cmd(database, source_id, species, assembly, manifest_name, yes):
    """Remove lookup source(s) and their marker rules."""
    if not yes:
        raise click.UsageError("Refusing to remove a source without --yes.")
    try:
        stats = remove_source(
            database,
            source_id=source_id,
            species=species,
            assembly=assembly,
            manifest_name=manifest_name,
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc
    click.echo(
        f"Removed {stats.sources_removed} lookup source(s) and "
        f"{stats.marker_rules_removed} marker rule(s)."
    )


@db_cmd.command("warnings")
@click.option("--database", required=True, type=click.Path(exists=True),
              help="SQLite database path")
@click.option("--source-id", required=False, type=int, help="Filter to one source id")
@click.option("--format", "output_format", default="table", show_default=True,
              type=click.Choice(["table", "csv", "json"]),
              help="Output format")
def db_warnings_cmd(database, source_id, output_format):
    """List import warnings recorded in the database."""
    rows = list_import_warnings(database, source_id=source_id)
    fieldnames = [
        "id", "source_id", "species", "assembly", "manifest_name",
        "warning_type", "marker_name", "row_count", "message",
    ]
    _echo_rows(rows, fieldnames, output_format)


@db_cmd.command("duplicates")
@click.option("--database", required=True, type=click.Path(exists=True),
              help="SQLite database path")
@click.option("--species", required=False, help="Species name")
@click.option("--assembly", required=False, help="Reference assembly name")
@click.option("--limit", required=False, type=int, help="Maximum rows to print")
@click.option("--format", "output_format", default="table", show_default=True,
              type=click.Choice(["table", "csv", "json"]),
              help="Output format")
def db_duplicates_cmd(database, species, assembly, limit, output_format):
    """Report marker names that occur in multiple lookup sources."""
    rows = duplicate_marker_report(
        database,
        species=species,
        assembly=assembly,
        limit=limit,
    )
    fieldnames = [
        "species", "assembly", "marker_name", "rule_count",
        "source_count", "manifest_names",
    ]
    _echo_rows(rows, fieldnames, output_format)


@db_cmd.command("unresolved")
@click.option("--database", required=True, type=click.Path(exists=True),
              help="SQLite database path")
@click.option("--species", required=False, help="Species name")
@click.option("--assembly", required=False, help="Reference assembly name")
@click.option("--manifest-name", required=False, help="Manifest or panel name")
@click.option("--limit", required=False, type=int, help="Maximum rows to print")
@click.option("--format", "output_format", default="table", show_default=True,
              type=click.Choice(["table", "csv", "json"]),
              help="Output format")
def db_unresolved_cmd(database, species, assembly, manifest_name, limit, output_format):
    """List marker rules without a positioned chromosome/base coordinate."""
    rows = unresolved_marker_report(
        database,
        species=species,
        assembly=assembly,
        manifest_name=manifest_name,
        limit=limit,
    )
    fieldnames = [
        "source_id", "species", "assembly", "manifest_name",
        "marker_name", "determination_type", "chromosome", "position",
    ]
    _echo_rows(rows, fieldnames, output_format)


@db_cmd.command("validate")
@click.option("--database", required=True, type=click.Path(exists=True),
              help="SQLite database path")
@click.option("--format", "output_format", default="table", show_default=True,
              type=click.Choice(["table", "csv", "json"]),
              help="Output format")
def db_validate_cmd(database, output_format):
    """Validate database integrity and common maintenance risks."""
    result = validate_database(database)
    fieldnames = ["check", "status", "count", "details"]
    _echo_rows(result.checks, fieldnames, output_format)
    if output_format == "table":
        click.echo(f"Overall status: {result.status}")
    if result.status == "fail":
        raise click.ClickException("Database validation failed.")


@db_cmd.command("vacuum")
@click.option("--database", required=True, type=click.Path(exists=True),
              help="SQLite database path")
def db_vacuum_cmd(database):
    """Reclaim free space in the SQLite database after deletes."""
    vacuum_database(database)
    click.echo(f"Vacuumed database: {database}")


@main.command("convert")
@click.option("--genotypes", required=False, type=click.Path(exists=True),
              help="Input genotype file (CSV)")
@click.option("--genotypes-dir", required=False, type=click.Path(exists=True, file_okay=False),
              help="Directory of genotype CSV files to convert")
@click.option("--pattern", default="*.csv", show_default=True,
              help="File pattern used with --genotypes-dir")
@click.option("--lookup", required=False, type=click.Path(exists=True),
              help="Lookup CSV produced by the build command (lookup.csv or conversion.csv)")
@click.option("--database", required=False, type=click.Path(exists=True),
              help="SQLite conversion database for CSV conversion")
@click.option("--species", required=False,
              help="Species name for --database CSV conversion")
@click.option("--assembly", required=False,
              help="Reference assembly name for --database CSV conversion")
@click.option("--manifest-name", required=False,
              help="Manifest name for --database CSV conversion; inferred from input markers if omitted")
@click.option("--resolve-mixed-manifests/--no-resolve-mixed-manifests",
              default=False, show_default=True,
              help="Resolve database rules per marker for inputs containing markers from multiple manifests")
@click.option("--on-ambiguous-marker", default="fail", show_default=True,
              type=click.Choice(["fail", "skip"]),
              help="How mixed-manifest mode handles conflicting rules that cannot be resolved")
@click.option("--resolution-report", required=False,
              help="CSV report path for mixed-manifest marker rule decisions")
@click.option("--from-format", "from_fmt", required=True,
              type=click.Choice(["AB", "TOP", "FORWARD", "DESIGN", "PLUS", "VCF"],
                                case_sensitive=False),
              help="Format encoding of the input genotypes")
@click.option("--to-format", "to_fmt", required=True,
              type=click.Choice(["AB", "TOP", "FORWARD", "DESIGN", "PLUS", "VCF"],
                                case_sensitive=False),
              help="Target format encoding for the output")
@click.option("--output", required=False,
              help="Output file path for single-file conversion")
@click.option("--outdir", required=False,
              help="Output directory for --genotypes-dir")
@click.option("--suffix", default=".converted.csv", show_default=True,
              help="Filename suffix for batch outputs")
@click.option("--overwrite/--no-overwrite", default=False, show_default=True,
              help="Allow batch conversion to replace existing output files")
@click.option("--layout", default="wide", show_default=True,
              type=click.Choice(
                  ["wide", "long", "illumina-matrix", "illumina-long", "affymetrix-matrix"],
                  case_sensitive=False,
              ),
              help="Input file layout")
@click.option("--in-sep", default=None,
              help="Allele separator in input (default: auto-detect from /, space, or adjacent)")
@click.option("--out-sep", default="/", show_default=True,
              help="Allele separator in output")
@click.option("--sample-col", default="sample_id", show_default=True,
              help="Column name identifying the sample (wide layout: first column)")
@click.option("--marker-col", default="marker_name", show_default=True,
              help="Column name for marker (long layout only)")
@click.option("--genotype-col", default="genotype", show_default=True,
              help="Column name for genotype (long layout only)")
@click.option("--on-unconvertible-marker", default="exclude", show_default=True,
              type=click.Choice(["exclude", "fail", "keep"]),
              help="How to handle markers missing from lookup or lacking target alleles")
def convert_cmd(genotypes, genotypes_dir, pattern, lookup, database, species,
                assembly, manifest_name, resolve_mixed_manifests,
                on_ambiguous_marker, resolution_report, from_fmt, to_fmt,
                output, outdir, suffix, overwrite, layout, in_sep, out_sep,
                sample_col, marker_col, genotype_col, on_unconvertible_marker):
    """Convert genotypes between format encodings (e.g. TOP → PLUS).

    Examples:

      genotype-converter convert --genotypes mydata.csv --lookup lookup.csv
      --from-format TOP --to-format PLUS --output converted.csv

      genotype-converter convert --genotypes-dir genotypes/ --lookup lookup.csv
      --from-format TOP --to-format PLUS --outdir converted/
    """
    from_fmt = from_fmt.upper()
    to_fmt = to_fmt.upper()
    _require_one_input(genotypes, genotypes_dir, "--genotypes", "--genotypes-dir")
    marker_names = None
    marker_sources = None
    if database and (not manifest_name or resolve_mixed_manifests):
        marker_sources = _csv_marker_source_map(
            _csv_input_paths(genotypes, genotypes_dir, pattern),
            layout=layout,
            sample_col=sample_col,
            marker_col=marker_col,
        )
        marker_names = list(marker_sources)
    if database and resolve_mixed_manifests and not resolution_report:
        resolution_report = _default_resolution_report_path(
            output=output,
            outdir=outdir,
            output_prefix=None,
            genotypes=genotypes,
            bfile=None,
            pfile=None,
            command="convert",
        )
    table = _load_conversion_table(
        lookup=lookup,
        database=database,
        species=species,
        assembly=assembly,
        manifest_name=manifest_name,
        context="CSV",
        marker_names=marker_names,
        marker_input_paths=marker_sources,
        resolve_mixed_manifests=resolve_mixed_manifests,
        on_ambiguous_marker=on_ambiguous_marker,
        resolution_report=resolution_report,
    )

    if genotypes:
        _require_output_for_mode(output, "--output", "--genotypes")
        stats = convert_file(
            input_path=genotypes,
            output_path=output,
            table=table,
            from_fmt=from_fmt,
            to_fmt=to_fmt,
            layout=layout,
            in_sep=in_sep,
            out_sep=out_sep,
            sample_col=sample_col,
            marker_col=marker_col,
            genotype_col=genotype_col,
            on_unconvertible_marker=on_unconvertible_marker,
        )
        _echo_csv_single_summary(stats, output)
        return

    _require_output_for_mode(outdir, "--outdir", "--genotypes-dir")
    try:
        stats = convert_batch(
            input_dir=genotypes_dir,
            output_dir=outdir,
            pattern=pattern,
            suffix=suffix,
            overwrite=overwrite,
            table=table,
            from_fmt=from_fmt,
            to_fmt=to_fmt,
            layout=layout,
            in_sep=in_sep,
            out_sep=out_sep,
            sample_col=sample_col,
            marker_col=marker_col,
            genotype_col=genotype_col,
            on_unconvertible_marker=on_unconvertible_marker,
        )
    except FileExistsError as exc:
        raise click.ClickException(str(exc)) from exc

    _echo_csv_batch_summary(stats, outdir)


@main.command("convert-plink")
@click.option("--bfile", required=False,
              help="Input PLINK binary prefix (expects .bed, .bim, and .fam)")
@click.option("--bfile-dir", required=False, type=click.Path(exists=True, file_okay=False),
              help="Directory of PLINK binary filesets to convert")
@click.option("--pattern", default="*.bed", show_default=True,
              help="File pattern used with --bfile-dir")
@click.option("--lookup", required=False, type=click.Path(exists=True),
              help="Lookup CSV produced by the build command")
@click.option("--database", required=False, type=click.Path(exists=True),
              help="SQLite conversion database")
@click.option("--species", required=False,
              help="Species name for --database conversion")
@click.option("--assembly", required=False,
              help="Reference assembly name for --database conversion")
@click.option("--manifest-name", required=False,
              help="Manifest name for --database conversion; inferred from input markers if omitted")
@click.option("--resolve-mixed-manifests/--no-resolve-mixed-manifests",
              default=False, show_default=True,
              help="Resolve database rules per marker for inputs containing markers from multiple manifests")
@click.option("--on-ambiguous-marker", default="fail", show_default=True,
              type=click.Choice(["fail", "skip"]),
              help="How mixed-manifest mode handles conflicting rules that cannot be resolved")
@click.option("--resolution-report", required=False,
              help="CSV report path for mixed-manifest marker rule decisions")
@click.option("--from-format", "from_fmt", required=True,
              type=click.Choice(["AB", "TOP", "FORWARD", "DESIGN", "PLUS"],
                                case_sensitive=False),
              help="Format encoding of the input .bim allele labels")
@click.option("--to-format", "to_fmt", required=True,
              type=click.Choice(["AB", "TOP", "FORWARD", "DESIGN", "PLUS"],
                                case_sensitive=False),
              help="Target format encoding for the output .bim allele labels")
@click.option("--out", "output_prefix", required=False,
              help="Output PLINK binary prefix for single-fileset conversion")
@click.option("--outdir", required=False,
              help="Output directory for --bfile-dir")
@click.option("--suffix", default=".converted", show_default=True,
              help="Filename suffix for batch output prefixes")
@click.option("--overwrite/--no-overwrite", default=False, show_default=True,
              help="Allow batch conversion to replace existing output files")
@click.option("--update-position/--keep-position", default=False, show_default=True,
              help="Also replace .bim chromosome/base-pair columns with lookup positions")
@click.option("--on-unconvertible-marker", default="exclude", show_default=True,
              type=click.Choice(["exclude", "fail", "keep"]),
              help="How to handle .bim markers missing from lookup or lacking target alleles")
@click.option("--plink", "plink_command", default="plink", show_default=True,
              help="PLINK executable used when --on-unconvertible-marker=exclude")
def convert_plink_cmd(bfile, bfile_dir, pattern, lookup, database, species,
                      assembly, manifest_name, resolve_mixed_manifests,
                      on_ambiguous_marker, resolution_report, from_fmt, to_fmt,
                      output_prefix, outdir, suffix, overwrite, update_position,
                      on_unconvertible_marker, plink_command):
    """Convert allele labels in a PLINK bed/bim/fam fileset."""
    _require_one_input(bfile, bfile_dir, "--bfile", "--bfile-dir")
    marker_names = None
    marker_sources = None
    if database and (not manifest_name or resolve_mixed_manifests):
        marker_sources = _plink_bfile_marker_source_map(
            _plink_bfile_input_prefixes(bfile, bfile_dir, pattern)
        )
        marker_names = list(marker_sources)
    if database and resolve_mixed_manifests and not resolution_report:
        resolution_report = _default_resolution_report_path(
            output=None,
            outdir=outdir,
            output_prefix=output_prefix,
            genotypes=None,
            bfile=bfile,
            pfile=None,
            command="convert-plink",
        )
    table = _load_conversion_table(
        lookup=lookup,
        database=database,
        species=species,
        assembly=assembly,
        manifest_name=manifest_name,
        context="PLINK",
        marker_names=marker_names,
        marker_input_paths=marker_sources,
        resolve_mixed_manifests=resolve_mixed_manifests,
        on_ambiguous_marker=on_ambiguous_marker,
        resolution_report=resolution_report,
    )

    if bfile:
        _require_output_for_mode(output_prefix, "--out", "--bfile")
        stats = convert_plink_bfile(
            input_prefix=bfile,
            output_prefix=output_prefix,
            table=table,
            from_fmt=from_fmt.upper(),
            to_fmt=to_fmt.upper(),
            update_position=update_position,
            on_unconvertible_marker=on_unconvertible_marker,
            plink_command=plink_command,
        )
        _echo_plink_single_summary(stats, "PLINK fileset")
        return

    _require_output_for_mode(outdir, "--outdir", "--bfile-dir")
    try:
        batch_stats = convert_plink_bfile_batch(
            input_dir=bfile_dir,
            output_dir=outdir,
            pattern=pattern,
            suffix=suffix,
            overwrite=overwrite,
            table=table,
            from_fmt=from_fmt.upper(),
            to_fmt=to_fmt.upper(),
            update_position=update_position,
            on_unconvertible_marker=on_unconvertible_marker,
            plink_command=plink_command,
        )
    except FileExistsError as exc:
        raise click.ClickException(str(exc)) from exc
    _echo_plink_batch_summary(batch_stats, outdir, "PLINK")


@main.command("convert-pfile")
@click.option("--pfile", required=False,
              help="Input PLINK 2 prefix (expects .pgen, .pvar, and .psam)")
@click.option("--pfile-dir", required=False, type=click.Path(exists=True, file_okay=False),
              help="Directory of PLINK 2 filesets to convert")
@click.option("--pattern", default="*.pgen", show_default=True,
              help="File pattern used with --pfile-dir")
@click.option("--lookup", required=False, type=click.Path(exists=True),
              help="Lookup CSV produced by the build command")
@click.option("--database", required=False, type=click.Path(exists=True),
              help="SQLite conversion database")
@click.option("--species", required=False,
              help="Species name for --database conversion")
@click.option("--assembly", required=False,
              help="Reference assembly name for --database conversion")
@click.option("--manifest-name", required=False,
              help="Manifest name for --database conversion; inferred from input markers if omitted")
@click.option("--resolve-mixed-manifests/--no-resolve-mixed-manifests",
              default=False, show_default=True,
              help="Resolve database rules per marker for inputs containing markers from multiple manifests")
@click.option("--on-ambiguous-marker", default="fail", show_default=True,
              type=click.Choice(["fail", "skip"]),
              help="How mixed-manifest mode handles conflicting rules that cannot be resolved")
@click.option("--resolution-report", required=False,
              help="CSV report path for mixed-manifest marker rule decisions")
@click.option("--from-format", "from_fmt", required=True,
              type=click.Choice(["AB", "TOP", "FORWARD", "DESIGN", "PLUS"],
                                case_sensitive=False),
              help="Format encoding of the input .pvar allele labels")
@click.option("--to-format", "to_fmt", required=True,
              type=click.Choice(["AB", "TOP", "FORWARD", "DESIGN", "PLUS"],
                                case_sensitive=False),
              help="Target format encoding for the output .pvar allele labels")
@click.option("--out", "output_prefix", required=False,
              help="Output PLINK 2 prefix for single-fileset conversion")
@click.option("--outdir", required=False,
              help="Output directory for --pfile-dir")
@click.option("--suffix", default=".converted", show_default=True,
              help="Filename suffix for batch output prefixes")
@click.option("--overwrite/--no-overwrite", default=False, show_default=True,
              help="Allow batch conversion to replace existing output files")
@click.option("--update-position/--keep-position", default=False, show_default=True,
              help="Also replace .pvar chromosome/base-pair columns with lookup positions")
@click.option("--on-unconvertible-marker", default="exclude", show_default=True,
              type=click.Choice(["exclude", "fail", "keep"]),
              help="How to handle .pvar markers missing from lookup or lacking target alleles")
@click.option("--plink2", "plink2_command", default="plink2", show_default=True,
              help="PLINK2 executable used when --on-unconvertible-marker=exclude")
def convert_pfile_cmd(pfile, pfile_dir, pattern, lookup, database, species,
                      assembly, manifest_name, resolve_mixed_manifests,
                      on_ambiguous_marker, resolution_report, from_fmt, to_fmt,
                      output_prefix, outdir, suffix, overwrite, update_position,
                      on_unconvertible_marker, plink2_command):
    """Convert allele labels in a PLINK 2 pgen/pvar/psam fileset."""
    _require_one_input(pfile, pfile_dir, "--pfile", "--pfile-dir")
    marker_names = None
    marker_sources = None
    if database and (not manifest_name or resolve_mixed_manifests):
        marker_sources = _plink_pfile_marker_source_map(
            _plink_pfile_input_prefixes(pfile, pfile_dir, pattern)
        )
        marker_names = list(marker_sources)
    if database and resolve_mixed_manifests and not resolution_report:
        resolution_report = _default_resolution_report_path(
            output=None,
            outdir=outdir,
            output_prefix=output_prefix,
            genotypes=None,
            bfile=None,
            pfile=pfile,
            command="convert-pfile",
        )
    table = _load_conversion_table(
        lookup=lookup,
        database=database,
        species=species,
        assembly=assembly,
        manifest_name=manifest_name,
        context="PLINK 2",
        marker_names=marker_names,
        marker_input_paths=marker_sources,
        resolve_mixed_manifests=resolve_mixed_manifests,
        on_ambiguous_marker=on_ambiguous_marker,
        resolution_report=resolution_report,
    )

    if pfile:
        _require_output_for_mode(output_prefix, "--out", "--pfile")
        stats = convert_plink_pfile(
            input_prefix=pfile,
            output_prefix=output_prefix,
            table=table,
            from_fmt=from_fmt.upper(),
            to_fmt=to_fmt.upper(),
            update_position=update_position,
            on_unconvertible_marker=on_unconvertible_marker,
            plink2_command=plink2_command,
        )
        _echo_plink_single_summary(stats, "PLINK 2 fileset")
        return

    _require_output_for_mode(outdir, "--outdir", "--pfile-dir")
    try:
        batch_stats = convert_plink_pfile_batch(
            input_dir=pfile_dir,
            output_dir=outdir,
            pattern=pattern,
            suffix=suffix,
            overwrite=overwrite,
            table=table,
            from_fmt=from_fmt.upper(),
            to_fmt=to_fmt.upper(),
            update_position=update_position,
            on_unconvertible_marker=on_unconvertible_marker,
            plink2_command=plink2_command,
        )
    except FileExistsError as exc:
        raise click.ClickException(str(exc)) from exc
    _echo_plink_batch_summary(batch_stats, outdir, "PLINK 2")
