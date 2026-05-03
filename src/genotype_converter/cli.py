from __future__ import annotations

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
    discover_source_folders,
    import_lookup,
    init_database,
    list_assemblies,
    list_manifests,
    list_species,
    load_lookup_table_from_database,
    query_marker,
    rows_to_csv,
    rows_to_json,
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


def _echo_plink_single_summary(stats, label: str) -> None:
    click.echo(
        f"Converted {label}: {stats.variants_converted}/{stats.variants_total} "
        f"variants, {stats.alleles_changed} allele labels changed."
    )
    if stats.variants_missing_lookup:
        click.echo(f"Markers missing from lookup: {stats.variants_missing_lookup}")
    click.echo(f"Wrote {stats.genotype_path}")
    click.echo(f"Wrote {stats.variant_path}")
    click.echo(f"Wrote {stats.sample_path}")


def _echo_plink_batch_summary(stats, outdir, label: str) -> None:
    click.echo(
        f"Converted {stats.filesets_converted}/{stats.filesets_total} "
        f"{label} filesets to {Path(outdir)}"
    )
    for item in stats.stats:
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
        click.echo("\t".join(str(row.get(field, "") or "") for field in fieldnames))


def _manifest_name_from_path(path: str) -> str:
    return Path(path).stem.replace(".", "_")


def _load_conversion_table(lookup, database, species, assembly, manifest_name, context: str):
    _require_one_input(lookup, database, "--lookup", "--database")
    if lookup:
        return load_lookup_table(lookup)
    missing = [
        flag for flag, value in [
            ("--species", species),
            ("--assembly", assembly),
            ("--manifest-name", manifest_name),
        ] if not value
    ]
    if missing:
        raise click.UsageError(
            f"{', '.join(missing)} required with --database for {context} conversion."
        )
    try:
        return load_lookup_table_from_database(
            database_path=database,
            species=species,
            assembly=assembly,
            manifest_name=manifest_name,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc


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
def db_import_lookup_cmd(database, lookup, species, assembly, manifest_name,
                         manifest_path, reference_name, reference_path,
                         tool_version, notes, replace):
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
    )
    click.echo(
        f"Imported {stats.rows_imported} marker rules from {lookup} "
        f"as source {stats.source_id}."
    )
    if stats.rows_replaced:
        click.echo("Replaced existing source with the same lookup checksum.")


@db_cmd.command("build")
@click.option("--source-root", required=True, type=click.Path(exists=True, file_okay=False),
              help="Root folder organized as species/assembly/manifests and references")
@click.option("--database", required=True, help="SQLite database path")
@click.option("--build-outdir", default="database_build", show_default=True,
              help="Directory for generated build outputs")
@click.option("--workers", default=1, show_default=True,
              help="Worker processes for each build. Use 1 for large references unless memory is available.")
@click.option("--replace/--no-replace", default=False, show_default=True,
              help="Replace existing imports with the same species, assembly, manifest, and lookup checksum")
@click.option("--progress/--no-progress", default=True, show_default=True,
              help="Show build progress")
def db_build_cmd(source_root, database, build_outdir, workers, replace, progress):
    """Build lookup files from source folders and import them into SQLite."""
    folders = discover_source_folders(source_root)
    if not folders:
        raise click.ClickException(f"No source folders found under {source_root}")

    total_imported = 0
    total_rules = 0
    for folder in folders:
        if not folder.manifest_paths:
            raise click.ClickException(
                f"No manifest files found in {Path(folder.root_path) / 'manifests'}"
            )
        if len(folder.reference_paths) != 1:
            raise click.ClickException(
                f"Expected exactly one reference file in "
                f"{Path(folder.root_path) / 'references'}, found {len(folder.reference_paths)}"
            )
        reference_path = folder.reference_paths[0]
        for manifest_path in folder.manifest_paths:
            manifest_name = _manifest_name_from_path(manifest_path)
            click.echo(
                f"Building {folder.species}/{folder.assembly}/{manifest_name}"
            )
            build_stats = run(
                manifest_path=manifest_path,
                reference_path=reference_path,
                outdir=build_outdir,
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
            )
            total_imported += 1
            total_rules += import_stats.rows_imported
            click.echo(
                f"Imported {import_stats.rows_imported} marker rules from {lookup_path}"
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
              help="Root folder organized as species/assembly/manifests and references")
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
            "manifest_count": len(folder.manifest_paths),
            "reference_count": len(folder.reference_paths),
            "manifest_paths": ";".join(folder.manifest_paths),
            "reference_paths": ";".join(folder.reference_paths),
        }
        for folder in folders
    ]
    fieldnames = [
        "species", "assembly", "root_path", "manifest_count", "reference_count",
        "manifest_paths", "reference_paths",
    ]
    _echo_rows(rows, fieldnames, output_format)


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
              help="Manifest name for --database CSV conversion")
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
              type=click.Choice(["wide", "long"], case_sensitive=False),
              help="Input file layout: wide (samples × markers) or long (one row per sample×marker)")
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
def convert_cmd(genotypes, genotypes_dir, pattern, lookup, database, species,
                assembly, manifest_name, from_fmt, to_fmt, output, outdir,
                suffix, overwrite, layout, in_sep, out_sep, sample_col,
                marker_col, genotype_col):
    """Convert genotypes between format encodings (e.g. TOP → PLUS).

    Examples:

      genotype-converter convert --genotypes mydata.csv --lookup lookup.csv
      --from-format TOP --to-format PLUS --output converted.csv

      genotype-converter convert --genotypes-dir genotypes/ --lookup lookup.csv
      --from-format TOP --to-format PLUS --outdir converted/
    """
    from_fmt = from_fmt.upper()
    to_fmt = to_fmt.upper()
    table = _load_conversion_table(
        lookup=lookup,
        database=database,
        species=species,
        assembly=assembly,
        manifest_name=manifest_name,
        context="CSV",
    )

    _require_one_input(genotypes, genotypes_dir, "--genotypes", "--genotypes-dir")

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
        )
        click.echo(
            f"Converted genotypes written to {output}: "
            f"{stats.genotypes_changed}/{stats.genotype_cells_total} genotype cells changed, "
            f"{stats.alleles_changed} allele labels changed."
        )
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
              help="Manifest name for --database conversion")
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
@click.option("--require-all-markers/--allow-missing-markers", default=True, show_default=True,
              help="Fail if any .bim marker is absent from the lookup table")
def convert_plink_cmd(bfile, bfile_dir, pattern, lookup, database, species,
                      assembly, manifest_name, from_fmt, to_fmt, output_prefix,
                      outdir, suffix, overwrite, update_position,
                      require_all_markers):
    """Convert allele labels in a PLINK bed/bim/fam fileset."""
    table = _load_conversion_table(
        lookup=lookup,
        database=database,
        species=species,
        assembly=assembly,
        manifest_name=manifest_name,
        context="PLINK",
    )
    _require_one_input(bfile, bfile_dir, "--bfile", "--bfile-dir")

    if bfile:
        _require_output_for_mode(output_prefix, "--out", "--bfile")
        stats = convert_plink_bfile(
            input_prefix=bfile,
            output_prefix=output_prefix,
            table=table,
            from_fmt=from_fmt.upper(),
            to_fmt=to_fmt.upper(),
            update_position=update_position,
            require_all_markers=require_all_markers,
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
            require_all_markers=require_all_markers,
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
              help="Manifest name for --database conversion")
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
@click.option("--require-all-markers/--allow-missing-markers", default=True, show_default=True,
              help="Fail if any .pvar marker is absent from the lookup table")
def convert_pfile_cmd(pfile, pfile_dir, pattern, lookup, database, species,
                      assembly, manifest_name, from_fmt, to_fmt, output_prefix,
                      outdir, suffix, overwrite, update_position,
                      require_all_markers):
    """Convert allele labels in a PLINK 2 pgen/pvar/psam fileset."""
    table = _load_conversion_table(
        lookup=lookup,
        database=database,
        species=species,
        assembly=assembly,
        manifest_name=manifest_name,
        context="PLINK 2",
    )
    _require_one_input(pfile, pfile_dir, "--pfile", "--pfile-dir")

    if pfile:
        _require_output_for_mode(output_prefix, "--out", "--pfile")
        stats = convert_plink_pfile(
            input_prefix=pfile,
            output_prefix=output_prefix,
            table=table,
            from_fmt=from_fmt.upper(),
            to_fmt=to_fmt.upper(),
            update_position=update_position,
            require_all_markers=require_all_markers,
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
            require_all_markers=require_all_markers,
        )
    except FileExistsError as exc:
        raise click.ClickException(str(exc)) from exc
    _echo_plink_batch_summary(batch_stats, outdir, "PLINK 2")
