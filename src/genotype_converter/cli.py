from __future__ import annotations

import click

from .pipeline import run
from .convert_genotypes import load_lookup_table, convert_wide, convert_long


@click.group()
def main():
    """Build genotype conversion files and convert genotype data between formats."""


@main.command("build")
@click.option("--manifest", required=True, type=click.Path(exists=True),
              help="Illumina or Affymetrix manifest CSV")
@click.option("--reference", required=True, type=click.Path(exists=True),
              help="Reference genome FASTA")
@click.option("--outdir", default="output", show_default=True,
              help="Output directory")
@click.option("--species", default="all", show_default=True,
              help="Species name for output subdirectory")
@click.option("--workers", default=0, show_default=True,
              help="Worker processes (0 = use all cores)")
@click.option("--align/--no-align", default=False, show_default=True,
              help="Write alignment.txt debug file")
@click.option("--parquet/--no-parquet", default=False, show_default=True,
              help="Also write lookup.parquet (requires pyarrow)")
def build_cmd(manifest, reference, outdir, species, workers, align, parquet):
    """Align manifest variants and build conversion/position files."""
    stats = run(
        manifest_path=manifest,
        reference_path=reference,
        outdir=outdir,
        species=species,
        workers=workers,
        save_alignment=align,
        save_parquet=parquet,
    )
    click.echo(
        f"Build complete: {stats.total_markers} markers "
        f"({stats.n_snp} SNPs, {stats.n_indel} indels), "
        f"{stats.n_positioned} positioned, {stats.n_not_positioned} not positioned."
    )
    summary_path = next((f for f in stats.output_files if f.endswith(".summary.txt")), None)
    if summary_path:
        click.echo(f"Full summary: {summary_path}")


@main.command("convert")
@click.option("--genotypes", required=True, type=click.Path(exists=True),
              help="Input genotype file (CSV)")
@click.option("--lookup", required=True, type=click.Path(exists=True),
              help="Lookup CSV produced by the build command (lookup.csv or conversion.csv)")
@click.option("--from-format", "from_fmt", required=True,
              type=click.Choice(["AB", "TOP", "FORWARD", "DESIGN", "PLUS", "VCF"],
                                case_sensitive=False),
              help="Format encoding of the input genotypes")
@click.option("--to-format", "to_fmt", required=True,
              type=click.Choice(["AB", "TOP", "FORWARD", "DESIGN", "PLUS", "VCF"],
                                case_sensitive=False),
              help="Target format encoding for the output")
@click.option("--output", required=True,
              help="Output file path")
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
def convert_cmd(genotypes, lookup, from_fmt, to_fmt, output, layout, in_sep,
                out_sep, sample_col, marker_col, genotype_col):
    """Convert genotypes between format encodings (e.g. TOP → PLUS).

    \b
    Examples:
      # Wide format: convert TOP-encoded genotypes to PLUS
      genotype-converter convert \\
        --genotypes mydata.csv --lookup manifest.reference.lookup.csv \\
        --from-format TOP --to-format PLUS --output converted.csv

      # Long format
      genotype-converter convert \\
        --genotypes mydata_long.csv --lookup manifest.reference.lookup.csv \\
        --from-format FORWARD --to-format VCF --layout long --output converted.csv
    """
    from_fmt = from_fmt.upper()
    to_fmt = to_fmt.upper()
    table = load_lookup_table(lookup)

    if layout == "wide":
        convert_wide(
            input_path=genotypes,
            output_path=output,
            table=table,
            from_fmt=from_fmt,
            to_fmt=to_fmt,
            in_sep=in_sep,
            out_sep=out_sep,
            sample_col=sample_col,
        )
    else:
        convert_long(
            input_path=genotypes,
            output_path=output,
            table=table,
            from_fmt=from_fmt,
            to_fmt=to_fmt,
            in_sep=in_sep,
            out_sep=out_sep,
            sample_col=sample_col,
            marker_col=marker_col,
            genotype_col=genotype_col,
        )
    click.echo(f"Converted genotypes written to {output}")
