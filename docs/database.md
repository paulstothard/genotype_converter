# SQLite Conversion Database

The SQLite database is optional. The normal lookup-file workflow remains
supported and does not require a database.

Current database support covers import, inspection, source-folder discovery, CSV
genotype conversion, PLINK 1 conversion, PLINK 2 conversion, and batch
conversion for those same genotype layouts.

## Install Notes

No extra dependency is required. The database layer uses Python's standard
`sqlite3` module.

## Schema

### `lookup_sources`

One row represents one imported lookup CSV in a specific context:

- species
- assembly
- manifest name
- manifest path and SHA-256, when provided
- reference name, path, and SHA-256, when provided
- lookup path and SHA-256
- import timestamp
- optional tool version and notes

The manifest name is required because the same marker name can appear in
multiple manifests.

### `marker_rules`

One row represents one marker rule from one imported lookup source. It stores
position, REF/ALT, determination type, and allele values for AB, TOP, FORWARD,
DESIGN, PLUS, and VCF notation.

Marker names are unique only within one `lookup_sources` row. Querying a marker
without a manifest filter can return multiple rules.

## Commands

```bash
genotype-converter db init --database genotype_converter.sqlite
```

```bash
genotype-converter db import-lookup \
  --database genotype_converter.sqlite \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --manifest-name bovinehd_manifest_b
```

```bash
genotype-converter db marker \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --marker SNP2 \
  --format table
```

`db marker` supports `--format table`, `--format csv`, and `--format json`.

## Conversion

Genotype conversion can use a database instead of `--lookup` when the command
supports database mode:

```bash
genotype-converter convert \
  --genotypes mydata.csv \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --manifest-name bovinehd_manifest_b \
  --from-format TOP \
  --to-format PLUS \
  --output converted.csv
```

```bash
genotype-converter convert-plink \
  --bfile mydata_top \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --manifest-name bovinehd_manifest_b \
  --from-format TOP \
  --to-format PLUS \
  --out mydata_plus
```

```bash
genotype-converter convert-pfile \
  --pfile mydata_top \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --manifest-name bovinehd_manifest_b \
  --from-format TOP \
  --to-format PLUS \
  --out mydata_plus
```

Folder conversion uses the same database context:

```bash
genotype-converter convert \
  --genotypes-dir genotype_files \
  --pattern "*.csv" \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --manifest-name bovinehd_manifest_b \
  --from-format TOP \
  --to-format PLUS \
  --layout wide \
  --outdir converted_genotype_files
```

```bash
genotype-converter convert-plink \
  --bfile-dir plink_files \
  --pattern "*.bed" \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --manifest-name bovinehd_manifest_b \
  --from-format TOP \
  --to-format PLUS \
  --outdir converted_plink_files
```

```bash
genotype-converter convert-pfile \
  --pfile-dir pfiles \
  --pattern "*.pgen" \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --manifest-name bovinehd_manifest_b \
  --from-format TOP \
  --to-format PLUS \
  --outdir converted_pfiles
```

If `--manifest-name` is omitted, the converter tries to infer the manifest from
the input marker IDs. It counts which imported lookup source for the requested
species and assembly matches the most input markers. Inference is conservative:
the command fails if no imported source matches or if two sources tie for the
best match. In those cases, rerun with an explicit `--manifest-name`.

Example inferred conversion:

```bash
genotype-converter convert \
  --genotypes mydata.csv \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --from-format TOP \
  --to-format PLUS \
  --output converted.csv
```

## Source Folder Discovery

Source folders can be inspected without running a build:

```bash
genotype-converter db discover-sources \
  --source-root database_sources \
  --format table
```

Expected folder shape:

```text
database_sources/
  bos_taurus/
    ARS_UCD_v2_0/
      manifests/
      references/
```

Discovery reports species, assembly, manifest files, and reference files. It
does not run alignment or import anything.

## Source Folder Build

The same folder layout can be built and imported:

```bash
genotype-converter db build \
  --source-root database_sources \
  --database genotype_converter.sqlite \
  --build-outdir database_build \
  --workers 1
```

Each assembly folder must contain one reference file in `references/` and one or
more manifest files in `manifests/`. The command runs `build` once for each
manifest, writes the normal build outputs under `--build-outdir`, and imports
the generated lookup CSV into SQLite.

Use `--workers 1` for large references unless the machine has enough memory for
additional minimap2 indexes.
