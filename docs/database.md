# SQLite Conversion Database

The SQLite database is optional. The normal lookup-file workflow remains
supported and does not require a database.

Current database support covers import, inspection, source-folder discovery, and
CSV genotype conversion. PLINK conversion still uses lookup CSV files.

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

## CSV Conversion

CSV genotype conversion can use a database instead of `--lookup`:

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

For now, `--manifest-name` is required. This avoids guessing when the same
marker name appears in more than one manifest. Manifest inference from
neighboring markers is planned but not implemented.

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
