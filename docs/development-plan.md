# Development Plan

This file records larger design ideas that should be preserved while the code
evolves. It is not a release promise; items here still need detailed design,
tests, and migration planning before implementation.

## SQLite Conversion Database

Add an optional SQLite-backed conversion database that can store SNP and indel
conversion information across many manifests, species, and assemblies.

### Goals

- Let users convert genotype files by marker name without manually selecting a
  single lookup CSV for every run.
- Support multiple species and reference assemblies.
- Support multiple manifests per species, with conversion records generated for
  each manifest/reference pair.
- Preserve marker provenance so the same marker name can appear in more than
  one manifest without losing which rule came from which source.
- Handle SNPs and indels consistently.
- Allow conversion of one genotype file or whole folders of genotype files.

### Proposed Workflow

Use a structured data folder such as:

```text
database_sources/
  bos_taurus/
    manifests/
      bovinehd-manifest-b.csv
      other-panel.csv
    references/
      ARS_UCD_v2_0/
        ARS_UCD_v2.0.fa
      ARS_UCD1_2/
        ARS-UCD1.2.fa
```

The database builder would scan this folder, run or reuse `build` outputs for
each manifest/reference pair within a species, and load the resulting
conversion records into an SQLite database.

Possible commands:

```bash
genotype-converter db init --database genotype_converter.sqlite

genotype-converter db build \
  --source-root database_sources \
  --database genotype_converter.sqlite \
  --build-outdir database_build \
  --workers 1

genotype-converter db import-lookup \
  --database genotype_converter.sqlite \
  --lookup manifest.reference.lookup.csv \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --manifest-name bovinehd_manifest_b

genotype-converter db list-species --database genotype_converter.sqlite
genotype-converter db list-assemblies --database genotype_converter.sqlite --species bos_taurus
genotype-converter db list-manifests --database genotype_converter.sqlite --species bos_taurus

genotype-converter convert \
  --genotypes input.csv \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --from-format TOP \
  --to-format PLUS \
  --output converted.csv
```

The exact command names can change, but the database should have clear
maintenance operations: initialize, build/update, inspect/list, validate, and
possibly remove stale manifest entries.

Stage 1 is implemented for existing lookup files:

- initialize a SQLite database
- import a built lookup CSV with species, assembly, and manifest provenance
- list species, assemblies, and imported manifests
- query rules for a marker
- export marker query results as table, CSV, or JSON
- discover source-folder contents without running build
- build and import source folders with species-level manifests and one reference
  FASTA per `references/<assembly>/` folder

Stage 2 has database conversion:

- `genotype-converter convert` can use `--database` instead of `--lookup` for
  single CSV files and CSV folders
- `genotype-converter convert-plink` can use `--database` instead of `--lookup`
  for single PLINK 1 filesets and PLINK 1 folders
- `genotype-converter convert-pfile` can use `--database` instead of `--lookup`
  for single PLINK 2 filesets and PLINK 2 folders
- `--species` and `--assembly` are required in database mode
- `--manifest-name` can be supplied explicitly, or omitted to use conservative
  marker-count manifest inference

Current inference uses all input marker IDs and chooses the imported manifest
with the most matches. It fails on no-match and tied-best cases. More nuanced
neighbor-window inference from input marker order is available through
`--resolve-mixed-manifests`.

### Duplicate Marker Names

The same marker name can appear in multiple manifests with different conversion
rules. The database should not collapse those records into one unqualified rule.
The current schema retains:

- marker name
- manifest/panel name
- species
- assembly/reference
- chromosome and position
- REF/ALT
- AB/TOP/FORWARD/DESIGN/PLUS/VCF conversion values
- determination type
- source file checksums
- build timestamp and tool version

Marker names are unique only within one imported lookup source. Querying a
marker without `--manifest-name` may return multiple rows, one per manifest.

When a genotype file contains duplicate-rule markers and the user has not
specified a manifest, mixed-manifest mode can infer the most likely manifest
context from nearby markers in the input order. The implemented heuristic:

1. For each ambiguous marker, look at a window of neighboring input markers.
2. Count which manifest provides rules for the most neighboring markers.
3. Prefer the rule from that manifest if it is clearly ahead.
4. If no manifest is clearly ahead, report the marker as ambiguous and leave it
   unconverted or require an explicit user choice.

This heuristic is conservative and explainable. The converter reports which rule
was selected and why in a marker-resolution CSV, especially when multiple
manifest rules exist.

### Should Existing Outputs Stay?

Yes. Keep the current CSV/parquet outputs.

Reasons:

- They are simple, inspectable, and easy to archive with publications or
  analyses.
- They are useful for debugging and validation against older pipelines.
- They let users run `convert` without adopting a database workflow.
- They provide a stable interchange format that can also be imported into the
  database.

The database should be optional, not required. The existing lookup-file workflow
should remain supported:

```bash
genotype-converter convert --lookup manifest.reference.lookup.csv ...
```

Database-backed conversion should be an additional path for users who manage
many species, assemblies, manifests, or genotype batches.

### Batch Conversion

Batch conversion is now supported for CSV wide, CSV long, PLINK 1 binary
filesets, and PLINK 2 filesets. Batch conversion can use either a lookup CSV or
an explicit database context.

Lookup-file command shape:

```bash
genotype-converter convert \
  --genotypes-dir genotypes/ \
  --pattern "*.csv" \
  --lookup manifest.reference.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --outdir converted/
```

```bash
genotype-converter convert-plink \
  --bfile-dir plink_files/ \
  --pattern "*.bed" \
  --lookup manifest.reference.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --outdir converted/
```

```bash
genotype-converter convert-pfile \
  --pfile-dir pfiles/ \
  --pattern "*.pgen" \
  --lookup manifest.reference.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --outdir converted/
```

Database-backed command shape:

```bash
genotype-converter convert \
  --genotypes-dir genotypes/ \
  --pattern "*.csv" \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --manifest-name bovinehd_manifest_b \
  --from-format TOP \
  --to-format PLUS \
  --outdir converted/
```

```bash
genotype-converter convert-plink \
  --bfile-dir plink_files/ \
  --pattern "*.bed" \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --manifest-name bovinehd_manifest_b \
  --from-format TOP \
  --to-format PLUS \
  --outdir converted/
```

```bash
genotype-converter convert-pfile \
  --pfile-dir pfiles/ \
  --pattern "*.pgen" \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --manifest-name bovinehd_manifest_b \
  --from-format TOP \
  --to-format PLUS \
  --outdir converted/
```

Batch conversion does:

- preserve filenames or use a predictable suffix
- avoid overwriting outputs unless explicitly requested
- write `conversion_summary.csv` for each batch

Batch conversion still needs to:

- report missing, ambiguous, and unconverted markers per file
- refine mixed-manifest reports so batch outputs can also include per-file
  resolution summaries

### Post-v0.1.1 Planning Notes

The v0.1.1 release covers the core database conversion workflow:

- database initialization, source-folder build/import, lookup import, and query
- CSV, PLINK 1, and PLINK 2 conversion from lookup CSVs or SQLite
- whole-input manifest inference when `--manifest-name` is omitted
- per-marker mixed-manifest resolution with marker-resolution reports
- ambiguity handling with fail or skip behavior
- a runnable mixed-manifest example under `examples/mixed_manifests/`

Future work to preserve:

- Add database maintenance commands for removing stale sources, replacing by
  manifest/source context, and validating database contents.
- Improve batch reports with per-file mixed-manifest summaries in addition to
  the global marker-resolution report.
- Record a short bovine validation status note after the next full validation
  run, including command, date, code version, and remaining discrepancy counts.
- Decide whether distribution remains GitHub/source-install only or should add
  package publishing.
- Refine nearby-marker/window scoring only if real mixed-manifest data shows
  the current conservative heuristic is too simple.

### Database Test Fixture

The repository includes a tiny source-folder fixture at
`tests/data/database_sources/`. It is intended for future SQLite database tests
and mirrors the proposed user-facing folder organization:

```text
tests/data/database_sources/
  bos_taurus/
    manifests/
    references/
      ARS_UCD_v2_0/
      ARS_UCD1_2/
    genotypes/
    expected/
```

Future database code should first prove that it can discover and build from this
small fixture. Do not use the full bovine validation panel as the first database
test target.

### Design Questions To Resolve

- Exact SQLite schema and indexes.
- Whether to store all old-style output rows directly or normalize into marker,
  manifest, assembly, and conversion-rule tables.
- How to version database records when the build algorithm changes.
- How users should choose species/assembly interactively versus by command-line
  flags.
- How to handle marker aliases and alternate marker names.
- How to expose ambiguity reports in a machine-readable form.
- Whether database building should require reference FASTA files every time or
  can import already-built lookup files.
