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
- Support multiple manifests per species/assembly.
- Preserve marker provenance so the same marker name can appear in more than
  one manifest without losing which rule came from which source.
- Handle SNPs and indels consistently.
- Allow conversion of one genotype file or whole folders of genotype files.

### Proposed Workflow

Use a structured data folder such as:

```text
database_sources/
  bos_taurus/
    ARS_UCD_v2_0/
      references/
        ARS_UCD_v2.0.fa
      manifests/
        bovinehd-manifest-b.csv
        other-panel.csv
    ARS_UCD1_2/
      references/
        ARS-UCD1.2.fa
      manifests/
        ...
```

The database builder would scan this folder, run or reuse `build` outputs for
each manifest/reference pair, and load the resulting conversion records into an
SQLite database.

Possible commands:

```bash
genotype-converter db init --database genotype_converter.sqlite

genotype-converter db build \
  --source-root database_sources \
  --database genotype_converter.sqlite

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

### Duplicate Marker Names

The same marker name can appear in multiple manifests with different conversion
rules. The database should not collapse those records into one unqualified rule.
It should retain at least:

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

When a genotype file contains duplicate-rule markers and the user has not
specified a manifest, the converter could infer the most likely manifest context
from nearby markers in the input order. A practical heuristic:

1. For each ambiguous marker, look at a window of neighboring input markers.
2. Count which manifest provides rules for the most neighboring markers.
3. Prefer the rule from that manifest if it is clearly ahead.
4. If no manifest is clearly ahead, report the marker as ambiguous and leave it
   unconverted or require an explicit user choice.

This heuristic must be conservative and explainable. The converter should report
which rule was selected and why, especially when multiple manifest rules exist.

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

Add support for converting whole folders of genotype files.

Possible command shape:

```bash
genotype-converter convert \
  --genotypes-dir genotypes/ \
  --pattern "*.csv" \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --from-format TOP \
  --to-format PLUS \
  --outdir converted/
```

Batch conversion should:

- preserve filenames or use a predictable suffix
- produce one summary report for the batch
- report missing, ambiguous, and unconverted markers per file
- support lookup-file mode and database mode
- avoid overwriting outputs unless explicitly requested

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

