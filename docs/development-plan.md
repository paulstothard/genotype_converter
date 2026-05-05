# Development Plan

This file records repo-level design status and remaining work. It is not a
release promise. User-facing usage details belong in `README.md`,
`docs/database.md`, and `docs/genotype-formats.md`.

## Current Status

The core workflows are implemented:

- Build lookup, position, conversion, and wide-format outputs from Illumina or
  Affymetrix manifests and a reference FASTA.
- Convert genotype datasets from lookup CSVs.
- Build and inspect SQLite conversion databases from species folders containing
  manifests and reference genomes.
- Convert genotype datasets from SQLite databases.
- Convert individual files and folders for CSV wide, CSV long,
  Illumina/GSGT, Affymetrix/Axiom, PLINK 1 binary, and PLINK 2 pfile inputs.
- Exclude unconvertible markers by default and write marker-level reports.
- Resolve mixed-manifest inputs with marker-resolution reports.

The database workflow is optional. Lookup CSV workflows remain supported and
should stay supported because they are simple, inspectable, easy to archive, and
useful for validation against older outputs.

## Database Source Layout

The supported source-folder layout is:

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
    genotypes/
      optional-local-examples
```

Manifests are independent of reference genomes. For each species, the database
builder pairs each manifest with each reference assembly under that species.
Each `references/<assembly>/` folder should contain exactly one FASTA.

## Duplicate Marker Names

The same marker name can appear in multiple manifests with different conversion
rules. The database must not collapse those records into one unqualified rule.

The schema retains:

- marker name and alternate marker name
- manifest/panel name
- species and assembly/reference
- chromosome and position
- REF/ALT
- AB/TOP/FORWARD/DESIGN/PLUS/VCF conversion values
- determination type
- source file checksums
- import/build metadata

Marker names are unique only within one imported lookup source. Querying a
marker without `--manifest-name` can return multiple rows, one per source.

## Manifest Resolution

Database-backed conversion supports three modes:

- Explicit source selection with `--manifest-name`.
- Whole-input manifest inference when `--manifest-name` is omitted.
- Per-marker mixed-manifest resolution with `--resolve-mixed-manifests`.

Whole-input inference uses all input marker IDs and chooses the imported
manifest with the most matches. It fails on no-match and tied-best cases.

Mixed-manifest mode resolves one marker at a time. For an ambiguous marker, it
looks at neighboring input markers, counts which manifest provides rules for
those neighbors, and selects the clearly supported rule. If no manifest is
clearly ahead, the marker is reported as ambiguous and handled according to
`--on-ambiguous-marker`.

This heuristic should remain conservative and explainable. The
marker-resolution CSV is part of the workflow and should stay machine-readable.

## Unconvertible Marker Policy

`convert`, `convert-plink`, and `convert-pfile` all use:

```text
--on-unconvertible-marker exclude|fail|keep
```

The default is `exclude`.

- Text formats remove unconvertible marker columns or rows.
- PLINK 1 uses PLINK to remove variants before rewriting `.bim`, so
  `.bed/.bim/.fam` stay synchronized.
- PLINK 2 uses PLINK2 to remove variants before rewriting `.pvar`, so
  `.pgen/.pvar/.psam` stay synchronized.

Do not manually drop variants from PLINK metadata files without also rewriting
the matching genotype file.

## Validation Workflows

Keep helper scripts in `validation/` and use them instead of reconstructing long
commands by hand.

Important validation paths:

- `validation/scripts/` for full BovineHD old/new validation.
- `validation/mixed_manifest/scripts/` for mixed-manifest database validation.
- `validation/plink_comparison/scripts/` for collaborator PLINK comparisons.

Large validation inputs, generated databases, references, reports, and converted
outputs should remain ignored by Git.

The mixed-manifest reference downloader intentionally keeps assembled molecules,
unlocalized scaffolds, and unplaced scaffolds from NCBI RefSeq assembly reports.
That keeps validation comparable with older outputs that may include unassigned
contig placements.

## Remaining Work

- Add a short validation status note after each substantial full-panel
  validation run, including command, date, code version, reference used, and
  remaining discrepancy counts.
- Improve batch mixed-manifest reports if real use shows that the global
  marker-resolution report is not enough for multi-file review.
- Refine nearby-marker/window scoring only if real mixed-manifest data shows the
  current conservative heuristic is too simple.
- Decide whether distribution should remain GitHub/source-install only or later
  add package publishing.
- Consider generating real tiny PLINK 2 pfiles for validation if PLINK2 becomes
  reliably available in the development environment.

## Agent Notes

Keep `AGENTS.md` in the repository. It is not user-facing documentation; it is a
repo-local operating guide for AI coding agents. It records validation safety
rules, large-file handling, and probe-positioning details that are easy to lose
between sessions.
