# Mixed-Manifest Validation Workspace

Use this folder to test database-backed conversion when genotype inputs contain
markers from multiple manifests, assemblies, or species.

Only this README and `.gitkeep` placeholder files should be committed. Put real
manifests, reference FASTA files, genotype inputs, SQLite databases, build
outputs, converted genotypes, and reports here locally; they are intentionally
ignored.

## Folder Layout

```text
validation/mixed_manifest/
  sources/
    bos_taurus/
      ARS_UCD1_2/
        manifests/
          <manifest-a>.csv
        references/          # required directory name; put one FASTA here
          <ARS-UCD1.2-reference>.fa
        genotypes/
          <genotype-file-or-folder-inputs>
      ARS_UCD_v2_0/
        manifests/
          <manifest-a>.csv
          <manifest-b>.csv
        references/          # required directory name; put one FASTA here
          <reference>.fa
        genotypes/
          <genotype-file-or-folder-inputs>
    sus_scrofa/
      Sscrofa11_1/
        manifests/
          <pig-manifest>.csv
        references/          # required directory name; put one FASTA here
          <Sscrofa11.1-reference>.fa
        genotypes/
          <genotype-file-or-folder-inputs>
  database_build/
  converted/
  reports/
```

Provenance-preserving filenames are preferred. Do not rename files just to make
them generic.

## Naming Conventions

Use stable scientific species labels for folder names:

- cattle or bovine: `bos_taurus`
- pig or porcine: `sus_scrofa`

Use the reference genome or assembly name for the assembly folder, for example
`ARS_UCD_v2_0`, `ARS_UCD1_2`, or `Sscrofa11_1`. The CLI `--species` and
`--assembly` values must exactly match these folder names.

Each `sources/<species>/<assembly>/` folder is one build target. It can contain
many manifests, but it should contain exactly one matching reference FASTA in
`references/`. A reference FASTA can contain many chromosomes or contigs. If the
same assembly is represented by substantially different FASTA files, such as
soft-masked and unmasked versions, create separate assembly folders with clear
names instead of mixing them in one folder.

## What To Add

- Put manifest CSVs for an assembly in
  `sources/<species>/<assembly>/manifests/`.
- Put exactly one matching reference FASTA in
  `sources/<species>/<assembly>/references/`. The directory name is plural
  because it is part of the database source-folder convention.
- Put genotype files to convert in
  `sources/<species>/<assembly>/genotypes/`.

If a reference is supplied as multiple FASTA files, make a single combined FASTA
for this validation workspace before building the database.

## Build The SQLite Database

From the repository root:

```bash
conda activate genotype-converter-env

genotype-converter db build \
  --source-root validation/mixed_manifest/sources \
  --database validation/mixed_manifest/mixed_manifest.sqlite \
  --build-outdir validation/mixed_manifest/database_build \
  --workers 1 \
  --progress
```

This discovers every complete `sources/<species>/<assembly>/` folder under the
source root. Use `--workers 1` for large livestock references unless the machine
has enough memory for one minimap2 reference index per worker.

## Inspect Imported Manifests

```bash
genotype-converter db list-manifests \
  --database validation/mixed_manifest/mixed_manifest.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0
```

Change `--species` and `--assembly` to inspect another folder, such as
`--species sus_scrofa --assembly Sscrofa11_1`.

## Convert A Folder Of CSV Genotypes

```bash
genotype-converter convert \
  --genotypes-dir validation/mixed_manifest/sources/bos_taurus/ARS_UCD_v2_0/genotypes \
  --pattern "*.csv" \
  --database validation/mixed_manifest/mixed_manifest.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --resolve-mixed-manifests \
  --on-ambiguous-marker skip \
  --from-format TOP \
  --to-format PLUS \
  --layout wide \
  --outdir validation/mixed_manifest/converted \
  --resolution-report validation/mixed_manifest/reports/manifest_resolution.csv
```

One conversion command targets one species and assembly. The SQLite database can
hold multiple species and assemblies, but genotype conversion should be run
separately for cattle and pig inputs unless a future workflow explicitly adds
cross-species dispatch.

Use `--on-ambiguous-marker fail` when you want the command to stop on unresolved
conflicting marker rules.

## Convert PLINK Files

For PLINK 1 binary filesets:

```bash
genotype-converter convert-plink \
  --bfile-dir validation/mixed_manifest/sources/bos_taurus/ARS_UCD_v2_0/genotypes \
  --pattern "*.bed" \
  --database validation/mixed_manifest/mixed_manifest.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --resolve-mixed-manifests \
  --on-ambiguous-marker skip \
  --from-format TOP \
  --to-format PLUS \
  --outdir validation/mixed_manifest/converted \
  --resolution-report validation/mixed_manifest/reports/manifest_resolution.csv
```

For PLINK 2 p-files:

```bash
genotype-converter convert-pfile \
  --pfile-dir validation/mixed_manifest/sources/bos_taurus/ARS_UCD_v2_0/genotypes \
  --pattern "*.pgen" \
  --database validation/mixed_manifest/mixed_manifest.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --resolve-mixed-manifests \
  --on-ambiguous-marker skip \
  --from-format TOP \
  --to-format PLUS \
  --outdir validation/mixed_manifest/converted \
  --resolution-report validation/mixed_manifest/reports/manifest_resolution.csv
```

## What To Review

Start with:

- `reports/manifest_resolution.csv`
- `converted/conversion_summary.csv`
- any markers reported as `ambiguous` or `missing`
- markers where `candidate_count` is greater than 1
- markers with unexpected `selected_manifest_name`

Record conclusions in `reports/` with the command, date, code version, input
filenames, and a short explanation of any unresolved markers.
