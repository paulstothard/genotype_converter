# Mixed-Manifest Validation Workspace

Use this folder to test database-backed conversion when a large genotype input
contains markers from multiple bovine manifests.

Only this README and `.gitkeep` placeholder files should be committed. Put real
manifests, reference FASTA files, genotype inputs, SQLite databases, build
outputs, converted genotypes, and reports here locally; they are intentionally
ignored.

## Folder Layout

```text
validation/mixed_manifest/
  sources/
    bos_taurus/
      ARS_UCD_v2_0/
        manifests/
          <manifest-a>.csv
          <manifest-b>.csv
        references/          # required directory name; put one FASTA here
          <reference>.fa
        genotypes/
          <genotype-file-or-folder-inputs>
  database_build/
  converted/
  reports/
```

Provenance-preserving filenames are preferred. Do not rename files just to make
them generic.

## What To Add

- Put all bovine manifest CSVs for the test in
  `sources/bos_taurus/ARS_UCD_v2_0/manifests/`.
- Put exactly one matching reference FASTA in
  `sources/bos_taurus/ARS_UCD_v2_0/references/`. The directory name is plural
  because it is part of the database source-folder convention, but each
  species/assembly folder should contain one reference file for a build.
- Put genotype files to convert in
  `sources/bos_taurus/ARS_UCD_v2_0/genotypes/`.

If you need a different assembly, create a sibling assembly folder under
`sources/bos_taurus/` using the same internal layout.

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

Use `--workers 1` for large bovine references unless the machine has enough
memory for one minimap2 reference index per worker.

## Inspect Imported Manifests

```bash
genotype-converter db list-manifests \
  --database validation/mixed_manifest/mixed_manifest.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0
```

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
