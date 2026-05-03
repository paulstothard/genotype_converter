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
      manifests/
        <manifest-a>.csv
        <manifest-b>.csv
      references/
        ARS_UCD1_2/
          <ARS-UCD1.2-reference>.fa
        ARS_UCD_v2_0/
          <ARS-UCD-v2.0-reference>.fa
        UMD3_1/
          <UMD3.1-reference>.fa
      genotypes/
        <cattle-genotype-file-or-folder-inputs>
    sus_scrofa/
      manifests/
        <pig-manifest>.csv
      references/
        Sscrofa11_1/
          <Sscrofa11.1-reference>.fa
      genotypes/
        <pig-genotype-file-or-folder-inputs>
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

Put manifests directly under the species, independent of reference genomes. Put
reference genomes under `references/<assembly>/`, for example `ARS_UCD1_2`,
`ARS_UCD_v2_0`, `UMD3_1`, or `Sscrofa11_1`. The CLI `--species` and
`--assembly` values must exactly match these folder names.

Each `sources/<species>/` folder can contain many manifests and many reference
genomes. The database builder generates conversion information for every
manifest/reference pair within that species. For example, if `bos_taurus` has
two manifests and two reference genome folders, the build creates four lookup
sources in the SQLite database.

Each `references/<assembly>/` folder should contain exactly one matching FASTA.
A reference FASTA can contain many chromosomes or contigs. If the same assembly
is represented by substantially different FASTA files, such as soft-masked and
unmasked versions, create separate assembly folders with clear names instead of
mixing them in one folder.

## What To Add

- Put manifest CSVs for a species in `sources/<species>/manifests/`.
- Put exactly one matching reference FASTA in
  `sources/<species>/references/<assembly>/`.
- Put genotype files to convert in `sources/<species>/genotypes/`.

If a reference is supplied as multiple FASTA files, make a single combined FASTA
for this validation workspace before building the database.

## Download Reference Genomes

The helper script downloads NCBI RefSeq genome FASTA files and filters them to
assembled molecules only. That keeps chromosomes and assembled sex chromosomes
where they are part of the assembly, while excluding unplaced and unlocalized
contigs.

```bash
conda activate genotype-converter-env

python validation/mixed_manifest/download_references.py
```

The current downloads are:

| Species folder | Reference folder | NCBI assembly | RefSeq accession |
| --- | --- | --- | --- |
| `bos_taurus` | `ARS_UCD1_2` | ARS-UCD1.2 | `GCF_002263795.1` |
| `bos_taurus` | `ARS_UCD_v2_0` | ARS-UCD2.0 | `GCF_002263795.3` |
| `bos_taurus` | `UMD3_1` | Bos_taurus_UMD_3.1 | `GCF_000003055.4` |
| `sus_scrofa` | `Sscrofa11_1` | Sscrofa11.1 | `GCF_000003025.6` |

The downloaded FASTA and assembly report files are ignored by Git.

## Build The SQLite Database

From the repository root:

```bash
validation/mixed_manifest/run_discover_sources.sh

validation/mixed_manifest/run_database_build.sh --yes
```

This discovers each species-level manifest set and pairs it with every reference
assembly under that species. Use `--workers 1` for large livestock references
unless the machine has enough memory for one minimap2 reference index per
worker. The build script requires `--yes` because the mixed-manifest workload
can be large. Use `validation/mixed_manifest/run_database_build.sh --help` for
overrides.

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
  --genotypes-dir validation/mixed_manifest/sources/bos_taurus/genotypes \
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
  --bfile-dir validation/mixed_manifest/sources/bos_taurus/genotypes \
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
  --pfile-dir validation/mixed_manifest/sources/bos_taurus/genotypes \
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
