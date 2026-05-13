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
  scripts/
    run_make_example_genotypes.sh
    run_discover_sources.sh
    run_database_build.sh
    run_example_conversions.sh
    run_check_example_conversions.sh
    make_example_genotypes.py
    download_references.py
    check_example_conversions.py
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

## Generate Small Genotype Examples

The helper script scans the local manifest folders and writes small AB-coded
genotype examples for each species:

```bash
validation/mixed_manifest/scripts/run_make_example_genotypes.sh
```

It writes files under:

```text
sources/<species>/genotypes/synthetic_mixed_manifest/
  marker_selection.csv
  mixed_manifest_wide_ab.csv
  mixed_manifest_long_ab.csv
  illumina_gsgt_matrix_ab.txt
  illumina_gsgt_long_multiformat.txt
  affymetrix_axiom_dual_call_matrix.txt
  mixed_manifest_plink1_ab.bed
  mixed_manifest_plink1_ab.bim
  mixed_manifest_plink1_ab.fam
  mixed_manifest_plink2_ab.pgen
  mixed_manifest_plink2_ab.pvar
  mixed_manifest_plink2_ab.psam
```

`marker_selection.csv` records why each marker was selected, including whether
it is unique to one manifest or shared across manifests. The PLINK examples are
for this converter's metadata rewrite tests: `.bim` and `.pvar` contain real
selected marker IDs. The `.bed` file is a valid small PLINK 1 genotype matrix
with missing calls, while `.pgen` is a placeholder used only for metadata
rewrite paths that do not require PLINK 2 to filter variants.

The Illumina/GSGT and Affymetrix/Axiom files are real vendor genotype formats
seen in the older `snp_conversion` test fixtures. They are generated here to
exercise direct parser support alongside the CSV, PLINK 1, and PLINK 2 inputs
above.

## Download Reference Genomes

The user-facing `genotype-converter reference download-ncbi` command can
download one NCBI RefSeq genome FASTA at a time into the source-folder layout.
This validation workspace also keeps a convenience script that downloads the
current cattle and pig references listed below. Both paths filter FASTA records
to the sequences listed in the NCBI assembly report as assembled molecules,
unlocalized scaffolds, or unplaced scaffolds. That keeps chromosomes, assembled
sex chromosomes where they are part of the assembly, and unassigned contigs that
can matter when comparing against older conversion outputs.

```bash
conda activate genotype-converter-env

python validation/mixed_manifest/scripts/download_references.py
```

Equivalent single-reference command:

```bash
genotype-converter reference download-ncbi \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --accession GCF_002263795.3 \
  --ncbi-name ARS-UCD2.0 \
  --source-root validation/mixed_manifest/sources
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
validation/mixed_manifest/scripts/run_discover_sources.sh

validation/mixed_manifest/scripts/run_database_build.sh --yes
```

This discovers each species-level manifest set and pairs it with every reference
assembly under that species. Use `--workers 1` for large livestock references
unless the machine has enough memory for one minimap2 reference index per
worker. The build script requires `--yes` because the mixed-manifest workload
can be large. Use
`validation/mixed_manifest/scripts/run_database_build.sh --help` for overrides.

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

After the SQLite database has been built, run all generated example conversions:

```bash
validation/mixed_manifest/scripts/run_example_conversions.sh --overwrite

validation/mixed_manifest/scripts/run_check_example_conversions.sh
```

The conversion script runs CSV wide, CSV long, Illumina/GSGT matrix,
Illumina/GSGT long, Affymetrix/Axiom matrix, PLINK 1, and PLINK 2 examples for
each species/reference assembly. The check script writes:

```text
reports/example_conversion_check.csv
reports/example_conversion_check.md
```

The PLINK 2 synthetic example uses `--on-unconvertible-marker keep` because its
`.pgen` file is a lightweight metadata placeholder. Real PLINK 2 files should
use the normal default `exclude` policy when PLINK 2 is installed.

For one manual CSV conversion:

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

Some real source lookups contain the same marker ID more than once within one
manifest-derived source. Identical duplicate rows are collapsed on import.
Conflicting duplicates in the same source are skipped and reported because the
converter cannot safely select one rule by marker name alone. Duplicate marker
names across different manifests are still retained and handled by explicit
`--manifest-name` or `--resolve-mixed-manifests`.

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
