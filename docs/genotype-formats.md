# Genotype Input and Output Formats

This guide describes the genotype files accepted by `genotype-converter convert`,
`genotype-converter convert-plink`, and `genotype-converter convert-pfile`.

## Allele Encodings

| Name | Meaning |
|---|---|
| `AB` | A/B allele coding (`A` = allele A, `B` = allele B). |
| `TOP` | Illumina TOP strand. |
| `FORWARD` | Illumina/dbSNP forward-strand encoding from manifest F/R annotations. |
| `DESIGN` | Illumina probe design strand. |
| `PLUS` | Genomic plus-strand / reference-strand encoding. |
| `VCF` | REF/ALT notation, reported as `REF` or `ALT`. |

In Illumina genotyping manifests, `FORWARD` and `PLUS` are not synonyms.
Forward/reverse is a source or dbSNP-oriented convention from the manifest,
whereas plus/minus is the genomic reference-strand convention.

## Supported Genotype File Layouts

| Command | Input | Output | What Changes |
|---|---|---|---|
| `convert --layout wide` | CSV, one sample per row and one marker per column | CSV | Genotype cells are rewritten. Unconvertible marker columns are excluded by default. Single files and folders are supported. |
| `convert --layout long` | CSV, one sample-marker genotype per row | CSV | The genotype column is rewritten. Unconvertible marker rows are excluded by default. Single files and folders are supported. |
| `convert --layout illumina-matrix` | Illumina GenomeStudio/GSGT matrix report | Illumina/GSGT matrix report | Matrix genotype calls are rewritten. Unconvertible marker rows are excluded by default. |
| `convert --layout illumina-long` | Illumina GenomeStudio/GSGT long report | Illumina/GSGT long report | Target allele columns are filled or rewritten. Unconvertible marker rows are excluded by default. |
| `convert --layout affymetrix-matrix` | Affymetrix/Axiom paired-call matrix | Affymetrix/Axiom paired-call matrix | The AB or native nucleotide call in each sample pair is rewritten. Unconvertible marker rows are excluded by default. |
| `convert-plink` | PLINK 1 binary fileset: `.bed`, `.bim`, `.fam` | PLINK 1 binary fileset | Allele labels in `.bim` are rewritten. Unconvertible variants are excluded with PLINK by default so `.bed/.bim/.fam` stay synchronized. Single filesets and folders are supported. |
| `convert-pfile` | PLINK 2 fileset: `.pgen`, `.pvar`, `.psam` | PLINK 2 fileset | Biallelic allele labels in `.pvar` are rewritten. Unconvertible variants are excluded with PLINK2 by default so `.pgen/.pvar/.psam` stay synchronized. Single filesets and folders are supported. |

PLINK text formats such as `.ped/.map` are not converted directly. Convert them
to PLINK binary with PLINK first, then use `convert-plink` or `convert-pfile`.

## CSV Wide

Wide CSV has one row per sample and one column per marker.

```text
sample_id,SNP1,SNP2,SNP3
SAMPLE001,A/A,A/C,G/G
SAMPLE002,A/G,C/C,A/G
```

Example:

```bash
genotype-converter convert \
  --genotypes mydata_top.csv \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --layout wide \
  --output mydata_plus.csv
```

Options:

| Flag | Default | Description |
|---|---|---|
| `--genotypes` | required | Input genotype CSV. |
| `--genotypes-dir` | unset | Directory of genotype CSV files. Use this instead of `--genotypes` for batch conversion. |
| `--pattern` | `*.csv` | File pattern used with `--genotypes-dir`. |
| `--lookup` | required unless `--database` is used | Lookup CSV from `build`. |
| `--database` | unset | SQLite conversion database. |
| `--species` | required with `--database` | Species name for database-backed conversion. |
| `--assembly` | required with `--database` | Reference assembly name for database-backed conversion. |
| `--manifest-name` | optional with `--database` | Manifest or panel name for database-backed conversion. If omitted, the command tries conservative marker-based inference. |
| `--resolve-mixed-manifests` | off | Resolve database rules per marker for inputs containing markers from multiple manifests. |
| `--on-ambiguous-marker` | `fail` | In mixed-manifest mode, either fail on unresolved conflicting rules or skip them unchanged. |
| `--resolution-report` | auto | CSV report path for mixed-manifest marker rule decisions. |
| `--from-format` | required | Input encoding: `AB`, `TOP`, `FORWARD`, `DESIGN`, `PLUS`, or `VCF`. |
| `--to-format` | required | Output encoding. |
| `--output` | required for single-file mode | Output file path. |
| `--outdir` | required for batch mode | Output directory for converted batch files. |
| `--suffix` | `.converted.csv` | Filename suffix for batch outputs. |
| `--overwrite` | off | Allow batch mode to replace existing outputs. |
| `--layout` | `wide` | Input layout. |
| `--in-sep` | auto | Allele separator in input. Auto-detects `/`, space, tab, or adjacent single-character alleles. |
| `--out-sep` | `/` | Allele separator in output. |
| `--sample-col` | `sample_id` | Column name identifying the sample. |
| `--on-unconvertible-marker` | `exclude` | How to handle markers that are missing from the lookup table or lack target allele labels: `exclude`, `fail`, or `keep`. |

With the default `exclude` policy, unconvertible markers are removed from text
outputs before conversion: wide CSV marker columns are removed, long CSV marker
rows are removed, and Illumina/Affymetrix marker rows are removed. The command
writes `<output>.marker_conversion_report.csv` and, when markers are excluded,
`<output>.exclude_markers.txt`. Use `--on-unconvertible-marker fail` for strict
audits, or `--on-unconvertible-marker keep` to preserve the original file shape
and leave unresolved allele values unchanged.

Batch example:

```bash
genotype-converter convert \
  --genotypes-dir genotype_files/ \
  --pattern "*.csv" \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --layout wide \
  --outdir converted_genotype_files/
```

For an input named `sample.csv`, the default output is
`sample.converted.csv`. Existing batch outputs are protected unless
`--overwrite` is supplied. Batch mode also writes `conversion_summary.csv` in
the output directory.

CSV conversion can also use an SQLite database instead of `--lookup`:

```bash
genotype-converter convert \
  --genotypes mydata_top.csv \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --manifest-name bovinehd_manifest_b \
  --from-format TOP \
  --to-format PLUS \
  --layout wide \
  --output mydata_plus.csv
```

If `--manifest-name` is omitted, the command chooses the imported lookup source
that matches the most input markers for the requested species and assembly. It
fails if no source matches or if the best match is tied.

That default inference chooses one manifest for the whole input. Use
`--resolve-mixed-manifests` for files that intentionally contain markers from
more than one manifest. Mixed-manifest mode writes a marker-resolution report
with the input path, selected manifest, candidate manifests, candidate count,
and selection reason for each marker.

The CSV batch summary includes:

| Column | Meaning |
|---|---|
| `input_path` | Input genotype CSV. |
| `output_path` | Converted genotype CSV. |
| `rows_total` | Input data rows processed. |
| `markers_total` | Marker columns in wide layout, or distinct marker IDs in long layout. |
| `genotype_cells_total` | Genotype cells considered for conversion. |
| `genotypes_parsed` | Genotypes recognized as two allele labels. |
| `genotypes_changed` | Genotype cells whose output value differs from input. |
| `missing_or_unparsed_genotypes` | Missing or unrecognized genotype cells left unchanged. |
| `genotypes_excluded` | Genotype cells or long-format rows removed because their marker was excluded. |
| `alleles_changed` | Individual allele labels changed. |
| `unknown_alleles` | Allele labels left unchanged in `keep` mode because no complete conversion rule was available. |
| `markers_missing_lookup` | Distinct marker IDs not found in the lookup table. |
| `markers_incomplete_mapping` | Distinct marker IDs found in the lookup table but not fully convertible to the target format. |
| `markers_excluded` | Distinct marker IDs removed from the output. |
| `marker_report_path` | CSV report listing unconvertible markers and reasons. |
| `exclude_marker_path` | Marker ID list used when markers were excluded. |

For single-file CSV conversion, the command-line summary also reports missing
or unparsed genotype cells, unresolved allele labels, excluded markers, and
marker-report paths when those counts are non-zero.

## CSV Long

Long CSV has one row per sample-marker genotype.

```text
sample_id,marker_name,genotype
SAMPLE001,SNP1,A/A
SAMPLE001,SNP2,A/C
SAMPLE002,SNP1,A/G
SAMPLE002,SNP2,C/C
```

Example:

```bash
genotype-converter convert \
  --genotypes mydata_top_long.csv \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --layout long \
  --output mydata_plus_long.csv
```

Additional options:

| Flag | Default | Description |
|---|---|---|
| `--marker-col` | `marker_name` | Column name for marker IDs. |
| `--genotype-col` | `genotype` | Column name for genotype calls. |

## Illumina GenomeStudio Matrix

Illumina GenomeStudio/GSGT matrix reports have a `[Header]` section followed by
`[Data]`. Marker names are rows and sample IDs are columns.

```text
[Header]
GSGT Version	2.0.4
[Data]
	SAMPLE001	SAMPLE002
SNP1	AA	AG
SNP2	CC	AC
```

Use `--layout illumina-matrix`. The output preserves the header and matrix
shape, rewriting genotype calls in the sample columns.

```bash
genotype-converter convert \
  --genotypes gsgt_top.txt \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --layout illumina-matrix \
  --out-sep "" \
  --output gsgt_plus.txt
```

## Illumina GenomeStudio Long

Illumina/GSGT long reports have one row per marker/sample and can contain
multiple allele encodings in separate columns.

```text
[Header]
GSGT Version	2.0.4
[Data]
SNP Name	Sample ID	Allele1 - Top	Allele2 - Top	Allele1 - Forward	Allele2 - Forward	Allele1 - AB	Allele2 - AB	Allele1 - Design	Allele2 - Design	Allele1 - Plus	Allele2 - Plus	GC Score
SNP1	SAMPLE001	A	G	-	-	A	B	-	-	-	-	0.99
```

Use `--layout illumina-long`. The converter reads the two columns matching
`--from-format` and fills or rewrites the two columns matching `--to-format`.
`VCF` is not supported for this layout because GenomeStudio long reports do not
have VCF genotype columns.

```bash
genotype-converter convert \
  --genotypes gsgt_long_top.txt \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --layout illumina-long \
  --output gsgt_long_plus.txt
```

## Affymetrix/Axiom Matrix

Affymetrix/Axiom paired-call matrices repeat each sample column. The first
column in each pair is the AB call and the second is the native nucleotide call.

```text
probeset_id	SAMPLE001	SAMPLE001	SAMPLE002	SAMPLE002
AX-1	AA	TT	AB	TC
AX-2	NoCall	---	BB	GG
```

Use `--layout affymetrix-matrix`. If `--from-format AB`, the converter reads
the first column in each sample pair. Otherwise it reads the native nucleotide
column. If `--to-format AB`, it rewrites the first column; otherwise it rewrites
the native nucleotide column. `VCF` is not supported for this layout.

```bash
genotype-converter convert \
  --genotypes axiom_ab_and_native.txt \
  --lookup output/cattle/genome/axiom.genome.lookup.csv \
  --from-format AB \
  --to-format PLUS \
  --layout affymetrix-matrix \
  --output axiom_plus.txt
```

## CSV Missing Data

These missing genotype codes are passed through unchanged:

```text
0
00
NA
N/A
--
-
.
0/0
00/00
```

If one allele in a genotype is missing, the full genotype is left unchanged.

## CSV Marker Rules

Wide CSV marker columns must match `marker_name` values in the lookup table.
Long CSV marker-column values must also match `marker_name` values. Missing
required columns, unknown markers, unrecognized lookup files, and unrecognized
manifests raise errors instead of producing partial output.

For `convert` text layouts, marker lookup problems and allele conversion
problems are handled separately. The table below uses the same terms as the
PLINK tables later in this document.

| Problem type | Exact condition | Output action | Accounting/reporting |
|---|---|---|---|
| Marker lookup missing | A wide marker column, long `marker_name`, or report marker is not present in the lookup table. | Exclude the marker from the output. | `markers_missing_lookup` and `markers_excluded` increment once per marker. Marker report row has `reason=missing_lookup` and `action=exclude`. |
| Source allele missing from rule | Marker ID is present, but a genotype allele is not listed in that marker's `--from-format` rule. | Exclude the marker from the output. | `markers_incomplete_mapping` and `markers_excluded` increment once per marker. Marker report row has `reason=input_allele_not_in_<FORMAT>` and `action=exclude`. |
| Target allele missing from rule | Marker ID is present and the input allele is recognized, but the requested `--to-format` allele is blank. | Exclude the marker from the output. | `markers_incomplete_mapping` and `markers_excluded` increment once per marker. Marker report row has `reason=missing_<FORMAT>_allele` and `action=exclude`. |
| Genotype missing or unparsed | Genotype is a supported missing code, contains a missing allele, or cannot be parsed as two alleles. | Keep the full genotype cell unchanged. | `missing_or_unparsed_genotypes` increments once per genotype cell. |

This policy applies to `wide`, `long`, `illumina-matrix`, `illumina-long`, and
`affymetrix-matrix` layouts. In `keep` mode, text outputs keep the same
row/column shape and unresolved allele values increment `unknown_alleles`
instead of excluding the marker. In `fail` mode, any marker lookup or
allele-rule problem stops the conversion after writing the marker report.

## PLINK 1 Binary: `.bed/.bim/.fam`

PLINK 1 binary datasets are three-file sets:

```text
mydata_top.bed
mydata_top.bim
mydata_top.fam
```

`convert-plink` writes a complete output fileset:

```text
mydata_plus.bed
mydata_plus.bim
mydata_plus.fam
```

When all variants are convertible, the genotype bit matrix in `.bed` is copied
unchanged. If any variants are unconvertible under the default policy, PLINK is
used to remove them before allele columns 5 and 6 in `.bim` are rewritten. The
`.fam` file is copied unchanged.

Example `.bim` input in TOP encoding:

```text
1 SNP1 0 300 A G
1 SNP2 0 700 A C
```

After `--from-format TOP --to-format PLUS`:

```text
1 SNP1 0 300 A G
1 SNP2 0 700 T G
```

Example:

```bash
genotype-converter convert-plink \
  --bfile mydata_top \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --out mydata_plus
```

Options:

| Flag | Default | Description |
|---|---|---|
| `--bfile` | required for single-fileset mode | Input PLINK 1 binary prefix, without `.bed/.bim/.fam`. |
| `--bfile-dir` | unset | Directory of PLINK 1 binary filesets. Use this instead of `--bfile` for batch conversion. |
| `--pattern` | `*.bed` | File pattern used with `--bfile-dir`. The pattern should match `.bed` files. |
| `--lookup` | required unless `--database` is used | Lookup CSV from `build`. |
| `--database` | unset | SQLite conversion database. |
| `--species` | required with `--database` | Species name for database-backed conversion. |
| `--assembly` | required with `--database` | Reference assembly name for database-backed conversion. |
| `--manifest-name` | optional with `--database` | Manifest or panel name for database-backed conversion. If omitted, the command tries conservative marker-based inference. |
| `--resolve-mixed-manifests` | off | Resolve database rules per marker for inputs containing markers from multiple manifests. |
| `--on-ambiguous-marker` | `fail` | In mixed-manifest mode, either fail on unresolved conflicting rules or skip them unchanged. |
| `--resolution-report` | auto | CSV report path for mixed-manifest marker rule decisions. |
| `--from-format` | required | Input `.bim` allele encoding: `AB`, `TOP`, `FORWARD`, `DESIGN`, or `PLUS`. |
| `--to-format` | required | Output `.bim` allele encoding: `AB`, `TOP`, `FORWARD`, `DESIGN`, or `PLUS`. |
| `--out` | required for single-fileset mode | Output PLINK 1 binary prefix. |
| `--outdir` | required for batch mode | Output directory for converted PLINK 1 filesets. |
| `--suffix` | `.converted` | Filename suffix for batch output prefixes. |
| `--overwrite` | off | Allow batch mode to replace existing outputs. |
| `--update-position` / `--keep-position` | keep | Also replace `.bim` chromosome and base-pair columns from the lookup table. |
| `--on-unconvertible-marker` | `exclude` | How to handle `.bim` markers that are missing from the lookup table or lack target allele labels: `exclude`, `fail`, or `keep`. |
| `--plink` | `plink` | PLINK executable used when unconvertible markers are excluded. |

By default, unconvertible PLINK variants are removed with PLINK before allele
labels are rewritten. This keeps `.bed`, `.bim`, and `.fam` synchronized. The
command writes `<out>.marker_conversion_report.csv` and, when variants are
excluded, `<out>.exclude_markers.txt`. Use `--on-unconvertible-marker fail` for
strict audits, or `--on-unconvertible-marker keep` only when you intentionally
want unresolved variants left unchanged.

PLINK unconvertible-marker handling uses the same problem categories as text
conversion, but the default output action is different because the genotype
matrix has to stay synchronized with the variant list.

| Problem type | Exact condition | Default `exclude` output action | Accounting/reporting |
|---|---|---|---|
| Marker lookup missing | The `.bim` marker ID is not present in the lookup table. | Write the marker ID to `<out>.exclude_markers.txt` and remove the variant with PLINK before writing final `.bed/.bim/.fam`. | `variants_missing_lookup` increments once per variant. Marker report row has `reason=missing_lookup` and `action=exclude`. |
| Source allele missing from rule | Marker ID is present, but a `.bim` allele is not listed in that marker's `--from-format` rule. | Write the marker ID to `<out>.exclude_markers.txt` and remove the variant with PLINK. | `variants_incomplete_mapping` increments once per variant. Marker report row has `reason=input_allele_not_in_<FORMAT>` and `action=exclude`. |
| Target allele missing from rule | Marker ID is present and the `.bim` allele is recognized, but the requested `--to-format` allele is blank. | Write the marker ID to `<out>.exclude_markers.txt` and remove the variant with PLINK. | `variants_incomplete_mapping` increments once per variant. Marker report row has `reason=missing_<FORMAT>_allele` and `action=exclude`. |
| PLINK allele is missing | Allele value in `.bim` is PLINK missing allele `0`. | Keep `0` as missing. This does not by itself make the variant unconvertible. | No unconvertible-marker reason. The variant is counted normally if the remaining allele labels are convertible. |

Batch example:

```bash
genotype-converter convert-plink \
  --bfile-dir plink_files/ \
  --pattern "*.bed" \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --outdir converted_plink_files/
```

For an input fileset prefix `sample`, the default outputs are
`sample.converted.bed`, `sample.converted.bim`, and `sample.converted.fam`.
Existing batch outputs are protected unless `--overwrite` is supplied. Batch
mode also writes `conversion_summary.csv` in the output directory.

Database example:

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

The PLINK batch summary includes:

| Column | Meaning |
|---|---|
| `input_prefix` | Input fileset prefix. |
| `output_prefix` | Converted fileset prefix. |
| `variants_total` | Variant rows processed in `.bim` or `.pvar`. |
| `variants_converted` | Variant rows whose allele labels were converted or confirmed convertible. In `keep` mode, unconvertible rows may still be present in the output but are not counted here. |
| `variants_missing_lookup` | Variant rows whose marker IDs were not found. |
| `variants_incomplete_mapping` | Variant rows found in the lookup table but not fully convertible to the target format. |
| `variants_excluded` | Variant rows removed from the output fileset. |
| `alleles_changed` | Individual allele labels changed. |
| `genotype_path` | Output `.bed` or `.pgen` path. |
| `variant_path` | Output `.bim` or `.pvar` path. |
| `sample_path` | Output `.fam` or `.psam` path. |
| `marker_report_path` | CSV report listing unconvertible markers and reasons. |
| `exclude_marker_path` | Marker ID list passed to PLINK/PLINK2 when variants were excluded. |

`VCF` is not a `convert-plink` target because PLINK `.bim` allele columns should
contain allele labels such as `A`, `C`, `I`, or `D`, not `REF` or `ALT`
keywords.

## PLINK 2: `.pgen/.pvar/.psam`

PLINK 2 datasets are three-file sets:

```text
mydata_top.pgen
mydata_top.pvar
mydata_top.psam
```

`convert-pfile` writes a complete output fileset:

```text
mydata_plus.pgen
mydata_plus.pvar
mydata_plus.psam
```

When all variants are convertible, the genotype matrix in `.pgen` is copied
unchanged. If any variants are unconvertible under the default policy, PLINK2 is
used to remove them before biallelic allele labels in `.pvar` are rewritten. The
`.psam` file is copied unchanged. Multiallelic `.pvar` rows are rejected for now
because the lookup table is biallelic.

Example `.pvar` input in TOP encoding:

```text
#CHROM POS ID   REF ALT
1      300 SNP1 A   G
1      700 SNP2 A   C
```

After `--from-format TOP --to-format PLUS`:

```text
#CHROM POS ID   REF ALT
1      300 SNP1 A   G
1      700 SNP2 T   G
```

Example:

```bash
genotype-converter convert-pfile \
  --pfile mydata_top \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --out mydata_plus
```

Options:

| Flag | Default | Description |
|---|---|---|
| `--pfile` | required for single-fileset mode | Input PLINK 2 prefix, without `.pgen/.pvar/.psam`. |
| `--pfile-dir` | unset | Directory of PLINK 2 filesets. Use this instead of `--pfile` for batch conversion. |
| `--pattern` | `*.pgen` | File pattern used with `--pfile-dir`. The pattern should match `.pgen` files. |
| `--lookup` | required unless `--database` is used | Lookup CSV from `build`. |
| `--database` | unset | SQLite conversion database. |
| `--species` | required with `--database` | Species name for database-backed conversion. |
| `--assembly` | required with `--database` | Reference assembly name for database-backed conversion. |
| `--manifest-name` | optional with `--database` | Manifest or panel name for database-backed conversion. If omitted, the command tries conservative marker-based inference. |
| `--resolve-mixed-manifests` | off | Resolve database rules per marker for inputs containing markers from multiple manifests. |
| `--on-ambiguous-marker` | `fail` | In mixed-manifest mode, either fail on unresolved conflicting rules or skip them unchanged. |
| `--resolution-report` | auto | CSV report path for mixed-manifest marker rule decisions. |
| `--from-format` | required | Input `.pvar` allele encoding: `AB`, `TOP`, `FORWARD`, `DESIGN`, or `PLUS`. |
| `--to-format` | required | Output `.pvar` allele encoding: `AB`, `TOP`, `FORWARD`, `DESIGN`, or `PLUS`. |
| `--out` | required for single-fileset mode | Output PLINK 2 prefix. |
| `--outdir` | required for batch mode | Output directory for converted PLINK 2 filesets. |
| `--suffix` | `.converted` | Filename suffix for batch output prefixes. |
| `--overwrite` | off | Allow batch mode to replace existing outputs. |
| `--update-position` / `--keep-position` | keep | Also replace `.pvar` chromosome and base-pair columns from the lookup table. |
| `--on-unconvertible-marker` | `exclude` | How to handle `.pvar` markers that are missing from the lookup table or lack target allele labels: `exclude`, `fail`, or `keep`. |
| `--plink2` | `plink2` | PLINK2 executable used when unconvertible markers are excluded. |

By default, unconvertible PLINK 2 variants are removed with PLINK2 before allele
labels are rewritten. This keeps `.pgen`, `.pvar`, and `.psam` synchronized and
produces the same marker report and exclude-list files described for PLINK 1.
The same accounting is used: missing `.pvar` IDs increment
`variants_missing_lookup`; recognized IDs with incomplete source or target
allele rules increment `variants_incomplete_mapping`; excluded variants
increment `variants_excluded`.
If PLINK2 is not available on `PATH`, install it separately and pass the
executable with `--plink2 /path/to/plink2`, or use
`--on-unconvertible-marker fail` for a strict report without writing filtered
p-files.

Batch example:

```bash
genotype-converter convert-pfile \
  --pfile-dir pfiles/ \
  --pattern "*.pgen" \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --outdir converted_pfiles/
```

For an input fileset prefix `sample`, the default outputs are
`sample.converted.pgen`, `sample.converted.pvar`, and `sample.converted.psam`.
Existing batch outputs are protected unless `--overwrite` is supplied. Batch
mode also writes `conversion_summary.csv` in the output directory using the
PLINK batch summary columns described above.

Database example:

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

`VCF` is not a `convert-pfile` target. The command rewrites allele labels in the
existing `.pvar`; it does not transform PLINK genotype records into a VCF
representation.

## Unsupported Genotype Inputs

These formats are not directly converted:

| Format | Suggested Workflow |
|---|---|
| PLINK `.ped/.map` | Convert to PLINK binary with PLINK, then run `convert-plink`. |
| VCF input | Convert to a supported table or PLINK format first. |
| BGEN | Convert to PLINK 2 or another supported format first. |

Unsupported formats can still be useful for downstream export. For example,
after running `convert-plink`, PLINK or PLINK 2 can export the recoded fileset
to VCF or other formats.
