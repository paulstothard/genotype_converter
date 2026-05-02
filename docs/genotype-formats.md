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
| `convert --layout wide` | CSV, one sample per row and one marker per column | CSV | Genotype cells are rewritten. |
| `convert --layout long` | CSV, one sample-marker genotype per row | CSV | The genotype column is rewritten. |
| `convert-plink` | PLINK 1 binary fileset: `.bed`, `.bim`, `.fam` | PLINK 1 binary fileset | Allele labels in `.bim` are rewritten; `.bed` and `.fam` are copied unchanged. |
| `convert-pfile` | PLINK 2 fileset: `.pgen`, `.pvar`, `.psam` | PLINK 2 fileset | Biallelic allele labels in `.pvar` are rewritten; `.pgen` and `.psam` are copied unchanged. |

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
| `--lookup` | required | Lookup CSV from `build`. |
| `--from-format` | required | Input encoding: `AB`, `TOP`, `FORWARD`, `DESIGN`, `PLUS`, or `VCF`. |
| `--to-format` | required | Output encoding. |
| `--output` | required | Output file path. |
| `--layout` | `wide` | Input layout. |
| `--in-sep` | auto | Allele separator in input. Auto-detects `/`, space, tab, or adjacent single-character alleles. |
| `--out-sep` | `/` | Allele separator in output. |
| `--sample-col` | `sample_id` | Column name identifying the sample. |

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

## CSV Missing Data

These missing genotype codes are passed through unchanged:

```text
0
00
NA
N/A
--
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

Unknown allele values for a known marker are currently passed through unchanged.

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

The genotype bit matrix in `.bed` is not unpacked or rewritten. The `.fam` file
is copied unchanged. Only allele columns 5 and 6 in `.bim` are rewritten.

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
| `--bfile` | required | Input PLINK 1 binary prefix, without `.bed/.bim/.fam`. |
| `--lookup` | required | Lookup CSV from `build`. |
| `--from-format` | required | Input `.bim` allele encoding: `AB`, `TOP`, `FORWARD`, `DESIGN`, or `PLUS`. |
| `--to-format` | required | Output `.bim` allele encoding: `AB`, `TOP`, `FORWARD`, `DESIGN`, or `PLUS`. |
| `--out` | required | Output PLINK 1 binary prefix. |
| `--update-position` / `--keep-position` | keep | Also replace `.bim` chromosome and base-pair columns from the lookup table. |
| `--require-all-markers` / `--allow-missing-markers` | require | Fail if any `.bim` marker is absent from the lookup table. |

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

The genotype matrix in `.pgen` is not unpacked or rewritten. The `.psam` file is
copied unchanged. Biallelic allele labels in `.pvar` are rewritten. Multiallelic
`.pvar` rows are rejected for now because the lookup table is biallelic.

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
| `--pfile` | required | Input PLINK 2 prefix, without `.pgen/.pvar/.psam`. |
| `--lookup` | required | Lookup CSV from `build`. |
| `--from-format` | required | Input `.pvar` allele encoding: `AB`, `TOP`, `FORWARD`, `DESIGN`, or `PLUS`. |
| `--to-format` | required | Output `.pvar` allele encoding: `AB`, `TOP`, `FORWARD`, `DESIGN`, or `PLUS`. |
| `--out` | required | Output PLINK 2 prefix. |
| `--update-position` / `--keep-position` | keep | Also replace `.pvar` chromosome and base-pair columns from the lookup table. |
| `--require-all-markers` / `--allow-missing-markers` | require | Fail if any `.pvar` marker is absent from the lookup table. |

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
