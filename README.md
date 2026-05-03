# genotype_converter

Convert Illumina and Affymetrix SNP chip manifests to standardised format tables,
then re-encode genotype data files between any supported format encodings.

---

## Overview

**`build`** — takes a SNP chip manifest and a reference genome FASTA, aligns each
variant's flanking sequence to the reference using minimap2, and writes a set of
output files that record the chromosomal position and allele encoding for every
variant.

**`convert`** — takes a genotype data file and the lookup table produced by `build`,
and rewrites the allele calls from one encoding to another (for example, Illumina
TOP format to the genomic PLUS format used by GWAS pipelines).

---

## Supported format encodings

| Name | Description |
|---|---|
| **AB** | A/B allele coding (`A` = allele A, `B` = allele B) |
| **TOP** | Illumina TOP strand |
| **FORWARD** | Illumina/dbSNP forward-strand encoding from manifest F/R annotations |
| **DESIGN** | Illumina probe design strand |
| **PLUS** | Genomic plus-strand / reference-strand encoding |
| **VCF** | REF/ALT notation (reported as `REF` or `ALT`) |

In Illumina genotyping manifests, **FORWARD** and **PLUS** are not synonyms.
Forward/reverse is a source or dbSNP-oriented convention from the manifest,
whereas plus/minus is the genomic reference-strand convention.

---

## Installation

Requires Python ≥ 3.8. The `mappy` package (a Python binding for minimap2)
compiles a small C extension on install. On macOS, Xcode Command Line Tools are
required (`xcode-select --install`); on Linux, `gcc` and `zlib-dev`/`zlib-devel`.

### Option A — conda (recommended)

```bash
git clone https://github.com/paulstothard/genotype_converter.git
cd genotype_converter

conda env create -f environment.yml
conda activate genotype-converter-env

genotype-converter --help
```

### Option B — pip + venv

```bash
git clone https://github.com/paulstothard/genotype_converter.git
cd genotype_converter

python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e .
```

### Verify

```bash
genotype-converter --help
```

### Developer install

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

### Optional: Parquet output

```bash
pip install pyarrow
```

### Development plan

For planned features and design notes, see [Development Plan](docs/development-plan.md).

---

## build

Align a manifest to a reference genome and produce lookup / position / conversion
files.

```bash
genotype-converter build \
  --manifest path/to/manifest.csv \
  --reference path/to/genome.fa \
  --outdir output/ \
  --species cattle
```

| Flag | Default | Description |
|---|---|---|
| `--manifest` | (required) | Illumina or Affymetrix manifest CSV |
| `--reference` | (required) | Reference genome FASTA |
| `--outdir` | `output` | Root output directory |
| `--species` | `all` | Species name; used as a subdirectory label |
| `--workers` | `1` | Parallel alignment worker processes. Each worker loads the reference index, so increase carefully for large genomes. |
| `--align` / `--no-align` | off | Write a detailed alignment display file |
| `--parquet` / `--no-parquet` | off | Also write `lookup.parquet` (requires pyarrow) |
| `--progress` / `--no-progress` | on | Show alignment progress while building. |

For full mammalian genomes, start with `--workers 1`. Raising `--workers`
can speed up small references or machines with abundant memory, but each worker
loads its own minimap2 reference index. Values below 1 are treated as 1.

### Output files

All files land in `<outdir>/<species>/<ref_stem>/`:

| File | Description |
|---|---|
| `*.lookup.csv` | **Primary distributable.** One row per variant; columns named `A_in_TOP`, `B_in_PLUS`, etc., plus `determination_type` for QC. |
| `*.position.csv` | Chromosome, 1-based position, VCF REF/ALT, and `determination_type` per variant. |
| `*.conversion.csv` | Two rows per variant (one for allele A, one for allele B) with all format encodings. |
| `*.wide.csv` | All of the above in one row per variant. |
| `*.summary.txt` | Marker counts by type, alignment success rate, SHA-256 checksums, per-chromosome distribution. |
| `*.alignment.txt` | Alignment display per variant (`--align` flag). |
| `*.lookup.parquet` | Binary version of `lookup.csv` (`--parquet` flag). |

### Build summary

Example `*.summary.txt`:

```
genotype_converter build summary
==================================================
Generated:  Friday May 01 14:22:33 2026

Input files
-----------
Manifest:   /data/bovine_snp50.csv
  SHA256:   3a7f9c2d…
Reference:  /data/ARS-UCD1.2.fa
  SHA256:   c891a04f…

Run parameters
--------------
Species:    cattle
Workers:    8

Marker counts
-------------
Total:      54609
  SNPs:     54547  (99.9%)
  Indels:   62  (0.1%)
  By manifest type:
    Illumina: 54609  (100.0%)

Alignment results
-----------------
Positioned:     54423  (99.7%)
Not positioned: 186  (0.3%)

Per-chromosome marker counts
---------------------------
  1                    3541
  2                    3042
  ...

Output files
------------
  output/cattle/ARS-UCD1_2/bovine_snp50.ARS-UCD1_2.position.csv
  output/cattle/ARS-UCD1_2/bovine_snp50.ARS-UCD1_2.lookup.csv
  ...
```

---

## convert

Re-encode a genotype data file from one format to another using the lookup table
produced by `build`. For complete format details and before/after examples, see
[Genotype Input and Output Formats](docs/genotype-formats.md).

Supported genotype file layouts:

| Command | Input | Output | Notes |
|---|---|---|---|
| `convert --layout wide` | CSV, one sample per row and one marker per column | CSV | Best for small datasets, examples, and debugging. Can also process a folder of CSV files. |
| `convert --layout long` | CSV, one sample-marker genotype per row | CSV | Useful for database-style genotype tables. Can also process a folder of CSV files. |
| `convert-plink` | PLINK 1 binary fileset: `.bed`, `.bim`, `.fam` | PLINK 1 binary fileset | Rewrites allele labels in `.bim`; copies `.bed` and `.fam` unchanged. Can also process a folder of filesets. |
| `convert-pfile` | PLINK 2 fileset: `.pgen`, `.pvar`, `.psam` | PLINK 2 fileset | Rewrites biallelic `REF`/`ALT` labels in `.pvar`; copies `.pgen` and `.psam` unchanged. Can also process a folder of filesets. |

PLINK text formats such as `.ped/.map` are not converted directly. Convert them
to PLINK binary with PLINK first, then use `convert-plink` or `convert-pfile`.

```bash
genotype-converter convert \
  --genotypes mydata.csv \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --output converted.csv
```

To convert a folder of CSV genotype files with the same layout:

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

Batch conversion writes one output per input file using the default suffix
`.converted.csv`, for example `sample.csv` becomes `sample.converted.csv`.
It also writes `conversion_summary.csv` in the output directory with per-file
conversion counts. Existing batch outputs are not overwritten unless
`--overwrite` is supplied.

### PLINK 1 binary

```bash
genotype-converter convert-plink \
  --bfile mydata_top \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --out mydata_plus
```

To convert a folder of PLINK 1 binary filesets:

```bash
genotype-converter convert-plink \
  --bfile-dir plink_files/ \
  --pattern "*.bed" \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --outdir converted_plink_files/
```

An input prefix `herd_a` produces `herd_a.converted.bed`,
`herd_a.converted.bim`, and `herd_a.converted.fam` by default. Batch mode also
writes `conversion_summary.csv`.

### PLINK 2 p-files

```bash
genotype-converter convert-pfile \
  --pfile mydata_top \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --out mydata_plus
```

To convert a folder of PLINK 2 filesets:

```bash
genotype-converter convert-pfile \
  --pfile-dir pfiles/ \
  --pattern "*.pgen" \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --outdir converted_pfiles/
```

An input prefix `herd_a` produces `herd_a.converted.pgen`,
`herd_a.converted.pvar`, and `herd_a.converted.psam` by default. Batch mode also
writes `conversion_summary.csv`.

### Optional SQLite database

The lookup CSV workflow remains the main conversion path. For projects with many
species, assemblies, or manifests, lookup files can also be imported into an
optional SQLite database for inspection and database-backed conversion.
SQLite support uses Python's standard library; no extra package is required. See
[SQLite Conversion Database](docs/database.md) for schema details.

```bash
genotype-converter db init --database genotype_converter.sqlite

genotype-converter db import-lookup \
  --database genotype_converter.sqlite \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --manifest-name bovinehd_manifest_b

genotype-converter db marker \
  --database genotype_converter.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --marker SNP2 \
  --format table
```

The database tracks manifest identity with each imported lookup source, so the
same marker name can appear in multiple manifests without being merged into a
single rule.

CSV genotype conversion can read from the database when the manifest context is
specified explicitly:

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

PLINK commands can also use the database with the same explicit context:

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

The same database options work in batch mode. For example:

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
  --outdir converted_pfiles/
```

If `--manifest-name` is omitted, the converter infers the manifest from the
input marker IDs by choosing the imported source that matches the most markers
for the requested species and assembly. It fails instead of guessing if no
source matches or if two sources tie. Use `--manifest-name` when you know the
panel or when inference reports ambiguity.

For genotype files that intentionally combine markers from multiple manifests,
use `--resolve-mixed-manifests`. That mode resolves rules per marker and writes
a marker-resolution report. Conflicting duplicate marker rules fail by default
unless local marker context identifies one manifest clearly; use
`--on-ambiguous-marker skip` to leave unresolved markers unchanged.
See `examples/mixed_manifests/` for a small runnable mixed-manifest example.

To inspect a proposed source folder without running a build:

```bash
genotype-converter db discover-sources \
  --source-root database_sources \
  --format table
```

To build and import a source folder:

```bash
genotype-converter db build \
  --source-root database_sources \
  --database genotype_converter.sqlite \
  --build-outdir database_build \
  --workers 1
```

Each species/assembly folder should contain one reference file and one or more
manifest files.

---

## End-to-end example

This example uses the test data bundled in `tests/data/`.

### Step 1: build the lookup table

```bash
genotype-converter build \
  --manifest tests/data/manifest.csv \
  --reference tests/data/reference.fa \
  --outdir output/ \
  --species test
```

Output:

```
Build complete: 8 markers (5 SNPs, 3 indels), 8 positioned, 0 not positioned.
Full summary: output/test/reference/manifest.reference.summary.txt
```

Files written to `output/test/reference/`:

```
manifest.reference.lookup.csv
manifest.reference.position.csv
manifest.reference.conversion.csv
manifest.reference.wide.csv
manifest.reference.summary.txt
```

### Step 2: inspect the lookup table

Each row is one variant. Column names are self-describing:

```
marker_name,alt_marker_name,chromosome,position,ref_allele,alt_allele,
A_in_AB,B_in_AB,A_in_TOP,B_in_TOP,A_in_FORWARD,B_in_FORWARD,
A_in_DESIGN,B_in_DESIGN,A_in_PLUS,B_in_PLUS,A_vcf,B_vcf
SNP1,SNP1-0_T_F_1511658221,1,300,A,G,A,B,A,G,A,G,A,G,A,G,REF,ALT
SNP2,SNP2-0_B_F_2328966441,1,700,T,G,A,B,A,C,T,G,T,G,T,G,REF,ALT
```

### Step 3: re-encode a genotype file from TOP to PLUS

Create `mydata.csv`:

```
sample_id,SNP1,SNP2,SNP3,SNP4,SNP5
SAMPLE001,A/A,A/C,A/C,A/A,A/A
SAMPLE002,A/G,A/A,C/C,A/G,A/C
SAMPLE003,G/G,C/C,A/C,G/G,C/C
```

```bash
genotype-converter convert \
  --genotypes mydata.csv \
  --lookup output/test/reference/manifest.reference.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --output mydata_plus.csv
```

### Step 4: re-encode from AB notation

```bash
genotype-converter convert \
  --genotypes mydata_ab.csv \
  --lookup output/test/reference/manifest.reference.lookup.csv \
  --from-format AB \
  --to-format FORWARD \
  --output mydata_forward.csv
```

AB input can use separators or adjacent allele labels:

```
sample_id,SNP1,SNP2
SAMPLE001,A/B,A/A
SAMPLE002,BB,AB
```

### Batch CSV input

```bash
genotype-converter convert \
  --genotypes-dir genotype_files/ \
  --pattern "*.csv" \
  --lookup output/test/reference/manifest.reference.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --layout wide \
  --outdir converted_genotype_files/
```

For input files `herd_a.csv` and `herd_b.csv`, the default output names are
`herd_a.converted.csv` and `herd_b.converted.csv`. Batch mode also writes
`conversion_summary.csv` with row, genotype, and allele-count summaries for each
input file.

### Long-format input

```bash
genotype-converter convert \
  --genotypes mydata_long.csv \
  --lookup output/test/reference/manifest.reference.lookup.csv \
  --from-format TOP \
  --to-format PLUS \
  --layout long \
  --output mydata_long_plus.csv
```

---

## How it works

1. **Query preparation.** The flanking sequence from the manifest
   (e.g. `AAACCC[A/G]TTTGGG`) has the variant bracket replaced with a
   tracked synthetic `N`, and non-GATCN characters are removed. Existing
   flanking `N` bases are preserved because they are part of the manifest's
   sequence context. The result is the minimap2 query.

2. **Alignment.** The query is aligned to the reference with minimap2
   (`sr` preset, up to 5 hits). When multiple hits are returned, the hit with
   the highest mapping quality is used.

3. **Position recovery.** For Illumina SNPs with probe sequences, the probe is
   matched to the left or right side of the aligned flanking sequence. Some
   probes stop next to the assayed base; others include the assayed allele as
   the terminal probe base. The reported SNP position is the reference base
   implied by that probe placement: the base immediately after a left-side
   adjacent probe, immediately before a right-side adjacent probe, or the
   allele-including probe base itself. This keeps gap-adjacent SNPs tied to the
   assayed probe context. If probe placement cannot be determined, the position
   of `N` in the query is mapped by walking the CIGAR string. When the
   probe-derived base is found by a unique probe match in the selected
   reference neighborhood, that probe placement defines the assayed base. If
   only query-side probe placement is available and it differs from direct
   CIGAR mapping, the manifest alleles are used only if exactly one candidate
   base matches those alleles; otherwise the probe-derived position is kept and
   the site is counted as ambiguous in the build summary. For SNPs with a CIGAR
   gap near the assayed site, a short reference window around the minimap2 hit
   is locally realigned before applying the probe and allele rules.

4. **VCF REF/ALT.** For SNPs, the reference base at the variant position
   determines which allele is REF and which is ALT. For indels, an anchor-based
   VCF representation is computed: the anchor base (one before the indel site)
   is prepended to both REF and ALT alleles.

5. **Format encoding.** Strand relationships (TOP/BOT, FORWARD, DESIGN, PLUS)
   are derived from the manifest's `IlmnStrand`, `SourceStrand`, and `IlmnID`
   fields plus the observed reference alignment. Affymetrix alleles are treated
   as forward-strand alleles and then projected onto PLUS using the alignment
   orientation.

---

## Performance

Typical runtimes on a modern server (8+ cores):

| Panel size | Approx. time |
|---|---|
| 50K SNPs (Bovine SNP50) | ~2 min |
| 150K SNPs (Bovine HD) | ~5 min |
| 800K SNPs (high-density panels) | ~20–30 min |

---

## Citation

If you use genotype_converter in your research, please cite:

Grant JR, Herman EK, Barlow LD, Miglior F, Schenkel FS, Baes CF, Stothard P.
A large structural variant collection in Holstein cattle and associated database
for variant discovery, characterization, and application. BMC Genomics.
2024;25(1):903. doi: 10.1186/s12864-024-10812-2. PMID: 39350025.

<details>
<summary>BibTeX</summary>

```bibtex
@article{Grant2024LargeStructuralVariant,
  author  = {Grant, Jason R and Herman, Emily K and Barlow, Lael D and
             Miglior, Filippo and Schenkel, Flavio S and Baes, Christine F and
             Stothard, Paul},
  title   = {A large structural variant collection in {Holstein} cattle and
             associated database for variant discovery, characterization, and application},
  journal = {BMC Genomics},
  year    = {2024},
  volume  = {25},
  number  = {1},
  pages   = {903},
  doi     = {10.1186/s12864-024-10812-2},
  pmid    = {39350025}
}
```

</details>

---

## License

MIT
