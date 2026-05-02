# Full-Scale Validation Data

Use this folder for release validation against the original pipeline. Keep large
real data files local unless you explicitly decide they should be committed.

## Directory Layout

Place files under these directories. Provenance-preserving filenames are fine;
do not rename source files just to match the generic examples.

```text
validation/
  data/
    <panel-manifest>.csv
    <reference-genome>.fa
  old_pipeline/
    <panel>.<reference>.position.csv
    <panel>.<reference>.conversion.csv
    <panel>.<reference>.wide.csv        # optional
    <panel>.<reference>.alignment.txt   # optional
    <panel>.<reference>.blast.csv       # optional
  new_pipeline/
    generated files from genotype-converter build
  reports/
    comparison notes and diff summaries
```

## Required Files

For a single validation panel, include one full manifest, one matching reference
FASTA, and the old pipeline outputs for that same manifest/reference pair.

| File | Required | Description |
|---|---:|---|
| `data/*.csv` | yes | Full Illumina or Affymetrix manifest used by the old pipeline. |
| `data/*.fa` or `data/*.fasta` | yes | Exact reference genome FASTA used by the old pipeline. |
| `old_pipeline/*.position.csv` | yes | Position output from the old pipeline. |
| `old_pipeline/*.conversion.csv` | yes | Conversion output from the old pipeline. |
| `old_pipeline/*.wide.csv` | no | Old wide output, if available. |
| `old_pipeline/*.alignment.txt` | no | Old alignment display, useful for manual debugging. |
| `old_pipeline/*.blast.csv` | no | Old BLAST details, useful for manual debugging. |

Current validation set:

```text
data/bovinehd-manifest-b.csv
data/ARS_UCD_v2.0.fa
old_pipeline/bovinehd-manifest-b.ARS_UCD_v2_0.position.csv
old_pipeline/bovinehd-manifest-b.ARS_UCD_v2_0.conversion.csv
old_pipeline/bovinehd-manifest-b.ARS_UCD_v2_0.wide.csv
old_pipeline/bovinehd-manifest-b.ARS_UCD_v2_0.alignment.txt
old_pipeline/bovinehd-manifest-b.ARS_UCD_v2_0.blast.csv
```

If you want to keep multiple validation panels, create named subfolders such as
`validation/panels/bovine_snp50/` with the same internal layout.

## Run the New Pipeline

From the repository root:

```bash
genotype-converter build \
  --manifest validation/data/bovinehd-manifest-b.csv \
  --reference validation/data/ARS_UCD_v2.0.fa \
  --outdir validation/new_pipeline \
  --species bos_taurus \
  --workers 1
```

Use `--workers 1` explicitly for this full genome. Higher worker counts should
only be used when the machine has enough memory for one reference index per
worker. Add `--align` only when debugging a small subset; for the full bovine HD
panel it writes a very large alignment display file.

The generated files will be under:

```text
validation/new_pipeline/bos_taurus/ARS_UCD_v2_0/
```

## What to Compare

Start with these checks:

1. Marker counts in the new `*.summary.txt` versus the old pipeline run.
2. `chromosome`, `position`, `VCF_REF`, and `VCF_ALT` in position files.
3. `AB`, `TOP`, `FORWARD`, `DESIGN`, `PLUS`, and `VCF` in conversion files.
4. Any markers missing from either pipeline output.
5. Any markers with no new alignment or ambiguous-looking mappings.

Expected differences should be recorded in `reports/`, including the command
used, date, input file checksums, and a short explanation.

The helper script compares the current bovine HD validation files:

```bash
python validation/compare_outputs.py
```

It writes:

```text
validation/reports/comparison_summary.md
validation/reports/position_mismatches.csv
validation/reports/conversion_mismatches.csv
```

## Notes

- Use the same reference FASTA as the old pipeline. A different assembly or
  chromosome naming convention will produce noisy differences.
- The new pipeline writes AB output as literal `A` and `B` labels. TOP,
  FORWARD, DESIGN, PLUS, and VCF columns carry the converted allele values.
- Do not edit old pipeline outputs after copying them here. If cleanup is
  needed for comparison, write a separate normalized copy or document the exact
  normalization command in `reports/`.
