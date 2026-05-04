# PLINK Conversion Comparison Workspace

Use this workspace for collaborator-provided PLINK examples. Keep the real
files local; inputs, converted outputs, and reports are ignored by Git.

## Folder Layout

```text
validation/plink_comparison/
  input/
    <example>/
      original.bed
      original.bim
      original.fam
      expected_plus.bed
      expected_plus.bim
      expected_plus.fam
      notes.txt
  converted/
    outputs from this converter
  reports/
    comparison CSV and Markdown reports
  scripts/
    run_plink_comparison.sh
    compare_plink_bim.py
```

Record the original allele format and the target reference assembly in
`notes.txt` or in the report notes. The target is expected to be PLUS on a
specific reference assembly.

## Compare Existing BIM Files

```bash
python validation/plink_comparison/scripts/compare_plink_bim.py \
  --expected validation/plink_comparison/input/mike_example/expected_plus.bim \
  --actual validation/plink_comparison/converted/mike_example.converted.bim \
  --report-prefix validation/plink_comparison/reports/mike_example
```

The report classifies each marker as exact, allele-order swapped, complement,
swapped-complement, mismatch, missing, or extra.

## Run Conversion And Compare

With a database:

```bash
validation/plink_comparison/scripts/run_plink_comparison.sh \
  --original-bfile validation/plink_comparison/input/mike_example/original \
  --expected-bfile validation/plink_comparison/input/mike_example/expected_plus \
  --out-prefix validation/plink_comparison/converted/mike_example \
  --report-prefix validation/plink_comparison/reports/mike_example \
  --database validation/mixed_manifest/mixed_manifest.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --from-format TOP \
  --to-format PLUS
```

Add `--manifest-name <name>` when the source panel is known. Add
`--resolve-mixed-manifests --on-ambiguous-marker skip` for mixed-manifest
inputs.
