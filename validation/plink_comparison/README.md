# PLINK Conversion Comparison Workspace

Use this workspace for collaborator-provided PLINK examples. Keep the real
files local; inputs, converted outputs, and reports are ignored by Git.

## Folder Layout

```text
validation/plink_comparison/
  input/
    <example>/
      mike_top.bed
      mike_top.bim
      mike_top.fam
      mike_plus.bed
      mike_plus.bim
      mike_plus.fam
      notes.txt
  converted/
    genotype_converter PLUS outputs
  reports/
    comparison CSV and Markdown reports
  scripts/
    run_plink_comparison.sh
    compare_plink_bim.py
```

Record Mike's source allele format and the target reference assembly in
`notes.txt` or in the report notes.

## Compare Existing BIM Files

```bash
python validation/plink_comparison/scripts/compare_plink_bim.py \
  --mike-top validation/plink_comparison/input/mike_example/mike_top.bim \
  --mike-plus validation/plink_comparison/input/mike_example/mike_plus.bim \
  --genotype-converter-plus validation/plink_comparison/converted/mike_example.bim \
  --report-prefix validation/plink_comparison/reports/mike_example
```

The report classifies each marker as exact, allele-order swapped, complement,
swapped-complement, mismatch, missing from genotype_converter PLUS, or extra in
genotype_converter PLUS.

## Run Conversion And Compare

With a database:

```bash
validation/plink_comparison/scripts/run_plink_comparison.sh \
  --mike-top-bfile validation/plink_comparison/input/mike_example/mike_top \
  --mike-plus-bfile validation/plink_comparison/input/mike_example/mike_plus \
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
