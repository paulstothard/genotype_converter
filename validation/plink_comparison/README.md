# PLINK Conversion Comparison Workspace

Use this workspace for collaborator-provided PLINK examples. Keep the real
files local; inputs, converted outputs, and reports are ignored by Git.

## Folder Layout

```text
validation/plink_comparison/
  input/
    <example>/
      collaborator_top.bed
      collaborator_top.bim
      collaborator_top.fam
      collaborator_plus.bed
      collaborator_plus.bim
      collaborator_plus.fam
      notes.txt
  converted/
    genotype_converter PLUS outputs
  reports/
    comparison CSV and Markdown reports
  scripts/
    run_plink_comparison.sh
    compare_plink_bim.py
```

Record the collaborator source allele format and target reference assembly in
`notes.txt` or in the report notes.

## Compare Existing BIM Files

```bash
python validation/plink_comparison/scripts/compare_plink_bim.py \
  --collaborator-top validation/plink_comparison/input/collaborator_example/collaborator_top.bim \
  --collaborator-plus validation/plink_comparison/input/collaborator_example/collaborator_plus.bim \
  --genotype-converter-plus validation/plink_comparison/converted/collaborator_example.bim \
  --report-prefix validation/plink_comparison/reports/collaborator_example
```

The report classifies each marker as exact, allele-order swapped, complement,
swapped-complement, mismatch, missing from genotype_converter PLUS, or extra in
genotype_converter PLUS.

## Run Conversion And Compare

With a database:

```bash
validation/plink_comparison/scripts/run_plink_comparison.sh \
  --collaborator-top-bfile validation/plink_comparison/input/collaborator_example/collaborator_top \
  --collaborator-plus-bfile validation/plink_comparison/input/collaborator_example/collaborator_plus \
  --out-prefix validation/plink_comparison/converted/collaborator_example \
  --report-prefix validation/plink_comparison/reports/collaborator_example \
  --database validation/mixed_manifest/mixed_manifest.sqlite \
  --species bos_taurus \
  --assembly ARS_UCD_v2_0 \
  --from-format TOP \
  --to-format PLUS
```

Add `--manifest-name <name>` when the source panel is known. Add
`--resolve-mixed-manifests --on-ambiguous-marker skip` for mixed-manifest
inputs.
