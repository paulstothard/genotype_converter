# Mixed-Manifest Scripts

Runnable helpers for the mixed-manifest validation workspace live here. The
parent folder holds local manifests, references, generated examples, converted
outputs, and reports.

Typical order:

```bash
validation/mixed_manifest/scripts/run_make_example_genotypes.sh
validation/mixed_manifest/scripts/run_discover_sources.sh
validation/mixed_manifest/scripts/run_database_build.sh --yes
validation/mixed_manifest/scripts/run_example_conversions.sh --overwrite
validation/mixed_manifest/scripts/run_check_example_conversions.sh
```

The database build can be large. Run discovery first, and use `--workers 1`
unless the machine has memory for more.
