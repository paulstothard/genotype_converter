# Database Source Fixture

This folder is a small, safe fixture for future SQLite conversion-database
tests. It mirrors the proposed user-facing source-folder layout without using
large validation files.

Expected layout:

```text
database_sources/
  bos_taurus/
    manifests/
      tiny_bovine_manifest.csv
    references/
      ARS_UCD_v2_0/
        tiny_reference.fa
      ARS_UCD1_2/
        tiny_reference.fa
    genotypes/
      tiny_top_wide.csv
      tiny_top_long.csv
    expected/
      tiny_lookup_markers.txt
```

Future database tests should use this fixture to verify that a database builder
can discover species-level manifests and pair them with each reference assembly
before using any full validation data.

The fixture intentionally uses tiny synthetic test data. Do not replace these
files with real panel-scale manifests or generated validation output.
