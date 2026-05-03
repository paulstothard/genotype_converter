# Validation Scripts

Runnable full-panel validation helpers live here so the validation root can stay
focused on inputs, outputs, and reports.

```bash
validation/scripts/run_bovine_hd_build.sh
validation/scripts/run_bovine_hd_compare.sh
```

Do not run the full bovine HD build unless you intend to launch the large local
validation job. Use `--workers 1` unless the machine has memory for more.
