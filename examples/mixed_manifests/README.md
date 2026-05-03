# Mixed-Manifest Database Conversion Example

This tiny example shows how to convert one genotype file containing markers
from more than one manifest. It is intentionally small and synthetic.

Build a database from the two lookup files:

```bash
genotype-converter db init --database examples/mixed_manifests/example.sqlite

genotype-converter db import-lookup \
  --database examples/mixed_manifests/example.sqlite \
  --lookup examples/mixed_manifests/panel_a.lookup.csv \
  --species bos_taurus \
  --assembly demo_assembly \
  --manifest-name panel_a

genotype-converter db import-lookup \
  --database examples/mixed_manifests/example.sqlite \
  --lookup examples/mixed_manifests/panel_b.lookup.csv \
  --species bos_taurus \
  --assembly demo_assembly \
  --manifest-name panel_b
```

Convert with per-marker mixed-manifest resolution:

```bash
genotype-converter convert \
  --genotypes examples/mixed_manifests/mixed_genotypes.csv \
  --database examples/mixed_manifests/example.sqlite \
  --species bos_taurus \
  --assembly demo_assembly \
  --resolve-mixed-manifests \
  --on-ambiguous-marker skip \
  --from-format TOP \
  --to-format PLUS \
  --output examples/mixed_manifests/mixed_genotypes.plus.csv \
  --resolution-report examples/mixed_manifests/manifest_resolution.csv
```

`manifest_resolution.csv` records the input path, marker name, selected
manifest, candidate manifests, candidate count, and reason for each decision.
If a shared marker has conflicting rules and local marker context does not
clearly support one manifest, `--on-ambiguous-marker skip` leaves that marker
unchanged. The default is `--on-ambiguous-marker fail`.
