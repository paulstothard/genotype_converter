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
  scripts/
    run_bovine_hd_build.sh
    run_bovine_hd_compare.sh
    compare_outputs.py
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
validation/scripts/run_bovine_hd_build.sh
```

Use `--workers 1` explicitly for this full genome. Higher worker counts should
only be used when the machine has enough memory for one reference index per
worker. Add `--align` only when debugging a small subset; for the full bovine HD
panel it writes a very large alignment display file.

The helper script defaults to `--workers 1`, `--progress`, and the current local
bovine HD files. Use `validation/scripts/run_bovine_hd_build.sh --help` for
overrides.

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

For Illumina SNPs, pay special attention to gap-adjacent probe matches. Some
`AlleleA_ProbeSeq` and `AlleleB_ProbeSeq` values stop next to the assayed base;
others include the assayed allele as the terminal probe base. For adjacent
probes, a left-side probe reports the first reference base immediately after the
probe, and a right-side probe reports the first reference base immediately
before the probe. For allele-including probes, the reported SNP position is the
matching allele base in the probe itself. In local gaps or short repeats, this
probe-derived base can differ from the coordinate obtained by mapping the
synthetic `N` character directly through the CIGAR string.

The synthetic `N` that marks the assayed site is tracked separately from real
flanking `N` bases. Real flanking `N`s should remain in the minimap2 query so
the aligned sequence length and spacing reflect the manifest rather than an
artificially shortened sequence.

When the probe can be placed uniquely on the selected reference neighborhood,
that reference-side probe placement defines the assayed base. The direct
CIGAR-mapped position of the synthetic `N` is useful as a check, but it should
not overrule a clear probe placement. If only query-side probe placement is
available and it differs from direct CIGAR mapping, the pipeline uses the
manifest alleles as a conservative tie-breaker only when exactly one candidate
reference base matches those alleles. If both candidates match, or neither
candidate matches, the probe-derived position is retained and the site is
classified as ambiguous. Build summaries include alignment determination counts
such as `SNP_PROBE_ADJACENT`, `SNP_ALLELE_ASSISTED`, and `SNP_AMBIGUOUS`.

For SNPs with a CIGAR gap near the tracked assayed site, the pipeline performs a
short local realignment against the selected reference window before applying
the probe-adjacent and allele-assisted rules. This refinement is intended for
small gap/repeat placement differences and should not be expected to resolve
different-chromosome or far-distance alignment disagreements.

The build outputs include a `determination_type` column in `position.csv`,
`wide.csv`, `lookup.csv`, and optional `lookup.parquet`. Use this column to
filter markers for QC. In the bovine HD validation panel, the remaining nearby
old/new position discrepancies are usually gap-adjacent: a scan after the
probe-orientation fixes found 150 of 165 close discrepancies had a CIGAR
insertion/deletion near the tracked assayed site in the selected alignment.

The helper script compares the current bovine HD validation files:

```bash
validation/scripts/run_bovine_hd_compare.sh
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
