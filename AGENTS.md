# Repository Notes for Agents

## Validation Safety

- Do not run the full bovine HD validation panel unless the user explicitly asks.
  The local reference FASTA is large, and `--align` creates a very large
  alignment display file.
- Prefer `--workers 1` for full validation unless the user confirms the machine
  has memory for one minimap2 reference index per worker.
- Use the validation helper scripts instead of reconstructing long commands by
  hand. The main scripts are:
  - `validation/scripts/run_bovine_hd_build.sh`
  - `validation/scripts/run_bovine_hd_compare.sh`
  - `validation/mixed_manifest/scripts/run_make_example_genotypes.sh`
  - `validation/mixed_manifest/scripts/run_discover_sources.sh`
  - `validation/mixed_manifest/scripts/run_database_build.sh`
  - `validation/mixed_manifest/scripts/run_example_conversions.sh`
  - `validation/mixed_manifest/scripts/run_check_example_conversions.sh`
  - `validation/plink_comparison/scripts/run_plink_comparison.sh`
  - `validation/plink_comparison/scripts/compare_plink_bim.py`
- Do not run `validation/scripts/run_bovine_hd_build.sh` or
  `validation/mixed_manifest/scripts/run_database_build.sh` unless the user explicitly
  asks for the large build. The mixed-manifest database build script also
  requires `--yes` to make accidental starts harder.
- Do not stage or commit large validation inputs or generated outputs under
  `validation/data/`, `validation/old_pipeline/`, or `validation/new_pipeline/`.
  These are intentionally ignored.
- Do not stage or commit collaborator PLINK files or generated reports under
  `validation/plink_comparison/input/`,
  `validation/plink_comparison/converted/`, or
  `validation/plink_comparison/reports/`.
- Do not change old pipeline outputs. If a comparison needs normalization,
  create a separate report or helper output and document the command.
- Do not change test expectations just to make a failing test pass. If expected
  output changes, explain the biological or file-format reason.

## Genotype Conversion Policy

- A lookup CSV is the `*.lookup.csv` output from one `build` run for one
  manifest/reference pair. `--lookup` points directly to that single file.
- SQLite database mode is for many manifests, many assemblies, or
  mixed-manifest genotype inputs. Keep lookup CSV workflows supported; the
  database must remain optional.
- `build` writes a site-only `*.sites.vcf` by default. `export-vcf` can also
  write one from a lookup CSV or one SQLite database source. This VCF has
  nucleotide `REF`/`ALT` alleles and no `FORMAT` or sample genotype columns.
  Do not add `0/1` genotype indexes unless a future workflow explicitly has
  sample genotypes to encode.
- Database source folders use `database_sources/<species>/manifests/` and
  `database_sources/<species>/references/<assembly>/`. Manifests are
  independent of reference genomes; the database builder pairs every manifest
  with every reference assembly for that species.
- `reference download-ncbi` is the user-facing reference download helper. It
  writes into `database_sources/<species>/references/<assembly>/`, preserves the
  assembly report, and keeps assembled molecules plus unlocalized/unplaced
  scaffolds by default.
- The same marker name can appear in multiple manifests. Do not collapse rules
  across manifest/source contexts, and keep mixed-manifest resolution
  conservative and reportable.
- `convert`, `convert-plink`, and `convert-pfile` all use
  `--on-unconvertible-marker exclude|fail|keep`.
- The default is `exclude` for text formats and PLINK formats.
- For text `convert` layouts, `exclude` removes unconvertible marker columns
  from wide CSV and marker rows from long, Illumina, and Affymetrix-style
  outputs. It writes `<output>.marker_conversion_report.csv` and, when markers
  are excluded, `<output>.exclude_markers.txt`.
- For `convert-plink`, `exclude` must use PLINK to remove variants before
  rewriting `.bim`, so `.bed/.bim/.fam` remain synchronized. Do not manually
  filter `.bim` without rewriting `.bed`.
- For `convert-pfile`, `exclude` must use PLINK2 to remove variants before
  rewriting `.pvar`, so `.pgen/.pvar/.psam` remain synchronized.
- Missing marker IDs and incomplete allele mappings are tracked separately.
  Text summaries use `markers_missing_lookup`, `markers_incomplete_mapping`,
  `markers_excluded`, and `genotypes_excluded`. PLINK summaries use
  `variants_missing_lookup`, `variants_incomplete_mapping`, and
  `variants_excluded`.
- `keep` is for audit/debugging runs where unresolved allele values should be
  left unchanged. `fail` is for strict checking after writing the marker report.
- The conda environment installs PLINK 1, not PLINK 2. `convert-pfile` should
  fail clearly when PLINK2 is needed for filtering and cannot be found.

## Documentation Hygiene

- Keep `README.md`, `docs/database.md`, `docs/genotype-formats.md`, and
  `docs/development-plan.md` aligned when behavior changes.
- The README should keep two distinct end-to-end examples: one for the single
  lookup CSV workflow and one for the SQLite database workflow.
- User-facing docs should explain workflows in terms of user decisions and
  inputs. Do not include internal slogans, conversation artifacts, or
  collaborator-specific names/details in tracked docs or helper scripts.
- Do not add performance/runtime tables unless they come from a reproducible
  benchmark or recorded validation run. Use qualitative memory and worker
  guidance instead.
- `docs/development-plan.md` records repo-level status and future work. Remove
  or update items when they are implemented so it does not become a stale TODO
  list.

## Probe-Adjacent SNP Positioning

For Illumina SNP manifests, `AlleleA_ProbeSeq` and `AlleleB_ProbeSeq` describe
the assay probe context. A probe may stop next to the assayed base, or it may
include the assayed allele as its terminal base. When probe placement is
detectable:

- A left-side probe reports the first reference base immediately after the probe
  in the alignment when the probe is adjacent-only.
- A right-side probe reports the first reference base immediately before the
  probe in the alignment when the probe is adjacent-only.
- An allele-including probe reports the matching allele base in the probe
  itself.
- If `AlleleA_ProbeSeq` and `AlleleB_ProbeSeq` imply incompatible sides or
  incompatible assayed positions, classify the site as ambiguous unless
  oriented manifest alleles uniquely resolve the candidate position.
- Reverse-complement probe matches are valid because the selected reference hit
  may be on the minus strand.
- When the selected hit is on the minus strand, orient the probe placement as a
  unit: flip the side and reverse-complement the probe sequence. Do not flip
  only the side label; that causes the probe lookup to fail and silently falls
  back to direct CIGAR mapping near gaps.
- This probe-adjacent rule is especially important when the synthetic `N` in
  `SourceSeq` is next to a CIGAR insertion/deletion. In those cases, direct
  `N`-to-reference CIGAR mapping can choose the wrong neighboring base for the
  assay.
- The synthetic `N` replacing the bracketed allele site must be tracked by
  position. Do not strip real flanking `N` bases from the minimap2 query; they
  are part of the manifest sequence context and should not be confused with the
  assayed-site placeholder.
- If the probe can be placed uniquely on the selected reference neighborhood,
  that reference-side probe placement defines the assayed base. Direct CIGAR
  mapping of the synthetic `N` is a fallback/check and must not overrule a
  clear reference-side probe placement.
- If only query-side probe placement is available and the probe-derived
  coordinate differs from the direct CIGAR coordinate, manifest alleles may be
  used as a tie-breaker only when exactly one candidate reference base is
  compatible with the oriented SNP alleles. If zero or multiple candidates are
  compatible, keep the probe-derived coordinate and classify the site as
  ambiguous.
- For SNPs with a CIGAR gap near the tracked assayed site, local realignment of
  a short reference window is allowed before probe/allele resolution. Keep this
  refinement SNP-only; indel placement and VCF anchoring have separate logic.
- If probe side cannot be determined, fall back to direct CIGAR mapping of the
  synthetic `N`.
- Build outputs include `determination_type` in position, wide, lookup, and
  optional Parquet outputs. Preserve this QC column when changing output
  schemas.
- Remaining nearby old/new position discrepancies in the bovine HD validation
  set are usually gap-adjacent. After the probe-orientation fixes, 150 of 165
  close discrepancies had a CIGAR insertion/deletion near the tracked assayed
  site in the selected minimap2 hit.

User-facing docs should describe this as probe-based or probe-adjacent SNP
positioning. Avoid internal slogans or implementation jargon in user-facing
docs.
