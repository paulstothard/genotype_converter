# genotype_converter — Development Notes

## What this is
A Python rewrite of `genotype_conversion_file_builder`. Uses minimap2 (via `mappy`)
instead of BLAST, Python instead of Perl, and `multiprocessing` instead of Nextflow.

## Project layout

```
src/genotype_converter/
  manifest.py          — parse Illumina + Affymetrix manifests → ManifestRecord
  aligner.py           — minimap2 alignment via mappy → AlignmentResult
  conversion.py        — strand/format conversion logic → VariantResult
  output.py            — write position/conversion/wide/lookup/alignment files
  pipeline.py          — orchestrate build (parse → align → convert → write)
  convert_genotypes.py — load lookup table + convert genotype files
  cli.py               — Click CLI (two subcommands: build, convert)

tests/
  data/                — self-contained test data (copied from old repo + 3 new indels)
    manifest.csv       — 5 Illumina SNPs + 3 synthetic indels
    reference.fa       — 1000 bp artificial reference
  conftest.py          — shared fixtures (manifest/reference paths, pipeline run)
  test_manifest.py
  test_aligner.py
  test_conversion.py
  test_convert.py
  test_integration.py
```

## Terminology (use consistently in docs and code)

| Term | Meaning |
|---|---|
| **build** | The pipeline phase that aligns a manifest to a reference genome and produces lookup/conversion files |
| **convert** | The phase that takes an existing genotype data file and rewrites allele calls from one encoding to another using a lookup file |
| **lookup file** | `*.lookup.csv` — the primary distributable output; one row per variant, self-describing column names (`A_in_TOP`, `B_in_PLUS`, etc.) |
| **conversion file** | `*.conversion.csv` — legacy two-row-per-variant format kept for backward compatibility |
| **position file** | `*.position.csv` — chromosome, position, VCF REF/ALT per variant |

## How to run

```bash
# activate environment
conda activate python-env   # Python 3.12

# install editable
pip install -e ".[dev]"

# run tests
pytest tests/ -v

# build conversion files
genotype-converter build \
  --manifest path/to/manifest.csv \
  --reference path/to/genome.fa \
  --outdir output/ --species cattle

# convert a genotype file from TOP to PLUS format
genotype-converter convert \
  --genotypes mydata.csv \
  --lookup output/cattle/genome/manifest.genome.lookup.csv \
  --from-format TOP --to-format PLUS \
  --output converted.csv
```

## Output files produced by `build`

| File | Description |
|---|---|
| `*.position.csv` | Chromosome + 1-based position + VCF REF/ALT per variant |
| `*.conversion.csv` | Two rows per variant: AB/TOP/FORWARD/DESIGN/PLUS/VCF for allele A and B |
| `*.wide.csv` | One row per variant with all of the above |
| `*.lookup.csv` | **Primary distributable** — one row per variant, columns named `A_in_TOP`, `B_in_PLUS`, etc. |
| `*.alignment.txt` | Debug alignment display (`--align` flag) |
| `*.lookup.parquet` | Binary version of lookup.csv (`--parquet` flag, requires pyarrow) |

## Current status (as of 2026-04-30)

### Completed this session
- [x] **Proper indel VCF representation** — `_resolve_indel_vcf` in `aligner.py` now produces
  correct anchor-based REF/ALT rather than the placeholder `ref_base/ref_base` the original used.
  - Insertion `[-/SEQ]`: REF = anchor, ALT = anchor + SEQ, position = anchor position
  - Deletion `[SEQ/-]`: N aligns to the *last* deleted base (CIGAR: `left_len M (del_len-1) D 1M`);
    true anchor = ref_pos_0 - del_len; deleted bases read from reference for accuracy
  - `indel_ref_is_del` field added to `AlignmentResult` so `_set_vcf_ab` can correctly assign
    REF/ALT to allele A vs B
  - Fallback to old behaviour for symbolic `[I/D]` notation (no actual sequence in flanking)
- [x] **Test data now self-contained** — copied to `tests/data/`, no longer depends on sibling repo
- [x] **3 synthetic indels added** to test manifest and verified:
  - INDEL1: insertion `[-/ACGT]`, `SNP=[I/D]` → pos 449, REF=A, ALT=AACGT; A=ALT, B=REF
  - INDEL2: deletion `[TCGA/-]`, `SNP=[I/D]` → pos 325, REF=GTCGA, ALT=G; A=REF, B=ALT
  - INDEL3: insertion `[-/TTCC]`, `SNP=[D/I]` → pos 649, REF=C, ALT=CTTCC; A=REF, B=ALT
- [x] **40 tests passing** (unit + integration, including 7 new indel-specific tests)

### Remaining tasks (priority order)

#### 1. Documentation — README rewrite
Write a proper README covering:
- Brief explanation of how the program works end-to-end (manifest → flanking alignment →
  variant site location → strand/format conversion)
- Explain the minimap2 alignment approach: how CIGAR walking locates the variant within the
  flanking sequence, how gaps/mismatches in the flanking are handled, what happens when a
  variant aligns to multiple locations (currently: take highest mapq hit; mapq=0 means ambiguous)
- Explain both subcommands (`build` and `convert`) with concrete examples
- Document all supported genotype input formats for `convert` (wide/long layouts,
  allele separator variants). Clarify what column names are expected.
- Real-world performance note: minimap2 is much faster than BLAST; a 50K SNP panel
  takes a few minutes, an 800K HD panel ~15–30 min depending on CPU count

#### 2. Large-scale validation
Run `build` on one of Paul's full real manifests + the corresponding reference genome,
then diff the position/conversion output against the old Perl pipeline's output.
The SNPchimp comparison data in `../genotype_conversion_file_builder/test/comparing_results_with_snpchimp/`
is a useful reference. This will expose any systematic discrepancies.

#### 3. Multimapping and alignment failure handling
- Clarify and test what happens to variants that don't align (currently: all fields
  are None → they appear in output with empty position/VCF columns; consider skipping
  or flagging them explicitly)
- Clarify what happens when mapq=0 (ambiguous/multimapping): currently we take the
  highest-mapq hit regardless; consider warning or marking these as unreliable
- Add tests for these edge cases

#### 4. Batch `convert` — multiple files / directories
The `convert` subcommand currently processes one file at a time. Add support for:
- Passing multiple `--genotypes` files in one invocation
- Passing a directory and processing all matching files within it

#### 5. Error messages
Add informative error messages for common failure cases, e.g.:
- Manifest file not recognised (neither Illumina nor Affymetrix format)
- Required columns missing from a genotype file
- Lookup file doesn't contain the expected marker names
- Reference FASTA not indexed or empty

#### 6. "Parse, don't validate" — data structure improvements
Where the code currently does defensive checks on raw strings (e.g., checking whether
a field is "." or empty), prefer returning typed values or raising at parse time so
that downstream code works with guaranteed-valid objects. Use enums or dataclasses
with required fields rather than Optional everywhere.

#### 7. Dockerfile
Write a Dockerfile (the original pipeline had one). Should install `mappy`, `click`,
and optionally `pyarrow`, expose the `genotype-converter` CLI.

#### 8. Code quality review
The new code was written from scratch in Python, not translated from Perl. However,
worth a pass to confirm: no Perl-isms crept in, idiomatic use of dataclasses and typing,
no unnecessary globals beyond the multiprocessing worker state, clean separation of
concerns between modules.

## Key algorithm notes (for README)

### How variant site location works
1. The flanking sequence from the manifest (e.g. `AAACCC[A/G]TTTGGG`) is prepared as a
   minimap2 query by replacing the `[X/Y]` bracket with a single `N` and removing
   non-GATCN characters.
2. The query is aligned to the reference with minimap2 (`sr` preset, `best_n=5`).
3. The position of `N` in the query is mapped to a reference position by walking the
   CIGAR string (`_cigar_query_to_ref` in `aligner.py`). Gaps in the flanking relative
   to the reference are handled: insertions in the query report the left-adjacent
   reference position; deletions in the query advance the reference counter silently.
4. Mismatches in the flanking are fine — minimap2 handles them and the CIGAR walk still
   finds the right position.
5. For indels, the position is the VCF anchor base (one before the indel site).

### Multimapping
When multiple hits are returned, we take the one with the highest `mapq`. A `mapq` of 0
means the aligner considers multiple equally-good mappings possible — these are currently
kept but not flagged. A future improvement would be to report them as uncertain.

## Environment
- Developed/tested with: conda env `python-env` (Python 3.12)
- Release prep: create and validate a project-specific conda env named
  `genotype-converter`, then update docs/commands away from the generic
  `python-env` name.
- Runtime deps: `mappy>=2.24`, `click>=8.0`
- Optional: `pyarrow` for Parquet output
- Dev deps: `pytest>=7.0`, `pytest-cov`
