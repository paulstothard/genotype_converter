#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${CONDA_ENV:-genotype-converter-env}"
MANIFEST="${MANIFEST:-validation/data/bovinehd-manifest-b.csv}"
REFERENCE="${REFERENCE:-validation/data/ARS_UCD_v2.0.fa}"
OUTDIR="${OUTDIR:-validation/new_pipeline}"
SPECIES="${SPECIES:-bos_taurus}"
WORKERS="${WORKERS:-1}"
ALIGN_ARGS=()
PROGRESS_ARGS=(--progress)

usage() {
  cat <<'EOF'
Run the full bovine HD build with safe defaults.

Usage:
  validation/scripts/run_bovine_hd_build.sh [options]

Options:
  --align              Also write the large alignment display file.
  --no-progress        Disable build progress messages.
  --workers N          Worker processes. Default: 1.
  --manifest PATH      Manifest CSV. Default: validation/data/bovinehd-manifest-b.csv
  --reference PATH     Reference FASTA. Default: validation/data/ARS_UCD_v2.0.fa
  --outdir PATH        Output directory. Default: validation/new_pipeline
  --species NAME       Species key. Default: bos_taurus
  --env NAME           Conda environment. Default: genotype-converter-env
  -h, --help           Show this help.

Environment overrides:
  CONDA_ENV, MANIFEST, REFERENCE, OUTDIR, SPECIES, WORKERS
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --align)
      ALIGN_ARGS=(--align)
      shift
      ;;
    --no-progress)
      PROGRESS_ARGS=(--no-progress)
      shift
      ;;
    --workers)
      WORKERS="$2"
      shift 2
      ;;
    --manifest)
      MANIFEST="$2"
      shift 2
      ;;
    --reference)
      REFERENCE="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    --species)
      SPECIES="$2"
      shift 2
      ;;
    --env)
      ENV_NAME="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

echo "Running bovine HD build with workers=${WORKERS}"
conda run -n "${ENV_NAME}" genotype-converter build \
  --manifest "${MANIFEST}" \
  --reference "${REFERENCE}" \
  --outdir "${OUTDIR}" \
  --species "${SPECIES}" \
  --workers "${WORKERS}" \
  "${PROGRESS_ARGS[@]}" \
  "${ALIGN_ARGS[@]}"
