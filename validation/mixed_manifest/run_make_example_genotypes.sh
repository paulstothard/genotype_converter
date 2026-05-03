#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${CONDA_ENV:-genotype-converter-env}"

usage() {
  cat <<'EOF'
Generate small synthetic genotype examples from the local mixed-manifest folders.

The script scans each species manifest folder, selects shared and unique
markers, and writes AB-coded wide CSV, long CSV, PLINK 1, and PLINK 2 examples
under sources/<species>/genotypes/synthetic_mixed_manifest/.

Usage:
  validation/mixed_manifest/run_make_example_genotypes.sh [options]

Options:
  --env NAME    Conda environment. Default: genotype-converter-env
  -h, --help    Show this help.

Environment overrides:
  CONDA_ENV
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
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

conda run -n "${ENV_NAME}" python validation/mixed_manifest/make_example_genotypes.py
