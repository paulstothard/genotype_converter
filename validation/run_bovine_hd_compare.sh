#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${CONDA_ENV:-genotype-converter-env}"

usage() {
  cat <<'EOF'
Compare the current full bovine HD validation outputs.

Usage:
  validation/run_bovine_hd_compare.sh [options]

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

conda run -n "${ENV_NAME}" python validation/compare_outputs.py
