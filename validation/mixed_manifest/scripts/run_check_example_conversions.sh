#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${CONDA_ENV:-genotype-converter-env}"

usage() {
  cat <<'EOF'
Check mixed-manifest example conversion outputs.

Usage:
  validation/mixed_manifest/scripts/run_check_example_conversions.sh [options]

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

conda run -n "${ENV_NAME}" python validation/mixed_manifest/scripts/check_example_conversions.py
