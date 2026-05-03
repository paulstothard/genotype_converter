#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${CONDA_ENV:-genotype-converter-env}"
SOURCE_ROOT="${SOURCE_ROOT:-validation/mixed_manifest/sources}"
FORMAT="${FORMAT:-table}"

usage() {
  cat <<'EOF'
Inspect mixed-manifest database source discovery without running any builds.

Usage:
  validation/mixed_manifest/scripts/run_discover_sources.sh [options]

Options:
  --source-root PATH   Source root. Default: validation/mixed_manifest/sources
  --format FORMAT      table, csv, or json. Default: table.
  --env NAME           Conda environment. Default: genotype-converter-env
  -h, --help           Show this help.

Environment overrides:
  CONDA_ENV, SOURCE_ROOT, FORMAT
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --source-root)
      SOURCE_ROOT="$2"
      shift 2
      ;;
    --format)
      FORMAT="$2"
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

conda run -n "${ENV_NAME}" genotype-converter db discover-sources \
  --source-root "${SOURCE_ROOT}" \
  --format "${FORMAT}"
