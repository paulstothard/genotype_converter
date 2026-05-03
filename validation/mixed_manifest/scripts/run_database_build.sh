#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${CONDA_ENV:-genotype-converter-env}"
SOURCE_ROOT="${SOURCE_ROOT:-validation/mixed_manifest/sources}"
DATABASE="${DATABASE:-validation/mixed_manifest/mixed_manifest.sqlite}"
BUILD_OUTDIR="${BUILD_OUTDIR:-validation/mixed_manifest/database_build}"
WORKERS="${WORKERS:-1}"
REPLACE_ARGS=()
PROGRESS_ARGS=(--progress)
YES=0

usage() {
  cat <<'EOF'
Build the mixed-manifest SQLite database from local source folders.

This can be a very large run. It builds every manifest/reference pair found
under the source root. Use the discover script first to check the workload.

Usage:
  validation/mixed_manifest/scripts/run_database_build.sh --yes [options]

Options:
  --yes                Required confirmation to run the build.
  --replace            Replace existing identical lookup imports.
  --no-progress        Disable progress messages.
  --workers N          Worker processes. Default: 1.
  --source-root PATH   Source root. Default: validation/mixed_manifest/sources
  --database PATH      SQLite path. Default: validation/mixed_manifest/mixed_manifest.sqlite
  --build-outdir PATH  Build output directory. Default: validation/mixed_manifest/database_build
  --env NAME           Conda environment. Default: genotype-converter-env
  -h, --help           Show this help.

Environment overrides:
  CONDA_ENV, SOURCE_ROOT, DATABASE, BUILD_OUTDIR, WORKERS
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --yes)
      YES=1
      shift
      ;;
    --replace)
      REPLACE_ARGS=(--replace)
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
    --source-root)
      SOURCE_ROOT="$2"
      shift 2
      ;;
    --database)
      DATABASE="$2"
      shift 2
      ;;
    --build-outdir)
      BUILD_OUTDIR="$2"
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

if [[ "${YES}" -ne 1 ]]; then
  echo "Refusing to start a large database build without --yes." >&2
  echo "Run validation/mixed_manifest/scripts/run_discover_sources.sh first to inspect the workload." >&2
  exit 2
fi

echo "Running mixed-manifest database build with workers=${WORKERS}"
conda run -n "${ENV_NAME}" genotype-converter db build \
  --source-root "${SOURCE_ROOT}" \
  --database "${DATABASE}" \
  --build-outdir "${BUILD_OUTDIR}" \
  --workers "${WORKERS}" \
  "${PROGRESS_ARGS[@]}" \
  "${REPLACE_ARGS[@]}"
