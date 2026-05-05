#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${CONDA_ENV:-genotype-converter-env}"
MIKE_TOP_BFILE=""
MIKE_PLUS_BFILE=""
OUT_PREFIX=""
REPORT_PREFIX=""
LOOKUP=""
DATABASE=""
SPECIES=""
ASSEMBLY=""
MANIFEST_NAME=""
FROM_FORMAT=""
TO_FORMAT="PLUS"
RESOLVE_ARGS=()
ON_AMBIGUOUS="fail"
UPDATE_POSITION_ARGS=()
ON_UNCONVERTIBLE="exclude"
COMPARE_ARGS=()

usage() {
  cat <<'EOF'
Run our PLINK conversion and compare it with Mike's PLINK PLUS file.

Usage:
  validation/plink_comparison/scripts/run_plink_comparison.sh [options]

Required:
  --mike-top-bfile PREFIX   Mike's TOP PLINK fileset prefix.
  --mike-plus-bfile PREFIX  Mike's PLUS PLINK fileset prefix.
  --out-prefix PREFIX       Output prefix for this converter.
  --report-prefix PREFIX    Report prefix for comparison CSV/Markdown.
  --from-format FORMAT      Mike TOP allele format.

Lookup source:
  --lookup PATH             Lookup CSV from build.
  --database PATH           SQLite database.
  --species NAME            Species for database mode.
  --assembly NAME           Assembly for database mode.
  --manifest-name NAME      Optional manifest/panel name for database mode.

Options:
  --to-format FORMAT        Target allele format. Default: PLUS.
  --resolve-mixed-manifests Resolve database rules per marker.
  --on-ambiguous-marker MODE  fail or skip. Default: fail.
  --update-position         Update BIM chromosome/base-pair columns.
  --on-unconvertible-marker MODE
                          exclude, fail, or keep. Default: exclude.
  --fail-on-differences   Return nonzero if the comparison finds differences.
  --env NAME                Conda environment. Default: genotype-converter-env
  -h, --help                Show this help.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mike-top-bfile) MIKE_TOP_BFILE="$2"; shift 2 ;;
    --mike-plus-bfile) MIKE_PLUS_BFILE="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    --report-prefix) REPORT_PREFIX="$2"; shift 2 ;;
    --lookup) LOOKUP="$2"; shift 2 ;;
    --database) DATABASE="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --assembly) ASSEMBLY="$2"; shift 2 ;;
    --manifest-name) MANIFEST_NAME="$2"; shift 2 ;;
    --from-format) FROM_FORMAT="$2"; shift 2 ;;
    --to-format) TO_FORMAT="$2"; shift 2 ;;
    --resolve-mixed-manifests) RESOLVE_ARGS=(--resolve-mixed-manifests); shift ;;
    --on-ambiguous-marker) ON_AMBIGUOUS="$2"; shift 2 ;;
    --update-position) UPDATE_POSITION_ARGS=(--update-position); shift ;;
    --on-unconvertible-marker) ON_UNCONVERTIBLE="$2"; shift 2 ;;
    --fail-on-differences) COMPARE_ARGS=(--fail-on-differences); shift ;;
    --env) ENV_NAME="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for value_name in MIKE_TOP_BFILE MIKE_PLUS_BFILE OUT_PREFIX REPORT_PREFIX FROM_FORMAT; do
  if [[ -z "${!value_name}" ]]; then
    echo "Missing required option for ${value_name}" >&2
    usage >&2
    exit 2
  fi
done

lookup_args=()
if [[ -n "${LOOKUP}" ]]; then
  lookup_args=(--lookup "${LOOKUP}")
elif [[ -n "${DATABASE}" ]]; then
  lookup_args=(--database "${DATABASE}" --species "${SPECIES}" --assembly "${ASSEMBLY}")
  if [[ -n "${MANIFEST_NAME}" ]]; then
    lookup_args+=(--manifest-name "${MANIFEST_NAME}")
  fi
  lookup_args+=("${RESOLVE_ARGS[@]}" --on-ambiguous-marker "${ON_AMBIGUOUS}")
else
  echo "Provide --lookup or --database." >&2
  exit 2
fi

mkdir -p "$(dirname "${OUT_PREFIX}")" "$(dirname "${REPORT_PREFIX}")"

conda run -n "${ENV_NAME}" genotype-converter convert-plink \
  --bfile "${MIKE_TOP_BFILE}" \
  "${lookup_args[@]}" \
  --from-format "${FROM_FORMAT}" \
  --to-format "${TO_FORMAT}" \
  --out "${OUT_PREFIX}" \
  "${UPDATE_POSITION_ARGS[@]}" \
  --on-unconvertible-marker "${ON_UNCONVERTIBLE}" \
  --overwrite

conda run -n "${ENV_NAME}" python validation/plink_comparison/scripts/compare_plink_bim.py \
  --mike-top "${MIKE_TOP_BFILE}.bim" \
  --mike-plus "${MIKE_PLUS_BFILE}.bim" \
  --genotype-converter-plus "${OUT_PREFIX}.bim" \
  --report-prefix "${REPORT_PREFIX}" \
  "${COMPARE_ARGS[@]}"
