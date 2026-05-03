#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${CONDA_ENV:-genotype-converter-env}"
DATABASE="${DATABASE:-validation/mixed_manifest/mixed_manifest.sqlite}"
SOURCE_ROOT="${SOURCE_ROOT:-validation/mixed_manifest/sources}"
OUT_ROOT="${OUT_ROOT:-validation/mixed_manifest/converted/example_conversions}"
REPORT_ROOT="${REPORT_ROOT:-validation/mixed_manifest/reports/example_conversions}"
TO_FORMAT="${TO_FORMAT:-PLUS}"
ON_AMBIGUOUS="${ON_AMBIGUOUS:-skip}"
OVERWRITE=0

usage() {
  cat <<'EOF'
Run database-backed conversions for the generated mixed-manifest genotype examples.

Usage:
  validation/mixed_manifest/run_example_conversions.sh [options]

Options:
  --overwrite          Replace previous example conversion outputs.
  --database PATH      SQLite database. Default: validation/mixed_manifest/mixed_manifest.sqlite
  --source-root PATH   Source root. Default: validation/mixed_manifest/sources
  --out-root PATH      Output root. Default: validation/mixed_manifest/converted/example_conversions
  --report-root PATH   Report root. Default: validation/mixed_manifest/reports/example_conversions
  --to-format FORMAT   Target encoding. Default: PLUS.
  --on-ambiguous MODE  fail or skip. Default: skip.
  --env NAME           Conda environment. Default: genotype-converter-env
  -h, --help           Show this help.

Environment overrides:
  CONDA_ENV, DATABASE, SOURCE_ROOT, OUT_ROOT, REPORT_ROOT, TO_FORMAT, ON_AMBIGUOUS
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --overwrite)
      OVERWRITE=1
      shift
      ;;
    --database)
      DATABASE="$2"
      shift 2
      ;;
    --source-root)
      SOURCE_ROOT="$2"
      shift 2
      ;;
    --out-root)
      OUT_ROOT="$2"
      shift 2
      ;;
    --report-root)
      REPORT_ROOT="$2"
      shift 2
      ;;
    --to-format)
      TO_FORMAT="$2"
      shift 2
      ;;
    --on-ambiguous)
      ON_AMBIGUOUS="$2"
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

if [[ ! -f "${DATABASE}" ]]; then
  echo "Database not found: ${DATABASE}" >&2
  echo "Run validation/mixed_manifest/run_database_build.sh --yes first." >&2
  exit 2
fi

overwrite_args=()
if [[ "${OVERWRITE}" -eq 1 ]]; then
  overwrite_args=(--overwrite)
fi

convert_text() {
  local species="$1"
  local assembly="$2"
  local input="$3"
  local layout="$4"
  local output="$5"
  local report="$6"
  local out_sep="$7"

  conda run -n "${ENV_NAME}" genotype-converter convert \
    --genotypes "${input}" \
    --database "${DATABASE}" \
    --species "${species}" \
    --assembly "${assembly}" \
    --resolve-mixed-manifests \
    --on-ambiguous-marker "${ON_AMBIGUOUS}" \
    --resolution-report "${report}" \
    --from-format AB \
    --to-format "${TO_FORMAT}" \
    --layout "${layout}" \
    --out-sep "${out_sep}" \
    --output "${output}"
}

for species_dir in "${SOURCE_ROOT}"/*; do
  [[ -d "${species_dir}" ]] || continue
  species="$(basename "${species_dir}")"
  genotype_dir="${species_dir}/genotypes/synthetic_mixed_manifest"
  [[ -d "${genotype_dir}" ]] || continue

  for reference_dir in "${species_dir}/references"/*; do
    [[ -d "${reference_dir}" ]] || continue
    assembly="$(basename "${reference_dir}")"
    out_dir="${OUT_ROOT}/${species}/${assembly}"
    report_dir="${REPORT_ROOT}/${species}/${assembly}"
    mkdir -p "${out_dir}" "${report_dir}"

    echo "Converting ${species}/${assembly}"
    convert_text "${species}" "${assembly}" \
      "${genotype_dir}/mixed_manifest_wide_ab.csv" "wide" \
      "${out_dir}/mixed_manifest_wide_plus.csv" \
      "${report_dir}/mixed_manifest_wide_resolution.csv" "/"
    convert_text "${species}" "${assembly}" \
      "${genotype_dir}/mixed_manifest_long_ab.csv" "long" \
      "${out_dir}/mixed_manifest_long_plus.csv" \
      "${report_dir}/mixed_manifest_long_resolution.csv" "/"
    convert_text "${species}" "${assembly}" \
      "${genotype_dir}/illumina_gsgt_matrix_ab.txt" "illumina-matrix" \
      "${out_dir}/illumina_gsgt_matrix_plus.txt" \
      "${report_dir}/illumina_gsgt_matrix_resolution.csv" ""
    convert_text "${species}" "${assembly}" \
      "${genotype_dir}/illumina_gsgt_long_multiformat.txt" "illumina-long" \
      "${out_dir}/illumina_gsgt_long_plus.txt" \
      "${report_dir}/illumina_gsgt_long_resolution.csv" "/"
    convert_text "${species}" "${assembly}" \
      "${genotype_dir}/affymetrix_axiom_dual_call_matrix.txt" "affymetrix-matrix" \
      "${out_dir}/affymetrix_axiom_dual_call_matrix_plus.txt" \
      "${report_dir}/affymetrix_axiom_dual_call_matrix_resolution.csv" "/"

    conda run -n "${ENV_NAME}" genotype-converter convert-plink \
      --bfile "${genotype_dir}/mixed_manifest_plink1_ab" \
      --database "${DATABASE}" \
      --species "${species}" \
      --assembly "${assembly}" \
      --resolve-mixed-manifests \
      --on-ambiguous-marker "${ON_AMBIGUOUS}" \
      --resolution-report "${report_dir}/mixed_manifest_plink1_resolution.csv" \
      --from-format AB \
      --to-format "${TO_FORMAT}" \
      --out "${out_dir}/mixed_manifest_plink1_plus" \
      "${overwrite_args[@]}"

    conda run -n "${ENV_NAME}" genotype-converter convert-pfile \
      --pfile "${genotype_dir}/mixed_manifest_plink2_ab" \
      --database "${DATABASE}" \
      --species "${species}" \
      --assembly "${assembly}" \
      --resolve-mixed-manifests \
      --on-ambiguous-marker "${ON_AMBIGUOUS}" \
      --resolution-report "${report_dir}/mixed_manifest_plink2_resolution.csv" \
      --from-format AB \
      --to-format "${TO_FORMAT}" \
      --out "${out_dir}/mixed_manifest_plink2_plus" \
      "${overwrite_args[@]}"
  done
done
