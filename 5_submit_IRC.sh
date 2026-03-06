#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${SCRIPT_DIR}"
SYSTEM="${1:-}"
DRY_RUN="${DRY_RUN:-0}"
ORCA_CMD="${ORCA_CMD:-orca}"

discover_inputs() {
  local root="$1"
  local system="$2"
  if [[ -n "${system}" ]]; then
    find "${root}/${system}/IRC" -type f -name "*.inp" ! -name "*_IRC.inp" ! -name "*_TS_IRC.inp" 2>/dev/null || true
  else
    find "${root}" -type f -path "*/IRC/*.inp" ! -name "*_IRC.inp" ! -name "*_TS_IRC.inp" || true
  fi
}

mapfile -t INPUTS < <(discover_inputs "${ROOT}" "${SYSTEM}" | sort)

if [[ "${#INPUTS[@]}" -eq 0 ]]; then
  echo "No IRC input files found."
  exit 1
fi

echo "Note: ORCA IRC jobs may take several minutes per input."

for inp in "${INPUTS[@]}"; do
  inp_dir="$(dirname "${inp}")"
  inp_file="$(basename "${inp}")"
  out_file="${inp_file%.inp}.out"
  out="${inp_dir}/${out_file}"
  echo "Input: ${inp}"
  echo "Output: ${out}"
  echo "Command: (cd ${inp_dir} && module load orca && ${ORCA_CMD} ${inp_file} > ${out_file})"
  if [[ "${DRY_RUN}" == "1" ]]; then
    continue
  fi
  (
    cd "${inp_dir}"
    module load orca
    "${ORCA_CMD}" "${inp_file}" > "${out_file}"
  )
  echo "Run finished."
done
