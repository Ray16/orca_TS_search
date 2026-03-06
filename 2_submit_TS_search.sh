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
    find "${root}/${system}/TS_search" -type f -name "*.inp" 2>/dev/null || true
  else
    find "${root}" -type f -path "*/TS_search/*.inp" || true
  fi
}

mapfile -t INPUTS < <(discover_inputs "${ROOT}" "${SYSTEM}" | sort)

if [[ "${#INPUTS[@]}" -eq 0 ]]; then
  echo "No TS-search input files found."
  exit 1
fi

echo "Note: ORCA TS-search jobs may take several minutes per input."

for inp in "${INPUTS[@]}"; do
  out="${inp%.inp}.out"
  echo "Input: ${inp}"
  echo "Output: ${out}"
  echo "Command: module load orca && ${ORCA_CMD} ${inp} > ${out}"
  if [[ "${DRY_RUN}" == "1" ]]; then
    continue
  fi
  module load orca
  "${ORCA_CMD}" "${inp}" > "${out}"
  echo "Run finished."
done
