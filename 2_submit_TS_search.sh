#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${SCRIPT_DIR}"
SYSTEM="${1:-}"
DRY_RUN="${DRY_RUN:-0}"
# Use absolute path for ORCA as required by parallel runs
ORCA_CMD="/software/orca-6.1.0-el8-x86_64/bin/orca"

# Aggressively suppress OpenMPI OpenFabrics warnings and aggregation on login nodes
export OMPI_MCA_btl_base_warn_component_unused=0
export OMPI_MCA_btl="self,vader,tcp"
export OMPI_MCA_pml="ob1"
export OMPI_MCA_btl_openib_warn_no_device_params_found=0
export OMPI_MCA_orte_base_help_aggregate=0
export OMPI_MCA_btl_openib_allow_ib=0

discover_inputs() {
  local root="$1"
  local system="$2"
  local find_cmd
  if [[ -n "${system}" ]]; then
    find_cmd="find ${root}/${system}/TS_search -type f -name \"*.inp\""
  else
    find_cmd="find ${root} -type f -path \"*/TS_search/*.inp\""
  fi
  
  # Filter out internal ORCA atom-wise inputs (e.g. *_atom53.inp)
  eval "${find_cmd}" | grep -v "_atom" || true
}

is_completed() {
  local out_file="$1"
  if [[ ! -f "${out_file}" ]]; then
    return 1
  fi
  if command -v rg >/dev/null 2>&1; then
    if rg -q "ORCA TERMINATED NORMALLY" "${out_file}"; then
      return 0
    fi
  else
    if grep -q "ORCA TERMINATED NORMALLY" "${out_file}"; then
      return 0
    fi
  fi
  return 1
}

mapfile -t INPUTS < <(discover_inputs "${ROOT}" "${SYSTEM}" | sort)

if [[ "${#INPUTS[@]}" -eq 0 ]]; then
  echo "No TS-search input files found."
  exit 1
fi

echo "Note: ORCA TS-search jobs may take several minutes per input."
echo "Running jobs SEQUENTIALLY to avoid login node overload."

for inp in "${INPUTS[@]}"; do
  # Ensure we have absolute path for the input
  abs_inp=$(realpath "${inp}")
  inp_dir=$(dirname "${abs_inp}")
  inp_base=$(basename "${abs_inp}" .inp)
  out="${inp_dir}/${inp_base}.out"

  echo "---------------------------------------------------"
  echo "Input: ${abs_inp}"
  echo "Output: ${out}"
  
  if is_completed "${out}"; then
    echo "Status: already completed (ORCA TERMINATED NORMALLY). Skipping."
    continue
  fi

  if [[ "${DRY_RUN}" == "1" ]]; then
    echo "Command: cd ${inp_dir} && module load orca && ${ORCA_CMD} ${inp_base}.inp > ${inp_base}.out"
    continue
  fi

  # Clean start: remove old output and temporary files in the target directory
  echo "Cleaning and starting ORCA..."
  rm -f "${out}"
  rm -f "${inp_dir}/${inp_base}"*.tmp "${inp_dir}/${inp_base}"*.gbw "${inp_dir}/${inp_base}"*.hess

  (
    cd "${inp_dir}"
    module load orca
    "${ORCA_CMD}" "${inp_base}.inp" > "${inp_base}.out"
  )
  
  if is_completed "${out}"; then
    echo "Run finished successfully: ${out}"
  else
    echo "Run finished with errors or incomplete: ${out}"
    # Check if we should stop on error
    # exit 1 
  fi
done

echo "All jobs processed."
