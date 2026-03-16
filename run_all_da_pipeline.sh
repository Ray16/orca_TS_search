#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS_DIR="$SCRIPT_DIR/runs"
cd "$SCRIPT_DIR"

export OMPI_MCA_btl_base_warn_component_unused=0
export OMPI_MCA_btl="self,vader,tcp"
export OMPI_MCA_pml="ob1"
export OMPI_MCA_btl_openib_warn_no_device_params_found=0
export OMPI_MCA_orte_base_help_aggregate=0
export OMPI_MCA_btl_openib_allow_ib=0
module load orca 2>/dev/null || true

ORCA_CMD="/software/orca-6.1.0-el8-x86_64/bin/orca"

SYSTEMS=(
  da_f da_cl da_br da_i
)

FAILED=()
SUCCEEDED=()

# Step 0: Generate TS guess XYZ files
echo "=== Step 0: Generating Diels-Alder TS guess XYZ files ==="
python3 0_prep_initial_TS_structure.py --family da --overwrite

for sys in "${SYSTEMS[@]}"; do
  echo "========================================"
  echo "System: $sys"
  echo "========================================"

  ts_dir="$RUNS_DIR/$sys/TS_search"
  irc_dir="$RUNS_DIR/$sys/IRC"
  rp_dir="$RUNS_DIR/$sys/geo_opt_reactant_product"

  # Clean previous outputs
  rm -f "$ts_dir/$sys".out "$ts_dir/$sys"*.tmp "$ts_dir/$sys"*.gbw "$ts_dir/$sys"*.hess \
        "$ts_dir/$sys"*.xyz "$ts_dir/$sys"*.opt "$ts_dir/$sys"*.engrad \
        "$ts_dir/$sys"*.densities "$ts_dir/$sys"*.densitiesinfo \
        "$ts_dir/$sys"*.property.txt "$ts_dir/$sys"*.bibtex "$ts_dir/$sys"*.carthess \
        "$ts_dir/${sys}_atom"* 2>/dev/null || true
  rm -rf "$irc_dir" "$rp_dir" 2>/dev/null || true

  # Step 1: Generate ORCA input from TS guess (neutral, singlet)
  echo "  Step 1a: Preparing TS search input..."
  python3 1_prep_TS_search.py "$sys" --charge 0 --mult 1 \
    --root "$RUNS_DIR" \
    --xyz-file "$SCRIPT_DIR/TS_guess_xyz/${sys}_input.xyz"
  if [[ ! -f "$ts_dir/$sys.inp" ]]; then
    echo "  FAILED: TS input not generated for $sys"
    FAILED+=("$sys:no_inp")
    continue
  fi

  # Step 1b: OptTS + Freq
  echo "  Step 1b: OptTS + Freq..."
  (cd "$ts_dir" && "$ORCA_CMD" "$sys.inp" > "$sys.out" 2>&1)

  if ! grep -q "HURRAY" "$ts_dir/$sys.out" 2>/dev/null; then
    echo "  FAILED: TS optimization did not converge for $sys"
    FAILED+=("$sys:TS_no_converge")
    continue
  fi

  imag_freq=$(grep "imaginary mode" "$ts_dir/$sys.out" | tail -1 | grep -oP '[-]\d+\.\d+')
  echo "  TS converged. Imaginary freq: $imag_freq cm^-1"
  python3 3_frequency_analysis.py "$sys" 2>/dev/null || true

  # Step 2: IRC
  echo "  Step 2: IRC..."
  python3 4_prep_IRC.py "$sys" 2>/dev/null
  if [[ ! -f "$irc_dir/$sys.inp" ]]; then
    echo "  FAILED: IRC input not generated for $sys"
    FAILED+=("$sys:IRC_prep_fail")
    continue
  fi
  (cd "$irc_dir" && "$ORCA_CMD" "$sys.inp" > "$sys.out" 2>&1)

  if ! grep -q "ORCA TERMINATED NORMALLY" "$irc_dir/$sys.out" 2>/dev/null; then
    echo "  FAILED: IRC did not complete for $sys"
    FAILED+=("$sys:IRC_fail")
    continue
  fi
  echo "  IRC done."
  python3 6_plot_irc_energy_profile.py "$sys" 2>/dev/null || true

  # Step 3: Geo opt + Freq for reactant/product
  echo "  Step 3: Reactant/Product Opt+Freq..."
  python3 7_prep_geo_opt_reactant_product.py "$sys" 2>/dev/null

  if [[ ! -f "$rp_dir/reactant_optfreq.inp" ]] || [[ ! -f "$rp_dir/product_optfreq.inp" ]]; then
    echo "  FAILED: reactant/product inputs not generated for $sys"
    FAILED+=("$sys:RP_prep_fail")
    continue
  fi

  (cd "$rp_dir" && "$ORCA_CMD" reactant_optfreq.inp > reactant_optfreq.out 2>&1)
  (cd "$rp_dir" && "$ORCA_CMD" product_optfreq.inp > product_optfreq.out 2>&1)
  echo "  Reactant/Product Opt+Freq done."

  # Step 4: DeltaH
  echo "  Step 4: DeltaH..."
  python3 9_compute_deltaH.py "$sys" 2>&1 || true

  SUCCEEDED+=("$sys")
  echo "  $sys COMPLETE"
  echo ""
done

echo "========================================"
echo "SUMMARY"
echo "========================================"
echo "Succeeded (${#SUCCEEDED[@]}): ${SUCCEEDED[*]}"
if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo "Failed (${#FAILED[@]}): ${FAILED[*]}"
fi
echo ""

# Activate conda for matplotlib/numpy
source /software/python-anaconda-2021.05-el8-x86_64/etc/profile.d/conda.sh 2>/dev/null || true
conda activate base 2>/dev/null || true

# Step 5: Combined plots
echo "=== Generating combined frequency and IRC plots ==="
python3 3_frequency_analysis.py --all --combined 2>&1 || true
python3 6_plot_irc_energy_profile.py --combined 2>&1 || true

# Step 6: Polanyi analysis
if [[ ${#SUCCEEDED[@]} -ge 2 ]]; then
  echo "=== Running Polanyi Analysis ==="
  python3 10_polanyi_analysis.py --root "$RUNS_DIR" 2>&1
fi

echo "=== ALL DONE ==="
