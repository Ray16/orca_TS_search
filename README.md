# Evans-Polanyi Principle Validation

Computational validation of the [Evans-Polanyi principle](https://en.wikipedia.org/wiki/Evans%E2%80%93Polanyi_principle) — the linear relationship between activation energy (E_a) and reaction enthalpy (ΔH) within a family of analogous reactions — using quantum chemistry calculations at the GFN2-xTB level with ORCA.

The pipeline automates transition state (TS) search, intrinsic reaction coordinate (IRC) calculation, and thermochemistry extraction to produce Evans-Polanyi correlation plots (E_a vs ΔH) for each reaction family.

## Reaction Families

- **SN2**: 12 halide-exchange reactions (X⁻ + CH₃Y → CH₃X + Y⁻, where X, Y = F, Cl, Br, I). In the system name `sn2_{nuc}_{lg}`, the first halogen is the nucleophile (nuc, X) and the second is the leaving group (lg, Y).
- **Diels-Alder**: 4 [4+2] cycloadditions (1,3-butadiene + CH₂=CHX → 4-X-cyclohex-2-ene, where X = F, Cl, Br, I). In the system name `da_{X}`, the halogen is the substituent on the dienophile.

## Environment Setup

```bash
module load orca
pip install pandas matplotlib numpy
```

## Quick Start

Run the full pipeline for a reaction family — from TS search through Evans-Polanyi analysis:

```bash
# SN2 systems (~5 min total)
bash run_all_sn2_pipeline.sh

# Diels-Alder systems (~2 min total)
bash run_all_da_pipeline.sh
```

Each script handles: TS optimization → IRC → reactant/product geometry opt → ΔH computation → Evans-Polanyi fit. Outputs (correlation plots, fit reports) are written to `output/`.

## Workflow Details

The pipeline consists of 10 steps. Steps 1–9 compute E_a and ΔH for each individual reaction; step 10 fits the Evans-Polanyi relation across all systems in a family.

Replace `[system]` with the system name (e.g., `sn2_cl_br`, `da_i`):

```bash
# 1. Prepare TS search input
python 1_prep_TS_search.py [system] --charge <charge> --mult <multiplicity>

# 2. Run TS optimization + frequency calculation
bash 2_submit_TS_search.sh [system]

# 3. Analyze frequencies (checks for exactly one imaginary mode)
python 3_frequency_analysis.py [system]

# 4. Prepare IRC input from TS result
python 4_prep_IRC.py [system]

# 5. Run IRC
bash 5_submit_IRC.sh [system]

# 6. Plot IRC energy profile
python 6_plot_irc_energy_profile.py [system]

# 7. Prepare reactant/product geometry optimization inputs
python 7_prep_geo_opt_reactant_product.py [system]

# 8. Run reactant/product opt + frequency
bash 8_submit_geo_opt_reactant_product.sh [system]

# 9. Compute reaction enthalpy (ΔH)
python 9_compute_deltaH.py [system]

# 10. Evans-Polanyi analysis (run after multiple systems are complete)
python 10_polanyi_analysis.py
```

### Input Preparation (optional)

Initial TS guess structures are provided in `TS_guess_xyz/` as `<system>_input.xyz`. To regenerate them:

```bash
# SN2 systems
python 0_prep_initial_TS_structure.py --nuc <nucleophile> --lg <leaving_group>

# Diels-Alder systems
python 0_prep_initial_TS_structure.py --family da
```

## Notes

- **`[system]`** is optional for steps 1, 3, 4, and 6 if `TS_guess_xyz/` contains exactly one `*_input.xyz` file.
- **Charge and multiplicity** must be set correctly:
  - **SN2**: e.g., F⁻ + CH₃Cl → CH₃F + Cl⁻: charge = `-1`, multiplicity = `1`.
  - **Diels-Alder**: e.g., butadiene + CH₂=CHI: charge = `0`, multiplicity = `1`.
- **Step 3** confirms the TS has exactly one imaginary frequency — if not, the TS guess or optimization needs to be revisited before proceeding.
- **SN2 double-well potential**: In SN2 reactions involving F as nucleophile or leaving group, pre-reactive ion-dipole complexes (van der Waals wells) create a double-well potential energy surface. The IRC can pass through these shallow wells, causing a nearby point to appear slightly higher in energy than the optimized TS due to numerical noise in the IRC stepping. As long as the TS has exactly one imaginary frequency and the energy difference with the apparent maximum is small, the TS is valid.
- **Reactant/product assignment**: By convention, negative IRC steps correspond to the reactant side and positive steps to the product side. However, the actual direction is arbitrary — it depends on the sign of the imaginary mode eigenvector, which varies between systems. The pipeline determines the correct assignment by examining endpoint geometries:
  - **SN2**: the true reactant has the leaving group bonded to C and the nucleophile far from C.
  - **Diels-Alder**: the true product (cyclohexene ring) has more C-C bonds (< 2.0 Å) than the reactant (separated diene + dienophile).

  This geometry check is applied at step 7 (when writing `reactant.xyz`/`product.xyz`), so the files are correctly named from the start. Steps 6 and 10 apply the same check as a safety measure. The IRC plots preserve the raw step data as computed and label R/P at the correct endpoints accordingly, so R and P may appear on either side of the plot depending on the system.
- **IRC energy profile vs. true E_a**: The energy differences visible in the IRC plots are approximate — the IRC endpoints are snapshots along the reaction path that have not been fully geometry-optimized. The true E_a and ΔH used in the Evans-Polanyi analysis come from separately optimized reactant and product geometries (steps 7–9), which relax to proper energy minima.
- **Barrierless reactions**: For some highly exothermic SN2 reactions (e.g., F⁻ + CH₃I), the activation energy from the optimized geometries can be slightly negative due to numerical noise. These are clamped to E_a = 0 (barrierless).
