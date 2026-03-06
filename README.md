# ORCA TS search/IRC Workflow

Performing transition state search and intrinsic reaction calculation on UChicago RCC cluster.

## Environment Setup
- Make sure you have `pandas` and `matplotlib` Python packages installed via pip (e.g. `pip install pandas matplotlib`).
- Make sure orca is loaded via `module load orca`

## Input Preparetion
Place an initial guess for the transition state structure under `TS_guess_xyz` folder (xyz format), name it as `<sys_name>_input.xyz`, where `<sys_name>` is the name of the system.

## TS search and IRC

```bash
python 1_prep_TS_search.py [system] --charge [charge] --mult [multiplicity]
bash 2_submit_TS_search.sh [system]
python 3_frequency_analysis.py [system]
python 4_prep_IRC.py [system]
bash 5_submit_IRC.sh [system]
python 6_plot_irc_energy_profile.py [system]
```

## Notes
- [charge] and [multiplicity] should be manually set for the system. For example, for the following SN2 reaction:
F⁻ + CH₃F → [F···CH₃···F]⁻ → CH₃F + F⁻  
charge should be -1 and multiplicity should be 1 (singlet, all electrons paired). Multiplicity is set by counting the number of unpaired electrons in the system and using the formula multiplicity = 2S + 1, where S is the total spin (½ per unpaired electron), so a closed-shell system with no unpaired electrons has multiplicity = 1 (singlet).
- `system` is optional for steps `1/3/4/6` if `TS_guess_xyz` contains exactly one `*_input.xyz`.
- Steps `3` and `6` write outputs to `output/`.
