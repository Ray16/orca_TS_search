# ORCA TS/IRC Workflow

## Setup

- Python 3.8+
- Python packages:
  - `pandas`
  - `matplotlib`
- ORCA available in shell (`module load orca` on cluster)

Install Python packages:

```bash
pip install pandas matplotlib
```

## Run Order

```bash
python 1_prep_TS_search.py [system]
bash 2_submit_TS_search.sh [system]
python 3_frequency_analysis.py [system]
python 4_prep_IRC.py [system]
bash 5_submit_IRC.sh [system]
python 6_plot_irc_energy_profile.py [system]
```

Notes:
- `system` is optional for steps `1/3/4/6` if `TS_guess_xyz` contains exactly one `*_input.xyz`.
- Steps `3` and `6` write outputs to `output/`.
