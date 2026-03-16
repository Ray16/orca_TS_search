import argparse
import os
import re

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd

from workflow_utils import resolve_system_name, resolve_system_dir


HARTREE_TO_KCAL_MOL = 627.509474


def parse_xyz_energies(xyz_path):
    if not os.path.exists(xyz_path):
        return []

    energies = []
    with open(xyz_path, "r", encoding="utf-8") as f:
        for line in f:
            if "Coordinates from ORCA-job" not in line:
                continue
            m = re.search(r"\bE\s+(-?\d+\.\d+)", line)
            if m:
                energies.append(float(m.group(1)))
    return energies


def parse_first_final_sp_energy(out_path):
    if not os.path.exists(out_path):
        return None
    pattern = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
    with open(out_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            m = pattern.search(line)
            if m:
                return float(m.group(1))
    return None


def detect_irc_dir(target_dir):
    if os.path.basename(os.path.normpath(target_dir)) == "IRC":
        return target_dir
    irc_dir = os.path.join(target_dir, "IRC")
    if os.path.isdir(irc_dir):
        return irc_dir
    raise FileNotFoundError(f"IRC directory not found under {target_dir}")


def normalize_stem(value):
    return value[:-4] if value.endswith("_IRC") else value


def detect_stem(irc_dir, explicit_base=None):
    if explicit_base:
        return normalize_stem(explicit_base)

    candidates = sorted(
        normalize_stem(f[:-4])
        for f in os.listdir(irc_dir)
        if f.endswith(".inp") and not f.endswith("_IRC.inp") and not f.endswith("_TS_IRC.inp")
    )
    if not candidates:
        raise FileNotFoundError(
            f"No IRC input '.inp' found in {irc_dir}; pass --base explicitly."
        )
    if len(candidates) > 1:
        print(f"Multiple IRC bases found, using: {candidates[0]}")
    return candidates[0]


def build_profile(irc_dir, stem):
    f_trj = os.path.join(irc_dir, f"{stem}_IRC_F_trj.xyz")
    b_trj = os.path.join(irc_dir, f"{stem}_IRC_B_trj.xyz")
    if not os.path.exists(f_trj):
        legacy_f_trj = os.path.join(irc_dir, f"{stem}_IRC_IRC_F_trj.xyz")
        if os.path.exists(legacy_f_trj):
            f_trj = legacy_f_trj
    if not os.path.exists(b_trj):
        legacy_b_trj = os.path.join(irc_dir, f"{stem}_IRC_IRC_B_trj.xyz")
        if os.path.exists(legacy_b_trj):
            b_trj = legacy_b_trj

    ts_out = os.path.join(irc_dir, f"{stem}.out")
    ts_e = parse_first_final_sp_energy(ts_out)
    if ts_e is None:
        raise FileNotFoundError(
            f"Could not find 'FINAL SINGLE POINT ENERGY' in {ts_out}"
        )
    f_energies = parse_xyz_energies(f_trj)
    b_energies = parse_xyz_energies(b_trj)

    if not f_energies and not b_energies:
        raise FileNotFoundError(
            "No IRC trajectory energies found. Expected files like "
            f"{os.path.basename(f_trj)} and/or {os.path.basename(b_trj)}"
        )

    rows = []
    rows.append(
        {
            "irc_step": 0,
            "direction": "TS",
            "energy_hartree": ts_e,
        }
    )

    for i, e in enumerate(b_energies, start=1):
        rows.append(
            {
                "irc_step": -i,
                "direction": "B",
                "energy_hartree": e,
            }
        )

    for i, e in enumerate(f_energies, start=1):
        rows.append(
            {
                "irc_step": i,
                "direction": "F",
                "energy_hartree": e,
            }
        )

    df = pd.DataFrame(rows).sort_values("irc_step").reset_index(drop=True)
    e0 = df["energy_hartree"].min()
    df["rel_energy_kcal_mol"] = (df["energy_hartree"] - e0) * HARTREE_TO_KCAL_MOL
    return df


SN2_ORDER = [
    "sn2_f_cl", "sn2_f_br", "sn2_f_i",
    "sn2_cl_f", "sn2_cl_br", "sn2_cl_i",
    "sn2_br_f", "sn2_br_cl", "sn2_br_i",
    "sn2_i_f",  "sn2_i_cl", "sn2_i_br",
]


def _plot_irc_on_ax(ax, df, sysname):
    short = re.sub(r"^(sn2|da)_", "", sysname).replace("_", u"\u2192")
    ax.plot(df["irc_step"], df["rel_energy_kcal_mol"], marker="o", markersize=3, linewidth=1.5)
    ts_rows = df[df["direction"] == "TS"]
    if not ts_rows.empty:
        ts_step = float(ts_rows.iloc[0]["irc_step"])
        ts_rel = float(ts_rows.iloc[0]["rel_energy_kcal_mol"])
        ax.scatter([ts_step], [ts_rel], color="crimson", zorder=5, s=40)
        ax.axhline(ts_rel, color="crimson", linestyle="--", linewidth=0.8, alpha=0.8)
        ax.text(0.03, 0.96, f"{ts_rel:.1f} kcal/mol",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=7, color="crimson", weight="bold",
                bbox={"facecolor": "white", "edgecolor": "crimson",
                      "alpha": 0.85, "boxstyle": "round,pad=0.2"})
    ax.set_title(short, fontsize=9, fontweight="bold")
    ax.set_xlabel("IRC Step", fontsize=7)
    ax.set_ylabel("Rel. E (kcal/mol)", fontsize=7)
    ax.tick_params(labelsize=6)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(True, linestyle="--", alpha=0.5)


def plot_combined_irc(systems_data, output_path, nrows, ncols, title):
    """systems_data: list of (sysname, df) tuples."""
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows))
    axes_flat = [ax for row in (axes if nrows > 1 else [axes]) for ax in (row if ncols > 1 else [row])]

    for idx, (sysname, df) in enumerate(systems_data):
        _plot_irc_on_ax(axes_flat[idx], df, sysname)

    for idx in range(len(systems_data), nrows * ncols):
        axes_flat[idx].set_visible(False)

    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    print(f"Combined IRC plot saved: {output_path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Create IRC energy profile CSV and PNG from ORCA IRC trajectory files."
    )
    parser.add_argument(
        "system",
        nargs="?",
        default=None,
        help="System name, e.g., sn2 (optional if TS_guess_xyz has exactly one *_input.xyz).",
    )
    parser.add_argument("--base", default=None, help="System base name (e.g., sn2).")
    parser.add_argument("--csv", default=None, help="Optional output CSV filename.")
    parser.add_argument("--png", default=None, help="Optional output plot filename.")
    parser.add_argument(
        "--combined",
        action="store_true",
        help="Save one combined subplot figure per family instead of individual plots.",
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    # Determine list of systems to process
    if args.combined:
        runs_dir = os.path.join(script_dir, "runs")
        all_systems = sorted(
            d for d in os.listdir(runs_dir)
            if os.path.isdir(os.path.join(runs_dir, d, "IRC"))
        ) if os.path.isdir(runs_dir) else []
        if args.system:
            all_systems = [args.system]
    else:
        all_systems = [resolve_system_name(args.system, script_dir)]

    # Collect profiles
    collected = []
    for sysname in all_systems:
        target_dir = resolve_system_dir(sysname, script_dir)
        if not os.path.isdir(target_dir):
            print(f"Skipping {sysname}: directory not found")
            continue
        try:
            irc_dir = detect_irc_dir(target_dir)
            stem = detect_stem(irc_dir, args.base)
            df = build_profile(irc_dir, stem)
            real_sysname = os.path.basename(os.path.dirname(irc_dir))
            csv_name = args.csv if (args.csv and not args.combined) else f"{real_sysname}_irc_energy_profile.csv"
            df.to_csv(os.path.join(output_dir, csv_name), index=False)
            print(f"Wrote CSV: {os.path.join(output_dir, csv_name)}")
            collected.append((real_sysname, df))
        except FileNotFoundError as e:
            print(f"Skipping {sysname}: {e}")

    if args.combined:
        sn2_data = [(s, d) for s, d in collected if s in SN2_ORDER]
        sn2_data.sort(key=lambda x: SN2_ORDER.index(x[0]) if x[0] in SN2_ORDER else 999)
        if sn2_data:
            plot_combined_irc(sn2_data, os.path.join(output_dir, "sn2_irc_combined.png"),
                              nrows=4, ncols=3, title="SN2 IRC Energy Profiles")
    else:
        if collected:
            sysname, df = collected[0]
            stem = detect_stem(detect_irc_dir(resolve_system_dir(sysname, script_dir)), args.base)
            png_name = args.png if args.png else f"{sysname}_irc_energy_profile.png"
            png_path = os.path.join(output_dir, png_name)
            fig, ax = plt.subplots(figsize=(9, 5.5))
            _plot_irc_on_ax(ax, df, sysname)
            ax.set_title(f"IRC Energy Profile ({stem}_IRC)")
            fig.tight_layout()
            fig.savefig(png_path, dpi=300)
            print(f"Wrote plot: {png_path}")
            plt.close(fig)


if __name__ == "__main__":
    main()
