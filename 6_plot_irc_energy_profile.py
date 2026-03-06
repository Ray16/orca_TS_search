import argparse
import os
import re

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd


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
    parser.add_argument(
        "--base",
        default=None,
        help="System base name (e.g., sn2).",
    )
    parser.add_argument(
        "--csv",
        default=None,
        help="Optional output CSV filename (defaults to <sysname>_irc_energy_profile.csv in output/).",
    )
    parser.add_argument(
        "--png",
        default=None,
        help="Optional output plot filename (defaults to <sysname>_irc_energy_profile.png in output/).",
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    def resolve_system_name(system_arg, root_dir):
        if system_arg:
            return system_arg
        guess_dir = os.path.join(root_dir, "TS_guess_xyz")
        pattern = "_input.xyz"
        candidates = []
        if os.path.isdir(guess_dir):
            candidates = sorted(
                f[: -len(pattern)]
                for f in os.listdir(guess_dir)
                if f.endswith(pattern)
            )
        if len(candidates) == 1:
            print(f"Auto-detected system from TS_guess_xyz: {candidates[0]}")
            return candidates[0]
        if len(candidates) > 1:
            print("Error: multiple *_input.xyz files found in TS_guess_xyz; pass system explicitly.")
            print("Candidates: " + ", ".join(candidates))
            raise SystemExit(1)
        print("Error: no *_input.xyz found in TS_guess_xyz; pass system explicitly.")
        raise SystemExit(1)

    system_name = resolve_system_name(args.system, script_dir)
    target_dir = system_name if os.path.isabs(system_name) else os.path.join(script_dir, system_name)
    target_dir = os.path.abspath(target_dir)
    if not os.path.isdir(target_dir):
        print(f"Error: target directory not found: {target_dir}")
        raise SystemExit(1)

    try:
        irc_dir = detect_irc_dir(target_dir)
        stem = detect_stem(irc_dir, args.base)
        df = build_profile(irc_dir, stem)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        raise SystemExit(1)

    sysname = os.path.basename(os.path.dirname(irc_dir))
    output_dir = os.path.join(script_dir, "output")
    os.makedirs(output_dir, exist_ok=True)
    csv_name = args.csv if args.csv else f"{sysname}_irc_energy_profile.csv"
    png_name = args.png if args.png else f"{sysname}_irc_energy_profile.png"
    csv_path = os.path.join(output_dir, csv_name)
    png_path = os.path.join(output_dir, png_name)
    df.to_csv(csv_path, index=False)
    print(f"Wrote CSV: {csv_path}")

    plt.figure(figsize=(9, 5.5))
    plt.plot(df["irc_step"], df["rel_energy_kcal_mol"], marker="o", linewidth=1.5)
    ts_rows = df[df["direction"] == "TS"]
    if not ts_rows.empty:
        ts_step = float(ts_rows.iloc[0]["irc_step"])
        ts_rel = float(ts_rows.iloc[0]["rel_energy_kcal_mol"])
        plt.scatter([ts_step], [ts_rel], color="crimson", zorder=5, label="TS")
        plt.axhline(ts_rel, color="crimson", linestyle="--", linewidth=1.0, alpha=0.8)
        plt.text(
            0.02,
            0.96,
            f"Barrier: {ts_rel:.2f} kcal/mol",
            transform=plt.gca().transAxes,
            ha="left",
            va="top",
            fontsize=10,
            color="crimson",
            weight="bold",
            bbox={"facecolor": "white", "edgecolor": "crimson", "alpha": 0.9, "boxstyle": "round,pad=0.25"},
        )
    plt.xlabel("IRC Step")
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.ylabel("Relative Energy (kcal/mol)")
    plt.title(f"IRC Energy Profile ({stem}_IRC)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.savefig(png_path, dpi=300)
    print(f"Wrote plot: {png_path}")


if __name__ == "__main__":
    main()
