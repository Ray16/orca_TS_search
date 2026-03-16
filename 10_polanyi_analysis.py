import argparse
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np

HARTREE_TO_KCAL = 627.509474


def extract_total_enthalpy(out_path):
    """Extract the last 'Total enthalpy' (or 'TOTAL ENTHALPY') value in Eh from an ORCA output."""
    if not os.path.exists(out_path):
        return None
    with open(out_path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()
    pattern = re.compile(r"Total\s+enthalpy.*?(-?\d+\.\d+)", flags=re.IGNORECASE)
    matches = pattern.findall(text)
    if matches:
        return float(matches[-1])
    # XTB format: "TOTAL ENTHALPY            -12.374686 Eh"
    pattern2 = re.compile(r"TOTAL\s+ENTHALPY\s+(-?\d+\.\d+)", flags=re.IGNORECASE)
    matches = pattern2.findall(text)
    if matches:
        return float(matches[-1])
    return None


def extract_final_energy(out_path):
    """Fallback: extract FINAL SINGLE POINT ENERGY if enthalpy not available."""
    if not os.path.exists(out_path):
        return None
    with open(out_path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()
    pattern = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
    matches = pattern.findall(text)
    if matches:
        return float(matches[-1])
    return None


def get_enthalpy_or_energy(out_path, label=""):
    """Try enthalpy first, fall back to electronic energy."""
    h = extract_total_enthalpy(out_path)
    if h is not None:
        return h, "enthalpy"
    e = extract_final_energy(out_path)
    if e is not None:
        print(f"  Warning: using electronic energy (no enthalpy) for {label}: {out_path}")
        return e, "energy"
    return None, None


def list_systems(root_dir):
    """List all system directories for Polanyi analysis."""
    systems = sorted(
        d for d in os.listdir(root_dir)
        if (d.startswith("sn2_") or d.startswith("da_"))
        and os.path.isdir(os.path.join(root_dir, d))
    )
    return systems


def collect_data(root_dir, systems):
    """For each system, collect H(TS), H(reactant), H(product)."""
    rows = []
    for sysname in systems:
        sys_dir = os.path.join(root_dir, sysname)

        # TS enthalpy from OptTS+Freq output
        ts_out = os.path.join(sys_dir, "TS_search", f"{sysname}.out")
        h_ts, ts_type = get_enthalpy_or_energy(ts_out, f"{sysname} TS")

        # Reactant/product enthalpies from geo opt + freq outputs
        rp_dir = os.path.join(sys_dir, "geo_opt_reactant_product")
        react_out = os.path.join(rp_dir, "reactant_optfreq.out")
        prod_out = os.path.join(rp_dir, "product_optfreq.out")
        h_react, r_type = get_enthalpy_or_energy(react_out, f"{sysname} reactant")
        h_prod, p_type = get_enthalpy_or_energy(prod_out, f"{sysname} product")

        if h_ts is None or h_react is None or h_prod is None:
            missing = []
            if h_ts is None:
                missing.append("TS")
            if h_react is None:
                missing.append("reactant")
            if h_prod is None:
                missing.append("product")
            print(f"  Skipping {sysname}: missing {', '.join(missing)} enthalpy")
            continue

        # For Diels-Alder (always exothermic), product must have lower enthalpy.
        # IRC endpoint ordering can vary; swap if needed.
        if sysname.startswith("da_") and h_react < h_prod:
            h_react, h_prod = h_prod, h_react
            print(f"  Note: swapped reactant/product for {sysname} (reactant had lower H)")

        ea_kcal = (h_ts - h_react) * HARTREE_TO_KCAL
        dh_kcal = (h_prod - h_react) * HARTREE_TO_KCAL

        rows.append({
            "system": sysname,
            "H_TS": h_ts,
            "H_reactant": h_react,
            "H_product": h_prod,
            "Ea_kcal": ea_kcal,
            "dH_kcal": dh_kcal,
        })
        print(f"  {sysname}: Ea = {ea_kcal:.2f} kcal/mol, dH = {dh_kcal:.2f} kcal/mol")

    return rows


def fit_polanyi(rows):
    """Fit Evans-Polanyi: Ea = E0 + alpha * dH. Returns (E0, alpha, R^2)."""
    dh = np.array([r["dH_kcal"] for r in rows])
    ea = np.array([r["Ea_kcal"] for r in rows])
    # Linear fit: Ea = alpha * dH + E0
    coeffs = np.polyfit(dh, ea, 1)
    alpha, e0 = coeffs
    ea_pred = np.polyval(coeffs, dh)
    ss_res = np.sum((ea - ea_pred) ** 2)
    ss_tot = np.sum((ea - np.mean(ea)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return e0, alpha, r2


def _repel_labels(ax, fig, xs, ys, shorts, fontsize=7.5, pad_pts=4, iters=200, step=1.5):
    """
    Place annotation text boxes without overlap using iterative force-based repulsion.
    Returns list of (ox, oy) offset-in-points for each label.
    Works entirely in display (pixel) coordinates, then converts back to offset points.
    """
    renderer = fig.canvas.get_renderer()
    dpi = fig.dpi
    pts_to_px = dpi / 72.0

    # Initial offsets in points (slightly above each data point)
    n = len(xs)
    ox = np.full(n, 5.0)
    oy = np.full(n, 5.0)

    # Alternate initial quadrant to spread labels around dense clusters
    for i in range(n):
        ox[i] =  8.0 if i % 2 == 0 else -8.0
        oy[i] =  6.0 if (i // 2) % 2 == 0 else -14.0

    # Convert data points to display coords once
    pts_disp = ax.transData.transform(np.column_stack([xs, ys]))  # (n, 2) in pixels

    # Estimate label sizes in pixels using a dummy text artist
    widths = np.zeros(n)
    heights = np.zeros(n)
    dummy = ax.annotate("", (0, 0))
    for i, txt in enumerate(shorts):
        dummy.set_text(txt)
        dummy.set_fontsize(fontsize)
        try:
            bb = dummy.get_window_extent(renderer)
            widths[i] = bb.width + pad_pts * pts_to_px
            heights[i] = bb.height + pad_pts * pts_to_px
        except Exception:
            widths[i] = len(txt) * fontsize * pts_to_px * 0.6
            heights[i] = fontsize * pts_to_px * 1.4
    dummy.remove()

    for _ in range(iters):
        # Label centres in display coords
        lx = pts_disp[:, 0] + ox * pts_to_px
        ly = pts_disp[:, 1] + oy * pts_to_px

        moved = False
        for i in range(n):
            for j in range(i + 1, n):
                dx = lx[i] - lx[j]
                dy = ly[i] - ly[j]
                # Overlap occurs when both x and y extents intersect
                overlap_x = abs(dx) < (widths[i] + widths[j]) / 2
                overlap_y = abs(dy) < (heights[i] + heights[j]) / 2
                if overlap_x and overlap_y:
                    dist = max(np.hypot(dx, dy), 0.1)
                    ux, uy = dx / dist, dy / dist
                    push = step * pts_to_px
                    ox[i] += ux * push / pts_to_px
                    oy[i] += uy * push / pts_to_px
                    ox[j] -= ux * push / pts_to_px
                    oy[j] -= uy * push / pts_to_px
                    lx[i] = pts_disp[i, 0] + ox[i] * pts_to_px
                    ly[i] = pts_disp[i, 1] + oy[i] * pts_to_px
                    lx[j] = pts_disp[j, 0] + ox[j] * pts_to_px
                    ly[j] = pts_disp[j, 1] + oy[j] * pts_to_px
                    moved = True
        if not moved:
            break

    return list(zip(ox.tolist(), oy.tolist()))


def plot_polanyi(rows, e0, alpha, r2, png_path):
    dh = np.array([r["dH_kcal"] for r in rows])
    ea = np.array([r["Ea_kcal"] for r in rows])
    labels = [r["system"] for r in rows]
    shorts = [re.sub(r"^(sn2|da)_", "", lbl).replace("_", "/") for lbl in labels]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(dh, ea, s=60, zorder=5, color="steelblue", edgecolors="black", linewidth=0.5)

    dh_range = np.linspace(dh.min() - 2, dh.max() + 2, 100)
    ea_fit = e0 + alpha * dh_range
    ax.plot(dh_range, ea_fit, "r--", linewidth=1.5,
            label=f"Polanyi: $E_a$ = {e0:.2f} + {alpha:.3f} $\\Delta H$\n$R^2$ = {r2:.4f}")

    ax.set_xlabel("$\\Delta H_{rxn}$ (kcal/mol)", fontsize=12)
    ax.set_ylabel("$E_a$ (kcal/mol)", fontsize=12)
    ax.set_title("Evans-Polanyi Relation (GFN2-xTB / ORCA)", fontsize=13)
    ax.legend(fontsize=10, loc="best")
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()

    # Draw figure once so renderer + transform are available before repulsion
    fig.canvas.draw()
    offsets = _repel_labels(ax, fig, dh, ea, shorts)
    for i, short in enumerate(shorts):
        ox, oy = offsets[i]
        ax.annotate(short, (dh[i], ea[i]), textcoords="offset points",
                    xytext=(ox, oy), fontsize=7.5,
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))

    fig.tight_layout()
    fig.savefig(png_path, dpi=300)
    print(f"Wrote plot: {png_path}")


def write_report(rows, e0, alpha, r2, report_path):
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("Evans-Polanyi Analysis (GFN2-xTB / ORCA)\n")
        f.write("=" * 55 + "\n\n")
        f.write(f"{'System':<18} {'Ea (kcal/mol)':>14} {'dH (kcal/mol)':>14}\n")
        f.write("-" * 48 + "\n")
        for r in rows:
            f.write(f"{r['system']:<18} {r['Ea_kcal']:>14.2f} {r['dH_kcal']:>14.2f}\n")
        f.write("\n")
        f.write(f"Evans-Polanyi fit: Ea = {e0:.2f} + {alpha:.3f} * dH\n")
        f.write(f"R^2 = {r2:.4f}\n")
        f.write(f"Alpha (Bronsted coefficient) = {alpha:.3f}\n")
        f.write(f"E0 (intrinsic barrier) = {e0:.2f} kcal/mol\n")
        f.write("\n")
        if 0.0 < alpha < 1.0:
            f.write("Alpha is between 0 and 1 as expected by Evans-Polanyi theory.\n")
        else:
            f.write(f"Warning: alpha = {alpha:.3f} is outside the expected [0, 1] range.\n")
        if r2 > 0.8:
            f.write(f"Good linear correlation (R^2 = {r2:.4f}): Polanyi relation holds.\n")
        else:
            f.write(f"Weak correlation (R^2 = {r2:.4f}): Polanyi relation may not hold well at this level of theory.\n")
    print(f"Wrote report: {report_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Evans-Polanyi analysis: fit Ea vs dH."
    )
    parser.add_argument(
        "--root", default=None,
        help="Root directory containing system subdirectories (defaults to script directory).",
    )
    parser.add_argument(
        "--systems", nargs="*", default=None,
        help="Specific system names to include (default: auto-detect from root directory).",
    )
    parser.add_argument(
        "--outdir", default=None,
        help="Output directory (defaults to output/).",
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root = os.path.abspath(args.root) if args.root else script_dir

    all_systems = args.systems if args.systems else list_systems(root)
    if not all_systems:
        print("No system directories found.")
        return 1

    print(f"Collecting data for {len(all_systems)} systems...")
    all_rows = collect_data(root, all_systems)

    out_dir = args.outdir or os.path.join(script_dir, "output")
    os.makedirs(out_dir, exist_ok=True)

    # Save per-system Ea and dH to CSV
    import csv
    ea_dh_path = os.path.join(out_dir, "ea_dh_all_systems.csv")
    with open(ea_dh_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["system", "Ea_kcal", "dH_kcal"])
        writer.writeheader()
        for r in all_rows:
            writer.writerow({"system": r["system"],
                             "Ea_kcal": f"{r['Ea_kcal']:.4f}",
                             "dH_kcal": f"{r['dH_kcal']:.4f}"})
    print(f"Wrote Ea/dH table: {ea_dh_path}")

    # Run Polanyi analysis separately for sn2 and e2 families, then combined
    families = {
        "sn2": [r for r in all_rows if r["system"].startswith("sn2_")],
        "da":  [r for r in all_rows if r["system"].startswith("da_")],
        "all": all_rows,
    }

    for family, rows in families.items():
        if len(rows) < 2:
            print(f"Skipping {family}: need ≥2 data points, got {len(rows)}.")
            continue
        e0, alpha, r2 = fit_polanyi(rows)
        print(f"\n[{family}] Evans-Polanyi fit: Ea = {e0:.2f} + {alpha:.3f} * dH  (R²={r2:.4f})")
        png_path    = os.path.join(out_dir, f"{family}_polanyi_plot.png")
        report_path = os.path.join(out_dir, f"{family}_polanyi_report.txt")
        plot_polanyi(rows, e0, alpha, r2, png_path)
        write_report(rows, e0, alpha, r2, report_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
