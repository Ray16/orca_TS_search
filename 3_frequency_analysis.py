import argparse
import os
import re

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import MaxNLocator

from workflow_utils import resolve_system_dir

def extract_frequencies(file_path):
    frequencies = []
    if not os.path.exists(file_path):
        print(f"Error: {file_path} not found.")
        return None

    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Improved regex for "VIBRATIONAL FREQUENCIES" header
    # ORCA 6 uses "VIBRATIONAL FREQUENCIES" followed by a newline and dashes
    blocks = list(re.finditer(r"VIBRATIONAL FREQUENCIES\s*\n-+\s*\n", content))
    if not blocks:
        print("No 'VIBRATIONAL FREQUENCIES' header found in the output file.")
        return None
    
    last_block_start = blocks[-1].end()
    # Look for lines like "   6:   -1121.98 cm**-1"
    freq_pattern = re.compile(r"\s+(\d+):\s+([-]?\d+\.\d+)\s+cm\*\*-1")
    
    # Search within a reasonable range after the block start
    search_text = content[last_block_start:last_block_start + 5000]
    for match in freq_pattern.finditer(search_text):
        index = int(match.group(1))
        value = float(match.group(2))
        frequencies.append({"Index": index, "Frequency": value})
    
    if not frequencies:
        print("Header found but no frequency lines were parsed. Check the file format.")
        return None

    return pd.DataFrame(frequencies)

def analyze_and_plot(df, output_csv, output_plot):
    if df is None or df.empty:
        return

    # Save to CSV
    df.to_csv(output_csv, index=False)
    print(f"Frequencies saved to {output_csv}")

    # Identify imaginary frequencies (negative values in ORCA output)
    imaginary = df[df['Frequency'] < 0]
    real = df[df['Frequency'] >= 0]

    plt.figure(figsize=(10, 6))
    
    # Plot real frequencies in blue
    if not real.empty:
        plt.stem(real['Index'], real['Frequency'], linefmt='b-', markerfmt='bo', basefmt='k-', label='Real Modes')
    
    # Plot imaginary frequencies in red
    if not imaginary.empty:
        plt.stem(imaginary['Index'], imaginary['Frequency'], linefmt='r-', markerfmt='ro', basefmt='k-', label='Imaginary Modes')
        for i, row in imaginary.iterrows():
            plt.annotate(f"{row['Frequency']:.2f}", (row['Index'], row['Frequency']), 
                         textcoords="offset points", xytext=(0,-15), ha='center', color='red', weight='bold')
            print(f"Found imaginary frequency: {row['Frequency']} cm^-1 at index {row['Index']}")

    plt.axhline(0, color='black', linewidth=0.8)
    plt.ylim(bottom=-1000)
    plt.xlabel('Mode Index')
    plt.ylabel('Frequency (cm^-1)')
    plt.title('ORCA Vibrational Frequency Analysis')
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.savefig(output_plot, dpi=300)
    print(f"Plot saved to {output_plot}")

def resolve_input_file(target_dir, explicit_input=None):
    if explicit_input:
        if os.path.isabs(explicit_input):
            candidate = explicit_input
        else:
            # Prefer path relative to target_dir first, then target_dir/TS_search.
            candidate = os.path.join(target_dir, explicit_input)
            if not os.path.exists(candidate):
                candidate = os.path.join(target_dir, "TS_search", explicit_input)
        if os.path.exists(candidate):
            return candidate
        print(f"Error: input file not found: {candidate}")
        return None

    ts_search_dir = os.path.join(target_dir, "TS_search")
    ts_out_files = []
    if os.path.isdir(ts_search_dir):
        ts_out_files = [
            os.path.join(ts_search_dir, f)
            for f in os.listdir(ts_search_dir)
            if f.endswith(".out")
        ]
        ts_out_files.sort()
        if ts_out_files:
            if len(ts_out_files) > 1:
                print(
                    f"Multiple TS_search outputs found in {ts_search_dir}; using {os.path.basename(ts_out_files[0])}"
                )
            return ts_out_files[0]

    out_files = [
        os.path.join(target_dir, f)
        for f in os.listdir(target_dir)
        if f.endswith(".out")
    ]
    if not out_files:
        print(f"Error: no .out file found in {target_dir}")
        return None

    legacy_ts_out_files = [f for f in out_files if f.endswith("_ts.out")]
    candidates = legacy_ts_out_files if legacy_ts_out_files else out_files
    candidates.sort()
    if len(candidates) > 1:
        print(f"Multiple .out files found in {target_dir}; using {os.path.basename(candidates[0])}")
    return candidates[0]


def list_systems_from_ts_guess(root_dir):
    guess_dir = os.path.join(root_dir, "TS_guess_xyz")
    pattern = "_input.xyz"
    if not os.path.isdir(guess_dir):
        print("Error: TS_guess_xyz directory not found; pass system explicitly.")
        return None
    candidates = sorted(
        f[: -len(pattern)] for f in os.listdir(guess_dir) if f.endswith(pattern)
    )
    if not candidates:
        print("Error: no *_input.xyz found in TS_guess_xyz; pass system explicitly.")
        return None
    return candidates


SN2_ORDER = [
    "sn2_f_cl", "sn2_f_br", "sn2_f_i",
    "sn2_cl_f", "sn2_cl_br", "sn2_cl_i",
    "sn2_br_f", "sn2_br_cl", "sn2_br_i",
    "sn2_i_f",  "sn2_i_cl", "sn2_i_br",
]


def plot_combined_frequencies(systems_data, output_path, nrows, ncols, title):
    """systems_data: list of (sysname, df) tuples in row-major order."""
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 3.5 * nrows))
    axes_flat = [ax for row in (axes if nrows > 1 else [axes]) for ax in (row if ncols > 1 else [row])]

    for idx, (sysname, df) in enumerate(systems_data):
        ax = axes_flat[idx]
        if df is None or df.empty:
            ax.set_visible(False)
            continue

        imaginary = df[df["Frequency"] < 0]
        real = df[df["Frequency"] >= 0]

        if not real.empty:
            ax.stem(real["Index"].to_numpy(), real["Frequency"].to_numpy(),
                    linefmt="b-", markerfmt="bo", basefmt="k-", label="Real")
        if not imaginary.empty:
            ax.stem(imaginary["Index"].to_numpy(), imaginary["Frequency"].to_numpy(),
                    linefmt="r-", markerfmt="ro", basefmt="k-", label="Imag")
            for _, row in imaginary.iterrows():
                ax.annotate(f"{row['Frequency']:.0f}", (row["Index"], row["Frequency"]),
                            textcoords="offset points", xytext=(0, -12),
                            ha="center", color="red", fontsize=7, weight="bold")

        ax.axhline(0, color="black", linewidth=0.8)
        # Extend y-axis downward so imaginary frequency labels don't overlap x-axis
        if not imaginary.empty:
            ymin = imaginary["Frequency"].min()
            ax.set_ylim(bottom=min(ymin * 2.5, ymin - 300))
        short = re.sub(r"^(sn2|da)_", "", sysname).replace("_", u"\u2192")
        ax.set_title(short, fontsize=9, fontweight="bold")
        ax.set_xlabel("Mode", fontsize=7)
        ax.set_ylabel("cm\u207b\u00b9", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(True, linestyle="--", alpha=0.5)

    for idx in range(len(systems_data), nrows * ncols):
        axes_flat[idx].set_visible(False)

    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    print(f"Combined frequency plot saved: {output_path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Extract ORCA vibrational frequencies and generate CSV/plot outputs."
    )
    parser.add_argument(
        "system",
        nargs="?",
        default=None,
        help="System name, e.g., sn2 (optional if TS_guess_xyz has exactly one *_input.xyz).",
    )
    parser.add_argument(
        "--input",
        dest="input_file",
        default=None,
        help="Optional input .out filename (relative to system) or absolute path.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Process all systems listed by *_input.xyz in TS_guess_xyz.",
    )
    parser.add_argument(
        "--combined",
        action="store_true",
        help="Save one combined subplot figure per family instead of individual plots.",
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    if args.all and args.system:
        print("Error: use either --all or a specific system, not both.")
        return 1

    systems = None
    if args.all:
        systems = list_systems_from_ts_guess(script_dir)
        if not systems:
            return 1
    else:
        systems = [args.system]

    output_dir = os.path.join(script_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    # Collect (sysname, df) for all systems regardless of combined mode
    collected = []
    for system in systems:
        target_dir = resolve_system_dir(system, script_dir)
        sysname = os.path.basename(os.path.normpath(target_dir))
        input_file = resolve_input_file(target_dir, args.input_file)
        if input_file is None:
            continue
        df = extract_frequencies(input_file)
        freq_csv_dir = os.path.join(output_dir, "frequency")
        os.makedirs(freq_csv_dir, exist_ok=True)
        csv_file = os.path.join(freq_csv_dir, f"{sysname}_frequencies.csv")
        if df is not None and not df.empty:
            df.to_csv(csv_file, index=False)
            print(f"Frequencies saved to {csv_file}")
        collected.append((sysname, df))

    if args.combined:
        sn2_data = [(s, d) for s, d in collected if s in SN2_ORDER]
        sn2_data.sort(key=lambda x: SN2_ORDER.index(x[0]) if x[0] in SN2_ORDER else 999)
        if sn2_data:
            plot_combined_frequencies(sn2_data, os.path.join(output_dir, "sn2_frequencies_combined.png"),
                                      nrows=4, ncols=3, title="SN2 Vibrational Frequencies")

        da_data = [(s, d) for s, d in collected if s.startswith("da_")]
        da_data.sort(key=lambda x: x[0])
        if da_data:
            plot_combined_frequencies(da_data, os.path.join(output_dir, "da_frequencies_combined.png"),
                                      nrows=2, ncols=2, title="Diels-Alder Vibrational Frequencies")
    else:
        for sysname, df in collected:
            plot_file = os.path.join(output_dir, f"{sysname}_frequency_plot.png")
            if df is not None and not df.empty:
                analyze_and_plot(df, os.path.join(output_dir, f"{sysname}_frequencies.csv"), plot_file)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
