import argparse
import os
import re

import matplotlib.pyplot as plt
import pandas as pd

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
    plt.xlabel('Mode Index')
    plt.ylabel('Frequency (cm^-1)')
    plt.title('ORCA Vibrational Frequency Analysis')
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
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    target_dir = resolve_system_dir(args.system, script_dir)
    sysname = os.path.basename(os.path.normpath(target_dir))

    input_file = resolve_input_file(target_dir, args.input_file)
    if input_file is None:
        return 1

    output_dir = os.path.join(script_dir, "output")
    os.makedirs(output_dir, exist_ok=True)
    csv_file = os.path.join(output_dir, f"{sysname}_frequencies.csv")
    plot_file = os.path.join(output_dir, f"{sysname}_frequency_plot.png")

    df = extract_frequencies(input_file)
    analyze_and_plot(df, csv_file, plot_file)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
