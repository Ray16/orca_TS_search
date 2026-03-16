import argparse
import os
import re

from workflow_utils import resolve_system_dir

HARTREE_TO_KCAL = 627.509474


def find_out_files(rp_dir, reactant_out, product_out):
    if reactant_out and product_out:
        return reactant_out, product_out

    if not os.path.isdir(rp_dir):
        raise FileNotFoundError(f"Reactant/Product directory not found: {rp_dir}")

    outs = sorted(
        f for f in os.listdir(rp_dir) if f.endswith(".out")
    )
    if not outs:
        raise FileNotFoundError(f"No .out files found in {rp_dir}")

    def pick_by_keyword(keyword):
        matches = [f for f in outs if keyword in f.lower()]
        return matches[0] if matches else None

    reactant = reactant_out or pick_by_keyword("reactant")
    product = product_out or pick_by_keyword("product")

    if reactant and product:
        return os.path.join(rp_dir, reactant), os.path.join(rp_dir, product)

    # If not found, fall back to first two outputs (sorted) to avoid silent failure.
    if len(outs) >= 2:
        return os.path.join(rp_dir, outs[0]), os.path.join(rp_dir, outs[1])
    raise FileNotFoundError(
        "Could not identify reactant/product outputs. Pass --reactant-out and --product-out explicitly."
    )


def extract_total_enthalpy(out_path):
    if not os.path.exists(out_path):
        raise FileNotFoundError(f"Missing output: {out_path}")
    with open(out_path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()

    # Common ORCA labels:
    # "Total enthalpy" or "Total Enthalpy" with value in Eh.
    pattern = re.compile(r"Total\s+enthalpy.*?(-?\d+\.\d+)", flags=re.IGNORECASE)
    matches = pattern.findall(text)
    if matches:
        return float(matches[-1])

    # Fallback: thermochemistry block sometimes labels "Enthalpy" on its own line.
    # Example: "Enthalpy         -XXX.XXXXXX Eh"
    pattern2 = re.compile(r"^\s*Enthalpy\s+(-?\d+\.\d+)\s+Eh", flags=re.MULTILINE)
    matches = pattern2.findall(text)
    if matches:
        return float(matches[-1])

    raise ValueError(f"Could not find Total enthalpy in {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute DeltaH_rxn from ORCA opt+freq outputs (reactant/product)."
    )
    parser.add_argument(
        "system",
        nargs="?",
        default=None,
        help="System name, e.g., sn2 (optional if TS_guess_xyz has exactly one *_input.xyz).",
    )
    parser.add_argument("--root", default=None, help="Optional root directory (defaults to script directory).")
    parser.add_argument(
        "--reactant-out",
        default=None,
        help="Reactant .out path (absolute or relative to the reactant/product directory).",
    )
    parser.add_argument(
        "--product-out",
        default=None,
        help="Product .out path (absolute or relative to the reactant/product directory).",
    )
    parser.add_argument(
        "--rp-dir",
        default=None,
        help="Reactant/product directory (relative to system dir unless absolute).",
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root = os.path.abspath(args.root) if args.root else script_dir
    system_dir = resolve_system_dir(args.system, root)
    sysname = os.path.basename(os.path.normpath(system_dir))

    if args.rp_dir:
        rp_dir = args.rp_dir if os.path.isabs(args.rp_dir) else os.path.join(system_dir, args.rp_dir)
    else:
        rp_dir = os.path.join(system_dir, "geo_opt_reactant_product")

    def resolve_out(path):
        if not path:
            return None
        if os.path.isabs(path):
            return path
        return os.path.join(rp_dir, path)

    reactant_out = resolve_out(args.reactant_out)
    product_out = resolve_out(args.product_out)

    reactant_out, product_out = find_out_files(rp_dir, reactant_out, product_out)

    h_react = extract_total_enthalpy(reactant_out)
    h_prod = extract_total_enthalpy(product_out)
    dh_eh = h_prod - h_react
    dh_kcal = dh_eh * HARTREE_TO_KCAL

    print(f"DeltaH (Eh): {dh_eh:.10f}")
    print(f"DeltaH (kcal/mol): {dh_kcal:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
