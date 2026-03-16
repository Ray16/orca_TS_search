import argparse
import os
import re

from workflow_utils import resolve_system_dir


def parse_charge_mult(inp_path):
    if not inp_path or not os.path.exists(inp_path):
        return None, None
    with open(inp_path, "r", encoding="utf-8") as f:
        text = f.read()
    m = re.search(r"^\*\s+xyz(?:file)?\s+(-?\d+)\s+(\d+)", text, flags=re.MULTILINE)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))


def parse_nprocs(inp_path):
    if not inp_path or not os.path.exists(inp_path):
        return None
    with open(inp_path, "r", encoding="utf-8") as f:
        text = f.read()
    m = re.search(r"%pal\b.*?^\s*nprocs\s+(\d+)\s*$.*?^end\s*$", text, flags=re.MULTILINE | re.DOTALL)
    if not m:
        return None
    return int(m.group(1))


def find_ts_input(system_dir):
    ts_dir = os.path.join(system_dir, "TS_search")
    if not os.path.isdir(ts_dir):
        return None
    inputs = sorted(
        os.path.join(ts_dir, f) for f in os.listdir(ts_dir) if f.endswith(".inp")
    )
    return inputs[0] if inputs else None


def read_multiframe_xyz(xyz_path):
    if not os.path.exists(xyz_path):
        raise FileNotFoundError(f"Missing IRC trajectory: {xyz_path}")

    with open(xyz_path, "r", encoding="utf-8") as f:
        lines = [line.rstrip("\n") for line in f]

    frames = []
    i = 0
    total = len(lines)
    while i < total:
        while i < total and not lines[i].strip():
            i += 1
        if i >= total:
            break
        try:
            natoms = int(lines[i].strip())
        except ValueError as exc:
            raise ValueError(f"Invalid XYZ frame at line {i + 1}: {lines[i]}") from exc
        i += 1
        if i >= total:
            raise ValueError("Unexpected end of file after atom count.")
        comment = lines[i]
        i += 1
        if i + natoms > total:
            raise ValueError("XYZ frame is truncated; not enough atom lines.")
        atoms = lines[i : i + natoms]
        i += natoms
        frames.append((natoms, comment, atoms))

    if len(frames) < 2:
        raise ValueError("Expected at least two frames in IRC trajectory.")
    return frames


def build_orca_input(atoms, charge, mult, method, basis, nprocs):
    if method.upper() == "XTB2":
        return (
            f"! XTB2 Opt Freq PAL8\n\n"
            f"%geom\n"
            f"  MaxIter 500\n"
            f"end\n\n"
            f"* xyz {charge} {mult}\n"
            + "\n".join(atoms)
            + "\n*\n"
        )
    return (
        f"! {method} {basis} Opt Freq TightSCF\n\n"
        f"%pal\n"
        f"  nprocs {nprocs}\n"
        f"end\n\n"
        f"* xyz {charge} {mult}\n"
        + "\n".join(atoms)
        + "\n*\n"
    )


def write_xyz(path, natoms, comment, atoms):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{natoms}\n")
        f.write(f"{comment}\n")
        f.write("\n".join(atoms))
        f.write("\n")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare ORCA Opt+Freq inputs for reactant/product from IRC trajectory endpoints."
    )
    parser.add_argument(
        "system",
        nargs="?",
        default=None,
        help="System name, e.g., sn2 (optional if TS_guess_xyz has exactly one *_input.xyz).",
    )
    parser.add_argument("--root", default=None, help="Optional root directory (defaults to script directory).")
    parser.add_argument("--irc-xyz", default=None, help="Path to IRC Full trajectory XYZ.")
    parser.add_argument("--method", default="XTB2")
    parser.add_argument("--basis", default="")
    parser.add_argument("--charge", type=int, default=None)
    parser.add_argument("--mult", type=int, default=None)
    parser.add_argument("--nprocs", type=int, default=None)
    parser.add_argument("--outdir", default=None, help="Output directory (defaults to <system>/geo_opt_reactant_product).")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root = os.path.abspath(args.root) if args.root else script_dir
    system_dir = resolve_system_dir(args.system, root)
    sysname = os.path.basename(os.path.normpath(system_dir))

    irc_xyz = args.irc_xyz or os.path.join(system_dir, "IRC", f"{sysname}_IRC_Full_trj.xyz")
    frames = read_multiframe_xyz(irc_xyz)
    reactant = frames[0]
    product = frames[-1]

    ts_inp = find_ts_input(system_dir)
    charge, mult = parse_charge_mult(ts_inp)
    nprocs = parse_nprocs(ts_inp)

    charge = args.charge if args.charge is not None else charge
    mult = args.mult if args.mult is not None else mult
    nprocs = args.nprocs if args.nprocs is not None else (nprocs if nprocs is not None else 1)

    if charge is None or mult is None:
        raise SystemExit("Error: charge/multiplicity not found. Pass --charge and --mult.")

    out_dir = args.outdir or os.path.join(system_dir, "geo_opt_reactant_product")
    os.makedirs(out_dir, exist_ok=True)

    reactant_xyz = os.path.join(out_dir, "reactant.xyz")
    product_xyz = os.path.join(out_dir, "product.xyz")
    write_xyz(reactant_xyz, *reactant)
    write_xyz(product_xyz, *product)

    reactant_inp = os.path.join(out_dir, "reactant_optfreq.inp")
    product_inp = os.path.join(out_dir, "product_optfreq.inp")

    with open(reactant_inp, "w", encoding="utf-8") as f:
        f.write(build_orca_input(reactant[2], charge, mult, args.method, args.basis, nprocs))
    with open(product_inp, "w", encoding="utf-8") as f:
        f.write(build_orca_input(product[2], charge, mult, args.method, args.basis, nprocs))

    print(f"Wrote: {reactant_inp}")
    print(f"Wrote: {product_inp}")
    print(f"Wrote: {reactant_xyz}")
    print(f"Wrote: {product_xyz}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
