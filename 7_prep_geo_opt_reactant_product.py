import argparse
import math
import os
import re

from workflow_utils import resolve_system_dir

# Map element keys used in system names to full element symbols
_HALOGEN_ELEMENT = {"f": "F", "cl": "Cl", "br": "Br", "i": "I", "h": "H"}


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


def _parse_atom_coords(atom_lines):
    """Parse atom lines into list of (element, x, y, z)."""
    atoms = []
    for line in atom_lines:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms


def _dist(a, b):
    return math.sqrt((a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2)


def _find_atom(atoms, element):
    for a in atoms:
        if a[0].upper() == element.upper():
            return a
    return None


def _extract_energy_from_comment(comment):
    """Extract energy from IRC trajectory comment line."""
    m = re.search(r"\bE\s+(-?\d+\.\d+)", comment)
    return float(m.group(1)) if m else None


def _should_swap_endpoints(sysname, endpoint_a, endpoint_b):
    """Determine if IRC endpoints need swapping so that endpoint_a is the
    true reactant and endpoint_b is the true product.

    For SN2 (sn2_{nuc}_{lg}): the true reactant has the leaving group bonded
    to C and the nucleophile far from C.  If endpoint_a has the nucleophile
    closer to C, it's actually the product — swap needed.

    For DA: always exothermic, so the true reactant has higher energy.
    If endpoint_a has lower energy, it's the product — swap needed.

    Returns True if swap needed, False otherwise.
    """
    _, comment_a, atoms_a_lines = endpoint_a
    _, comment_b, atoms_b_lines = endpoint_b

    m = re.match(r"sn2_([a-z]+)_([a-z]+)", sysname)
    if m:
        nuc_elem = _HALOGEN_ELEMENT.get(m.group(1))
        lg_elem = _HALOGEN_ELEMENT.get(m.group(2))
        if nuc_elem and lg_elem and nuc_elem != lg_elem:
            atoms = _parse_atom_coords(atoms_a_lines)
            c_atom = _find_atom(atoms, "C")
            nuc_atom = _find_atom(atoms, nuc_elem)
            lg_atom = _find_atom(atoms, lg_elem)
            if all([c_atom, nuc_atom, lg_atom]):
                d_nuc = _dist(c_atom, nuc_atom)
                d_lg = _dist(c_atom, lg_atom)
                # True reactant: nuc far from C, lg close to C
                # If nuc is closer, this endpoint is the product
                return d_nuc < d_lg
        return False

    if sysname.startswith("da_"):
        # DA product (cyclohexene) has more C-C bonds than the reactant
        # (separated diene + dienophile).  Count C-C bonds < 2.0 Å in each
        # endpoint: the one with more is the product.
        atoms_a = _parse_atom_coords(atoms_a_lines)
        atoms_b = _parse_atom_coords(atoms_b_lines)

        def _count_cc_bonds(atoms):
            carbons = [a for a in atoms if a[0] == "C"]
            count = 0
            for i in range(len(carbons)):
                for j in range(i + 1, len(carbons)):
                    if _dist(carbons[i], carbons[j]) < 2.0:
                        count += 1
            return count

        cc_a = _count_cc_bonds(atoms_a)
        cc_b = _count_cc_bonds(atoms_b)
        # If endpoint_a has more C-C bonds, it's the product (ring) → swap
        return cc_a > cc_b

    return False


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
    endpoint_a = frames[0]
    endpoint_b = frames[-1]

    if _should_swap_endpoints(sysname, endpoint_a, endpoint_b):
        reactant = endpoint_b
        product = endpoint_a
        print(f"  Note: swapped IRC endpoints for {sysname} (geometry-based)")
    else:
        reactant = endpoint_a
        product = endpoint_b

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
