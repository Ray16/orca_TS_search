import argparse
import os
import sys

from workflow_utils import resolve_system_name


TEMPLATE_HEADER = """! B3LYP def2-SVP OptTS Freq

%pal
  nprocs 1
end

%geom
  Calc_Hess true
  Recalc_Hess 5
end

"""


def xyz_to_orca_block(xyz_text, charge, mult):
    lines = [line.rstrip() for line in xyz_text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("XYZ content is empty.")

    # Already an ORCA block
    if lines[0].lstrip().startswith("* xyz"):
        if lines[-1].strip() != "*":
            raise ValueError("XYZ block must end with a single '*' line.")
        return "\n".join(lines) + "\n"

    # Standard XYZ format: natoms, comment, atom lines...
    try:
        natoms = int(lines[0].strip())
        start_idx = 2
    except ValueError:
        natoms = None
        start_idx = 0

    atom_lines = lines[start_idx:]
    if natoms is not None and len(atom_lines) < natoms:
        raise ValueError(
            f"XYZ appears to declare {natoms} atoms but only {len(atom_lines)} coordinate lines were found."
        )
    if natoms is not None:
        atom_lines = atom_lines[:natoms]

    if not atom_lines:
        raise ValueError("No atom coordinate lines found in XYZ content.")

    return f"* xyz {charge} {mult}\n" + "\n".join(atom_lines) + "\n*\n"


def read_xyz_block(args, xyz_guess_path):
    if args.xyz_file:
        with open(args.xyz_file, "r", encoding="utf-8") as f:
            block = f.read().strip()
    elif args.xyz_block:
        block = args.xyz_block.strip()
    elif os.path.exists(xyz_guess_path):
        with open(xyz_guess_path, "r", encoding="utf-8") as f:
            block = f.read().strip()
    else:
        # Supports: cat block.txt | python prep_TS_search.py sn2
        block = sys.stdin.read().strip()

    if not block:
        raise ValueError(
            "No geometry provided. Expected TS_guess_xyz/<system>_input.xyz, or use --xyz-file/--xyz-block/stdin."
        )
    return xyz_to_orca_block(block, args.charge, args.mult)


def build_input_text(xyz_block):
    return TEMPLATE_HEADER + xyz_block


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate <system>/TS_search/<system>.inp for ORCA TS search."
        )
    )
    parser.add_argument(
        "system",
        nargs="?",
        default=None,
        help="System name, e.g., sn2 (optional if TS_guess_xyz has exactly one *_input.xyz).",
    )
    parser.add_argument(
        "--xyz-file",
        help="Path to a file that contains the full xyz block including '* xyz ...' and trailing '*'.",
    )
    parser.add_argument(
        "--xyz-block",
        help="Full xyz block as a single quoted string.",
    )
    parser.add_argument(
        "--root",
        default=None,
        help="Optional root directory (defaults to script directory).",
    )
    parser.add_argument("--charge", type=int, default=-1, help="Default charge for plain XYZ inputs.")
    parser.add_argument("--mult", type=int, default=1, help="Default multiplicity for plain XYZ inputs.")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root = args.root if args.root else script_dir
    root = os.path.abspath(root)
    system_name = resolve_system_name(args.system, root)

    system_dir = os.path.join(root, system_name)
    out_dir = os.path.join(system_dir, "TS_search")
    out_name = f"{system_name}.inp"
    out_path = os.path.join(out_dir, out_name)
    xyz_guess_path = os.path.join(root, "TS_guess_xyz", f"{system_name}_input.xyz")

    xyz_block = read_xyz_block(args, xyz_guess_path)
    text = build_input_text(xyz_block)

    os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(text)

    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
