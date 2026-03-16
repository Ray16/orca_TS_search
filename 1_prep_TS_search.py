import argparse
import os
import sys

from workflow_utils import resolve_system_name


# XTB2 (GFN2-xTB): Fast semi-empirical method, suitable for TS search + IRC
# within a couple of minutes even for systems with heavy atoms like Iodine.
# Freq is combined with OptTS so that enthalpy H(TS) is available for Polanyi analysis.
TEMPLATE_HEADER = """! XTB2 OptTS Freq PAL8

%geom
  Calc_Hess true
  Recalc_Hess 3
  MaxIter 200
  Trust -0.1
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


def list_systems_from_ts_guess(root_dir):
    guess_dir = os.path.join(root_dir, "TS_guess_xyz")
    pattern = "_input.xyz"
    if not os.path.isdir(guess_dir):
        raise ValueError("TS_guess_xyz directory not found.")
    candidates = sorted(
        f[: -len(pattern)] for f in os.listdir(guess_dir) if f.endswith(pattern)
    )
    if not candidates:
        raise ValueError("No *_input.xyz files found in TS_guess_xyz.")
    return candidates


def read_systems_file(path):
    if not os.path.exists(path):
        raise ValueError(f"Systems file not found: {path}")
    systems = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            systems.append(line)
    if not systems:
        raise ValueError("No systems found in systems file.")
    return systems


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
        "--all",
        action="store_true",
        help="Generate inputs for all *_input.xyz files in TS_guess_xyz.",
    )
    parser.add_argument(
        "--systems-file",
        default=None,
        help="Path to a file listing system names (one per line).",
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
    parser.add_argument("--charge", type=int, required=True, help="Charge for plain XYZ inputs.")
    parser.add_argument("--mult", type=int, required=True, help="Multiplicity for plain XYZ inputs.")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root = args.root if args.root else script_dir
    root = os.path.abspath(root)
    if args.system and (args.all or args.systems_file):
        raise ValueError("Use either a single system or --all/--systems-file, not both.")
    if (args.all or args.systems_file) and (args.xyz_file or args.xyz_block):
        raise ValueError("--xyz-file/--xyz-block cannot be used with --all/--systems-file.")

    if args.all:
        system_names = list_systems_from_ts_guess(root)
    elif args.systems_file:
        system_names = read_systems_file(args.systems_file)
    else:
        system_names = [resolve_system_name(args.system, root)]

    for system_name in system_names:
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
