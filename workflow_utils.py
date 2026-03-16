import math
import os
import re
import sys


def script_root(anchor_file):
    return os.path.dirname(os.path.abspath(anchor_file))


def detect_system_from_ts_guess(root_dir):
    guess_dir = os.path.join(root_dir, "TS_guess_xyz")
    pattern = "_input.xyz"
    if not os.path.isdir(guess_dir):
        raise SystemExit(
            "Error: TS_guess_xyz directory not found; pass system explicitly."
        )

    candidates = sorted(
        f[: -len(pattern)] for f in os.listdir(guess_dir) if f.endswith(pattern)
    )
    if len(candidates) == 1:
        print(f"Auto-detected system from TS_guess_xyz: {candidates[0]}")
        return candidates[0]
    if len(candidates) > 1:
        raise SystemExit(
            "Error: multiple *_input.xyz files found in TS_guess_xyz; pass system explicitly.\n"
            f"Candidates: {', '.join(candidates)}"
        )
    raise SystemExit("Error: no *_input.xyz found in TS_guess_xyz; pass system explicitly.")


def resolve_system_name(system_arg, root_dir):
    return system_arg if system_arg else detect_system_from_ts_guess(root_dir)


def resolve_system_dir(system_arg, root_dir):
    system_name = resolve_system_name(system_arg, root_dir)
    if os.path.isabs(system_name):
        target_dir = system_name
    else:
        target_dir = os.path.abspath(os.path.join(root_dir, system_name))
        if not os.path.isdir(target_dir):
            runs_candidate = os.path.abspath(os.path.join(root_dir, "runs", system_name))
            if os.path.isdir(runs_candidate):
                target_dir = runs_candidate
    if not os.path.isdir(target_dir):
        raise SystemExit(f"Error: target directory not found: {target_dir}")
    return target_dir


# ---------------------------------------------------------------------------
# Reactant/product swap detection from endpoint geometries
# ---------------------------------------------------------------------------

_HALOGEN_ELEMENT = {"f": "F", "cl": "Cl", "br": "Br", "i": "I", "h": "H"}


def _parse_xyz_atoms(xyz_path):
    """Return list of (element, x, y, z) from a single-frame XYZ file."""
    atoms = []
    with open(xyz_path, "r", encoding="utf-8") as fh:
        lines = fh.readlines()
    natoms = int(lines[0].strip())
    for line in lines[2 : 2 + natoms]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms


def _distance(a, b):
    return math.sqrt((a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2)


def _find_atom(atoms, element):
    """Find the first atom matching element (case-insensitive)."""
    for a in atoms:
        if a[0].upper() == element.upper():
            return a
    return None


def _sn2_should_swap(sysname, system_dir):
    """For sn2_{nuc}_{lg}, the true reactant has the leaving group (lg) bonded
    to C and the nucleophile (nuc) far from C.  If the backward-endpoint
    geometry (reactant.xyz) has nuc closer to C than lg, the labels are
    swapped."""
    m = re.match(r"sn2_([a-z]+)_([a-z]+)", sysname)
    if not m:
        return None
    nuc_key, lg_key = m.group(1), m.group(2)
    nuc_elem = _HALOGEN_ELEMENT.get(nuc_key)
    lg_elem = _HALOGEN_ELEMENT.get(lg_key)
    if not nuc_elem or not lg_elem:
        return None
    if nuc_elem == lg_elem:
        return False  # identity reaction

    react_xyz = os.path.join(system_dir, "geo_opt_reactant_product", "reactant.xyz")
    if not os.path.exists(react_xyz):
        return None

    atoms = _parse_xyz_atoms(react_xyz)
    c_atom = _find_atom(atoms, "C")
    nuc_atom = _find_atom(atoms, nuc_elem)
    lg_atom = _find_atom(atoms, lg_elem)
    if not all([c_atom, nuc_atom, lg_atom]):
        return None

    d_nuc = _distance(c_atom, nuc_atom)
    d_lg = _distance(c_atom, lg_atom)
    return d_nuc < d_lg


def _da_should_swap(sysname, system_dir):
    """Diels-Alder: the product (cyclohexene ring) has more C-C bonds than the
    reactant (separated diene + dienophile).  If reactant.xyz has more C-C
    bonds (< 2.0 Å), the labels are swapped."""
    react_xyz = os.path.join(system_dir, "geo_opt_reactant_product", "reactant.xyz")
    prod_xyz = os.path.join(system_dir, "geo_opt_reactant_product", "product.xyz")
    if not os.path.exists(react_xyz) or not os.path.exists(prod_xyz):
        return None

    def _count_cc_bonds(atoms):
        carbons = [a for a in atoms if a[0] == "C"]
        count = 0
        for i in range(len(carbons)):
            for j in range(i + 1, len(carbons)):
                if _distance(carbons[i], carbons[j]) < 2.0:
                    count += 1
        return count

    atoms_r = _parse_xyz_atoms(react_xyz)
    atoms_p = _parse_xyz_atoms(prod_xyz)
    cc_r = _count_cc_bonds(atoms_r)
    cc_p = _count_cc_bonds(atoms_p)
    # If "reactant" has more C-C bonds, it's actually the product (ring)
    return cc_r > cc_p


def should_swap_rp(sysname, system_dir):
    """Determine if the IRC backward/forward endpoint labels ('reactant'/'product')
    need to be swapped for a given system.

    Dispatches to reaction-family-specific geometry checks:
      - SN2:  nucleophile should be far from C in the true reactant
      - DA:   product (cyclohexene ring) has more C-C bonds than reactant

    Returns True if swap needed, False if not, None if unable to determine.
    """
    if sysname.startswith("sn2_"):
        return _sn2_should_swap(sysname, system_dir)
    if sysname.startswith("da_"):
        return _da_should_swap(sysname, system_dir)
    return None


def main_guard(main_fn):
    try:
        code = main_fn()
    except SystemExit:
        raise
    except Exception as exc:
        print(f"Error: {exc}")
        raise SystemExit(1)
    raise SystemExit(code if isinstance(code, int) else 0)
