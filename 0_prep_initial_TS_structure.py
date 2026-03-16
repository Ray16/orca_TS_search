import argparse
import os
import math

# Traditional values for comparison
LEGACY_X_DIST = {"F": 1.70, "Cl": 2.00, "Br": 2.20, "I": 2.50}
# Standard C-X single-bond lengths (Angstroms)
C_X = {"F": 1.38, "Cl": 1.78, "Br": 1.94, "I": 2.14}


# Improved SN2 TS parameters based on optimized geometries.
# Base scale for symmetric SN2 transition states.
SN2_TS_SCALE_DEFAULT = 1.31

# Asymmetric ratios { (StrongerBondAtom, WeakerBondAtom): (Ratio_Stronger, Ratio_Weaker) }
# Bond Strength Rank: F > Cl > Br > I
# Ratios tuned for XTB2 (GFN2-xTB) TS searches.
# Key constraint: absolute C-X distances must be distinct enough
# for the optimizer to find the correct saddle point.
SN2_ASYM_RATIOS = {
    ("F", "Cl"): (1.38, 1.20),  # C-F~1.90, C-Cl~2.14
    ("F", "Br"): (1.38, 1.20),  # C-F~1.90, C-Br~2.33
    ("F", "I"): (1.38, 1.30),   # C-F~1.90, C-I~2.78
    ("Cl", "Br"): (1.40, 1.23), # C-Cl~2.49, C-Br~2.39 (works)
    ("Cl", "I"): (1.40, 1.25),  # C-Cl~2.49, C-I~2.68
    ("Br", "I"): (1.45, 1.20),  # C-Br~2.81, C-I~2.57
}

SN2_ASYM_STEP_DEFAULT = 0.06
SN2_LG_RANK = {"F": 0, "Cl": 1, "Br": 2, "I": 3}


def write_xyz(path, atoms, comment, overwrite=False):
    if os.path.exists(path) and not overwrite:
        print(f"Skip (exists): {path}")
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{len(atoms)}\n")
        f.write(comment + "\n")
        for sym, x, y, z in atoms:
            f.write(f"{sym:<2}  {x: .6f}  {y: .6f}  {z: .6f}\n")
    print(f"Wrote: {path}")


def get_sn2_ratios(x1, x2, ts_scale=SN2_TS_SCALE_DEFAULT):
    if x1 == x2:
        return ts_scale, ts_scale

    # Sort by bond strength (F > Cl > Br > I)
    order = ["F", "Cl", "Br", "I"]
    idx1, idx2 = order.index(x1), order.index(x2)
    
    if idx1 < idx2: # x1 is stronger
        key = (x1, x2)
        if key in SN2_ASYM_RATIOS:
            return SN2_ASYM_RATIOS[key]
    else: # x2 is stronger
        key = (x2, x1)
        if key in SN2_ASYM_RATIOS:
            r2, r1 = SN2_ASYM_RATIOS[key]
            return r1, r2
            
    # Fallback if not in table
    return ts_scale, ts_scale


def sn2_axial_distances(
    x1,
    x2,
    ts_scale=SN2_TS_SCALE_DEFAULT,
    asym_step=SN2_ASYM_STEP_DEFAULT,
    symmetric=False,
    legacy=False,
):
    if legacy:
        return LEGACY_X_DIST[x1], LEGACY_X_DIST[x2]

    if not symmetric:
        r1, r2 = get_sn2_ratios(x1, x2, ts_scale)
        return C_X[x1] * r1, C_X[x2] * r2

    d1 = C_X[x1] * ts_scale
    d2 = C_X[x2] * ts_scale
    return d1, d2


def gen_sn2(
    x1,
    x2,
    ts_scale=SN2_TS_SCALE_DEFAULT,
    asym_step=SN2_ASYM_STEP_DEFAULT,
    symmetric=False,
    legacy=False,
):
    d1, d2 = sn2_axial_distances(
        x1,
        x2,
        ts_scale=ts_scale,
        asym_step=asym_step,
        symmetric=symmetric,
        legacy=legacy,
    )
    
    # Asymmetric TS has a pyramidalized methyl group.
    # The carbon shifts towards the closer (more bonded) atom.
    # In an endothermic reaction, it looks like products (C is closer to Nu).
    # d1 is Nu distance, d2 is LG distance.
    shift = (d2 - d1) * 0.4 # Empirical shift factor
    
    # methyl planar parameters (from optimized CH3F/F TS)
    # H distance ~1.085, angle ~120
    h_dist = 1.085
    h_y = h_dist * math.cos(math.radians(30)) # 0.940
    h_x_off = h_dist * math.sin(math.radians(30)) # 0.542
    
    return [
        ("C", 0.0, 0.0, shift),
        ("H", 0.0, h_dist, 0.0),
        ("H", h_y, -h_x_off, 0.0),
        ("H", -h_y, -h_x_off, 0.0),
        (x1, 0.0, 0.0, d1 + shift),
        (x2, 0.0, 0.0, -d2 + shift),
    ]


# ---------------------------------------------------------------------------
# Diels-Alder TS: 1,3-butadiene + CH₂=CHX → 4-X-cyclohex-2-ene
#   charge=0, mult=1  (neutral, closed-shell, concerted pericyclic)
#
# Geometry: s-cis diene (xy plane) + dienophile above, forming a boat-like
# 6-membered ring.  Two new C-C bonds (C1-C6, C4-C5) at ~2.2 Å.
# ---------------------------------------------------------------------------
DA_FORM_DIST = 2.20      # forming C-C bond distance at TS (Å)
DA_CC_TS = 1.40           # C-C bond lengths in diene at TS (between single/double)
DA_DIENOPHILE_CC = 1.38   # dienophile C=C at TS (slightly stretched)
DA_CCC_ANGLE = 120.0      # C-C-C angle in diene at TS
DA_CH = 1.08              # C-H bond length


def _normalize(v):
    d = math.sqrt(sum(c * c for c in v))
    return tuple(c / d for c in v) if d > 1e-12 else v


def _cross(a, b):
    return (a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])


def _add(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])


def _scale(v, s):
    return (v[0]*s, v[1]*s, v[2]*s)


def _sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])


def _place_tetrahedral_h(center, bond1_vec, bond2_vec, n_h):
    """Place 1 or 2 H atoms to complete roughly tetrahedral geometry at center."""
    b1 = _normalize(bond1_vec)
    b2 = _normalize(bond2_vec)
    # Bisector of the two existing bonds (points "into" bonds)
    bis = _normalize(_add(b1, b2))
    # Anti-bisector (points "away" from bonds) — H atoms go here
    anti = _scale(bis, -1.0)
    # Normal to the plane of the two bonds
    n = _normalize(_cross(bond1_vec, bond2_vec))
    if n_h == 1:
        return [_add(center, _scale(anti, DA_CH))]
    # 2 H atoms: spread along normal direction
    spread = 0.62  # controls angular spread
    d1 = _normalize(_add(anti, _scale(n, spread)))
    d2 = _normalize(_add(anti, _scale(n, -spread)))
    return [_add(center, _scale(d1, DA_CH)),
            _add(center, _scale(d2, DA_CH))]


def _place_trigonal_h(center, bond1_vec, bond2_vec):
    """Place 1 H atom in the plane of two bonds (sp2), pointing outward."""
    b1 = _normalize(bond1_vec)
    b2 = _normalize(bond2_vec)
    anti = _scale(_normalize(_add(b1, b2)), -1.0)
    return [_add(center, _scale(anti, DA_CH))]


def gen_da(x=None):
    """
    TS guess for Diels-Alder: 1,3-butadiene + CH₂=CHX → 4-X-cyclohex-2-ene.
    If x is None, use unsubstituted ethylene (X=H).
    Arrangement:
      Diene C1=C2-C3=C4 in s-cis conformation (xy plane, z=0).
      Dienophile C5=C6 above, forming bonds C4-C5 and C1-C6.
    """
    d = DA_CC_TS
    ang = math.radians(180.0 - DA_CCC_ANGLE)

    # Diene carbons (s-cis, xy plane)
    C1 = (0.0, 0.0, 0.0)
    C2 = (d, 0.0, 0.0)
    C3 = (C2[0] + d * math.cos(ang), d * math.sin(ang), 0.0)
    # C4: rotate C3→C2 direction by DA_CCC_ANGLE for s-cis
    dx = C2[0] - C3[0]
    dy = C2[1] - C3[1]
    dd = math.sqrt(dx*dx + dy*dy)
    ux, uy = dx/dd, dy/dd
    ca, sa = math.cos(math.radians(DA_CCC_ANGLE)), math.sin(math.radians(DA_CCC_ANGLE))
    vx = ux * ca - uy * sa
    vy = ux * sa + uy * ca
    C4 = (C3[0] + d * vx, C3[1] + d * vy, 0.0)

    # Dienophile above diene
    mx = (C1[0] + C4[0]) / 2
    my = (C1[1] + C4[1]) / 2
    d14 = math.sqrt((C4[0]-C1[0])**2 + (C4[1]-C1[1])**2)
    ux14 = (C4[0]-C1[0]) / d14
    uy14 = (C4[1]-C1[1]) / d14
    hd = DA_DIENOPHILE_CC / 2
    C5xy = (mx + hd * ux14, my + hd * uy14)
    C6xy = (mx - hd * ux14, my - hd * uy14)
    # Height so C4-C5 = DA_FORM_DIST
    dxy2 = (C4[0]-C5xy[0])**2 + (C4[1]-C5xy[1])**2
    h = math.sqrt(DA_FORM_DIST**2 - dxy2)
    C5 = (C5xy[0], C5xy[1], h)
    C6 = (C6xy[0], C6xy[1], h)

    atoms = []
    # Ring carbons
    ring = [C1, C2, C3, C4, C5, C6]
    for c in ring:
        atoms.append(("C", c[0], c[1], c[2]))

    # H on C1 (sp3-like, 2 H): bonds to C2 and C6
    for pos in _place_tetrahedral_h(C1, _sub(C2, C1), _sub(C6, C1), 2):
        atoms.append(("H", pos[0], pos[1], pos[2]))

    # H on C2 (sp2-like, 1 H): bonds to C1 and C3
    for pos in _place_trigonal_h(C2, _sub(C1, C2), _sub(C3, C2)):
        atoms.append(("H", pos[0], pos[1], pos[2]))

    # H on C3 (sp2-like, 1 H): bonds to C2 and C4
    for pos in _place_trigonal_h(C3, _sub(C2, C3), _sub(C4, C3)):
        atoms.append(("H", pos[0], pos[1], pos[2]))

    # H on C4 (sp3-like, 2 H): bonds to C3 and C5
    for pos in _place_tetrahedral_h(C4, _sub(C3, C4), _sub(C5, C4), 2):
        atoms.append(("H", pos[0], pos[1], pos[2]))

    # H on C6 (sp3-like, 2 H): bonds to C5 and C1
    for pos in _place_tetrahedral_h(C6, _sub(C5, C6), _sub(C1, C6), 2):
        atoms.append(("H", pos[0], pos[1], pos[2]))

    # H (or X) on C5: bonds to C6 and C4
    c5_h_positions = _place_tetrahedral_h(C5, _sub(C6, C5), _sub(C4, C5), 2)
    if x is None:
        # Unsubstituted: 2 H on C5
        for pos in c5_h_positions:
            atoms.append(("H", pos[0], pos[1], pos[2]))
    else:
        # Substituted: 1 H + 1 X on C5
        # First position: H (stays)
        atoms.append(("H", c5_h_positions[0][0], c5_h_positions[0][1], c5_h_positions[0][2]))
        # Second position: X at C-X bond length (instead of C-H)
        x_dir = _normalize(_sub(c5_h_positions[1], C5))
        x_pos = _add(C5, _scale(x_dir, C_X[x]))
        atoms.append((x, x_pos[0], x_pos[1], x_pos[2]))

    return atoms


def main():
    parser = argparse.ArgumentParser(
        description="Generate initial TS-guess XYZ files."
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory (defaults to <script_dir>/TS_guess_xyz).",
    )
    parser.add_argument(
        "--family",
        choices=["all", "sn2", "da"],
        default="all",
        help="Which family to generate.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing files.",
    )
    parser.add_argument(
        "--sn2-legacy",
        action="store_true",
        help="Use legacy, shorter C-X distances for SN2 (more reactant-like).",
    )
    parser.add_argument(
        "--sn2-symmetric",
        action="store_true",
        help="Force symmetric axial distances for mixed SN2 pairs.",
    )
    parser.add_argument(
        "--sn2-ts-scale",
        type=float,
        default=SN2_TS_SCALE_DEFAULT,
        help="Scale factor applied to C-X single-bond lengths for SN2 TS guesses.",
    )
    parser.add_argument(
        "--sn2-asym-step",
        type=float,
        default=SN2_ASYM_STEP_DEFAULT,
        help="Asymmetry (angstrom) per leaving-group rank difference for SN2.",
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    outdir = args.outdir if args.outdir else os.path.join(script_dir, "TS_guess_xyz")
    outdir = os.path.abspath(outdir)

    if args.family in ("all", "sn2"):
        # Generate all 12 asymmetric pairs and 4 symmetric pairs
        halogens = ["F", "Cl", "Br", "I"]
        sn2_pairs = []
        for x1 in halogens:
            for x2 in halogens:
                if x1 == x2 and not args.sn2_symmetric:
                    # Original script didn't explicitly list symmetric pairs in sn2_pairs,
                    # but we add them here for completeness if needed.
                    continue
                sn2_pairs.append((x1, x2))
        
        # If we want to strictly match the original list of 8 pairs:
        # sn2_pairs = [("F", "Cl"), ("F", "Br"), ("F", "I"), ("Cl", "F"), ("Cl", "Br"), ("Cl", "I"), ("Br", "F"), ("Br", "I")]
        
        for x1, x2 in sn2_pairs:
            name = f"sn2_{x1.lower()}_{x2.lower()}_input.xyz"
            comment = f"sn2 {x1}/{x2} TS guess"
            write_xyz(
                os.path.join(outdir, name),
                gen_sn2(
                    x1,
                    x2,
                    ts_scale=args.sn2_ts_scale,
                    asym_step=args.sn2_asym_step,
                    symmetric=args.sn2_symmetric,
                    legacy=args.sn2_legacy,
                ),
                comment,
                overwrite=args.overwrite,
            )

    if args.family in ("all", "da"):
        # Halogen-substituted ethylenes
        for x in ["F", "Cl", "Br"]:
            name = f"da_{x.lower()}_input.xyz"
            comment = f"Diels-Alder TS: butadiene + CH2=CH{x}"
            write_xyz(
                os.path.join(outdir, name),
                gen_da(x),
                comment,
                overwrite=args.overwrite,
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
