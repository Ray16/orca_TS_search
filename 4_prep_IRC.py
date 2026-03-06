import argparse
import os
import re
import shlex
import subprocess
import sys


def find_ts_base(data_dir, base_name=None):
    if base_name:
        base = base_name
        xyz = os.path.join(data_dir, f"{base}.xyz")
        if not os.path.exists(xyz):
            raise FileNotFoundError(f"Missing TS geometry file: {xyz}")
        return base

    xyz_candidates = sorted(
        f[:-4]
        for f in os.listdir(data_dir)
        if f.endswith(".xyz") and "_trj" not in f
    )
    ts_candidates = [b for b in xyz_candidates if b.endswith("_TS_search") or b.endswith("_TS")]
    candidates = ts_candidates if ts_candidates else xyz_candidates
    if not candidates:
        raise FileNotFoundError(f"No suitable .xyz file found in {data_dir}")
    if len(candidates) > 1:
        print(f"Multiple TS candidates found, using: {candidates[0]}")
    return candidates[0]


def build_irc_stem(sysname, base):
    if base.endswith("_TS"):
        return base[:-3]
    return sysname


def parse_charge_mult(inp_path):
    if not os.path.exists(inp_path):
        return -1, 1

    with open(inp_path, "r", encoding="utf-8") as f:
        text = f.read()

    m = re.search(r"^\*\s+xyz(?:file)?\s+(-?\d+)\s+(\d+)", text, flags=re.MULTILINE)
    if not m:
        return -1, 1
    return int(m.group(1)), int(m.group(2))


def parse_nprocs(inp_path):
    if not os.path.exists(inp_path):
        return 1

    with open(inp_path, "r", encoding="utf-8") as f:
        text = f.read()

    m = re.search(r"%pal\b.*?^\s*nprocs\s+(\d+)\s*$.*?^end\s*$", text, flags=re.MULTILINE | re.DOTALL)
    if not m:
        return 1
    return int(m.group(1))


def build_irc_input(base, charge, mult, nprocs, method, basis, maxiter, rel_data_dir):
    return (
        f"! {method} {basis} TightSCF IRC\n\n"
        f"%pal\n"
        f"  nprocs {nprocs}\n"
        f"end\n\n"
        f"%moinp \"{rel_data_dir}/{base}.gbw\"\n\n"
        f"%irc\n"
        f"  Direction both\n"
        f"  MaxIter {maxiter}\n"
        f"  InitHess read\n"
        f"  Hess_Filename \"{rel_data_dir}/{base}.hess\"\n"
        f"end\n\n"
        f"* xyzfile {charge} {mult} {rel_data_dir}/{base}.xyz\n"
    )


def run_orca(target_dir, input_file, output_file, orca_cmd):
    cmd = shlex.split(orca_cmd) + [input_file]
    output_path = output_file if os.path.isabs(output_file) else os.path.join(target_dir, output_file)
    with open(output_path, "w", encoding="utf-8") as out:
        try:
            proc = subprocess.run(cmd, cwd=target_dir, stdout=out, stderr=subprocess.STDOUT)
        except FileNotFoundError:
            return None
    return proc.returncode


def infer_orca_failure_hint(output_file):
    if not os.path.exists(output_file):
        return None
    try:
        with open(output_file, "r", encoding="utf-8", errors="ignore") as f:
            text = f.read(4000)
    except OSError:
        return None

    low = text.lower()
    if "command not found" in low or "not found" in low:
        return "ORCA is not available in this shell. Run `module load orca` first."
    if "error while loading shared libraries" in low:
        return "ORCA was found, but runtime libraries are missing. Load the ORCA module/environment."
    if "cannot execute" in low or "exec format error" in low or "permission denied" in low:
        return "The ORCA executable could not be launched. Check --orca-cmd and permissions."
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Create and optionally run an ORCA IRC job from a TS directory."
    )
    parser.add_argument("system", nargs="?", default=None, help="System name, e.g., sn2 (optional if TS_guess_xyz has exactly one *_input.xyz).")
    parser.add_argument(
        "--base",
        default=None,
        help="Base filename for TS files without extension",
    )
    parser.add_argument("--method", default="B3LYP")
    parser.add_argument("--basis", default="def2-SVP")
    parser.add_argument("--maxiter", type=int, default=70)
    parser.add_argument("--run", action="store_true", help="Run ORCA after writing the IRC input.")
    parser.add_argument(
        "--orca-cmd",
        default="orca",
        help='ORCA executable command, e.g. "orca" or "/path/to/orca".',
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    def resolve_system_name(system_arg, root_dir):
        if system_arg:
            return system_arg
        guess_dir = os.path.join(root_dir, "TS_guess_xyz")
        pattern = "_input.xyz"
        candidates = []
        if os.path.isdir(guess_dir):
            candidates = sorted(
                f[: -len(pattern)]
                for f in os.listdir(guess_dir)
                if f.endswith(pattern)
            )
        if len(candidates) == 1:
            print(f"Auto-detected system from TS_guess_xyz: {candidates[0]}")
            return candidates[0]
        if len(candidates) > 1:
            print("Error: multiple *_input.xyz files found in TS_guess_xyz; pass system explicitly.")
            print("Candidates: " + ", ".join(candidates))
            return 1
        print("Error: no *_input.xyz found in TS_guess_xyz; pass system explicitly.")
        return 1

    system_name = resolve_system_name(args.system, script_dir)
    if system_name == 1:
        return 1
    target_dir = system_name if os.path.isabs(system_name) else os.path.join(script_dir, system_name)
    target_dir = os.path.abspath(target_dir)
    if not os.path.isdir(target_dir):
        print(f"Error: target directory not found: {target_dir}")
        return 1

    ts_search_dir = os.path.join(target_dir, "TS_search")
    data_dir = ts_search_dir if os.path.isdir(ts_search_dir) else target_dir

    try:
        base = find_ts_base(data_dir, args.base)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 1

    ts_inp = os.path.join(data_dir, f"{base}.inp")
    charge, mult = parse_charge_mult(ts_inp)
    nprocs = parse_nprocs(ts_inp)

    required = [f"{base}.xyz", f"{base}.hess", f"{base}.gbw"]
    missing = [f for f in required if not os.path.exists(os.path.join(data_dir, f))]
    if missing:
        print(f"Error: missing required TS files in {data_dir}: {', '.join(missing)}")
        return 1

    irc_dir = os.path.join(target_dir, "IRC")
    os.makedirs(irc_dir, exist_ok=True)

    sysname = os.path.basename(os.path.normpath(target_dir))
    irc_stem = build_irc_stem(sysname, base)
    irc_inp_name = f"{irc_stem}.inp"
    irc_out_name = f"{irc_stem}.out"
    irc_inp_path = os.path.join(irc_dir, irc_inp_name)
    rel_data_dir = os.path.relpath(data_dir, irc_dir).replace(os.sep, "/")

    irc_text = build_irc_input(
        base=base,
        charge=charge,
        mult=mult,
        nprocs=nprocs,
        method=args.method,
        basis=args.basis,
        maxiter=args.maxiter,
        rel_data_dir=rel_data_dir,
    )
    with open(irc_inp_path, "w", encoding="utf-8") as f:
        f.write(irc_text)
    print(f"Wrote IRC input: {irc_inp_path}")

    if not args.run:
        return 0

    code = run_orca(irc_dir, irc_inp_name, irc_out_name, args.orca_cmd)
    if code is None:
        print(
            "Error: ORCA executable was not found in PATH.\n"
            "Load ORCA first (e.g., `module load orca`) or pass --orca-cmd /full/path/to/orca."
        )
        print(f"See {os.path.join(irc_dir, irc_out_name)} for any captured output.")
        return 127
    if code != 0:
        hint = infer_orca_failure_hint(os.path.join(irc_dir, irc_out_name))
        if hint:
            print(f"Hint: {hint}")
        print(f"ORCA IRC run failed with exit code {code}. See {os.path.join(irc_dir, irc_out_name)}")
        return code
    print(f"ORCA IRC run completed. Output: {os.path.join(irc_dir, irc_out_name)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
