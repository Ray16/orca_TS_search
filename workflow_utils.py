import os
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


def main_guard(main_fn):
    try:
        code = main_fn()
    except SystemExit:
        raise
    except Exception as exc:
        print(f"Error: {exc}")
        raise SystemExit(1)
    raise SystemExit(code if isinstance(code, int) else 0)
