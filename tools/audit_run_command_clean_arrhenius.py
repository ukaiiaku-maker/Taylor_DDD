#!/usr/bin/env python3
from pathlib import Path
import argparse
import shlex
import sys

def get_flag(toks, flag):
    if flag not in toks:
        return None
    i = toks.index(flag)
    if i + 1 >= len(toks):
        return ""
    return toks[i+1]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--cross-scale", required=True)
    ap.add_argument("--peierls-scale", required=True)
    ap.add_argument("--pk-cap", default="0")
    ap.add_argument("--pileup-mode", default="same_obstacle_global")
    ap.add_argument("--pileup-max", default="16")
    args = ap.parse_args()

    root = Path(args.root)
    files = sorted(root.glob("T*/run_command.txt"))
    if not files:
        raise SystemExit(f"No run_command.txt files found under {root}")

    failures = []
    for f in files:
        toks = shlex.split(f.read_text(errors="ignore"))

        expected = {
            "--expfit-cross-scale": args.cross_scale,
            "--expfit-peierls-scale": args.peierls_scale,
            "--pk-reaction-cap-leff-to-x-factor": args.pk_cap,
            "--pk-reaction-pileup-force-mode": args.pileup_mode,
            "--pk-reaction-max-pileup-contributors": args.pileup_max,
            "--max-pinned-nodes-per-obstacle": "1",
            "--occupied-junction-mode": "queue_pin",
            "--saturated-obstacle-mode": "share_pin",
            "--plastic-strain-mode": "swept_area",
        }

        for flag, val in expected.items():
            got = get_flag(toks, flag)
            if got != val:
                failures.append((str(f), flag, val, got))

        if "--pk-reaction-cap-Leff-to-X-factor" in toks:
            failures.append((str(f), "--pk-reaction-cap-Leff-to-X-factor", "NOT PRESENT", get_flag(toks, "--pk-reaction-cap-Leff-to-X-factor")))

        if "--no-taylor-crossing-prefactor" not in toks:
            failures.append((str(f), "--no-taylor-crossing-prefactor", "PRESENT", "MISSING"))

    if failures:
        print("RUN COMMAND AUDIT FAILED")
        for f, flag, exp, got in failures:
            print(f"{f}: {flag}: expected {exp}, got {got}")
        sys.exit(2)

    print("RUN COMMAND AUDIT PASSED")
    print(f"checked {len(files)} case command(s) under {root}")

if __name__ == "__main__":
    main()
