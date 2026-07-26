#!/usr/bin/env python
"""Execute the workshop notebook's code cells headlessly, like a Jupyter kernel.

All code cells run in one shared namespace (so ``cfg`` from cell 0 is visible to
later cells), with the notebook's own directory placed on ``sys.path`` exactly as
ipykernel does. Per-cell status is printed so you can see how far the notebook gets
without opening JupyterLab.

Usage:
    python tools/run_cells.py prismatic_workshop.ipynb [--to N] [--cell-timeout S] [--keep-going]

    --to N            run only code cells 0..N (inclusive); default: all
    --cell-timeout S  abort a cell after S seconds and report it as STARTED-OK
                      (use to confirm a heavy cell begins real work without
                      waiting out multi-hour compute); 0 = no timeout
    --keep-going      continue after a failing cell instead of stopping
"""
import argparse
import json
import os
import signal
import sys
import traceback


class CellTimeout(Exception):
    pass


def _alarm(signum, frame):
    raise CellTimeout()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("notebook")
    ap.add_argument("--to", type=int, default=10**9)
    ap.add_argument("--cell-timeout", type=float, default=0.0)
    ap.add_argument("--keep-going", action="store_true")
    args = ap.parse_args()

    nb_path = os.path.abspath(args.notebook)
    nb = json.load(open(nb_path))
    codes = [c for c in nb["cells"] if c["cell_type"] == "code"]

    # Mimic ipykernel: the notebook's directory is importable (finds `initialize`, `utils`).
    sys.path.insert(0, os.path.dirname(nb_path))
    os.chdir(os.path.dirname(nb_path))

    if args.cell_timeout > 0:
        signal.signal(signal.SIGALRM, _alarm)

    ns = {"__name__": "__main__"}
    failures = 0
    for i, c in enumerate(codes):
        if i > args.to:
            break
        src = "".join(c["source"])
        first = next((ln for ln in src.splitlines()
                      if ln.strip() and not ln.strip().startswith("#")), "")
        if args.cell_timeout > 0:
            signal.setitimer(signal.ITIMER_REAL, args.cell_timeout)
        try:
            exec(compile(src, f"<cell {i}>", "exec"), ns)
            print(f"cell {i}: PASS                | {first[:72]}", flush=True)
        except CellTimeout:
            print(f"cell {i}: STARTED-OK (still running after {args.cell_timeout:.0f}s) | {first[:72]}", flush=True)
        except BaseException as e:
            tb = traceback.extract_tb(sys.exc_info()[2])
            loc = ""
            for frame in reversed(tb):
                if frame.filename.startswith("<cell") or "/PRISMATIC_tutorial/" in frame.filename:
                    loc = f"{frame.filename}:{frame.lineno} in {frame.name}"
                    break
            msg = (str(e).splitlines() or [""])[0]
            print(f"cell {i}: FAIL {type(e).__name__}: {msg[:120]}", flush=True)
            print(f"          at {loc}\n          | {first[:72]}", flush=True)
            failures += 1
            if not args.keep_going:
                break
        finally:
            if args.cell_timeout > 0:
                signal.setitimer(signal.ITIMER_REAL, 0)

    print(f"\n=== {failures} failed cell(s) ===", flush=True)
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
