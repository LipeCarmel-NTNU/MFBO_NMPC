"""Run one case from end to end.

This is the only Python command you type. It runs the design phase, waits, and
then runs the optimization phase.

    python run_pipeline.py --case case1

Start MATLAB on main_initialization first. MATLAB stays on that one command for
the whole run. When the design phase finishes, the driver sends its first
optimization request. main_initialization sees the phase code on that request,
stops serving, and calls main_BO, which serves the rest of the run.

You can also run one phase alone:

    python run_pipeline.py --case case1 --phase init
    python run_pipeline.py --case case1 --phase bo

Both phases resume. If a run stops, start MATLAB again and repeat the same
command.
"""

from __future__ import annotations

import argparse
import sys
import time

from pathlib import Path

from pipeline.console_log import ConsoleLog
from pipeline.driver import main as run_phase
from run_config import CASES

# Both phases copy their console output here, one file for each start.
LOG_ROOT = Path(__file__).resolve().parent / "results"

# The pause covers the handover. MATLAB has to notice the optimization request,
# leave the design loop, clear its workspace and start main_BO. The driver
# would wait for the result anyway, so this only keeps the console readable.
HANDOVER_PAUSE_S = 10.0


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", default="case1", choices=sorted(CASES),
                   help="which search space to optimize (default: case1)")
    p.add_argument("--phase", default="both", choices=("both", "init", "bo"),
                   help="run one phase instead of the whole case")
    p.add_argument("--pause", type=float, default=HANDOVER_PAUSE_S,
                   help="seconds to wait between the two phases")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse(sys.argv[1:] if argv is None else argv)

    if args.phase in ("both", "init"):
        with ConsoleLog(LOG_ROOT, "init"):
            print("=" * 70)
            print(f"PHASE 1 of 2: design of experiments, {args.case}")
            print("=" * 70)
            code = run_phase(["init", "--case", args.case])
            if code != 0:
                print(f"\n[run_pipeline] the design phase stopped with code {code}. "
                      f"The optimization phase does not start.")
                return code
            if args.phase == "init":
                return 0
            print(f"\n[run_pipeline] design complete. Waiting {args.pause:.0f} s for "
                  f"MATLAB to hand over to main_BO.")
        time.sleep(args.pause)

    if args.phase in ("both", "bo"):
        with ConsoleLog(LOG_ROOT, "bo"):
            print("=" * 70)
            print(f"PHASE 2 of 2: optimization, {args.case}")
            print("=" * 70)
            code = run_phase(["bo", "--case", args.case])
            if code != 0:
                print(f"\n[run_pipeline] the optimization phase stopped with code {code}.")
                return code
            print("\n[run_pipeline] done. You can stop MATLAB now.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
