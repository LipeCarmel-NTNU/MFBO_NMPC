"""Run one case with Python owning the MATLAB server.

    python run_supervised.py --case baseline

This is run_pipeline.py with the MATLAB half started, watched and relaunched by
the driver's own process. You type one command instead of two, and you do not
start MATLAB yourself.

The campaign writes to its own folder, results/<arm>_<case>_s<sobol_seed>, which
pipeline/matlab_interface.py derives from --case unless MFBO_RESULTS_DIR is set.
MATLAB inherits the variable, so both halves use the same folder. One process
runs one case: two campaigns started from one checkout would share
inbox/theta.txt and matlab.lock, so a second campaign needs its own clone.

MATLAB runs headless, so its command window does not exist. Everything it would
have printed goes to <results folder>/logs/matlab_console.log. To watch it as it
is written, in a second shell:

    python watch_matlab_log.py --case baseline

The relaunch count is unbounded. Both halves resume from the records they wrote,
so a relaunch continues the run rather than restarting it.

To resume one phase, as with run_pipeline.py:

    python run_supervised.py --case case1 --phase bo

run_pipeline.py still works. Use it when you would rather drive MATLAB yourself,
in the desktop for instance while looking at a candidate that failed.
"""

from __future__ import annotations

import argparse
import signal
import sys

import run_pipeline
from pipeline import driver
from pipeline.matlab_interface import RESULTS_DIR
from pipeline.matlab_supervisor import CONSOLE_LOG, MatlabSupervisor
from run_config import CASES


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", default="case1", choices=sorted(CASES),
                   help="which search space to optimize (default: case1)")
    p.add_argument("--phase", default="both", choices=("both", "init", "bo"),
                   help="run one phase instead of the whole case")
    p.add_argument("--pause", type=float, default=run_pipeline.HANDOVER_PAUSE_S,
                   help="seconds to wait between the two phases")
    p.add_argument("--matlab", default="matlab",
                   help="the MATLAB executable to launch (default: matlab)")
    p.add_argument("--ready-timeout", type=float, default=900.0,
                   help="seconds allowed for MATLAB startup and the parallel pool")
    p.add_argument("--wedge-timeout", type=float, default=600.0,
                   help="seconds a request may be pending with matlab.lock absent "
                        "before MATLAB is killed and relaunched")
    p.add_argument("--diary", action="store_true",
                   help="also write the MATLAB console through diary(), which "
                        "flushes each line, in case the redirected output lags")
    return p.parse_args(argv)


def _raise_on_sigterm() -> None:
    """Turn SIGTERM into an exception so that MATLAB is stopped on the way out.

    Slurm sends SIGTERM to the whole job step at a time limit or an scancel. The
    default action ends this process without running the teardown, which would
    leave MATLAB and its pool to the cgroup cleanup.
    """
    def handler(signum, frame):  # noqa: ARG001
        raise KeyboardInterrupt

    try:
        signal.signal(signal.SIGTERM, handler)
    except (ValueError, OSError, AttributeError):
        pass


def run_case(case: str, args) -> int:
    """Run one case start to finish under a MATLAB of its own."""
    print("#" * 70)
    print(f"CASE {case}")
    print("#" * 70)

    # The optimisation phase has its own entry point, so resuming it does not go
    # through the design server.
    entry = "main_BO" if args.phase == "bo" else "main_initialization"

    supervisor = MatlabSupervisor(
        entry=entry,
        matlab=args.matlab,
        ready_timeout_s=args.ready_timeout,
        wedge_timeout_s=args.wedge_timeout,
        diary=args.diary,
    )

    # start() returns once MATLAB is serving. Nothing may publish a request
    # before then: main_initialization deletes the inbox on startup.
    supervisor.start()
    supervisor.watch()
    # The driver calls this when an evaluation exceeds its timeout: it kills the
    # wedged MATLAB and relaunches it so the run continues. It is set only while a
    # supervisor is attached, so a manual run_pipeline leaves it None and stops
    # with instructions instead.
    driver.TIMEOUT_RESTART_HOOK = supervisor.abandon
    try:
        code = run_pipeline.main(["--case", case, "--phase", args.phase,
                                  "--pause", str(args.pause)])
    finally:
        driver.TIMEOUT_RESTART_HOOK = None
        supervisor.stop()

    if supervisor.restarts:
        print(f"[run_supervised] MATLAB was relaunched {supervisor.restarts} time(s) "
              f"during {case}. See {supervisor.restart_log}.")
    return code


def main(argv=None) -> int:
    args = parse(sys.argv[1:] if argv is None else argv)
    _raise_on_sigterm()

    print(f"[run_supervised] case: {args.case}")
    print(f"[run_supervised] results: {RESULTS_DIR}")
    print(f"[run_supervised] MATLAB console: {CONSOLE_LOG}")
    print(f"[run_supervised] follow it with: python watch_matlab_log.py --case {args.case}")

    try:
        return run_case(args.case, args)
    except KeyboardInterrupt:
        print(f"\n[run_supervised] interrupted during {args.case}. MATLAB is stopped.")
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
