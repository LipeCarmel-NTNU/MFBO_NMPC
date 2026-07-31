"""Run one case, or several in sequence, with Python owning the MATLAB server.

    python run_supervised.py --case case1
    python run_supervised.py --case case1 case2

This is run_pipeline.py with the MATLAB half started, watched and relaunched by
the driver's own process. You type one command instead of two, and you do not
start MATLAB yourself.

MATLAB runs headless, so its command window does not exist. Everything it would
have printed goes to results/logs/matlab_console.log. To watch it as it is
written, in a second shell:

    python watch_matlab_log.py

Several cases in one command
---------------------------
Both cases write the same paths, and the driver refuses to resume under a
manifest that declares a different case, so a sequence has to move the finished
tree aside between cases. That is what RERUN.md describes doing by hand, and
pipeline/case_archive.py does it here: results/ becomes
results_archive/<case>_<timestamp>/ and an empty results/ is left behind. The
move is a rename, so the earlier run is recoverable by moving it back.

Each case starts a MATLAB of its own at main_initialization. A case that ends in
a non-zero status stops the sequence, so the state that produced it stays where
it is and you can resume that case with the same command. Pass --keep-going to
run the remaining cases anyway.

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
from pipeline import case_archive
from pipeline.matlab_supervisor import CONSOLE_LOG, MatlabSupervisor
from run_config import CASES


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", nargs="+", default=["case1"], choices=sorted(CASES),
                   metavar="CASE",
                   help="one or more search spaces to optimize, in order "
                        f"(choices: {', '.join(sorted(CASES))}; default: case1)")
    p.add_argument("--phase", default="both", choices=("both", "init", "bo"),
                   help="run one phase instead of the whole case")
    p.add_argument("--pause", type=float, default=run_pipeline.HANDOVER_PAUSE_S,
                   help="seconds to wait between the two phases")
    p.add_argument("--keep-going", action="store_true",
                   help="start the next case even when this one failed")
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
    args = p.parse_args(argv)

    if len(args.case) > 1 and args.phase != "both":
        p.error("--phase applies to a single case. Name one case, or drop --phase.")
    seen = set()
    for case in args.case:
        if case in seen:
            p.error(f"case {case!r} is named twice, and the second run would archive "
                    f"the first")
        seen.add(case)
    return args


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

    case_archive.prepare_for(case)

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
    try:
        code = run_pipeline.main(["--case", case, "--phase", args.phase,
                                  "--pause", str(args.pause)])
    finally:
        supervisor.stop()

    if supervisor.restarts:
        print(f"[run_supervised] MATLAB was relaunched {supervisor.restarts} time(s) "
              f"during {case}. See {supervisor.restart_log}.")
    return code


def main(argv=None) -> int:
    args = parse(sys.argv[1:] if argv is None else argv)
    _raise_on_sigterm()

    print(f"[run_supervised] cases: {', '.join(args.case)}")
    print(f"[run_supervised] MATLAB console: {CONSOLE_LOG}")
    print("[run_supervised] follow it with: python watch_matlab_log.py")

    # A sequence restarted after an interruption must not begin again at the
    # first case. results/ names the case it holds, and the cases before that one
    # in this sequence have already run and been archived. Beginning at the first
    # case again would archive a tree that is still being filled and then rerun
    # an earlier case from nothing.
    cases = list(args.case)
    declared = case_archive.declared_case()
    if declared in cases and cases.index(declared) > 0:
        skipped = cases[:cases.index(declared)]
        cases = cases[cases.index(declared):]
        print(f"[run_supervised] results/ holds {declared}, so {', '.join(skipped)} "
              f"ran already and is under {case_archive.ARCHIVE_ROOT.name}/. "
              f"Continuing with {', '.join(cases)}. To run {skipped[0]} again, "
              f"name it on its own.")

    status = 0
    for index, case in enumerate(cases):
        try:
            code = run_case(case, args)
        except KeyboardInterrupt:
            print(f"\n[run_supervised] interrupted during {case}. MATLAB is stopped.")
            return 130
        if code == 0:
            continue

        status = code
        remaining = cases[index + 1:]
        if remaining and not args.keep_going:
            print(f"\n[run_supervised] {case} stopped with code {code}. "
                  f"{', '.join(remaining)} will not start, so the state that "
                  f"produced it stays in results/. Resume with the same command, "
                  f"or pass --keep-going to move on regardless.")
            return status
        if remaining:
            print(f"\n[run_supervised] {case} stopped with code {code}. "
                  f"Continuing to {remaining[0]} because --keep-going was given.")

    if len(cases) > 1 and status == 0:
        print(f"[run_supervised] all cases finished: {', '.join(cases)}. "
              f"Earlier cases are under {case_archive.ARCHIVE_ROOT.name}/, and the "
              f"last one is still in results/.")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
