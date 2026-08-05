"""File-based exchange between the optimisation driver and MATLAB.

The driver writes one evaluation request to inbox/theta.txt and waits for the
matching row to appear in the results CSV. MATLAB serves one request at a time
and holds matlab.lock while it is busy.

Each request carries an eval_id. Pairing a result with its request by that index
rather than by comparing theta values matters for two reasons: the acquisition
maximiser can propose the same point twice, so theta does not identify a row;
and a float that has been through text is not reliably equal to the one that was
sent, so equality comparison can fail on a row that did arrive.
"""

from __future__ import annotations

import csv
import os
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


BASE_DIR = Path(__file__).resolve().parents[1]
LOCK_FILE = BASE_DIR / "matlab.lock"
THETA_FILE = BASE_DIR / "inbox" / "theta.txt"

# The results tree may live in a subfolder instead of results/, so several cases
# can sit side by side without moving files between runs. MFBO_RESULTS_DIR sets
# it (relative to the project root, e.g. "results/running"); it defaults to
# "results". The MATLAB entry scripts read the same variable with getenv. The
# exchange (inbox/theta.txt, matlab.lock) stays at the project root regardless.
RESULTS_DIR = BASE_DIR / os.environ.get("MFBO_RESULTS_DIR", "results")

# The initialisation phase and the optimisation phase write separate files, so
# that the costs measured at the simulated fidelity stay distinguishable from
# the surrogate-scaled costs of the optimisation runs.
INIT_RESULTS_FILE = RESULTS_DIR / "init" / "results.csv"
BO_RESULTS_FILE = RESULTS_DIR / "results.csv"
# One failures file for both phases. The MATLAB serve loop writes it, and that
# loop does not know which phase a request belongs to.
FAILURES_FILE = RESULTS_DIR / "failures.csv"
INIT_OUT_DIR = RESULTS_DIR / "init"
BO_OUT_DIR = RESULTS_DIR

PHI_COEFFS_FILE = RESULTS_DIR / "surrogate" / "phi_coeffs.mat"

THETA_LEN = 12

# Phase codes sent with every request. The MATLAB server reads the code and
# routes the request. It therefore needs no state of its own to know which
# phase a run is in.
PHASE_DOE = 0    # report the cost measured at the simulated fidelity
PHASE_OPT = 1    # scale the measured cost by phi(z)
PHASE_CODES = {"init": PHASE_DOE, "bo": PHASE_OPT}

# Columns written by dependencies/io/results_csv_header.m. The reader checks
# them, so a schema change on either side is caught at the first read instead of
# silently shifting the meaning of a column.
RESULTS_COLUMNS = [
    "eval_id", "timestamp", "phase", "phi_vintage", "z",
    "SSE_measured", "SSdU_measured", "phi_SSE", "phi_SSdU",
    "SSE", "SSdU", "J", "runtime_s", "n_flag_not_one", "phi_floored",
    "wall_total_s", "wall_cases_s", "wall_phi_s", "wall_build_s", "wall_save_s",
] + [f"theta_{i}" for i in range(1, THETA_LEN + 1)]

_FLOAT_COLUMNS = {
    "phi_vintage", "z", "SSE_measured", "SSdU_measured", "phi_SSE", "phi_SSdU",
    "SSE", "SSdU", "J", "runtime_s", "n_flag_not_one", "phi_floored",
    "wall_total_s", "wall_cases_s", "wall_phi_s", "wall_build_s", "wall_save_s",
}


def results_file(phase: str) -> Path:
    if phase == "init":
        return INIT_RESULTS_FILE
    if phase == "bo":
        return BO_RESULTS_FILE
    raise ValueError(f"unknown phase {phase!r}, choose 'init' or 'bo'")


def failures_file(phase: str) -> Path:
    """Both phases share one failures file. phase is accepted for symmetry."""
    if phase not in ("init", "bo"):
        raise ValueError(f"unknown phase {phase!r}, choose 'init' or 'bo'")
    return FAILURES_FILE


def out_dir(phase: str) -> Path:
    if phase == "init":
        return INIT_OUT_DIR
    if phase == "bo":
        return BO_OUT_DIR
    raise ValueError(f"unknown phase {phase!r}, choose 'init' or 'bo'")


# ---------------------------------------------------------------------------
# Requests
# ---------------------------------------------------------------------------

def _validate_theta(theta: Sequence[float] | Iterable[float]) -> List[float]:
    """Check a proposal before it leaves the driver.

    The bounds are enforced here as well as in the driver because this is the
    last point at which a malformed vector can be caught with a message naming
    the component. Past it, MATLAB sees only numbers.
    """
    values = [float(v) for v in theta]
    if len(values) != THETA_LEN:
        raise ValueError(f"theta must contain {THETA_LEN} elements, received {len(values)}")

    z = values[0]
    if not (0.0 <= z <= 1.0):
        raise ValueError(f"theta[0]: z must lie in [0, 1], got {z}")

    for i, name in ((1, "theta_p"), (2, "theta_m")):
        v = values[i]
        if v != round(v) or v < 0:
            raise ValueError(f"theta[{i}]: {name} must be a non-negative integer, got {v}")
        values[i] = float(round(v))

    for i in range(3, 6):
        if not (-3.0 <= values[i] <= 3.0):
            raise ValueError(f"theta[{i}]: q[{i - 3}] must lie in [-3, 3], got {values[i]}")

    # The input-magnitude exponents are allowed outside [-3, 3]: case2 sets them
    # to a large negative value so that 10^r underflows to exactly zero, which
    # is how the penalty is removed from the NMPC objective.
    for i in range(9, 12):
        if not (-3.0 <= values[i] <= 3.0):
            raise ValueError(
                f"theta[{i}]: r_du[{i - 9}] must lie in [-3, 3], got {values[i]}")

    for i, v in enumerate(values):
        if v != v or v in (float("inf"), float("-inf")):
            raise ValueError(f"theta[{i}] is not finite: {v}")

    return values


def send_request(eval_id: int, phase: str, theta: Sequence[float], *,
                 lock_stale_s: float = 6 * 3600.0) -> None:
    """Publish one evaluation request.

    The line is

        eval_id phase_code theta_1 ... theta_12

    The phase code tells the MATLAB server whether to scale the measured cost
    by phi. The first request carrying the optimization code is also what makes
    main_initialization hand over to main_BO.

    The function writes a temporary file and renames it. MATLAB polling the
    file therefore sees the previous request or this one. The reader also
    checks the value count, so a partial line is rejected rather than parsed as
    a short theta.
    """
    if phase not in PHASE_CODES:
        raise ValueError(f"unknown phase {phase!r}, choose 'init' or 'bo'")
    values = _validate_theta(theta)
    wait_for_lock(lock_stale_s)

    THETA_FILE.parent.mkdir(parents=True, exist_ok=True)
    payload = " ".join([str(int(eval_id)), str(PHASE_CODES[phase])]
                       + [repr(v) for v in values]) + "\n"

    tmp = THETA_FILE.with_suffix(".tmp")
    with tmp.open("w", encoding="ascii") as fh:
        fh.write(payload)
        fh.flush()
        os.fsync(fh.fileno())

    last: Optional[Exception] = None
    for _ in range(10):
        try:
            tmp.replace(THETA_FILE)
            return
        except OSError as exc:
            last = exc
            time.sleep(0.15)
    raise RuntimeError(f"could not publish the request to {THETA_FILE}") from last


def wait_for_lock(stale_s: float = 6 * 3600.0, poll_s: float = 1.0) -> None:
    """Block while MATLAB is busy, ignoring a lock left behind by a crash.

    Without the staleness check, a MATLAB process killed mid-evaluation leaves
    the driver waiting forever on a file nobody will remove.
    """
    warned = False
    while LOCK_FILE.exists():
        try:
            age = time.time() - LOCK_FILE.stat().st_mtime
        except OSError:
            return
        if age > stale_s:
            print(f"[warn] lock file is {age / 3600:.1f} h old and is treated as "
                  f"abandoned: {LOCK_FILE}")
            try:
                LOCK_FILE.unlink()
            except OSError:
                pass
            return
        if not warned:
            print(f"[wait] MATLAB is busy ({LOCK_FILE.name} present) ...")
            warned = True
        time.sleep(poll_s)


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

def read_results(path: Path | str) -> List[Dict]:
    """Read a results CSV into a list of row dictionaries.

    A row still being written is dropped rather than parsed. MATLAB appends each
    row with one buffered write, so the window is small, but a short final line
    is still possible and must not raise: the caller polls and will see the
    complete row on the next pass.
    """
    results_path = Path(path)
    if not results_path.exists():
        return []

    with results_path.open("r", encoding="ascii", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            return []

        missing = [c for c in RESULTS_COLUMNS if c not in reader.fieldnames]
        if missing:
            raise ValueError(
                f"{results_path} is missing columns {missing}. It was written under a "
                f"different schema; move it aside or point the run at a new results "
                f"directory."
            )

        rows: List[Dict] = []
        for raw in reader:
            if any(raw.get(c) is None for c in RESULTS_COLUMNS):
                continue                      # truncated trailing row
            try:
                rows.append(_parse_row(raw))
            except (TypeError, ValueError):
                continue
    return rows


def _parse_row(raw: Dict[str, str]) -> Dict:
    row: Dict = {
        "eval_id": int(float(raw["eval_id"])),
        "timestamp": raw["timestamp"],
        "phase": raw["phase"],
    }
    for col in _FLOAT_COLUMNS:
        row[col] = float(raw[col])
    row["theta"] = [float(raw[f"theta_{i}"]) for i in range(1, THETA_LEN + 1)]
    return row


def read_failures(path: Path | str) -> List[Dict]:
    """Read the failure log MATLAB writes when an evaluation raises."""
    failures_path = Path(path)
    if not failures_path.exists():
        return []

    out: List[Dict] = []
    with failures_path.open("r", encoding="ascii", newline="") as fh:
        for raw in csv.DictReader(fh):
            if raw.get("eval_id") is None:
                continue
            try:
                out.append({
                    "eval_id": int(float(raw["eval_id"])),
                    "timestamp": raw.get("timestamp", ""),
                    "identifier": raw.get("identifier", ""),
                    "message": raw.get("message", ""),
                })
            except (TypeError, ValueError):
                continue
    return out


def wait_for_result(eval_id: int, results_path: Path, failures_path: Path, *,
                    poll_s: float = 1.0, timeout_s: float = 6 * 3600.0,
                    heartbeat_s: float = 600.0, idle_warn_s: float = 120.0) -> Dict:
    """Block until the given evaluation reports success or failure.

    Returns the results row on success. Raises EvaluationFailed when MATLAB
    recorded a failure for this eval_id, and TimeoutError when neither appears,
    which is the signature of a MATLAB process that died without writing.
    """
    t0 = time.time()
    next_heartbeat = t0 + heartbeat_s
    idle_since = None

    while True:
        for row in read_results(results_path):
            if row["eval_id"] == eval_id:
                return row

        for failure in read_failures(failures_path):
            if failure["eval_id"] == eval_id:
                raise EvaluationFailed(eval_id, failure)

        now = time.time()
        if now - t0 > timeout_s:
            raise TimeoutError(
                f"evaluation {eval_id} produced neither a results row nor a failure "
                f"row within {timeout_s / 3600:.1f} h. Check that MATLAB is still "
                f"running main_initialization.m or main_BO.m."
            )

        # MATLAB holds the lock for as long as it works on a request. A missing
        # lock while this function waits therefore means that MATLAB is not
        # working on the request. That is the signature of a MATLAB process that
        # died without writing a failure row, which otherwise looks exactly like
        # a slow evaluation and costs the whole timeout to discover.
        busy = LOCK_FILE.exists()
        if busy:
            idle_since = None
        elif idle_since is None:
            idle_since = now

        if now >= next_heartbeat:
            mins = (now - t0) / 60.0
            if idle_since is not None and now - idle_since > idle_warn_s:
                print(f"[warn] evaluation {eval_id} has been pending {mins:.0f} min and "
                      f"{LOCK_FILE.name} has been absent for "
                      f"{(now - idle_since) / 60:.0f} min. MATLAB holds that lock while "
                      f"it works, so it is probably no longer running. Check the MATLAB "
                      f"console, then restart the server and this command to resume.",
                      flush=True)
            else:
                print(f"[wait] evaluation {eval_id} running for {mins:.0f} min",
                      flush=True)
            next_heartbeat = now + heartbeat_s

        time.sleep(poll_s)


class EvaluationFailed(RuntimeError):
    """MATLAB recorded a failure for one evaluation."""

    def __init__(self, eval_id: int, failure: Dict):
        self.eval_id = eval_id
        self.failure = failure
        super().__init__(
            f"evaluation {eval_id} failed in MATLAB: "
            f"{failure.get('identifier', '')} {failure.get('message', '')}"
        )


def wait_for_matlab_ready(results_path: Path, max_wait_s: float = 120.0,
                          poll_s: float = 2.0) -> None:
    """Wait for MATLAB to create the results file at startup."""
    t0 = time.time()
    while not Path(results_path).exists():
        if time.time() - t0 > max_wait_s:
            print(f"[warn] {results_path} still absent after {max_wait_s:.0f} s. "
                  f"Continuing; the first request will wait for it.")
            return
        print(f"[wait] waiting for MATLAB to create {results_path} ...")
        time.sleep(poll_s)
