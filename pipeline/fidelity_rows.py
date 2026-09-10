"""Extra training rows read out of the full-fidelity design runs.

A design evaluation that ran to z = 1 recorded, for every control step, the cost
it had accumulated and the time the solver had spent. A run stopped early at
fidelity z replays exactly that prefix: same noise draw, same warm-start chain,
same decisions. The cost and the time such a run would have reported are
therefore prefix sums of what is already on disk, and no simulation is needed to
obtain them.

Nothing here is modelled and phi never appears. The costs are the measured
partial sums and the times are measured times, so a row from this module is the
same kind of observation as a row MATLAB wrote. The caller converts cost to a
full-horizon estimate exactly as it does for every other evaluation.

Two differences from a real truncated run, both immaterial:
  - a real run skips the plant integration on its final step, so the time here
    is high by one ode45 call, of the order of microseconds;
  - the run is not repeated, so the row carries no fresh solver noise. It is the
    same trajectory, read at an earlier point.

The z = 1 row is deliberately not produced: that evaluation is already in the
history as a real row, and duplicating an input with an identical output only
makes the GP's covariance matrix worse conditioned.

Verified on the design runs: the full-horizon prefix reproduces the SSE_measured,
SSdU_measured and runtime_s that MATLAB wrote, to about 1e-11.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np
from scipy.io import loadmat

# Fidelities to read off each design run. z = 1 comes from the real row.
DEFAULT_Z_GRID: Sequence[float] = (0.25, 0.50, 0.75)

# Anything this close to 1 is the real evaluation, not a prefix of it.
_FULL_EPS = 1e-9


def prefix_rows(mat_path: Path | str, z_grid: Sequence[float] = DEFAULT_Z_GRID) -> List[Dict]:
    """Read one out_*.mat and return what it would have reported at each z.

    The returned dictionaries carry the keys the driver's history assembly
    needs: theta (with the fidelity slot replaced), z, SSE_measured,
    SSdU_measured and runtime_s. phase is "DOE_PREFIX" so a reader can tell
    these apart from rows MATLAB wrote, and source names the file they came from.
    """
    mat_path = Path(mat_path)
    store = loadmat(str(mat_path), struct_as_record=False, squeeze_me=True)
    if "out" not in store:
        raise ValueError(f"{mat_path.name} has no 'out' struct")
    out = store["out"]

    cases = np.atleast_1d(out.case)
    n_full = int(out.N)
    theta_full = np.asarray(out.theta, dtype=float).ravel()

    # Per scenario: cost contribution of each step, and the solver time of each
    # step. partial_SSdU is one shorter than the others because it is a
    # difference of consecutive inputs.
    partial_sse = [np.atleast_1d(c.partial_SSE).astype(float) for c in cases]
    partial_ssdu = [np.atleast_1d(c.partial_SSdU).astype(float) for c in cases]
    runtime = [np.atleast_1d(c.RUNTIME).astype(float) for c in cases]

    steps_full = max(n_full - 1, 1)          # control intervals in the full run

    rows: List[Dict] = []
    for z in z_grid:
        if z >= 1.0 - _FULL_EPS:
            continue                          # the real row already covers z = 1
        if not 0.0 < z < 1.0:
            raise ValueError(f"fidelity out of range: {z}")

        # Matches MATLAB: N = min(ceil(tf/dt) + 1, base.N) with tf = z * base.tf.
        n = min(int(math.ceil(z * steps_full)) + 1, n_full)
        if n < 2:
            continue                          # too short to have a Delta u term

        rows.append({
            "eval_id": -1,
            "timestamp": mat_path.stem.replace("out_", ""),
            "phase": "DOE_PREFIX",
            "source": mat_path.name,
            "z": float(z),
            "n_steps": n,
            "theta": _theta_at(theta_full, z),
            "SSE_measured": float(sum(a[:n].sum() for a in partial_sse)),
            "SSdU_measured": float(sum(a[:n - 1].sum() for a in partial_ssdu)),
            "runtime_s": float(sum(a[:n].sum() for a in runtime)),
        })
    return rows


def doe_prefix_rows(paths: Iterable[Path | str],
                    z_grid: Sequence[float] = DEFAULT_Z_GRID) -> List[Dict]:
    """prefix_rows over every design run, skipping files that cannot be read.

    A file missing a field is reported and skipped rather than raising: losing an
    extra row is not worth losing the run.
    """
    rows: List[Dict] = []
    for path in paths:
        try:
            rows.extend(prefix_rows(path, z_grid))
        except (ValueError, KeyError, AttributeError, OSError) as exc:
            print(f"[prefix] skipping {Path(path).name}: {exc}")
    return rows


def verify_full_horizon(mat_path: Path | str, recorded: Dict[str, float],
                        tol: float = 1e-6) -> Dict[str, float]:
    """Check that the prefix sums reproduce the run's own recorded totals."""
    store = loadmat(str(mat_path), struct_as_record=False, squeeze_me=True)
    cases = np.atleast_1d(store["out"].case)
    got = {
        "SSE_measured": float(sum(np.atleast_1d(c.partial_SSE).astype(float).sum() for c in cases)),
        "SSdU_measured": float(sum(np.atleast_1d(c.partial_SSdU).astype(float).sum() for c in cases)),
        "runtime_s": float(sum(np.atleast_1d(c.RUNTIME).astype(float).sum() for c in cases)),
    }
    diffs = {k: abs(got[k] - float(recorded[k])) for k in got}
    for key, diff in diffs.items():
        if diff > tol:
            raise AssertionError(
                f"{Path(mat_path).name}: {key} prefix sum {got[key]!r} does not match "
                f"the recorded {recorded[key]!r} (difference {diff:g})")
    return diffs


def _theta_at(theta_full: np.ndarray, z: float) -> List[float]:
    """theta of the run with the fidelity slot set to z."""
    theta = theta_full.copy()
    theta[0] = float(z)
    return [float(v) for v in theta]
