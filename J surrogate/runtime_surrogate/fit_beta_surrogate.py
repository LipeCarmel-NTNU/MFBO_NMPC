"""Fit the fidelity surrogate phi(z) from NMPC simulation outputs.

The surrogate is the regularized incomplete beta function:

    phi(z) = I_z(a, b),   a > 0,   b > 0.

phi(z) is the share of the full-horizon cost that a run accumulates by fidelity
z. The optimizer divides a measured partial cost by phi(z) to estimate the
full-horizon cost. phi(0) = 0 and phi(1) = 1 hold exactly for every positive a
and b. A full-fidelity evaluation therefore divides by exactly one.

MODEL AND FIT

The design draws z from the Sobol sequence. A training run is truncated, so its
total J_i(1) is unknown. Write log J_i(z) = c_i + log phi(z) with c_i free. The
least-squares c_i is the mean of u_i(z) = log J_i(z) - log phi(z) over the
samples of that run. Substitute it and the cost becomes:

    cost(a, b) = mean_iz ( u_i(z) - mean_z u_i(z) )^2.

This is profile least squares. The fit never estimates c_i. The cost uses a mean
and not a sum, so one L2 strength has the same effect on a small training set
and a large one.

REGULARIZATION

The fit pulls the parameters toward a = b = 1, which is the identity phi(z) = z:

    cost_lambda(a, b) = cost(a, b) + lambda * ( (a - 1)^2 + (b - 1)^2 ).

Cross-validation over runs selects lambda. It never splits over samples, because
the profiled constant c_i makes one run the smallest unit that the fit can hold
out. The score is the unregularized cost on the held-out runs.

The selected lambda is the plain argmin of the mean held-out loss. The fit
applies no one-standard-error rule. Each fit records the full loss grid, so you
can apply another rule later without a refit.

The fit uses one fold partition. It does not redraw the partition. The
optimization uses one surrogate at a time, and the record keeps the loss grid
that produced it.

USE

The optimization driver imports this module and refits during a run. You can
also run it as a script to fit one directory of outputs:

    python fit_beta_surrogate.py results/init
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import scipy.io
from scipy.optimize import minimize
from scipy.special import betainc

# The project root holds run_config, which declares both fidelity limits.
# Z_MIN_PHI is the floor of this fit. Z_MIN_BO is the lower bound of the search
# space and is a separate limit on a separate step. Import the value so that
# this module and the run configuration cannot disagree.
_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from run_config import Z_MIN_PHI  # noqa: E402


# ---------------------------------------------------------------------------
# Fit settings. Each fit copies these values into its record. You can therefore
# reproduce a fit from its own record without this file.
# ---------------------------------------------------------------------------

HORIZON_HOURS = 10.0        # T in z = t / T
LAMBDA_GRID = (0.0, 1e-2, 1e0, 1e2)   # [0, logspace(-2, 2, 3)]
K_FOLD = 5                  # folds of the cross-validation, split over runs
RNG_SEED = 1                # fixes the fold partition
AB0 = (1.0, 1.0)            # start point, phi(z) = z
AB_LB = (0.05, 0.05)
AB_UB = (25.0, 25.0)
MIN_SAMPLES_PER_RUN = 2     # a run needs two samples after the fit profiles out c_i
AB_CLIP = (1e-6, 1e4)       # domain guard for betainc at trial points
PHI_FLOOR = 1e-12           # keeps the log defined where phi underflows

# Optimizer settings for the (a, b) fit. Each fit record keeps them.
OPT_METHOD = "L-BFGS-B"
OPT_MAX_ITER = 400
OPT_FTOL = 1e-12
OPT_GTOL = 1e-10
OPT_FD_STEP = 1e-6          # central-difference step, same as the MATLAB fit


@dataclass(frozen=True)
class TargetSpec:
    """Locates one cost target in an out_*.mat file.

    partial_field is the per-step cost of one case. first_time_index is the
    index in out.T of the first entry of partial_field. The input-change cost
    needs two samples, so its first entry belongs to out.T[1] and its index
    is 1.
    """

    partial_field: str
    first_time_index: int


TARGETS: Dict[str, TargetSpec] = {
    "SSE": TargetSpec("partial_SSE", 0),
    "SSdU": TargetSpec("partial_SSdU", 1),
}


@dataclass
class RunCurve:
    """One training run and its cumulative cost against fidelity."""

    source: str
    z: np.ndarray
    log_j: np.ndarray

    @property
    def n_samples(self) -> int:
        return int(self.z.size)


@dataclass
class FitResult:
    """The result of one fit, and the data needed to audit it."""

    target: str
    a: float
    b: float
    lam: float
    loss: float
    lambda_grid: List[float]
    cv_loss: List[float]
    cv_folds_scored: int
    n_runs_total: int
    n_runs_used: int
    n_samples: int
    sources: List[str] = field(default_factory=list)
    dropped_sources: List[str] = field(default_factory=list)
    optimizer: Dict = field(default_factory=dict)
    settings: Dict = field(default_factory=dict)
    wall_s: Dict = field(default_factory=dict)

    def to_dict(self) -> Dict:
        d = asdict(self)
        d["a"] = float(self.a)
        d["b"] = float(self.b)
        d["lam"] = float(self.lam)
        d["loss"] = float(self.loss)
        return d


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

def phi(z: np.ndarray, ab: Sequence[float]) -> np.ndarray:
    """Return phi(z) = I_z(a, b), the regularized incomplete beta function.

    phi increases with z. phi(0) = 0 and phi(1) = 1 for every positive a and b.

    This function clips its arguments to the domain that betainc accepts. The
    bounds hold at the solution but not at every trial point that a line search
    proposes. The clip limits on (a, b) are much wider than the fit bounds, so
    they act only on infeasible trial points.

    Note the argument order. scipy.special.betainc takes (a, b, x). MATLAB
    betainc takes (x, a, b). The two agree after you match the order.
    """
    a = float(np.clip(ab[0], AB_CLIP[0], AB_CLIP[1]))
    b = float(np.clip(ab[1], AB_CLIP[0], AB_CLIP[1]))
    x = np.clip(np.asarray(z, dtype=float), 0.0, 1.0)
    return betainc(a, b, x)


def _profiled_residuals(ab: Sequence[float], z: np.ndarray, log_j: np.ndarray,
                        group: np.ndarray, n_groups: int) -> np.ndarray:
    """Return u - mean_z u for each run. This is the residual after profiling."""
    u = log_j - np.log(np.maximum(phi(z, ab), PHI_FLOOR))
    counts = np.bincount(group, minlength=n_groups).astype(float)
    sums = np.bincount(group, weights=u, minlength=n_groups)
    u_bar = sums / np.maximum(counts, 1.0)
    return u - u_bar[group]


def phi_loss(ab: Sequence[float], z: np.ndarray, log_j: np.ndarray,
             group: np.ndarray, n_groups: int) -> float:
    """Return the mean squared deviation of u from its per-run mean.

    The centering on the per-run mean profiles out the unknown log J_i(1). The
    loss therefore holds no run total. The mean over samples makes the value
    comparable between training sets of different size.
    """
    r = _profiled_residuals(ab, z, log_j, group, n_groups)
    return float(np.mean(r ** 2))


def phi_cost(ab: Sequence[float], z: np.ndarray, log_j: np.ndarray,
             group: np.ndarray, n_groups: int, lam: float) -> float:
    """Return the data loss plus an L2 pull of (a, b) toward (1, 1).

    The penalty center is the identity phi(z) = z. Shrinkage therefore moves the
    surrogate toward a cost that accumulates at a constant rate.
    """
    ab = np.asarray(ab, dtype=float)
    return phi_loss(ab, z, log_j, group, n_groups) + lam * float(np.sum((ab - 1.0) ** 2))


def _central_gradient(fun, ab: np.ndarray, step: float = OPT_FD_STEP) -> np.ndarray:
    """Return the central-difference gradient in two dimensions.

    The profiled loss is smooth, but betainc makes its forward difference noisy
    at the step size that L-BFGS-B selects by default. The model has two
    parameters, so the central difference is cheap.
    """
    g = np.empty(2, dtype=float)
    for i in range(2):
        h = step * max(1.0, abs(float(ab[i])))
        up = np.array(ab, dtype=float)
        dn = np.array(ab, dtype=float)
        up[i] += h
        dn[i] -= h
        g[i] = (fun(up) - fun(dn)) / (2.0 * h)
    return g


def fit_phi(z: np.ndarray, log_j: np.ndarray, group: np.ndarray, n_groups: int,
            lam: float, ab0: Sequence[float] = AB0) -> Tuple[np.ndarray, Dict]:
    """Minimize phi_cost over the two shape parameters under box bounds."""
    def obj(p):
        return phi_cost(p, z, log_j, group, n_groups, lam)

    t0 = time.perf_counter()
    res = minimize(
        obj,
        np.asarray(ab0, dtype=float),
        method=OPT_METHOD,
        jac=lambda p: _central_gradient(obj, p),
        bounds=[(AB_LB[0], AB_UB[0]), (AB_LB[1], AB_UB[1])],
        options={"maxiter": OPT_MAX_ITER, "ftol": OPT_FTOL, "gtol": OPT_GTOL},
    )

    info = {
        "method": OPT_METHOD,
        "success": bool(res.success),
        "status": int(getattr(res, "status", -1)),
        "message": str(res.message),
        "nit": int(getattr(res, "nit", -1)),
        "nfev": int(getattr(res, "nfev", -1)),
        "fun": float(res.fun),
        "x0": [float(v) for v in ab0],
        "bounds": [list(AB_LB), list(AB_UB)],
        "maxiter": OPT_MAX_ITER,
        "ftol": OPT_FTOL,
        "gtol": OPT_GTOL,
        "fd_step": OPT_FD_STEP,
        "fd_type": "central",
        "wall_s": time.perf_counter() - t0,
    }
    return np.asarray(res.x, dtype=float), info


def select_lambda(z: np.ndarray, log_j: np.ndarray, group: np.ndarray, n_groups: int,
                  lambda_grid: Sequence[float] = LAMBDA_GRID,
                  k_fold: int = K_FOLD,
                  seed: int = RNG_SEED) -> Tuple[float, np.ndarray, int]:
    """Select the L2 strength by k-fold cross-validation over runs.

    The function splits the runs into k folds. It scores each fold in turn with
    parameters fitted on the other folds. The score is the unregularized
    phi_loss. It uses one partition, which the seed fixes.

    The function walks the grid from weak to strong shrinkage. It starts each
    fit from the previous solution, which keeps the optimizer iteration count
    small along the path.

    It returns the plain argmin of the mean held-out loss, the full loss grid,
    and the number of folds it scored. A tie takes the first grid entry, which
    is the weakest shrinkage among the tied values.
    """
    lambda_grid = list(lambda_grid)
    n_lam = len(lambda_grid)
    k = min(k_fold, n_groups)
    cv_loss = np.zeros(n_lam, dtype=float)
    n_scored = 0

    if k < 2:
        # One run cannot be split, so no lambda can be scored. Return the
        # weakest shrinkage and report zero folds. The caller records this.
        return float(lambda_grid[0]), np.full(n_lam, np.nan), 0

    rng = np.random.default_rng(seed)
    fold = np.mod(rng.permutation(n_groups), k)

    for f in range(k):
        held_groups = np.flatnonzero(fold == f)
        train_groups = np.flatnonzero(fold != f)
        if held_groups.size == 0 or train_groups.size == 0:
            continue

        tr = np.isin(group, train_groups)
        he = np.isin(group, held_groups)
        z_tr, lj_tr, g_tr, n_tr = _compact(z[tr], log_j[tr], group[tr])
        z_he, lj_he, g_he, n_he = _compact(z[he], log_j[he], group[he])

        ab = np.asarray(AB0, dtype=float)
        for i_lam, lam in enumerate(lambda_grid):
            ab, _ = fit_phi(z_tr, lj_tr, g_tr, n_tr, lam, ab0=ab)
            cv_loss[i_lam] += phi_loss(ab, z_he, lj_he, g_he, n_he)
        n_scored += 1

    if n_scored == 0:
        return float(lambda_grid[0]), np.full(n_lam, np.nan), 0

    cv_loss /= n_scored
    lam_best = float(lambda_grid[int(np.argmin(cv_loss))])
    return lam_best, cv_loss, n_scored


def _compact(z: np.ndarray, log_j: np.ndarray, group: np.ndarray):
    """Rebuild a contiguous group index over a subset of the samples."""
    _, g = np.unique(group, return_inverse=True)
    return z, log_j, g, int(g.max()) + 1 if g.size else 0


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def load_curve(path: Path, target: str, horizon_hours: float = HORIZON_HOURS) -> RunCurve:
    """Read one out_*.mat file into a cumulative cost curve.

    The function adds the two initial-condition cases together. This matches the
    definition of the truncated objective. first_time_index is the index in
    out.T of the first entry of the partial array.
    """
    spec = TARGETS[target]
    try:
        mat = scipy.io.loadmat(str(path), struct_as_record=False, squeeze_me=True)
    except NotImplementedError as exc:
        # scipy reads MATLAB v7 and below. A v7.3 file uses HDF5 and raises
        # here. Without this message the failure appears in the middle of a run
        # with no cause.
        raise RuntimeError(
            f"{Path(path).name} is a MATLAB v7.3 file, which this fit cannot read. "
            f"The simulation drivers must call save() without the '-v7.3' flag."
        ) from exc

    if "out" not in mat:
        raise ValueError(f"{path} has no 'out' struct")
    out = mat["out"]

    time_h = np.asarray(out.T, dtype=float).ravel()
    cases = out.case
    if np.ndim(cases) == 0:
        cases = [cases]
    if len(cases) < 2:
        raise ValueError(f"{path} holds {len(cases)} case(s), expected 2")

    p1 = np.asarray(getattr(cases[0], spec.partial_field), dtype=float).ravel()
    p2 = np.asarray(getattr(cases[1], spec.partial_field), dtype=float).ravel()

    n = int(min(p1.size, p2.size, time_h.size - spec.first_time_index))
    if n <= 0:
        raise ValueError(f"{path} holds no usable samples for {target}")

    z = time_h[spec.first_time_index: spec.first_time_index + n] / horizon_hours
    j = np.cumsum(p1[:n]) + np.cumsum(p2[:n])
    with np.errstate(divide="ignore"):
        log_j = np.log(j)
    return RunCurve(source=str(path), z=z, log_j=log_j)


def pool_runs(curves: Sequence[RunCurve], z_min: float = Z_MIN_PHI):
    """Return flat sample arrays over the usable runs, tagged by run.

    One flat vector lets the loss call betainc once for each evaluation instead
    of once for each run. A run needs two samples after the fit profiles out its
    constant. The function drops shorter runs and reports them.
    """
    z_parts: List[np.ndarray] = []
    lj_parts: List[np.ndarray] = []
    g_parts: List[np.ndarray] = []
    used: List[str] = []
    dropped: List[str] = []

    g = 0
    for c in curves:
        keep = (c.z >= z_min) & np.isfinite(c.log_j)
        n_keep = int(np.count_nonzero(keep))
        if n_keep < MIN_SAMPLES_PER_RUN:
            dropped.append(c.source)
            continue
        z_parts.append(c.z[keep])
        lj_parts.append(c.log_j[keep])
        g_parts.append(np.full(n_keep, g, dtype=int))
        used.append(c.source)
        g += 1

    if g == 0:
        return (np.empty(0), np.empty(0), np.empty(0, dtype=int), 0, used, dropped)

    return (np.concatenate(z_parts), np.concatenate(lj_parts),
            np.concatenate(g_parts), g, used, dropped)


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def fit_target(paths: Sequence[Path], target: str, *, z_min: float = Z_MIN_PHI,
               lambda_grid: Sequence[float] = LAMBDA_GRID, k_fold: int = K_FOLD,
               seed: int = RNG_SEED,
               horizon_hours: float = HORIZON_HOURS) -> FitResult:
    """Fit phi(z) for one target from an explicit list of output files.

    The caller passes the file list and not a directory. A fit made during a run
    therefore records exactly which evaluations informed it.
    """
    t_start = time.perf_counter()

    t0 = time.perf_counter()
    curves = [load_curve(Path(p), target, horizon_hours) for p in paths]
    z, log_j, group, n_groups, used, dropped = pool_runs(curves, z_min)
    t_load = time.perf_counter() - t0

    if n_groups == 0:
        raise ValueError(
            f"no usable run for target {target}: all {len(curves)} candidates hold "
            f"fewer than {MIN_SAMPLES_PER_RUN} samples at z >= {z_min}"
        )

    t0 = time.perf_counter()
    lam, cv_loss, n_scored = select_lambda(z, log_j, group, n_groups,
                                           lambda_grid, k_fold, seed)
    t_cv = time.perf_counter() - t0

    ab, opt_info = fit_phi(z, log_j, group, n_groups, lam)
    loss = phi_loss(ab, z, log_j, group, n_groups)

    return FitResult(
        target=target,
        a=float(ab[0]),
        b=float(ab[1]),
        lam=float(lam),
        loss=float(loss),
        lambda_grid=[float(v) for v in lambda_grid],
        cv_loss=[float(v) for v in cv_loss],
        cv_folds_scored=int(n_scored),
        n_runs_total=len(curves),
        n_runs_used=int(n_groups),
        n_samples=int(z.size),
        sources=[Path(s).name for s in used],
        dropped_sources=[Path(s).name for s in dropped],
        optimizer=opt_info,
        settings={
            "model": "phi(z) = I_z(a, b)",
            "horizon_hours": float(horizon_hours),
            "z_min": float(z_min),
            "k_fold": int(k_fold),
            "cv_partitions": 1,
            "cv_seed": int(seed),
            "lambda_selection": "argmin of the mean held-out loss, no one-standard-error rule",
            "min_samples_per_run": int(MIN_SAMPLES_PER_RUN),
            "phi_floor": PHI_FLOOR,
            "ab_clip": list(AB_CLIP),
        },
        wall_s={
            "load_and_pool": t_load,
            "cross_validation": t_cv,
            "final_fit": opt_info["wall_s"],
            "total": time.perf_counter() - t_start,
        },
    )


def fit_all_targets(paths: Sequence[Path], **kwargs) -> Dict[str, FitResult]:
    """Fit every target on the same file list."""
    return {name: fit_target(paths, name, **kwargs) for name in TARGETS}


def write_coefficients(results: Dict[str, FitResult], mat_path: Path, vintage: int,
                       created_at: str, extra: Dict | None = None) -> float:
    """Write the coefficient file that MATLAB loads. Return the wall time.

    The write is atomic. The payload goes to a temporary name in the same
    directory. A rename then puts it in place. A reader therefore sees the
    previous vintage or the new one, and never a partial file.
    """
    t0 = time.perf_counter()
    mat_path = Path(mat_path)
    mat_path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "vintage": np.array([[float(vintage)]]),
        "created_at": str(created_at),
        "model": "phi(z) = I_z(a, b), regularized incomplete beta",
    }
    for name, r in results.items():
        payload[f"a_{name}"] = np.array([[float(r.a)]])
        payload[f"b_{name}"] = np.array([[float(r.b)]])
        payload[f"lambda_{name}"] = np.array([[float(r.lam)]])
        payload[f"n_runs_{name}"] = np.array([[float(r.n_runs_used)]])
    if extra:
        payload.update(extra)

    tmp = mat_path.with_suffix(mat_path.suffix + ".tmp")
    scipy.io.savemat(str(tmp), payload, do_compression=False)

    last: Exception | None = None
    for _ in range(10):
        try:
            tmp.replace(mat_path)
            return time.perf_counter() - t0
        except OSError as exc:
            last = exc
            time.sleep(0.15)
    raise RuntimeError(f"could not publish the coefficients to {mat_path}") from last


def write_vintage_record(results: Dict[str, FitResult], json_path: Path, vintage: int,
                         created_at: str, extra: Dict | None = None) -> float:
    """Write the full record of one vintage. Return the wall time."""
    t0 = time.perf_counter()
    json_path = Path(json_path)
    json_path.parent.mkdir(parents=True, exist_ok=True)

    record = {
        "vintage": int(vintage),
        "created_at": str(created_at),
        "model": "phi(z) = I_z(a, b), regularized incomplete beta",
        "targets": {name: r.to_dict() for name, r in results.items()},
    }
    if extra:
        record["context"] = extra

    tmp = json_path.with_suffix(json_path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2, sort_keys=True)
        fh.flush()
    tmp.replace(json_path)
    return time.perf_counter() - t0


def _main(argv: Sequence[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        print("usage: fit_beta_surrogate.py <directory of out_*.mat> [more directories]")
        return 2

    paths: List[Path] = []
    for d in argv[1:]:
        paths.extend(sorted(Path(d).glob("out_*.mat")))
    if not paths:
        print(f"no out_*.mat found under {argv[1:]}")
        return 1

    print(f"fitting phi(z) on {len(paths)} files, z_min = {Z_MIN_PHI}")
    results = fit_all_targets(paths)
    for name, r in results.items():
        print(f"\n{name}: a = {r.a:.6f}, b = {r.b:.6f}, lambda = {r.lam:g}, "
              f"loss = {r.loss:.6e}")
        print(f"  runs used {r.n_runs_used}/{r.n_runs_total}, samples {r.n_samples}, "
              f"folds scored {r.cv_folds_scored}")
        for lam, cv in zip(r.lambda_grid, r.cv_loss):
            print(f"    lambda {lam:<8g} cv loss {cv:.6e}")
        print(f"  optimizer: {r.optimizer['message']} "
              f"(nit={r.optimizer['nit']}, nfev={r.optimizer['nfev']})")
        print(f"  wall: load {r.wall_s['load_and_pool']:.2f} s, "
              f"cv {r.wall_s['cross_validation']:.2f} s, "
              f"fit {r.wall_s['final_fit']:.2f} s, "
              f"total {r.wall_s['total']:.2f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv))
