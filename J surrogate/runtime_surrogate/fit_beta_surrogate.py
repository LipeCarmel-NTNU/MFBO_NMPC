"""Fit the fidelity surrogate f(z) = I_z(a, b) from NMPC simulation outputs.

Model
    J(z) = J(1) * f(z), with f the regularised incomplete beta function
        f(z) = I_z(a, b),   a > 0,   b > 0,
    fitted by profile least squares on log J. The design draws z from the Sobol
    sequence, so a training run is truncated and its total J_i(1) is unknown.
    Writing log J_i(z) = c_i + log f(z) with c_i free, the least-squares c_i is
    the mean of u_i(z) = log J_i(z) - log f(z) over that run's samples, and
    substituting it leaves

        cost(a, b) = mean_iz ( u_i(z) - mean_z u_i(z) )^2.

    The mean rather than the sum keeps the data term independent of the number
    of samples, so one L2 strength is comparable across training sizes and
    across the successive refits of a run.

Regularisation
    The fitted parameters are pulled towards a = b = 1, the identity f(z) = z:

        cost_lambda(a, b) = cost(a, b) + lambda * ( (a - 1)^2 + (b - 1)^2 ).

    lambda is selected by k-fold cross-validation over runs, never over samples,
    because the profiled per-run constant makes a run the indivisible unit. The
    partition is redrawn N_CV_REPEATS times and the held-out score is the
    unregularised cost on the held-out runs, averaged over folds and redraws.
    The selected lambda is the plain argmin of that average. No
    one-standard-error rule is applied: the returned lambda is the grid value
    with the lowest mean held-out loss, and the full loss grid is recorded so
    that any other selection rule can be applied after the fact.

Endpoints
    I_z(a, b) satisfies f(0) = 0 and f(1) = 1 exactly for every admissible
    a and b, so no endpoint constraint is imposed and no residual is left at
    z = 1. A full-fidelity evaluation is therefore divided by exactly one.

This module is imported by the optimisation driver, which refits during the run,
and it can be run as a script to fit one directory of outputs offline.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import scipy.io
from scipy.optimize import minimize
from scipy.special import betainc

# ---------------------------------------------------------------------------
# Fit settings. Every value here is copied into the vintage record, so a fit can
# be reproduced from its own record without reading this file.
# ---------------------------------------------------------------------------

HORIZON_HOURS = 10.0        # T in z = t / T
Z_MIN = 0.1                 # samples below this fidelity are dropped
LAMBDA_GRID = (0.0, 1e-2, 1e0, 1e2)   # [0, logspace(-2, 2, 3)]
K_FOLD = 5                  # folds of the cross-validation, split over runs
N_CV_REPEATS = 5            # fold partitions redrawn per fit
RNG_SEED = 1                # fixes the fold partitions
AB0 = (1.0, 1.0)            # warm start, f(z) = z
AB_LB = (0.05, 0.05)
AB_UB = (25.0, 25.0)
MIN_SAMPLES_PER_RUN = 2     # a run needs two samples once its constant is profiled out
AB_CLIP = (1e-6, 1e4)       # domain guard for betainc at trial points
F_FLOOR = 1e-12             # keeps log defined where f underflows

# Optimiser settings for the (a, b) fit, recorded alongside the result.
OPT_METHOD = "L-BFGS-B"
OPT_MAX_ITER = 400
OPT_FTOL = 1e-12
OPT_GTOL = 1e-10
OPT_FD_STEP = 1e-6          # central-difference step, matching the MATLAB fit


@dataclass(frozen=True)
class TargetSpec:
    """Locates one cost target inside an out_*.mat file.

    partial_field is the per-step cost of one case. first_time_index is the
    index in out.T that belongs to the first entry of partial_field. The
    input-change cost needs two samples, so its first entry belongs to
    out.T[1] and its offset is 1.
    """

    partial_field: str
    first_time_index: int


TARGETS: Dict[str, TargetSpec] = {
    "SSE": TargetSpec("partial_SSE", 0),
    "SSdU": TargetSpec("partial_SSdU", 1),
}


@dataclass
class RunCurve:
    """One training run: its cumulative cost against fidelity."""

    source: str
    z: np.ndarray
    log_j: np.ndarray

    @property
    def n_samples(self) -> int:
        return int(self.z.size)


@dataclass
class FitResult:
    """Everything needed to reproduce and audit one fit of one target."""

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

def beta_ratio(z: np.ndarray, ab: Sequence[float]) -> np.ndarray:
    """f(z) = I_z(a, b), the regularised incomplete beta function.

    Increasing in z with f(0) = 0 and f(1) = 1 for every positive a and b.
    Arguments are clipped to the domain betainc accepts before the call: the
    bounds passed to the optimiser hold at the solution but not at every trial
    point a line search proposes. The clip limits on (a, b) are far wider than
    any parameter box this module uses, so they act only on infeasible trial
    points.

    Note the argument order. scipy.special.betainc takes (a, b, x) whereas
    MATLAB's betainc takes (x, a, b); the two agree once the order is matched.
    """
    a = float(np.clip(ab[0], AB_CLIP[0], AB_CLIP[1]))
    b = float(np.clip(ab[1], AB_CLIP[0], AB_CLIP[1]))
    x = np.clip(np.asarray(z, dtype=float), 0.0, 1.0)
    return betainc(a, b, x)


def _profiled_residuals(ab: Sequence[float], z: np.ndarray, log_j: np.ndarray,
                        group: np.ndarray, n_groups: int) -> np.ndarray:
    """u - mean_z u per run, the residual left after profiling out log J_i(1)."""
    u = log_j - np.log(np.maximum(beta_ratio(z, ab), F_FLOOR))
    counts = np.bincount(group, minlength=n_groups).astype(float)
    sums = np.bincount(group, weights=u, minlength=n_groups)
    u_bar = sums / np.maximum(counts, 1.0)
    return u - u_bar[group]


def beta_loss(ab: Sequence[float], z: np.ndarray, log_j: np.ndarray,
              group: np.ndarray, n_groups: int) -> float:
    """Mean squared deviation of u = log J - log f from its per-run mean.

    Centring on the per-run mean is what profiling out the unknown log J_i(1)
    amounts to, so the loss holds no run total. Averaging over samples instead
    of summing makes the value comparable between training sets of different
    size, which one shared grid of L2 strengths requires.
    """
    r = _profiled_residuals(ab, z, log_j, group, n_groups)
    return float(np.mean(r ** 2))


def beta_cost(ab: Sequence[float], z: np.ndarray, log_j: np.ndarray,
              group: np.ndarray, n_groups: int, lam: float) -> float:
    """Data loss plus an L2 pull of (a, b) towards (1, 1).

    The penalty centre is the identity f(z) = z, so shrinkage moves the
    surrogate towards a cost that accumulates at a constant rate.
    """
    ab = np.asarray(ab, dtype=float)
    return beta_loss(ab, z, log_j, group, n_groups) + lam * float(np.sum((ab - 1.0) ** 2))


def _central_gradient(fun, ab: np.ndarray, step: float = OPT_FD_STEP) -> np.ndarray:
    """Central-difference gradient in two dimensions.

    The profiled loss is smooth but is evaluated through betainc, whose forward
    difference is noisy at the step sizes L-BFGS-B would choose by default. Two
    parameters make the central difference cheap.
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


def fit_beta(z: np.ndarray, log_j: np.ndarray, group: np.ndarray, n_groups: int,
             lam: float, ab0: Sequence[float] = AB0) -> Tuple[np.ndarray, Dict]:
    """Minimise beta_cost over the two shape parameters under box bounds."""
    def obj(p):
        return beta_cost(p, z, log_j, group, n_groups, lam)

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
    }
    return np.asarray(res.x, dtype=float), info


def select_lambda(z: np.ndarray, log_j: np.ndarray, group: np.ndarray, n_groups: int,
                  lambda_grid: Sequence[float] = LAMBDA_GRID,
                  k_fold: int = K_FOLD, n_repeats: int = N_CV_REPEATS,
                  seed: int = RNG_SEED) -> Tuple[float, np.ndarray, int]:
    """Pick the L2 strength by repeated k-fold cross-validation over runs.

    The runs of the training set are split into k folds, each fold in turn is
    scored by the unregularised beta_loss under parameters fitted on the other
    folds, and the partition is redrawn n_repeats times. The grid is walked from
    weak to strong shrinkage with the previous solution as the start point,
    which keeps the number of optimiser iterations small along the path.

    Returns the plain argmin of the averaged held-out loss, the full loss grid,
    and the number of folds that were scored. Ties take the first grid entry,
    which is the weakest shrinkage among the tied values.
    """
    lambda_grid = list(lambda_grid)
    n_lam = len(lambda_grid)
    k = min(k_fold, n_groups)
    cv_loss = np.zeros(n_lam, dtype=float)
    n_scored = 0

    if k < 2:
        # A single run cannot be split, so no lambda can be scored. Fall back to
        # the weakest shrinkage and report zero folds, which the caller records.
        return float(lambda_grid[0]), np.full(n_lam, np.nan), 0

    rng = np.random.default_rng(seed)
    for _ in range(int(n_repeats)):
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
                ab, _ = fit_beta(z_tr, lj_tr, g_tr, n_tr, lam, ab0=ab)
                cv_loss[i_lam] += beta_loss(ab, z_he, lj_he, g_he, n_he)
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

    The two initial-condition cases of the file are added together, matching the
    definition of the truncated objective. offset is the index in out.T that
    holds the first entry of the partial array.
    """
    spec = TARGETS[target]
    try:
        mat = scipy.io.loadmat(str(path), struct_as_record=False, squeeze_me=True)
    except NotImplementedError as exc:
        # scipy reads MATLAB v7 and below. A v7.3 file is HDF5 and raises here,
        # which would otherwise surface as an unexplained refit failure halfway
        # through a run.
        raise RuntimeError(
            f"{path.name} is a MATLAB v7.3 file, which this fit cannot read. The "
            f"simulation drivers must call save() without the '-v7.3' flag."
        ) from exc
    if "out" not in mat:
        raise ValueError(f"{path} has no 'out' struct")
    out = mat["out"]

    time = np.asarray(out.T, dtype=float).ravel()
    cases = out.case
    if np.ndim(cases) == 0:
        cases = [cases]
    if len(cases) < 2:
        raise ValueError(f"{path} holds {len(cases)} case(s), expected 2")

    p1 = np.asarray(getattr(cases[0], spec.partial_field), dtype=float).ravel()
    p2 = np.asarray(getattr(cases[1], spec.partial_field), dtype=float).ravel()

    n = int(min(p1.size, p2.size, time.size - spec.first_time_index))
    if n <= 0:
        raise ValueError(f"{path} holds no usable samples for {target}")

    z = time[spec.first_time_index: spec.first_time_index + n] / horizon_hours
    j = np.cumsum(p1[:n]) + np.cumsum(p2[:n])
    with np.errstate(divide="ignore"):
        log_j = np.log(j)
    return RunCurve(source=str(path), z=z, log_j=log_j)


def pool_runs(curves: Sequence[RunCurve], z_min: float = Z_MIN):
    """Flat sample arrays over the usable runs, tagged by run.

    One flat vector lets the loss call betainc once per evaluation instead of
    once per run. A run needs two samples to say anything after its own constant
    is profiled out, so shorter ones are dropped and reported.
    """
    z_parts: List[np.ndarray] = []
    lj_parts: List[np.ndarray] = []
    g_parts: List[np.ndarray] = []
    used: List[str] = []
    dropped: List[str] = []

    g = 0
    for c in curves:
        keep = (c.z >= z_min) & np.isfinite(c.log_j)
        if int(np.count_nonzero(keep)) < MIN_SAMPLES_PER_RUN:
            dropped.append(c.source)
            continue
        z_parts.append(c.z[keep])
        lj_parts.append(c.log_j[keep])
        g_parts.append(np.full(int(np.count_nonzero(keep)), g, dtype=int))
        used.append(c.source)
        g += 1

    if g == 0:
        return (np.empty(0), np.empty(0), np.empty(0, dtype=int), 0, used, dropped)

    return (np.concatenate(z_parts), np.concatenate(lj_parts),
            np.concatenate(g_parts), g, used, dropped)


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def fit_target(paths: Sequence[Path], target: str, *, z_min: float = Z_MIN,
               lambda_grid: Sequence[float] = LAMBDA_GRID, k_fold: int = K_FOLD,
               n_repeats: int = N_CV_REPEATS, seed: int = RNG_SEED,
               horizon_hours: float = HORIZON_HOURS) -> FitResult:
    """Fit f(z) for one target from an explicit list of output files.

    The caller passes the file list rather than a directory so that a vintage
    fitted mid-run records exactly which evaluations informed it.
    """
    curves = [load_curve(Path(p), target, horizon_hours) for p in paths]
    z, log_j, group, n_groups, used, dropped = pool_runs(curves, z_min)

    if n_groups == 0:
        raise ValueError(
            f"no usable run for target {target}: all {len(curves)} candidates fell "
            f"below {MIN_SAMPLES_PER_RUN} samples at z >= {z_min}"
        )

    lam, cv_loss, n_scored = select_lambda(
        z, log_j, group, n_groups, lambda_grid, k_fold, n_repeats, seed)
    ab, opt_info = fit_beta(z, log_j, group, n_groups, lam)
    loss = beta_loss(ab, z, log_j, group, n_groups)

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
            "horizon_hours": float(horizon_hours),
            "z_min": float(z_min),
            "k_fold": int(k_fold),
            "n_cv_repeats": int(n_repeats),
            "cv_seed": int(seed),
            "lambda_selection": "argmin of mean held-out loss, no one-standard-error rule",
            "min_samples_per_run": int(MIN_SAMPLES_PER_RUN),
            "f_floor": F_FLOOR,
            "ab_clip": list(AB_CLIP),
        },
    )


def fit_all_targets(paths: Sequence[Path], **kwargs) -> Dict[str, FitResult]:
    """Fit every target on the same file list."""
    return {name: fit_target(paths, name, **kwargs) for name in TARGETS}


def write_coefficients(results: Dict[str, FitResult], mat_path: Path, vintage: int,
                       created_at: str, extra: Dict | None = None) -> None:
    """Publish the coefficient handoff file that MATLAB loads.

    The write is atomic: the payload goes to a temporary name in the same
    directory and is then renamed over the target, so a reader either sees the
    previous vintage or the new one and never a partial file.
    """
    mat_path = Path(mat_path)
    mat_path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "vintage": np.array([[float(vintage)]]),
        "created_at": str(created_at),
        "model": "regularised_incomplete_beta",
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
    tmp.replace(mat_path)


def write_vintage_record(results: Dict[str, FitResult], json_path: Path, vintage: int,
                         created_at: str, extra: Dict | None = None) -> None:
    """Write the full audit record of one vintage."""
    json_path = Path(json_path)
    json_path.parent.mkdir(parents=True, exist_ok=True)

    record = {
        "vintage": int(vintage),
        "created_at": str(created_at),
        "model": "f(z) = I_z(a, b), regularised incomplete beta",
        "targets": {name: r.to_dict() for name, r in results.items()},
    }
    if extra:
        record["context"] = extra

    tmp = json_path.with_suffix(json_path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2, sort_keys=True)
        fh.flush()
    tmp.replace(json_path)


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

    print(f"fitting on {len(paths)} files")
    results = fit_all_targets(paths)
    for name, r in results.items():
        print(f"\n{name}: a = {r.a:.6f}, b = {r.b:.6f}, lambda = {r.lam:g}, "
              f"loss = {r.loss:.6e}")
        print(f"  runs used {r.n_runs_used}/{r.n_runs_total}, samples {r.n_samples}, "
              f"folds scored {r.cv_folds_scored}")
        for lam, cv in zip(r.lambda_grid, r.cv_loss):
            print(f"    lambda {lam:<8g} cv loss {cv:.6e}")
        print(f"  optimiser: {r.optimizer['message']} "
              f"(nit={r.optimizer['nit']}, nfev={r.optimizer['nfev']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv))
