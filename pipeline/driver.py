"""Runtime-aware multi-fidelity Bayesian optimization driver, in two phases.

Phase "init" evaluates the Sobol design and stops. These evaluations use no
surrogate. MATLAB reports the cost that it measured at the simulated fidelity
and stores the per-step cost trends. The driver then fits the fidelity surrogate
phi(z) = I_z(a, b) from those runs and calls it vintage 0.

Phase "bo" runs the acquisition loop against phi. It refits phi every
refit_every iterations.

A refit does not rescale an earlier row. An objective value keeps the estimate
that the vintage in force produced when MATLAB measured it. Every row records
that vintage. You can therefore recompute the whole history under one fit later.

    python run_pipeline.py --case case1        # runs both phases
    python run_pipeline.py --case case1 --phase init

Start MATLAB on main_initialization first. It hands over to main_BO by itself.

Both phases resume after an interruption. The driver rebuilds its state from the
results CSV that MATLAB writes. It compares that state against its own ledger.
It then continues from the first evaluation that neither file records.

The driver seeds each proposal from the iteration index. A resumed run therefore
proposes what an uninterrupted run would propose at the same index.

The driver measures the wall time of each step and stores it in the ledger.
"""

from __future__ import annotations

import sys
import time
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import torch
from botorch.acquisition.acquisition import AcquisitionFunction
from botorch.acquisition.multi_objective.logei import (
    qLogNoisyExpectedHypervolumeImprovement,
)
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.models.transforms import Normalize, Standardize
from botorch.optim import optimize_acqf
from botorch.sampling.normal import SobolQMCNormalSampler
from gpytorch.mlls import ExactMarginalLogLikelihood

from pipeline.matlab_interface import (
    PHI_COEFFS_FILE,
    EvaluationFailed,
    failures_file,
    out_dir,
    read_results,
    results_file,
    send_request,
    wait_for_matlab_ready,
    wait_for_result,
)
from pipeline.provenance import Registry, summarise_gp
from run_config import INTEGER_IDXS, THETA_D, RunConfig, parse_args

from pipeline.phi_surrogate import (
    fit_all_targets,
    write_coefficients,
    write_vintage_record,
)

# The project root is one level above this package. Every results path hangs
# off it, so moving this file deeper would need this line changed and nothing
# else.
BASE_DIR = Path(__file__).resolve().parents[1]
DEVICE = torch.device("cpu")
DTYPE = torch.double


# ---------------------------------------------------------------------------
# Search space helpers
# ---------------------------------------------------------------------------

class Space:
    """Bounds, fixed components and the projection between full and free theta."""

    def __init__(self, cfg: RunConfig):
        lb, ub = cfg.bounds()
        self.cfg = cfg
        self.spec = cfg.spec()
        self.lb = torch.tensor(lb, dtype=DTYPE, device=DEVICE)
        self.ub = torch.tensor(ub, dtype=DTYPE, device=DEVICE)
        self.opt_idxs = self.spec.opt_idxs
        self.lb_opt = self.lb[self.opt_idxs]
        self.ub_opt = self.ub[self.opt_idxs]

    def snap(self, X: torch.Tensor) -> torch.Tensor:
        """Clamp to bounds, round the integer components, reimpose the fixed ones."""
        X = X.clone()
        X = torch.max(torch.min(X, self.ub), self.lb)
        for idx in INTEGER_IDXS:
            X[..., idx] = torch.round(X[..., idx])
        for idx, val in self.spec.fixed.items():
            X[..., idx] = val
        return X

    def to_full(self, X_opt: torch.Tensor) -> torch.Tensor:
        shape = X_opt.shape[:-1] + (THETA_D,)
        X = torch.empty(shape, dtype=DTYPE, device=DEVICE)
        X[..., self.opt_idxs] = X_opt
        for idx, val in self.spec.fixed.items():
            X[..., idx] = val
        return self.snap(X)

    def to_opt(self, X_full: torch.Tensor) -> torch.Tensor:
        return self.snap(X_full)[..., self.opt_idxs]

    def acq_bounds(self) -> torch.Tensor:
        return torch.stack([self.lb_opt, self.ub_opt])


def sobol_design(space: Space, n: int, seed: int) -> torch.Tensor:
    """The full initialisation design, drawn once and indexed by position.

    Every point of the design is generated regardless of how many have already
    been evaluated, and the resumed run takes the ones it still needs by index.
    Drawing only the missing points would advance the Sobol stream by a
    different amount after each interruption, so the design would depend on when
    the run was interrupted.
    """
    d = len(space.opt_idxs)
    engine = torch.quasirandom.SobolEngine(dimension=d, scramble=True, seed=seed)
    raw = engine.draw(n * 4).to(device=DEVICE, dtype=DTYPE)
    candidates = space.lb_opt + (space.ub_opt - space.lb_opt) * raw

    chosen: List[torch.Tensor] = []
    for i in range(candidates.shape[0]):
        point = space.to_opt(space.to_full(candidates[i]))
        if any(torch.allclose(point, p, atol=1e-12) for p in chosen):
            continue                       # duplicate after rounding the integers
        chosen.append(point)
        if len(chosen) == n:
            break

    if len(chosen) < n:
        raise RuntimeError(
            f"the Sobol design produced only {len(chosen)} distinct points of the "
            f"{n} requested; widen the integer bounds or lower n_init")
    return torch.stack(chosen)


# ---------------------------------------------------------------------------
# Models and acquisition
# ---------------------------------------------------------------------------

def build_models(space: Space, X_opt: torch.Tensor, Y_obj: torch.Tensor,
                 Y_time: torch.Tensor, eps: float):
    """Fit the objective GP and the log-runtime GP on the data collected so far.

    The function returns the two models and the wall time of each fit.
    """
    bounds = space.acq_bounds().to(device=DEVICE, dtype=DTYPE)
    Y_time_log = torch.log(Y_time.clamp_min(eps))

    model_obj = SingleTaskGP(
        X_opt, Y_obj,
        input_transform=Normalize(d=X_opt.shape[-1], bounds=bounds),
        outcome_transform=Standardize(m=Y_obj.shape[-1]),
    )
    model_time = SingleTaskGP(
        X_opt, Y_time_log,
        input_transform=Normalize(d=X_opt.shape[-1], bounds=bounds),
        outcome_transform=Standardize(m=1),
    )

    t0 = time.perf_counter()
    fit_gpytorch_mll(ExactMarginalLogLikelihood(model_obj.likelihood, model_obj))
    wall_obj = time.perf_counter() - t0

    t0 = time.perf_counter()
    fit_gpytorch_mll(ExactMarginalLogLikelihood(model_time.likelihood, model_time))
    wall_time = time.perf_counter() - t0

    return model_obj, model_time, {"gp_objective": wall_obj, "gp_log_runtime": wall_time}


class RuntimeAwareLogAcq(AcquisitionFunction):
    """log alpha = log qLogNEHVI + w_z log(a_z + eps) - w_t gamma log(E[t] + eps).

    The fidelity bias a_z has its center at z = 1. It rewards a candidate that is
    closer to full fidelity. The runtime term penalizes a candidate whose
    predicted evaluation cost is large. The runtime GP uses log t, so the
    function takes the expectation in lognormal form.

    The two terms do not share units. qLogNEHVI carries the units of the
    hypervolume in standardized objective space. Its magnitude drifts as the
    frontier fills in. The effective weight of the runtime penalty is therefore
    not constant over a run.
    """

    def __init__(self, qlognehvi, model_time_log, *, ell_z: float, gamma: float,
                 eps: float, w_z: float, w_t: float, t_floor: float,
                 max_log_et: float = 50.0):
        super().__init__(model=qlognehvi.model)
        self.qlognehvi = qlognehvi
        self.model_time_log = model_time_log
        self.ell_z = ell_z
        self.gamma = gamma
        self.eps = eps
        self.w_z = w_z
        self.w_t = w_t
        self.t_floor = float(t_floor)
        self.max_log_et = float(max_log_et)

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        log_hv = torch.nan_to_num(self.qlognehvi(X), neginf=-1e6, posinf=1e6)

        z = X[..., 0]
        a_z = torch.exp(-((1.0 - z) ** 2) / (2.0 * self.ell_z ** 2)).mean(dim=-1)
        log_az = torch.nan_to_num(torch.log(a_z + self.eps), neginf=-1e6, posinf=1e6)

        post = self.model_time_log.posterior(X)
        log_et = (post.mean + 0.5 * post.variance).mean(dim=-2).squeeze(-1)
        log_et = torch.nan_to_num(log_et, neginf=-1e6, posinf=1e6).clamp_max(self.max_log_et)
        et = torch.exp(log_et).clamp_min(self.t_floor)
        log_pen = torch.nan_to_num(-self.gamma * torch.log(et + self.eps),
                                   neginf=-1e6, posinf=1e6)

        return torch.nan_to_num(log_hv + self.w_z * log_az + self.w_t * log_pen,
                                neginf=-1e6, posinf=1e6)


def propose(space: Space, cfg: RunConfig, X: torch.Tensor, Y_obj: torch.Tensor,
            Y_time: torch.Tensor, seed: int) -> Tuple[torch.Tensor, Dict]:
    """Fit the surrogates, maximise the acquisition, return the candidate.

    The seed is set from the iteration index rather than left to accumulate, so
    a proposal depends on the data and the index alone. A resumed run therefore
    proposes what an uninterrupted one would have proposed at the same index.
    """
    t_propose = time.perf_counter()
    torch.manual_seed(seed)

    X_opt = space.to_opt(X)
    model_obj, model_time, wall_fit = build_models(space, X_opt, Y_obj, Y_time, cfg.eps)

    y_min = Y_obj.min(dim=0).values
    y_range = (Y_obj.max(dim=0).values - y_min).clamp_min(1e-6)
    ref_point = (y_min - cfg.ref_point_backoff * y_range).tolist()

    acq_inner = qLogNoisyExpectedHypervolumeImprovement(
        model=model_obj,
        ref_point=ref_point,
        X_baseline=X_opt,
        sampler=SobolQMCNormalSampler(sample_shape=torch.Size([cfg.mc_samples])),
        prune_baseline=cfg.prune_baseline,
        cache_pending=True,
        max_iep=cfg.max_iep,
    )

    t_floor = float(Y_time.min().clamp_min(cfg.eps).item())
    acq = RuntimeAwareLogAcq(
        acq_inner, model_time,
        ell_z=cfg.ell_z, gamma=cfg.gamma_time, eps=cfg.eps,
        w_z=cfg.w_z, w_t=cfg.w_t, t_floor=t_floor,
    )

    t0 = time.perf_counter()
    cand_opt, acq_at_optimum = optimize_acqf(
        acq_function=acq,
        bounds=space.acq_bounds(),
        q=cfg.q_batch,
        num_restarts=cfg.num_restarts,
        raw_samples=cfg.raw_samples,
        options={"batch_limit": cfg.acq_batch_limit, "maxiter": cfg.acq_maxiter},
    )
    wall_acq = time.perf_counter() - t0

    cand_full = space.to_full(cand_opt)
    with torch.no_grad():
        acq_at_snapped = float(acq(space.to_opt(cand_full)).detach().cpu().view(-1)[0])

    diagnostics = {
        "acq_value_at_optimum": float(acq_at_optimum.detach().cpu().view(-1)[0]),
        "acq_value_at_snapped": acq_at_snapped,
        "acq_seed": int(seed),
        "ref_point": [float(v) for v in ref_point],
        "t_floor": t_floor,
        "n_train": int(X.shape[0]),
        "gp_objective": summarise_gp(model_obj, "objective"),
        "gp_log_runtime": summarise_gp(model_time, "log_runtime"),
        "wall_s": {
            "gp_objective_fit": wall_fit["gp_objective"],
            "gp_log_runtime_fit": wall_fit["gp_log_runtime"],
            "acquisition_maximisation": wall_acq,
            "total": time.perf_counter() - t_propose,
        },
    }
    return cand_full.view(-1), diagnostics


# ---------------------------------------------------------------------------
# History
# ---------------------------------------------------------------------------

def load_history(paths: List[Path]) -> List[Dict]:
    """Read every results row of the given files, in file order then eval_id."""
    rows: List[Dict] = []
    for path in paths:
        rows.extend(read_results(path))
    rows.sort(key=lambda r: r["eval_id"])
    return rows


def history_tensors(rows: List[Dict]) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Assemble the optimiser's view of the history.

    Y_obj is the negated objective vector as it was recorded. The values scaled
    by an earlier surrogate vintage are used unchanged: a refit does not revise
    them. The measured costs and the divisors are in the CSV, so the whole
    history can be recomputed under one vintage afterwards.
    """
    if not rows:
        return (torch.empty((0, THETA_D), dtype=DTYPE, device=DEVICE),
                torch.empty((0, 2), dtype=DTYPE, device=DEVICE),
                torch.empty((0, 1), dtype=DTYPE, device=DEVICE))

    X = torch.tensor([r["theta"] for r in rows], dtype=DTYPE, device=DEVICE)
    Y_obj = torch.tensor([[-r["SSE"], -r["SSdU"]] for r in rows],
                         dtype=DTYPE, device=DEVICE)
    Y_time = torch.tensor([[r["runtime_s"]] for r in rows],
                          dtype=DTYPE, device=DEVICE).clamp_min(1e-9)
    return X, Y_obj, Y_time


# ---------------------------------------------------------------------------
# Surrogate refit
# ---------------------------------------------------------------------------

def fit_vintage(cfg: RunConfig, registry: Registry, vintage: int,
                init_rows: List[Dict], bo_rows: List[Dict]) -> Dict:
    """Fit, publish and record one surrogate vintage.

    Vintage 0 is fitted on the initialisation runs alone. Vintage v is fitted on
    those plus the first v * refit_every optimisation runs. The file list is
    resolved from the ledger, so the record names exactly which evaluations
    informed the fit.
    """
    n_bo = vintage * cfg.refit_every
    used_rows = list(init_rows) + list(bo_rows[:n_bo])

    paths: List[Path] = []
    for row in used_rows:
        phase_dir = out_dir("init") if row["phase"] == "DOE" else out_dir("bo")
        candidate = phase_dir / f"out_{row['timestamp']}.mat"
        if candidate.exists():
            paths.append(candidate)
        else:
            print(f"[fit] missing trends file for eval {row['eval_id']}: {candidate.name}")

    if not paths:
        raise RuntimeError(f"vintage {vintage}: no out_*.mat file available to fit on")

    print(f"[fit] vintage {vintage}: fitting on {len(paths)} runs "
          f"({len(init_rows)} DOE + {n_bo} OPT requested)")
    t_fit = time.perf_counter()
    t0 = time.perf_counter()
    results = fit_all_targets(
        paths,
        lambda_grid=cfg.fit_lambda_grid,
        k_fold=cfg.fit_k_fold,
        seed=cfg.cv_seed,
        horizon_hours=cfg.horizon_hours,
    )
    elapsed = time.perf_counter() - t0

    created_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    context = {
        "case": cfg.case,
        "n_doe_rows": len(init_rows),
        "n_opt_rows_used": n_bo,
        "eval_ids_used": [r["eval_id"] for r in used_rows],
        "files_used": [p.name for p in paths],
        "files_missing": len(used_rows) - len(paths),
        "fit_wall_s": elapsed,
        "governs_iterations": [vintage * cfg.refit_every + 1,
                               (vintage + 1) * cfg.refit_every],
        "rescales_past_rows": False,
    }

    wall_record = write_vintage_record(results, registry.vintage_path(vintage),
                                       vintage, created_at, extra=context)
    wall_publish = write_coefficients(results, PHI_COEFFS_FILE, vintage, created_at)
    context["wall_s"] = {
        "fit_total": elapsed,
        "per_target": {n: r.wall_s for n, r in results.items()},
        "write_record": wall_record,
        "publish_coefficients": wall_publish,
        "total": time.perf_counter() - t_fit,
    }
    write_vintage_record(results, registry.vintage_path(vintage), vintage,
                         created_at, extra=context)

    for name, r in results.items():
        print(f"[fit]   {name}: a = {r.a:.6f}, b = {r.b:.6f}, lambda = {r.lam:g}, "
              f"loss = {r.loss:.4e} (runs {r.n_runs_used}/{r.n_runs_total})")
    print(f"[fit] vintage {vintage} published to {PHI_COEFFS_FILE.name} in {elapsed:.1f} s")

    return {"vintage": vintage, "created_at": created_at, "context": context,
            "targets": {n: r.to_dict() for n, r in results.items()}}


def ensure_vintage(cfg: RunConfig, registry: Registry, vintage: int,
                   init_rows: List[Dict], bo_rows: List[Dict]) -> None:
    """Publish the vintage that should govern the next iteration.

    Called before every proposal rather than only on refit iterations, so an
    interrupted run republishes the correct coefficients on restart instead of
    continuing against whatever file happened to survive.
    """
    record = registry.load_vintage(vintage)
    if record is not None and PHI_COEFFS_FILE.exists():
        import scipy.io
        published = float(scipy.io.loadmat(str(PHI_COEFFS_FILE))["vintage"].ravel()[0])
        if int(published) == vintage:
            return
        print(f"[fit] published coefficients are vintage {int(published)}, "
              f"iteration needs {vintage}; republishing from the record")
        write_coefficients(
            _results_from_record(record), PHI_COEFFS_FILE, vintage,
            record.get("created_at", ""))
        return

    fit_vintage(cfg, registry, vintage, init_rows, bo_rows)


def _results_from_record(record: Dict):
    """Rebuild the coefficient payload from a stored vintage record."""
    from pipeline.phi_surrogate import FitResult

    out = {}
    for name, t in record["targets"].items():
        out[name] = FitResult(
            target=name, a=t["a"], b=t["b"], lam=t["lam"], loss=t["loss"],
            lambda_grid=t["lambda_grid"], cv_loss=t["cv_loss"],
            cv_folds_scored=t["cv_folds_scored"], n_runs_total=t["n_runs_total"],
            n_runs_used=t["n_runs_used"], n_samples=t["n_samples"],
        )
    return out


# ---------------------------------------------------------------------------
# Phases
# ---------------------------------------------------------------------------

def run_initialization(cfg: RunConfig) -> None:
    space = Space(cfg)
    registry = Registry(BASE_DIR / "results", "init")
    registry.write_manifest(cfg.to_dict(), BASE_DIR)

    path_results = results_file("init")
    path_failures = failures_file("init")
    wait_for_matlab_ready(path_results, cfg.max_wait_matlab_s, cfg.matlab_wait_s)

    rows = load_history([path_results])
    registry.reconcile(rows)
    done = {r["eval_id"] for r in rows}

    design = sobol_design(space, cfg.n_init, cfg.sobol_seed)
    print(f"[init] {len(done)} of {cfg.n_init} design points already evaluated")

    # Iteration wall times of this session, used for the estimate of time left.
    # Rows from an earlier session are excluded: they may come from a different
    # machine or a different fmincon budget.
    iter_walls: List[float] = []

    for i in range(cfg.n_init):
        eval_id = i + 1
        if eval_id in done:
            continue

        theta = space.to_full(design[i]).view(-1)
        theta_list = [float(v) for v in theta.tolist()]
        sent_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")

        t_wait = time.perf_counter()
        send_request(eval_id, "init", theta_list, lock_stale_s=cfg.lock_stale_s)
        try:
            row = wait_for_result(eval_id, path_results, path_failures,
                                  poll_s=cfg.poll_s, timeout_s=cfg.eval_timeout_s)
        except EvaluationFailed as exc:
            registry.append_evaluation({
                "eval_id": eval_id, "phase": "DOE", "case": cfg.case,
                "theta": theta_list, "sent_at": sent_at,
                "failed": True, "failure": exc.failure,
            })
            print(f"[init] design point {eval_id} failed; it is recorded and skipped")
            continue

        wall_matlab = time.perf_counter() - t_wait
        registry.append_evaluation({
            "eval_id": eval_id, "phase": "DOE", "case": cfg.case,
            "doe_index": eval_id, "theta": theta_list,
            "sobol_seed": cfg.sobol_seed, "sent_at": sent_at,
            "recorded_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "phi_vintage_applied": None,
            "result": row,
            "wall_s": {"matlab_round_trip": wall_matlab},
        })
        iter_walls.append(wall_matlab)
        eta = _eta(iter_walls, cfg.n_init - eval_id)
        print(f"[init] {eval_id:3d}/{cfg.n_init}  z={row['z']:.4f}  "
              f"SSE={row['SSE']:.6g}  SSdU={row['SSdU']:.6g}  "
              f"solver={row['runtime_s']:.1f}s  matlab={wall_matlab:.1f}s"
              + (f"\n       {eta}" if eta else ""))

    rows = load_history([path_results])
    print(f"\n[init] initialisation complete: {len(rows)} evaluations")

    if not registry.vintage_exists(0):
        bo_registry = Registry(BASE_DIR / "results", "bo")
        fit_vintage(cfg, bo_registry, 0, rows, [])
        print("[init] vintage 0 fitted and published")
    print("[init] MATLAB hands over to main_BO on the first optimization request.")


def run_bo(cfg: RunConfig) -> None:
    space = Space(cfg)
    registry = Registry(BASE_DIR / "results", "bo")
    registry.write_manifest(cfg.to_dict(), BASE_DIR)

    path_init = results_file("init")
    path_results = results_file("bo")
    path_failures = failures_file("bo")

    # No readiness wait here. main_BO creates results.csv only after
    # main_initialization hands over, and the handover needs the first
    # optimization request. Waiting for the file before sending that request
    # would block for the whole timeout every time.

    init_rows = load_history([path_init])
    if len(init_rows) < cfg.n_init:
        raise RuntimeError(
            f"only {len(init_rows)} of {cfg.n_init} initialisation evaluations found "
            f"in {path_init}. Finish phase 'init' first.")

    bo_rows = load_history([path_results])
    registry.reconcile(bo_rows)
    print(f"[bo] resuming with {len(init_rows)} DOE and {len(bo_rows)} OPT evaluations")

    done = {r["eval_id"] for r in bo_rows}

    # See the note in run_initialization: only this session's timings feed the
    # estimate of time left.
    iter_walls: List[float] = []

    for k in range(len(bo_rows), cfg.n_iter):
        bo_idx = k + 1
        eval_id = cfg.n_init + bo_idx
        if eval_id in done:
            continue

        t_iter = time.perf_counter()
        vintage = (bo_idx - 1) // cfg.refit_every

        t0 = time.perf_counter()
        ensure_vintage(cfg, registry, vintage, init_rows, bo_rows)
        wall_vintage = time.perf_counter() - t0

        X, Y_obj, Y_time = history_tensors(init_rows + bo_rows)
        theta, diagnostics = propose(space, cfg, X, Y_obj, Y_time,
                                     seed=cfg.torch_seed + eval_id)
        theta_list = [float(v) for v in theta.tolist()]
        sent_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")

        t_wait = time.perf_counter()
        send_request(eval_id, "bo", theta_list, lock_stale_s=cfg.lock_stale_s)
        try:
            row = wait_for_result(eval_id, path_results, path_failures,
                                  poll_s=cfg.poll_s, timeout_s=cfg.eval_timeout_s)
        except EvaluationFailed as exc:
            registry.append_evaluation({
                "eval_id": eval_id, "phase": "OPT", "bo_index": bo_idx,
                "case": cfg.case, "theta": theta_list, "sent_at": sent_at,
                "phi_vintage_expected": vintage, "proposal": diagnostics,
                "acquisition_settings": _acq_settings(cfg),
                "failed": True, "failure": exc.failure,
            })
            print(f"[bo] iteration {bo_idx} failed in MATLAB; recorded and skipped")
            continue
        wall_matlab = time.perf_counter() - t_wait

        applied = int(row["phi_vintage"]) if row["phi_vintage"] == row["phi_vintage"] else None
        if applied != vintage:
            raise RuntimeError(
                f"evaluation {eval_id} was scaled by surrogate vintage {applied} but "
                f"iteration {bo_idx} required vintage {vintage}. MATLAB read a stale "
                f"coefficient file; the row is on disk but the run is stopped rather "
                f"than continuing with mixed vintages."
            )

        t0 = time.perf_counter()
        registry.append_evaluation({
            "eval_id": eval_id, "phase": "OPT", "bo_index": bo_idx,
            "case": cfg.case, "theta": theta_list, "sent_at": sent_at,
            "recorded_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "phi_vintage_expected": vintage,
            "phi_vintage_applied": applied,
            "proposal": diagnostics,
            "acquisition_settings": _acq_settings(cfg),
            "result": row,
            "wall_s": {
                "ensure_vintage": wall_vintage,
                "propose": diagnostics["wall_s"]["total"],
                "matlab_round_trip": wall_matlab,
                "iteration": time.perf_counter() - t_iter,
            },
        })
        wall_ledger = time.perf_counter() - t0
        bo_rows.append(row)

        print(f"[bo] {bo_idx:3d}/{cfg.n_iter} [phi v{vintage}]  z={row['z']:.4f}  "
              f"SSE={row['SSE']:.6g}  SSdU={row['SSdU']:.6g}  "
              f"acq={diagnostics['acq_value_at_snapped']:.4g}")
        wall_iter = time.perf_counter() - t_iter
        print(f"       wall: propose {diagnostics['wall_s']['total']:.1f}s "
              f"(gp {diagnostics['wall_s']['gp_objective_fit']:.1f}+"
              f"{diagnostics['wall_s']['gp_log_runtime_fit']:.1f}s, "
              f"acq {diagnostics['wall_s']['acquisition_maximisation']:.1f}s), "
              f"matlab {wall_matlab:.1f}s, ledger {wall_ledger:.3f}s, "
              f"iteration {wall_iter:.1f}s")

        iter_walls.append(wall_iter)
        eta = _eta(iter_walls, cfg.n_iter - bo_idx)
        if eta:
            print(f"       {eta}")

        if bo_idx % cfg.refit_every == 0:
            next_vintage = bo_idx // cfg.refit_every
            is_terminal = bo_idx >= cfg.n_iter
            if not is_terminal or cfg.refit_after_last:
                if not registry.vintage_exists(next_vintage):
                    fit_vintage(cfg, registry, next_vintage, init_rows, bo_rows)
                    if is_terminal:
                        print(f"[fit] vintage {next_vintage} is terminal: it is "
                              f"recorded for post-hoc use and governs no evaluation")

    print(f"\n[bo] budget complete: {len(bo_rows)} optimisation evaluations")
    _print_summary(init_rows + bo_rows, registry)


ETA_WINDOW = 5


def _eta(iter_walls: List[float], n_remaining: int) -> str:
    """Return an estimate of the time left, from the recent iteration times.

    The estimate uses the median of the last ETA_WINDOW iterations and not the
    mean over the whole run. Two reasons. The cost of one evaluation grows with
    the fidelity z, and the fidelity bias in the acquisition pushes z toward 1
    as the run proceeds, so early iterations are systematically cheaper than
    late ones and a full-run mean is biased low. The median also absorbs a
    single slow evaluation without dragging the estimate.

    The range comes from the fastest and the slowest iteration in the same
    window. It states how uncertain the estimate is instead of implying a
    precision the data does not support.

    An estimate needs at least two finished iterations. Before that the function
    returns an empty string and the caller prints nothing.
    """
    if n_remaining <= 0 or len(iter_walls) < 2:
        return ""

    window = iter_walls[-ETA_WINDOW:]
    ordered = sorted(window)
    mid = len(ordered) // 2
    typical = ordered[mid] if len(ordered) % 2 else 0.5 * (ordered[mid - 1] + ordered[mid])

    remaining_s = typical * n_remaining
    low_s = ordered[0] * n_remaining
    high_s = ordered[-1] * n_remaining
    finish = time.strftime("%a %H:%M", time.localtime(time.time() + remaining_s))

    return (f"{n_remaining} left, eta {_hours(remaining_s)} "
            f"[{_hours(low_s)}-{_hours(high_s)}], finish ~{finish}")


def _hours(seconds: float) -> str:
    """Format a duration as minutes below one hour, otherwise as hours."""
    if seconds < 3600:
        return f"{seconds / 60:.0f} min"
    return f"{seconds / 3600:.1f} h"


def _acq_settings(cfg: RunConfig) -> Dict:
    return {
        "ell_z": cfg.ell_z, "gamma_time": cfg.gamma_time,
        "w_z": cfg.w_z, "w_t": cfg.w_t, "eps": cfg.eps,
        "ref_point_backoff": cfg.ref_point_backoff,
        "num_restarts": cfg.num_restarts, "raw_samples": cfg.raw_samples,
        "acq_batch_limit": cfg.acq_batch_limit, "acq_maxiter": cfg.acq_maxiter,
        "mc_samples": cfg.mc_samples, "prune_baseline": cfg.prune_baseline,
        "max_iep": cfg.max_iep, "q_batch": cfg.q_batch,
    }


def _print_summary(rows: List[Dict], registry: Registry | None = None) -> None:
    """Print the totals of a finished run, including where the time went."""
    if not rows:
        print("[done] no evaluations recorded")
        return

    best_sse = min(rows, key=lambda r: r["SSE"])
    best_ssdu = min(rows, key=lambda r: r["SSdU"])
    solver_s = sum(r["runtime_s"] for r in rows)
    sim_wall_s = sum(r.get("wall_total_s", 0.0) or 0.0 for r in rows)

    print("\n[done] summary")
    print(f"  evaluations         {len(rows)}")
    print(f"  solver time         {solver_s / 3600:.2f} h")
    print(f"  simulation wall     {sim_wall_s / 3600:.2f} h")
    print(f"  lowest SSE          {best_sse['SSE']:.6g}  (eval {best_sse['eval_id']})")
    print(f"  lowest SSdU         {best_ssdu['SSdU']:.6g}  (eval {best_ssdu['eval_id']})")

    if registry is None:
        return

    gp_s = acq_s = fit_s = 0.0
    for rec in registry.load_evaluations():
        w = (rec.get("proposal") or {}).get("wall_s") or {}
        gp_s += w.get("gp_objective_fit", 0.0) + w.get("gp_log_runtime_fit", 0.0)
        acq_s += w.get("acquisition_maximisation", 0.0)
    for v in registry.list_vintages():
        rec = registry.load_vintage(v) or {}
        fit_s += ((rec.get("context") or {}).get("wall_s") or {}).get("total", 0.0)

    overhead = gp_s + acq_s + fit_s
    print(f"  GP fits             {gp_s / 60:.1f} min")
    print(f"  acquisition         {acq_s / 60:.1f} min")
    print(f"  phi fits            {fit_s / 60:.1f} min")
    if sim_wall_s > 0:
        print(f"  driver overhead     {100 * overhead / (sim_wall_s + overhead):.1f}% "
              f"of the total wall time")


def main(argv: List[str]) -> int:
    phase, cfg = parse_args(argv)
    torch.set_default_dtype(DTYPE)
    torch.manual_seed(cfg.torch_seed)

    print(f"[run] phase={phase} case={cfg.case} "
          f"({cfg.spec().dimension} free dimensions) "
          f"budget={cfg.n_init}+{cfg.n_iter}")

    try:
        if phase == "init":
            run_initialization(cfg)
        else:
            run_bo(cfg)
    except KeyboardInterrupt:
        print("\n[run] interrupted. Every completed evaluation is on disk; "
              "restart the same command to resume.")
        return 130
    except Exception:
        traceback.print_exc()
        print("\n[run] stopped on an error. Every completed evaluation is on disk; "
              "fix the cause and restart the same command to resume.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
