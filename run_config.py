"""Declared configuration of one optimization run.

This file names every setting that changes what a run does. The driver copies
it into the manifest before the first evaluation. No setting that affects a
result needs a hand edit of the driver. You select the case by name. The budget
is a number and not the point at which somebody stopped the process. The seeds
are explicit.

Select a case on the command line:

    python run_pipeline.py --case case2
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Dict, List, Tuple

# ---------------------------------------------------------------------------
# Fidelity limits
# ---------------------------------------------------------------------------
# Two different limits act on the fidelity z. They are not the same quantity and
# they must keep separate names.
#
# Z_MIN_BO is the lower bound of the search space. The acquisition maximizer
# cannot propose a fidelity below it. It is the limit that decides how short an
# evaluation the optimizer can buy.
#
# The simulation grid uses Ts = 1 min over a 10 h horizon. A run of fidelity z
# therefore covers 600 * z minutes and steps 600 * z times:
#
#     z = 0.01  ->  6.0 min  ->  6 steps, 7 grid points
#     z = 0.10  ->   60 min  ->  60 steps, 61 grid points
#     z = 1.00  ->  600 min  ->  600 steps, 601 grid points
#
# decode_theta.m clamps z into [0, 1]. That clamp guards against a malformed
# value and is not a fidelity policy. A fidelity below Z_MIN_BO covers less than
# one sampling interval, so the clamp never decides anything in practice and
# Z_MIN_BO is the floor that MATLAB sees.
Z_MIN_BO: float = 0.01

# Z_MIN_PHI is the floor of the surrogate fit. A sample below it leaves the fit.
# It is larger than Z_MIN_BO on purpose. Near z = 0 the cumulative cost of a run
# starts close to zero, and the log of that ratio is then dominated by the first
# few steps. A very short run therefore says little about the shape of phi.
#
# A run below Z_MIN_PHI still receives a phi correction. It contributes no
# sample to the fit of phi.
Z_MIN_PHI: float = 0.1

# ---------------------------------------------------------------------------
# Search space
# ---------------------------------------------------------------------------
# theta = [z, theta_p, theta_m, q1..q3, r_u1..r_u3, r_du1..r_du3]
THETA_NAMES: Tuple[str, ...] = (
    "z", "theta_p", "theta_m",
    "q1", "q2", "q3",
    "r_u1", "r_u2", "r_u3",
    "r_du1", "r_du2", "r_du3",
)
THETA_D = len(THETA_NAMES)

BASE_LB: Tuple[float, ...] = (Z_MIN_BO, 0.0, 0.0) + (-3.0,) * 9
BASE_UB: Tuple[float, ...] = (1.00, 15.0, 7.0) + (3.0,) * 9

INTEGER_IDXS: Tuple[int, ...] = (1, 2)   # theta_p, theta_m are rounded after proposal

# The exponent that stands for a disabled weight. decode_theta forms 10^r, so
# this underflows to exactly zero in MATLAB. It is not inside the [-3, 3] box:
# the encoding cannot switch a weight off from within the search space, because
# the smallest representable entry there is 1e-3.
DISABLED_EXPONENT = -1000.0


@dataclass(frozen=True)
class CaseSpec:
    """Which components of theta the optimiser is allowed to move."""

    name: str
    description: str
    fixed: Dict[int, float] = field(default_factory=dict)

    @property
    def opt_idxs(self) -> List[int]:
        return [i for i in range(THETA_D) if i not in self.fixed]

    @property
    def dimension(self) -> int:
        return len(self.opt_idxs)


CASES: Dict[str, CaseSpec] = {
    "case1": CaseSpec(
        name="case1",
        description=(
            "Full 12-dimensional search. Every component of theta is free in its "
            "bounds, so Q and R_u may take any diagonal value in [1e-3, 1e3]^3 and "
            "the input-magnitude penalty contributes with a magnitude chosen by BO."
        ),
        fixed={},
    ),
    "case2": CaseSpec(
        name="case2",
        description=(
            "Reduced 8-dimensional search. q1 is held at 0, so Q11 = 1, and the "
            "three input-magnitude exponents are set to the disabled value, so R_u "
            "vanishes from the NMPC objective. The remaining eight components keep "
            "the bounds of case1."
        ),
        fixed={3: 0.0, 6: DISABLED_EXPONENT, 7: DISABLED_EXPONENT, 8: DISABLED_EXPONENT},
    ),
    "baseline": CaseSpec(
        name="baseline",
        description=(
            "Single-fidelity, runtime-unaware counterpart of case2. z is held at "
            "1, so every evaluation, in the design and in the optimisation, runs "
            "the full 10 h horizon and phi(1) = I_1(a, b) = 1 scales it by exactly "
            "one. It carries case2's reduced structure as well: q1 is held at 0, so "
            "Q11 = 1, and the three input-magnitude exponents are disabled, so R_u "
            "vanishes from the NMPC objective. The remaining seven components, "
            "theta_p, theta_m, q2, q3 and the three R_du exponents, keep the "
            "bounds of case1. It is the plain qLogNEHVI counterpart of case2 and "
            "shares the budget and seeds."
        ),
        fixed={0: 1.0, 3: 0.0, 6: DISABLED_EXPONENT, 7: DISABLED_EXPONENT, 8: DISABLED_EXPONENT},
    ),
}


@dataclass(frozen=True)
class RunConfig:
    """One optimisation run, start to finish."""

    case: str = "case1"

    # Budget. The run makes n_init evaluations of the Sobol design, then n_iter
    # acquisition-driven evaluations. The total is n_init + n_iter. The loop
    # ends when it reaches n_iter. No other criterion stops it.
    n_init: int = 20
    n_iter: int = 100
    q_batch: int = 1

    # Surrogate refit. Vintage 0 uses the initialization runs alone. It governs
    # optimization iterations 1 to refit_every. The driver fits vintage v after
    # iteration v * refit_every completes. Vintage v then governs the next
    # refit_every iterations. The driver never rescales an earlier row. A row
    # keeps the estimate that the vintage in force produced when MATLAB measured
    # it.
    refit_every: int = 10
    refit_after_last: bool = True   # fit a last vintage for later analysis

    # Seeds.
    sobol_seed: int = 1234          # scrambled Sobol stream of the DOE
    torch_seed: int = 0             # GP fitting and acquisition multistarts
    cv_seed: int = 1                # fold partitions of the surrogate fit

    # Surrogate fit. lambda is the plain argmin of the mean held-out loss over
    # the grid. The fit applies no one-standard-error rule. Each vintage record
    # keeps the full loss grid, so you can apply another rule later without a
    # refit. The fit uses one fold partition, which cv_seed fixes.
    horizon_hours: float = 10.0
    fit_lambda_grid: Tuple[float, ...] = (0.0, 1e-2, 1e0, 1e2)
    fit_k_fold: int = 5             # folds, split over runs and never over samples

    # Acquisition, matching the log-domain form in the paper:
    #   log alpha = log qLogNEHVI + w_z log(a_z(z) + eps) - w_t gamma log(E[t] + eps)
    ell_z: float = 0.25
    gamma_time: float = 1.0
    w_z: float = 1.0
    w_t: float = 1.0
    eps: float = 1e-6
    ref_point_backoff: float = 0.2  # r_j = y_min - backoff * (y_max - y_min)

    # Acquisition maximiser.
    num_restarts: int = 20
    raw_samples: int = 1024
    acq_batch_limit: int = 5
    acq_maxiter: int = 250
    mc_samples: int = 128
    prune_baseline: bool = True
    max_iep: int = 256

    # Handshake with MATLAB.
    poll_s: float = 1.0
    matlab_wait_s: float = 2.0
    max_wait_matlab_s: float = 120.0
    # A single-fidelity evaluation runs the full 10 h horizon over two cases and
    # can take many hours. The cutoff is the wall time a driver waits for one row
    # before it gives the evaluation up. lock_stale_s must not sit below it: a
    # legitimately busy MATLAB holds matlab.lock for the whole evaluation, and a
    # shorter staleness would treat that lock as abandoned and unlink it while the
    # work is still running.
    eval_timeout_s: float = 10 * 3600.0
    lock_stale_s: float = 10 * 3600.0
    max_consecutive_failures: int = 5

    # Timeout imputation. An evaluation that exceeds eval_timeout_s is recorded as
    # a dominated point rather than skipped, so the optimiser does not propose it
    # again. The penalty is per objective, just beyond the worst measured so far,
    # so the imputed point is strictly dominated without introducing an outlier
    # that would distort the objective GP or the reference point:
    #     penalty_j = (1 + margin) * max_i observed_j
    # The fallback applies only before any objective has been measured, which can
    # happen if the very first design point times out.
    timeout_penalty_margin: float = 0.1
    timeout_penalty_fallback: float = 1e6

    def spec(self) -> CaseSpec:
        if self.case not in CASES:
            raise ValueError(
                f"unknown case {self.case!r}, choose one of {sorted(CASES)}")
        return CASES[self.case]

    @property
    def is_baseline(self) -> bool:
        """The single-fidelity, runtime-unaware variant.

        It fixes z at 1, so the fidelity surrogate scales by exactly one and the
        driver publishes a static phi rather than fitting one. The acquisition
        drops the runtime term and reduces to plain qLogNEHVI. These two
        behaviours follow from the case name alone, so no separate flag can drift
        out of step with it.
        """
        return self.case == "baseline"

    @property
    def runtime_aware(self) -> bool:
        return not self.is_baseline

    @property
    def refit_phi(self) -> bool:
        return not self.is_baseline

    @property
    def n_total(self) -> int:
        return self.n_init + self.n_iter

    def bounds(self) -> Tuple[List[float], List[float]]:
        """Lower and upper bounds over the full 12-dimensional theta."""
        lb = list(BASE_LB)
        ub = list(BASE_UB)
        for idx, val in self.spec().fixed.items():
            lb[idx] = val
            ub[idx] = val
        return lb, ub

    def to_dict(self) -> Dict:
        spec = self.spec()
        d = asdict(self)
        lb, ub = self.bounds()
        d.update({
            "z_min_bo": Z_MIN_BO,
            "z_min_phi": Z_MIN_PHI,
            "z_min_note": (
                "z_min_bo bounds the search space below. z_min_phi is the floor of "
                "the surrogate fit. They are separate limits on separate steps."
            ),
            "case_description": spec.description,
            "is_baseline": self.is_baseline,
            "runtime_aware": self.runtime_aware,
            "refit_phi": self.refit_phi,
            "fixed_components": {THETA_NAMES[i]: v for i, v in sorted(spec.fixed.items())},
            "optimised_components": [THETA_NAMES[i] for i in spec.opt_idxs],
            "search_dimension": spec.dimension,
            "theta_names": list(THETA_NAMES),
            "theta_dimension": THETA_D,
            "lower_bounds": lb,
            "upper_bounds": ub,
            "base_lower_bounds": list(BASE_LB),
            "base_upper_bounds": list(BASE_UB),
            "integer_components": [THETA_NAMES[i] for i in INTEGER_IDXS],
            "n_total_evaluations": self.n_total,
            "stopping_rule": (
                f"fixed budget: {self.n_init} initialisation evaluations followed by "
                f"{self.n_iter} acquisition-driven evaluations, declared before the "
                f"first evaluation"
            ),
            "refit_rule": (
                "the fidelity surrogate is not fitted: z is held at 1, so phi(1) = 1 "
                "scales every cost by exactly one. A static phi (a = b = 1, "
                "vintage 0) is published once and governs the whole run"
                if self.is_baseline else
                f"vintage 0 is fitted on the {self.n_init} initialisation runs and "
                f"governs optimisation iterations 1-{self.refit_every}; a refit "
                f"follows every {self.refit_every} optimisation iterations and "
                f"governs the next {self.refit_every}; past rows are never rescaled"
            ),
        })
        return d


def parse_args(argv: List[str]) -> Tuple[str, RunConfig]:
    """Read the phase and the case from the command line.

    Deliberately small: anything else that would change a result belongs in
    RunConfig, where the manifest picks it up, rather than in a flag that leaves
    no trace.
    """
    phase = "bo"
    case = "case1"

    rest = list(argv)
    if rest and not rest[0].startswith("-"):
        phase = rest.pop(0)

    while rest:
        token = rest.pop(0)
        if token == "--case":
            if not rest:
                raise ValueError("--case needs a value")
            case = rest.pop(0)
        elif token.startswith("--case="):
            case = token.split("=", 1)[1]
        else:
            raise ValueError(f"unknown argument {token!r}")

    if phase not in ("init", "bo"):
        raise ValueError(f"unknown phase {phase!r}, choose 'init' or 'bo'")

    cfg = RunConfig(case=case)
    cfg.spec()   # validates the case name
    return phase, cfg
