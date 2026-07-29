"""Declared configuration of one optimisation run.

Every setting that changes what the run does is named here and copied verbatim
into the manifest before the first evaluation. Nothing that affects a result is
left to a hand edit of the driver: the case is selected by name, the budget is
a number rather than a point at which somebody stopped the process, and the
seeds are explicit.

Select a case on the command line:

    python main.py init --case case2
    python main.py bo   --case case2
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Dict, List, Tuple

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

BASE_LB: Tuple[float, ...] = (0.01, 0.0, 0.0) + (-3.0,) * 9
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
}


@dataclass(frozen=True)
class RunConfig:
    """One optimisation run, start to finish."""

    case: str = "case1"

    # Budget. n_init evaluations of the Sobol design, then n_iter
    # acquisition-driven evaluations, for n_init + n_iter in total. The loop
    # ends when n_iter is reached; no other criterion stops it.
    n_init: int = 20
    n_iter: int = 100
    q_batch: int = 1

    # Surrogate refit. Vintage 0 is fitted on the initialisation runs alone and
    # governs optimisation iterations 1 to refit_every. Vintage v is fitted once
    # iteration v * refit_every completes and governs the next refit_every
    # iterations. Past evaluations are never rescaled: a row keeps the estimate
    # produced by the vintage in force when it was measured.
    refit_every: int = 10
    refit_after_last: bool = True   # fit a terminal vintage for post-hoc analysis

    # Seeds.
    sobol_seed: int = 1234          # scrambled Sobol stream of the DOE
    torch_seed: int = 0             # GP fitting and acquisition multistarts
    cv_seed: int = 1                # fold partitions of the surrogate fit

    # Surrogate fit. lambda is the plain argmin of the mean held-out loss over
    # the grid; no one-standard-error rule is applied. The full loss grid is
    # kept in each vintage record, so another rule can be applied afterwards
    # without refitting.
    horizon_hours: float = 10.0
    fit_z_min: float = 0.1          # samples below this fidelity leave the fit
    fit_lambda_grid: Tuple[float, ...] = (0.0, 1e-2, 1e0, 1e2)
    fit_k_fold: int = 5             # folds, split over runs and never over samples
    fit_cv_repeats: int = 5         # fold partitions redrawn per fit

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
    eval_timeout_s: float = 6 * 3600.0
    lock_stale_s: float = 6 * 3600.0
    max_consecutive_failures: int = 5

    def spec(self) -> CaseSpec:
        if self.case not in CASES:
            raise ValueError(
                f"unknown case {self.case!r}, choose one of {sorted(CASES)}")
        return CASES[self.case]

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
            "case_description": spec.description,
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
