"""Fit the runtime cost surrogate f(z) from NMPC simulation output files.

The surrogate maps the fidelity z = t / T_full to the fraction of the
full-horizon cost that the simulation has accumulated at time t. One out_*.mat
file holds one simulation pair. Each pair gives one curve of cost ratio
against z. The curves of many pairs enter one fit.

The model is a Chebyshev polynomial of order 4 in x = 2 z - 1. The fit is ridge
regression under two hard endpoint constraints, f(0) = 0 and f(1) = 1. The
ridge weight lambda comes from 10-fold cross-validation with 5 repeats. Folds
group by simulation pair, so all points of one pair leave the fit together.
The data term of the loss is a mean, not a sum, so one lambda value has the
same effect on a 20-pair fit and on a 120-pair fit.

Import the module and send a list of paths:

    from fit_runtime_surrogate import fit_from_paths
    result = fit_from_paths(paths, "J_TV")
    print(result.coefficients)

Run the module as a script to fit the first 20 files of results/run1. Those 20
files are the DOE phase of that run. Give a different run directory as the
first command-line argument.
"""

from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np
import scipy.io

FULL_HORIZON_HOURS = 10.0
Z_STEP = 0.01
CHEB_ORDER = 4
LAMBDA_GRID = np.concatenate(([0.0], np.logspace(-6.0, 0.0, 13)))
N_FOLDS = 10
N_CV_REPEATS = 5
RNG_SEED = 1
N_DOE_FILES = 20
Z_TOLERANCE = 1e-12


@dataclass(frozen=True)
class TargetSpec:
    """Names of the .mat fields that hold one cost target.

    partial_field is the per-step cost of one case. total_field is the
    full-horizon cost of one case. first_time_index is the index in out.T that
    belongs to the first entry of partial_field. The input-change cost needs
    two samples, so its first entry belongs to out.T[1], not to out.T[0].
    """

    partial_field: str
    total_field: str
    first_time_index: int


TARGETS = {
    "J_TV": TargetSpec("partial_SSdU", "SSdU", 1),
    "J_track": TargetSpec("partial_SSE", "SSE", 0),
}


@dataclass
class Pair:
    """One simulation pair, resampled on the shared fidelity grid.

    file_name is the source file. z holds the grid points that the simulation
    reached. ratio maps each target name to the cost ratio at those points.
    z_end is the last fidelity that the simulation reached.
    """

    file_name: str
    z: np.ndarray
    ratio: dict
    z_end: float


@dataclass
class FitResult:
    """One fitted surrogate and the numbers that produced it.

    coefficients holds c_0 to c_4 of the Chebyshev sum. lambda_used is the
    selected ridge weight. cv_loss holds the mean held-out squared error of
    every lambda in lambda_grid, in the same order. lambda_index is the
    position of lambda_used in lambda_grid. on_grid_edge is True when that
    position is the first or the last one. A True value at the last position
    means the grid is too narrow, so lambda was limited and not selected.
    n_pairs and n_points report the size of the fitted data.
    """

    target: str
    coefficients: np.ndarray
    lambda_used: float
    lambda_grid: np.ndarray
    cv_loss: np.ndarray
    lambda_index: int
    on_grid_edge: bool
    n_pairs: int
    n_points: int


def list_out_files(run_dir):
    """Return the out_*.mat paths of one run directory, sorted by name.

    The file names carry a timestamp, so the sorted order is the order of the
    optimization. The first files are therefore the DOE files.

    Raises FileNotFoundError when the directory is missing or holds no such
    file.
    """
    run_dir = Path(run_dir)
    if not run_dir.is_dir():
        raise FileNotFoundError(f"run directory does not exist: {run_dir}")
    paths = sorted(run_dir.glob("out_*.mat"))
    if not paths:
        raise FileNotFoundError(f"no out_*.mat files in {run_dir}")
    return paths


def read_out_struct(path):
    """Read one .mat file and return its out variable.

    Raises ValueError when the file holds no variable named out.
    """
    contents = scipy.io.loadmat(path, struct_as_record=False, squeeze_me=True)
    if "out" not in contents:
        raise ValueError(f"{Path(path).name} holds no variable named out")
    return contents["out"]


def fidelity_grid(step=Z_STEP):
    """Return the shared fidelity grid, from 0 to 1 with the given step."""
    n_steps = int(round(1.0 / step))
    return np.linspace(0.0, 1.0, n_steps + 1)


def interpolate_linear(x_source, y_source, x_query):
    """Interpolate y_source at x_query, and extrapolate outside the source.

    Inside the source range the function interpolates between neighbors.
    Outside the range it continues the straight line through the two nearest
    source points. x_source must increase.

    numpy.interp alone holds the end value outside the range. The two extra
    branches replace that flat behavior with a straight line, which matches
    interp1(..., "linear", "extrap") in MATLAB.
    """
    x_source = np.asarray(x_source, dtype=float)
    y_source = np.asarray(y_source, dtype=float)
    x_query = np.asarray(x_query, dtype=float)
    y_query = np.interp(x_query, x_source, y_source)

    below = x_query < x_source[0]
    if np.any(below):
        slope = (y_source[1] - y_source[0]) / (x_source[1] - x_source[0])
        y_query[below] = y_source[0] + slope * (x_query[below] - x_source[0])

    above = x_query > x_source[-1]
    if np.any(above):
        slope = (y_source[-1] - y_source[-2]) / (x_source[-1] - x_source[-2])
        y_query[above] = y_source[-1] + slope * (x_query[above] - x_source[-1])

    return y_query


def target_curve(out, target, horizon_hours=FULL_HORIZON_HOURS):
    """Build the raw cost ratio curve of one target from one out struct.

    The function adds the per-step cost of both cases, accumulates it, and
    divides by the sum of the two full-horizon costs. It returns the fidelity
    of every sample and the ratio at that sample. The ratio reaches 1 only
    when the simulation ran the full horizon.

    Raises ValueError when the struct has fewer than two cases, when the
    full-horizon cost is not a positive finite number, or when fewer than
    three samples remain.
    """
    spec = TARGETS[target]
    cases = np.atleast_1d(out.case)
    if cases.size < 2:
        raise ValueError(f"out.case holds {cases.size} case, two are needed")

    time_hours = np.asarray(out.T, dtype=float).ravel()
    partial_1 = np.asarray(getattr(cases[0], spec.partial_field), dtype=float).ravel()
    partial_2 = np.asarray(getattr(cases[1], spec.partial_field), dtype=float).ravel()

    total = float(np.asarray(getattr(out, spec.total_field)).ravel()[0])
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError(f"out.{spec.total_field} is {total}, a positive value is needed")

    start = spec.first_time_index
    n_samples = min(partial_1.size, partial_2.size, time_hours.size - start)
    if n_samples < 3:
        raise ValueError(f"{target} has {n_samples} samples, three are needed")

    cumulative = np.cumsum(partial_1[:n_samples]) + np.cumsum(partial_2[:n_samples])
    z_source = time_hours[start:start + n_samples] / horizon_hours
    return z_source, cumulative / total


def load_pair(path, grid, horizon_hours=FULL_HORIZON_HOURS):
    """Load one .mat file and resample every target on the shared grid.

    The pair keeps only the grid points up to the fidelity that both targets
    reached. All returned values are finite.

    Raises ValueError when the file cannot give at least three usable grid
    points. The message names the reason.
    """
    out = read_out_struct(path)

    z_source = {}
    ratio_source = {}
    for target in TARGETS:
        z_source[target], ratio_source[target] = target_curve(out, target, horizon_hours)

    z_end = min(float(z_source[target][-1]) for target in TARGETS)
    if not np.isfinite(z_end) or z_end <= grid[1]:
        raise ValueError(f"the run stops at z = {z_end:.4f}, which is at or below the first grid step")

    z = grid[grid <= z_end + Z_TOLERANCE]
    ratio = {}
    for target in TARGETS:
        ratio[target] = interpolate_linear(z_source[target], ratio_source[target], z)

    finite = np.ones(z.shape, dtype=bool)
    for target in TARGETS:
        finite &= np.isfinite(ratio[target])
    z = z[finite]
    for target in TARGETS:
        ratio[target] = ratio[target][finite]

    if z.size < 3:
        raise ValueError(f"{z.size} finite grid points remain, three are needed")

    return Pair(file_name=Path(path).name, z=z, ratio=ratio, z_end=z_end)


def load_pairs(paths, grid=None, horizon_hours=FULL_HORIZON_HOURS, verbose=True):
    """Load every path and return the usable pairs.

    The function skips a file that load_pair rejects. With verbose set it
    prints the file name and the reason for every skip, so a short return list
    is easy to explain.
    """
    if grid is None:
        grid = fidelity_grid()
    pairs = []
    for path in paths:
        try:
            pairs.append(load_pair(path, grid, horizon_hours))
        except (ValueError, AttributeError, KeyError) as reason:
            if verbose:
                print(f"skipped {Path(path).name}: {reason}")
    return pairs


def select_groups(pairs, target):
    """Return one (z, y) tuple per pair for the given target.

    A group is the unit that cross-validation holds out. All points of one
    pair stay in the same group.

    Raises KeyError when the target name is unknown.
    """
    if target not in TARGETS:
        raise KeyError(f"unknown target {target!r}, choose from {sorted(TARGETS)}")
    return [(pair.z, pair.ratio[target]) for pair in pairs]


def cheb_features(z, order=CHEB_ORDER):
    """Return the Chebyshev design matrix of the given order.

    Column k holds T_k(x) with x = 2 z - 1, so column 0 is all ones. The
    function uses the recurrence T_0 = 1, T_1 = x, T_k+1 = 2 x T_k - T_k-1.
    The matrix has one row per element of z and order + 1 columns.
    """
    z = np.asarray(z, dtype=float).ravel()
    x = 2.0 * z - 1.0
    features = np.zeros((z.size, order + 1))
    features[:, 0] = 1.0
    if order >= 1:
        features[:, 1] = x
    for k in range(2, order + 1):
        features[:, k] = 2.0 * x * features[:, k - 1] - features[:, k - 2]
    return features


def endpoint_constraints(order=CHEB_ORDER):
    """Return the matrix A and the vector b of the endpoint constraints.

    The two constraints are f(0) = 0 and f(1) = 1. Both are linear in the
    coefficients, so A c = b states them exactly. Row 0 evaluates the
    Chebyshev basis at z = 0 and row 1 evaluates it at z = 1.
    """
    A = np.vstack((cheb_features([0.0], order), cheb_features([1.0], order)))
    b = np.array([0.0, 1.0])
    return A, b


def stack_groups(groups, order=CHEB_ORDER):
    """Join the groups into one design matrix and one target vector.

    Raises ValueError when the list is empty.
    """
    if not groups:
        raise ValueError("no groups to stack")
    z = np.concatenate([np.asarray(z, dtype=float).ravel() for z, _ in groups])
    y = np.concatenate([np.asarray(y, dtype=float).ravel() for _, y in groups])
    return cheb_features(z, order), y


def solve_constrained_ridge(X, y, ridge_weight, order=CHEB_ORDER):
    """Return the coefficients of the ridge fit under the endpoint constraints.

    The problem is

        minimize    mean((y - X c) ** 2) + ridge_weight * sum(c ** 2)
        subject to  f(0) = 0 and f(1) = 1.

    The function solves the KKT system

        [2 (X'X / n + w I)   A'] [c ]   [2 X'y / n]
        [A                    0] [mu] = [b        ]

    with n the number of rows of X and w the ridge weight. It returns c and
    drops the multipliers mu.

    A numpy.linalg.LinAlgError means the system is singular. The usual cause is
    a training set with fewer distinct z values than coefficients.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    n_rows = X.shape[0]
    n_coefficients = order + 1

    A, b = endpoint_constraints(order)
    hessian = 2.0 * (X.T @ X / n_rows + ridge_weight * np.eye(n_coefficients))

    system = np.zeros((n_coefficients + 2, n_coefficients + 2))
    system[:n_coefficients, :n_coefficients] = hessian
    system[:n_coefficients, n_coefficients:] = A.T
    system[n_coefficients:, :n_coefficients] = A

    right_side = np.concatenate((2.0 * (X.T @ y) / n_rows, b))
    solution = np.linalg.solve(system, right_side)
    return solution[:n_coefficients]


def mean_squared_error(X, y, coefficients):
    """Return the mean squared residual of the given coefficients."""
    residual = np.asarray(y, dtype=float).ravel() - np.asarray(X, dtype=float) @ coefficients
    return float(np.mean(residual ** 2))


def make_fold_labels(n_groups, n_folds, rng):
    """Assign every group to one fold, in random order and in equal shares.

    The function returns one integer label per group, from 0 to n_folds - 1.
    The first folds take one extra group when the division leaves a remainder.
    n_folds drops to n_groups when there are fewer groups than folds.
    """
    n_folds = min(n_folds, n_groups)
    order = rng.permutation(n_groups)
    labels = np.zeros(n_groups, dtype=int)
    base_size = n_groups // n_folds
    remainder = n_groups - base_size * n_folds
    cursor = 0
    for fold in range(n_folds):
        size = base_size + (1 if fold < remainder else 0)
        labels[order[cursor:cursor + size]] = fold
        cursor += size
    return labels


def make_fold_plan(n_groups, n_folds=N_FOLDS, n_repeats=N_CV_REPEATS, seed=RNG_SEED):
    """Return one fold assignment per repeat.

    Every repeat is a complete and independent K-fold split of the same
    groups. Averaging over repeats lowers the effect of one lucky or unlucky
    split on the selected lambda. The seed fixes the whole plan, so two runs
    of the same data select the same lambda.

    All lambda values later share this plan, so the comparison between them
    uses the same splits.
    """
    rng = np.random.default_rng(seed)
    return [make_fold_labels(n_groups, n_folds, rng) for _ in range(n_repeats)]


def split_groups(groups, labels, fold):
    """Split the groups into the training part and the held-out part.

    The held-out part holds the groups whose label equals fold. The training
    part holds the rest. No group appears in both parts.
    """
    train = [group for group, label in zip(groups, labels) if label != fold]
    held_out = [group for group, label in zip(groups, labels) if label == fold]
    return train, held_out


def fold_loss(train_groups, held_out_groups, ridge_weight, order=CHEB_ORDER):
    """Fit on the training groups and return the error on the held-out groups.

    The returned number is a mean squared error over the held-out points. The
    fit never sees those points.

    Raises ValueError when either part is empty.
    """
    X_train, y_train = stack_groups(train_groups, order)
    X_held_out, y_held_out = stack_groups(held_out_groups, order)
    coefficients = solve_constrained_ridge(X_train, y_train, ridge_weight, order)
    return mean_squared_error(X_held_out, y_held_out, coefficients)


def cross_validation_loss(groups, plan, ridge_weight, order=CHEB_ORDER):
    """Return the mean held-out error of one lambda over the whole plan.

    The function walks every fold of every repeat, collects one held-out mean
    squared error each time, and averages them.
    """
    losses = []
    for labels in plan:
        for fold in np.unique(labels):
            train, held_out = split_groups(groups, labels, fold)
            losses.append(fold_loss(train, held_out, ridge_weight, order))
    return float(np.mean(losses))


def cross_validation_curve(groups, plan, lambda_grid=LAMBDA_GRID, order=CHEB_ORDER):
    """Return the mean held-out error of every lambda in the grid.

    The output has the same length and the same order as lambda_grid.
    """
    return np.array([
        cross_validation_loss(groups, plan, ridge_weight, order)
        for ridge_weight in lambda_grid
    ])


def pick_lambda(lambda_grid, cv_loss):
    """Return the lambda with the lowest mean held-out error, and its index.

    The first index wins a tie, so the rule prefers less shrinkage.
    """
    index = int(np.argmin(cv_loss))
    return float(lambda_grid[index]), index


def fit_surrogate(groups, target, lambda_grid=LAMBDA_GRID, order=CHEB_ORDER,
                  n_folds=N_FOLDS, n_repeats=N_CV_REPEATS, seed=RNG_SEED):
    """Select lambda by repeated cross-validation, then fit on all groups.

    The steps are: build the fold plan, score every lambda on that plan, take
    the lambda with the lowest error, and refit once on every group.

    Raises ValueError when there are fewer than two groups, because
    cross-validation then has nothing to hold out.
    """
    if len(groups) < 2:
        raise ValueError(f"{len(groups)} group given, two or more are needed")

    plan = make_fold_plan(len(groups), n_folds, n_repeats, seed)
    cv_loss = cross_validation_curve(groups, plan, lambda_grid, order)
    lambda_used, lambda_index = pick_lambda(lambda_grid, cv_loss)

    X, y = stack_groups(groups, order)
    coefficients = solve_constrained_ridge(X, y, lambda_used, order)

    # The last grid point means the search stopped at the edge of the grid.
    # The first grid point is lambda = 0, which is a real outcome.
    on_grid_edge = lambda_index in (0, len(lambda_grid) - 1)

    return FitResult(
        target=target,
        coefficients=coefficients,
        lambda_used=lambda_used,
        lambda_grid=np.asarray(lambda_grid, dtype=float),
        cv_loss=cv_loss,
        lambda_index=lambda_index,
        on_grid_edge=on_grid_edge,
        n_pairs=len(groups),
        n_points=y.size,
    )


def evaluate_surrogate(z, coefficients):
    """Return the surrogate value f(z) for the given coefficients."""
    order = len(coefficients) - 1
    return cheb_features(z, order) @ np.asarray(coefficients, dtype=float)


def fit_from_paths(paths, target, horizon_hours=FULL_HORIZON_HOURS,
                   lambda_grid=LAMBDA_GRID, order=CHEB_ORDER, n_folds=N_FOLDS,
                   n_repeats=N_CV_REPEATS, seed=RNG_SEED, verbose=True):
    """Load the files, take one target, and fit the surrogate.

    This is the entry point for another module. It runs load_pairs,
    select_groups, and fit_surrogate in that order. Call those three in turn
    when you need the pairs or the groups as well.
    """
    pairs = load_pairs(paths, fidelity_grid(), horizon_hours, verbose)
    groups = select_groups(pairs, target)
    return fit_surrogate(groups, target, lambda_grid, order, n_folds, n_repeats, seed)


def print_fit(result):
    """Print one fit in a fixed layout."""
    print(f"target            {result.target}")
    print(f"pairs, points     {result.n_pairs}, {result.n_points}")
    print(f"lambda            {result.lambda_used:.4e} (grid index {result.lambda_index})")
    print(f"lambda on edge    {result.on_grid_edge}")
    print(f"held-out mse      {result.cv_loss[result.lambda_index]:.6e}")
    print("coefficients      " + " ".join(f"{c: .6f}" for c in result.coefficients))
    z_report = np.array([0.2, 0.3, 0.4, 0.5])
    values = evaluate_surrogate(z_report, result.coefficients)
    for z, value in zip(z_report, values):
        print(f"f(z = {z:.2f})        {value: .6f}")


def main(run_dir=None):
    """Fit both targets on the first 20 files of one run directory.

    Those 20 files are the DOE phase of the run. The default directory is
    results/run1 of this project.
    """
    if run_dir is None:
        project_root = Path(__file__).resolve().parents[2]
        run_dir = project_root / "results" / "run1"

    paths = list_out_files(run_dir)[:N_DOE_FILES]
    print(f"run directory     {run_dir}")
    print(f"files             {len(paths)}")
    print(f"cross-validation  {N_FOLDS} folds, {N_CV_REPEATS} repeats, seed {RNG_SEED}")

    for target in TARGETS:
        print("")
        print_fit(fit_from_paths(paths, target))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
