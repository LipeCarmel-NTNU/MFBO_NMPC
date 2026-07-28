# tangent/

Generic linearization-based design suite for a nonlinear system `xdot = f(x,u)` with optional output map `y = h(x)` (default `y = x`). It generalizes the workflow in `TerminalCost/` (steady state, linearization, incremental LQR, closed-loop weight tuning, Lyapunov terminal cost) so it can be reused on any model.

## Pipeline

| Step | Function | Purpose |
|---|---|---|
| 1 | `trim_point` | Given `yss`, solve `[f(x,u); h(x)-yss] = 0` for `(xss, uss)`. NaN entries mark free quantities: `yss(i) = NaN` leaves output i unconstrained, `uss(i) = NaN` leaves input i free (solved for); numeric `uss(i)` fixes it. Non-square systems handled by Levenberg–Marquardt. |
| 2 | `linearize_fd` | Numeric (central-difference) Jacobians `A = df/dx`, `B = df/du`, `C = dh/dx` at the trim point. No Symbolic Toolbox required. |
| 3 | `build_plant` | ZOH-discretize and package as a `plant` struct, `mode = 'positional'` or `'incremental'` (delta-u form with augmented state `z = [x - xss; u_{k-1} - uss]`). |
| 4 | `lqr_synth` | `dlqr` gain. Positional: `x'Qx + u'Ru`. Incremental: `x'Qx + u'R1u + du'R2du` expanded to the cross-term form `(Qz, R, N)`. |
| 5 | `tuning_cost` + `lqr_tune` | `tuning_cost(var, ...)` scores a tuning vector by nonlinear closed-loop performance over initial-condition scenarios (default map `log10_weights`: diagonal weights in log10, `Q(1,1) = 1` fixed for scale), with input-move and weight-magnitude regularization. `lqr_tune(objective_func, var0)` minimizes any `@(var) scalar` objective. |
| 6 | `terminal_cost` | Given the fixed gain `K` and *evaluation* weights (e.g. NMPC stage-cost weights), solve `P = dlyap(Acl', Qbar)`. Then `z'Pz` is the exact infinite-horizon cost of the LQR policy under those weights — an upper bound on the optimal cost-to-go, hence a valid NMPC terminal cost. Includes a truncated-sum validation check. |

Support: `lqr_sim` (nonlinear closed-loop simulation with ZOH input, saturation, positional or incremental law), `log10_weights` (tuning-vector-to-weights map).

## Minimal usage

```matlab
f = @(x,u) my_model(x, u);          % xdot, column vector
h = @(x) x(2);                      % optional; default y = x

% NaN = free: here u1 solved for, u2 fixed at 0, u3 solved for
[xss, uss] = trim_point(f, yss, x_guess, [NaN; 0; NaN], h = h);
[A, B]     = linearize_fd(f, xss, uss);
plant      = build_plant(A, B, Ts, 'incremental');

% Either hand-picked weights ...
K = lqr_synth(plant, struct('Q', Q, 'R1', R1, 'R2', R2));

% ... or tuned against the nonlinear closed loop
obj = @(var) tuning_cost(var, f, plant, xss, uss, {x0_a, x0_b}, ...
    Tf = 20, y_weight = [10 1 1], reg_du = 1e4, reg_var = 1e2, ...
    u_min = zeros(nu,1), u_max = 0.4*ones(nu,1));
[var_opt, J] = lqr_tune(obj, zeros(1, plant.nx-1 + 2*plant.nu));
K = lqr_synth(plant, log10_weights(var_opt, plant));

% Terminal cost under the NMPC stage-cost weights
[P, info] = terminal_cost(plant, K, struct('Q', Q_nmpc, 'R1', R1_nmpc, 'R2', R2_nmpc));
```

`demo_tangent.m` runs the full pipeline on a self-contained pendulum model.

## Conventions

- All vectors are treated as columns internally; `Y`, `U` trajectories are returned row-per-sample (as in `LQR_simulation`).
- Incremental augmentation matches `TerminalCost/LQR_functions/incremental.m`: `Ai = [Ad Bd; 0 I]`, `Bi = [Bd; I]`, control `du = -K z`, absolute input `u_k = u_{k-1} + du_k`.
- The synthesis weights (step 4/5) and the evaluation weights (step 6) are deliberately separate: `P` certifies the *evaluation* cost of the *synthesized* policy.
- Function names avoid collisions with `TerminalCost/LQR_functions/` so both can be on the MATLAB path.

## Requirements

MATLAB R2019b+ (`arguments` blocks; R2021a+ for the `name = value` call syntax — older releases can use `'name', value`), Control System Toolbox (`c2d`, `dlqr`, `dlyap`), Optimization Toolbox (`fsolve`, `fminunc`; `fminsearch` works without it via `solver = 'fminsearch'`).
