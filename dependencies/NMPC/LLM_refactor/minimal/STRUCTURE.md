# NMPC Refactor — Minimal Version

A single-class, single-file NMPC. Optimizes for *reading the controller end-to-end in one sitting*. Everything optional is gated by a small set of properties with safe defaults; behavior is selected by `if` checks, not by polymorphism or composition.

## Fixed assumptions (not configurable)

- **Continuous model** `xdot = f(x, u)`, supplied as a function handle.
- **Multiple shooting** with one shooting node per sampling interval. Continuity is enforced as nonlinear equality constraints (already the structure used by `NMPC_abstract.confun`).
- **Scaling is always present** in the decision vector. Default scale = 1 ⇒ identity, but the code path is unconditional (no `if scaled`).
- **fmincon** as the solver.
- **RK4 / RKF45** as the integrator (one routine, fixed).
- **Single tracking + terminal cost form**: `Σ ‖y - y_sp‖_Q² + Σ ‖u - u_sp‖_R² + ‖x_N - x_sp‖_P²`. `P = 0` ⇒ no terminal contribution; same code path.

## Optional features (toggled by data, handled inline with `if`)

| Feature              | Toggle                            | Default     |
| -------------------- | --------------------------------- | ----------- |
| Terminal cost        | `P` non-empty / non-zero          | `P = []`    |
| State scaling        | `x_scale ≠ 1` (any entry)         | `ones(1,nx)`|
| Input scaling        | `u_scale ≠ 1` (any entry)         | `ones(1,nu)`|
| Soft state bounds    | `soft_mask` non-empty             | `[]`        |
| L1 slack penalty     | `rho_L1` (used iff slacks active) | `0`         |
| L2 slack penalty     | `rho_L2` (used iff slacks active) | `0`         |
| Δu penalty / bound   | `S` / `dumax` non-empty           | `[]`        |
| History logging      | `log_enabled` flag                | `false`     |

The class never branches on model type, solver type, or shooting type — those are fixed.

## File layout

```
minimal/
├── STRUCTURE.md     (this file)
├── NMPC.m           Single class. ~300–400 lines incl. comments.
└── tests/
    ├── test_parity_mono.m       Closed-loop vs NMPC_mono_dilution.
    ├── test_parity_terminal.m   Closed-loop vs NMPC_terminal.
    └── test_slacks.m            Soft state bounds, L1 exact-penalty.
```

That is the entire module. No subfolders for `model/`, `cost/`, `constraints/`, etc.

## Class skeleton

```matlab
classdef NMPC < handle
    properties
        % ----- model & horizon (required) -----
        f                         % @(x,u) -> xdot
        nx, nu, ny
        Ts, p, m

        % ----- tracking targets & weights (required) -----
        y_sp, u_sp                % setpoints (vectors or per-step matrices)
        Q, R                      % weights
        h_y                       % @(x) -> y, default identity (eye(ny,nx))

        % ----- bounds (required) -----
        Ymin, Ymax
        umin, umax

        % ----- optional terminal cost -----
        P    = []                 % empty ⇒ skipped
        x_sp = []

        % ----- optional scaling (default = 1) -----
        x_scale = []              % filled to ones(1,nx) on init
        u_scale = []              % filled to ones(1,nu) on init

        % ----- optional soft state bounds -----
        soft_mask = []            % logical 1×nx; [] ⇒ no slacks
        rho_L1    = 0
        rho_L2    = 0

        % ----- optional input movement -----
        S      = []
        dumax  = []

        % ----- optional history log -----
        log_enabled = false
        log         = struct([])

        % ----- solver -----
        optimizer_options

        % ----- runtime state (fallback machinery) -----
        latest_wopt = []
        latest_flag = NaN
    end

    methods
        function obj = NMPC(args)   % name-value constructor; fills defaults

        function [uk, x, u, info] = solve(obj, x_init, u_init)
            % Identical fallback flow to NMPC_abstract.solve:
            %   1. cold solve if no warm start or previous failure
            %   2. else moving-horizon warm start
            %   3. flag == -2 → continuity-restoring retry
            %   4. flag  < 0 → fallback to u(k+1|k-1) or zeros
            % Logs at the end if log_enabled.
        end

        % ----- internals (all private-by-convention) -----
        function J         = objfun(obj, w, u_prev)
        function [c, ceq]  = confun(obj, w, x_init)
        function [w0]      = guess_from_initial(obj, x_init, u_init)
        function [x, u, s] = unpack(obj, w)         % s = [] when no slacks
        function w         = pack(obj, x, u, s)
        function [A, b, wL, wU] = build_linear(obj) % bounds + soft-bound rows
        function x_next    = step(obj, x, u)        % RK4 inside
        function record    = make_log_record(obj, ...)
    end
end
```

## How the `if` checks are organized

There are exactly four optional branches, all in obvious places:

1. **`pack` / `unpack` / `build_linear`** — `if ~isempty(obj.soft_mask)` appends a slack block to `w`, adds linear rows `x − s ≤ Ymax`, `−x − s ≤ −Ymin` for the soft entries, and sets `s ≥ 0` in `wL`. Hard rows for non-soft entries stay in `wL`/`wU` as before.
2. **`objfun`** — terminal: `if ~isempty(obj.P), J = J + (xN-x_sp)'*P*(xN-x_sp); end`. Slack penalty: `if ~isempty(s), J = J + obj.rho_L1*sum(s(:)) + s(:)'*diag(obj.rho_L2)*s(:); end`. Δu: `if ~isempty(obj.S), J = J + sum_k norm_S(du_k); end`.
3. **Scaling** — applied *unconditionally* in `pack`/`unpack` (`x_phys = x_dec .* x_scale`). With default `x_scale = ones`, this is a no-op multiplication; no branch needed.
4. **Logging** — `if obj.log_enabled, obj.log(end+1) = make_log_record(...); end` at the tail of `solve`.

Everything else (continuity, tracking cost, bounds, fallback) is unconditional.

## Decision vector layout (always the same shape, slack block optional)

```
w = [ x(1,:)  x(2,:)  …  x(p+1,:)        % (p+1)·nx, scaled
      u(1,:)  u(2,:)  …  u(m,:)          %    m  ·nu, scaled
      s(1,:)  s(2,:)  …  s(p+1,:) ]      % (p+1)·n_soft, present iff soft_mask non-empty
```

`unpack` returns `s = []` when no soft mask is configured, so call sites can use `if ~isempty(s)`.

## Fallback mechanism

Copied verbatim in spirit from `NMPC_abstract.m:78-168`:

- Try a warm-started solve from the shifted previous solution.
- If `exitflag == -2`, rebuild `w0` to satisfy continuity from current `x_init` and retry.
- If still negative, return `u(k+1|k-1)` (or zeros on first step), warn.
- Persist `latest_wopt` only when `exitflag >= 0`.

## History log schema

Each `solve` appends one struct with these fields when `log_enabled = true`:

```
t_wall, x_init, u_init,
y_sp_k, u_sp_k, x_sp_k,
Q, R, P, rho_L1, rho_L2, S,
w0, w_opt, uk,
flag, fval, iters, first_order_opt
```

Sufficient to replay any single step by feeding `w0`, bounds, weights back into `fmincon` with the same options.

## Trade-offs vs the flexible version

| Aspect                | Minimal                                     | Flexible (sibling folder)                    |
| --------------------- | ------------------------------------------- | -------------------------------------------- |
| Files to read         | 1                                           | ~20 across model/cost/constraints/solver     |
| Adding a new model    | Replace `f` handle. Must be continuous.     | New `Model` class, plug in.                  |
| Discrete model        | Not supported (would require editing class) | `DiscreteModel` drops in.                    |
| New cost term         | Edit `objfun`                               | Add a `Cost` object via builder              |
| Different solver      | Edit `solve_optimization`                   | Implement `SolverBackend`                    |
| Custom constraint     | Edit `confun` or `build_linear`             | Add a `Constraint` object                    |
| Test surface          | Small, monolithic                           | Per-component unit tests possible            |
| Best when             | Scope is fixed (current + next model)       | Many model/cost variants are anticipated     |

## Forward compatibility note

The "next model" mentioned (same cost + economic term + state slacks) fits the minimal class with **two edits**:

1. Replace `f`.
2. In `objfun`, add `J = J + obj.w_econ * sum(phi_econ(x, u))` behind `if ~isempty(obj.w_econ)`, plus a property `w_econ` and `phi_econ` handle.

State slacks are already supported via `soft_mask` + `rho_L1`. So the upcoming change is roughly 10 lines, not a refactor.

## Migration plan

1. Land `minimal/NMPC.m`.
2. `test_parity_mono.m` and `test_parity_terminal.m` must match the legacy classes within tight tolerance.
3. Add `test_slacks.m` for the new soft-bound feature.
4. Switch downstream scripts when ready; legacy classes remain untouched until parity is signed off.
