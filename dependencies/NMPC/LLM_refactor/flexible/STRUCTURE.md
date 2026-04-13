# NMPC Refactor — Flexible Version

Composition-based NMPC. One concrete `NMPC` orchestrator built from interchangeable pieces. No inheritance tree; contracts are duck-typed (MATLAB `ismethod` dispatch). The high-level rationale is in `LLM_refactor/STRUCTURE.md`; this file is the implementation map.

## Folder layout

```
flexible/
├── STRUCTURE.md
├── NMPC.m                       Orchestrator. solve(), fallback, logging hook.
├── NMPCBuilder.m                Fluent assembler producing a validated NMPC.
├── DecisionLayout.m             Maps w ⇄ named blocks (x, u, slacks, aux).
├── Scaling.m                    Per-variable scale factors (default identity).
├── HistoryLogger.m              Optional per-call record store.
│
├── model/
│   ├── ContinuousModel.m        xdot = f(x,u); composes an Integrator.
│   └── DiscreteModel.m          x_next = f(x,u); no integrator.
├── integrator/
│   └── RK4.m                    Fixed-step RK4.
│
├── constraints/
│   ├── StateBounds.m            Hard or per-state-soft bounds. Registers slacks.
│   ├── InputBounds.m            Hard control bounds.
│   ├── InputMovement.m          |Δu| ≤ dumax. Per-call linear rows (need u_prev).
│   └── Continuity.m             Auto-added; nonlinear equality continuity.
│
├── cost/
│   ├── QuadraticTracking.m      Σ ‖y - y_sp‖_Q² + Σ ‖u - u_sp‖_R².
│   ├── TerminalCost.m           ‖x_N - x_sp‖_P².
│   ├── InputMovementCost.m      Σ ‖Δu‖_S².
│   ├── SlackPenalty.m           L1 (exact) + optional L2 on registered slacks.
│   └── EconomicCost.m           Σ w · ϕ(x, u). Forward-compat hook.
│
├── solver/
│   └── FminconBackend.m         Wraps fmincon. Single problem-struct interface.
│
└── tests/
    ├── build_basic_nmpc.m       Shared factory (mirror of minimal/tests).
    ├── test_smoke.m             End-to-end smoke checks.
    └── runtests_flexible.m      Runner: adds paths, runs all tests.
```

## Contracts

### Model (duck-typed)

```
m.nx                 % # states
m.nu                 % # inputs
x_next = m.step(x, u)   % one Ts ahead (continuous: integrator inside)
```

### Constraint (any subset of these methods)

```
register(c, layout)                                 % add slack/aux blocks before freeze
[A, b]      = linear_static(c, layout)              % constant linear rows in scaled w
[A, b]      = linear_per_call(c, layout, ctx)       % rebuilt each solve (e.g. Δu)
[c_, ceq_]  = nonlinear(c, parts, ctx)              % nonlinear part. parts = unpacked named blocks (physical)
overrides   = bound_overrides(c, layout)            % struct: block name → struct(lo, hi)
```

`ctx` is a struct: `{x_init, u_prev, model, scaling, layout, p, m, nx, nu, ny}`. Constraints only implement what they need; the orchestrator dispatches via `ismethod`.

### Cost (any subset)

```
register(c, layout)
J = evaluate(c, parts, ctx)
```

`parts` is a struct of physical arrays (`parts.x`, `parts.u`, `parts.slack_<name>`, …). Costs return a scalar; the orchestrator sums them.

### Solver backend

```
problem = struct('objfun', @, 'w0', ws0, 'wL', wL, 'wU', wU, ...
                 'A', A, 'b', b, 'Aeq', Aeq, 'beq', beq, ...
                 'nonlcon', @, 'options', opts);
[w_opt, fval, flag] = backend.solve(problem)
```

Only `FminconBackend` ships day-1; the seam exists for IPOPT/CasADi.

## Layout

`DecisionLayout` owns the `w_scaled` vector structure. Built incrementally:

```matlab
layout = DecisionLayout();
layout.add('x', p+1, nx);
layout.add('u', m,   nu);
% Constraints register their own blocks during build:
layout.add('slack_state', p+1, n_soft);   % from StateBounds(soft mode)
layout.freeze();
```

After freeze, `layout.pack(parts_phys, scaling)` and `layout.unpack(w_scaled, scaling)` are the only conversion API; cost/constraint code always sees physical values.

## Default behaviour

The Builder requires a model, horizon, bounds, and at least one cost. Everything else (scaling, terminal cost, slacks, Δu, economic, logging) is opt-in. Minimal opt-ins reproduce the legacy unscaled tracking controller.

## Forward-looking example

New model + economic term + state slacks:

```matlab
nmpc = NMPCBuilder() ...
    .withModel(ContinuousModel(@f_new, nx, nu, RK4(Ts))) ...
    .withHorizon(p, m) ...
    .addConstraint(StateBounds(Ymin, Ymax, soft=soft_mask)) ...
    .addConstraint(InputBounds(umin, umax)) ...
    .addCost(QuadraticTracking(Q, R, y_sp, u_sp)) ...
    .addCost(EconomicCost(@phi_econ, weight=w_econ)) ...
    .addCost(SlackPenalty(rho_L1=rho)) ...
    .build();
```

No new class; `EconomicCost` is already in the box.
