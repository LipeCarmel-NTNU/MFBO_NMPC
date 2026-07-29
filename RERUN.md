# Rerun procedure

Two processes talk through files: MATLAB simulates, Python proposes. Start
MATLAB first in both phases, because the driver waits for the results file that
MATLAB creates.

## Sequence

Phase 1, the design of experiments. 20 Sobol points, each simulated at the
fidelity its own `z` asks for, with no surrogate applied.

```
MATLAB> main_initialization
shell>  python main.py init --case case1
```

The driver fits and publishes surrogate vintage 0 when the design finishes, then
stops.

Phase 2, the optimisation. 100 acquisition-driven evaluations.

```
MATLAB> main_BO
shell>  python main.py bo --case case1
```

Repeat both phases with `--case case2` in a separate checkout or after moving
`results/` aside. The two cases write the same paths.

## Surrogate schedule

`f(z) = I_z(a, b)`, the regularised incomplete beta function, replaces the
Chebyshev polynomial. It satisfies `f(0) = 0` and `f(1) = 1` exactly, so a
full-fidelity evaluation is divided by exactly one.

| Vintage | Fitted on | Governs optimisation iterations |
|---|---|---|
| 0 | 20 DOE runs | 1–10 |
| 1 | 20 DOE + 10 OPT | 11–20 |
| … | … | … |
| 9 | 20 DOE + 90 OPT | 91–100 |
| 10 | 20 DOE + 100 OPT | none; recorded for post-hoc use |

Selection of the L2 strength is 5-fold cross-validation over runs, 5 redrawn
partitions, plain argmin of the mean held-out loss over `[0, 1e-2, 1, 1e2]`. No
one-standard-error rule.

**Past rows are never rescaled.** An evaluation keeps the estimate produced by
the vintage in force when it was measured, and `results.csv` records which
vintage that was next to the measured costs and the divisors. Recomputing the
whole history under any single fit afterwards needs no simulation:

```
SSE_under_new_fit = SSE_measured / I_z(a_new, b_new)
```

## What is recorded

`results/results.csv` and `results/init/results.csv` carry, per evaluation:
`eval_id`, timestamp, phase, `beta_vintage`, `z`, the measured costs before
scaling, the two divisors, the scaled costs, `J`, runtime, the count of fmincon
exit flags other than 1, whether the 0.01 divisor floor engaged, and all twelve
theta components at `%.17g`.

`results/registry/`:

- `manifest_<phase>.json` — declared configuration, every optimiser flag, the
  seeds, package versions, git commit and whether the working tree was dirty. On
  resume the configuration is compared against it and any difference stops the
  run.
- `working_tree_<phase>.diff` — the source diff, when the tree was dirty.
- `evaluations_<phase>.jsonl` — one line per evaluation: acquisition value at the
  optimum and at the snapped candidate, reference point, runtime floor, fitted GP
  hyperparameters (lengthscales, outputscale, noise, mean), the acquisition seed,
  the vintage expected and the vintage applied.
- `vintages/vintage_NN.json` — shape parameters, selected lambda, the whole
  cross-validation loss grid, optimiser exit status and iteration counts, the
  eval_ids and filenames that informed the fit, and which iterations it governs.

`results/failures.csv` — evaluations that raised in MATLAB, with the identifier
and message.

## Interruption

Both processes can be killed at any point. Restart the same two commands.

- Completed evaluations are recovered from `results.csv`, which MATLAB writes.
- A ledger line truncated by a kill is dropped on read.
- A crash between MATLAB's row and the driver's line is repaired: the row is
  folded back into the ledger and marked `recovered`, with its proposal metadata
  recorded as lost rather than invented.
- Proposals are seeded per iteration, so a resumed run proposes what an
  uninterrupted one would have proposed at the same index.
- The Sobol design is generated in full every time and indexed by position, so an
  interruption does not shift the stream.
- A `matlab.lock` left by a killed process is treated as abandoned after
  `lock_stale_s` rather than waited on forever.
- A MATLAB process that dies silently ends in a timeout after `eval_timeout_s`,
  not an indefinite wait.
- Five consecutive MATLAB failures stop the serve loop, on the grounds that the
  fault is no longer candidate-specific.

Two conditions stop the run rather than continue:

- The configuration differs from the manifest the run started under.
- A row comes back scaled by a different vintage than the iteration required,
  which would mean MATLAB read a stale coefficient file.

## Changing a setting

Edit `run_config.py`, not the drivers. Everything there reaches the manifest.
Changing a value mid-run is refused on the next start; move `results/` aside to
begin a new run.

## Files

| Path | Role |
|---|---|
| `run_config.py` | declared configuration, case definitions, budget |
| `main.py` | driver: design, proposals, refits, ledger |
| `provenance.py` | manifest, ledger, environment capture, reconciliation |
| `matlab_interface.py` | request and result exchange |
| `J surrogate/runtime_surrogate/fit_beta_surrogate.py` | the beta fit |
| `main_initialization.m`, `main_BO.m` | MATLAB entry points |
| `dependencies/io/serve_requests.m` | shared serve loop and failure handling |
| `dependencies/simulation/beta_eval.m` | `I_z(a, b)` |
| `dependencies/simulation/load_beta_coeffs.m` | coefficient load with validation |

`cheb_eval.m` and `load_cheb_coeffs.m` remain on disk but are off the production
path; the analysis scripts under `Result analysis/` still reference the old
coefficients for the published run.
