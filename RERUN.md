# Rerun procedure

Two processes exchange files. MATLAB runs the simulations. Python proposes the
candidates. Start MATLAB first in both phases. The driver waits for the results
file that MATLAB creates.

## Sequence

Phase 1 is the design of experiments. It evaluates 20 Sobol points. Each point
runs at the fidelity that its own z asks for. No surrogate acts in this phase.

1. In MATLAB, run `main_initialization`.
2. In a shell, run `python main.py init --case case1`.

The driver fits phi and writes vintage 0 when the design finishes. It then
stops.

Phase 2 is the optimization. It makes 100 acquisition-driven evaluations.

1. In MATLAB, run `main_BO`.
2. In a shell, run `python main.py bo --case case1`.

To run the second case, repeat both phases with `--case case2`. Move `results/`
aside first. The two cases write the same paths.

## The surrogate phi

phi(z) is the share of the full-horizon cost that a run accumulates by fidelity
z. The driver divides a measured partial cost by phi(z).

The model is the regularized incomplete beta function:

    phi(z) = I_z(a, b)

It replaces the Chebyshev polynomial. phi(0) = 0 and phi(1) = 1 hold exactly.
A full-fidelity evaluation therefore divides by exactly one.

Two limits act on the fidelity. They are different quantities.

- `z_min_bo` is 0.01. It bounds the search space below. At Ts = 1 min over a
  10 h horizon, z = 0.01 covers 6 min and 6 simulation steps.
- `z_min_phi` is 0.1. It is the floor of the fit. A sample below it leaves the
  fit. A run below it still receives a phi correction.

`decode_theta.m` clamps z into [0, 1]. That clamp guards against a malformed
value. It is not a fidelity policy.

## Refit schedule

| Vintage | Fitted on | Governs optimization iterations |
|---|---|---|
| 0 | 20 DOE runs | 1 to 10 |
| 1 | 20 DOE + 10 OPT | 11 to 20 |
| ... | ... | ... |
| 9 | 20 DOE + 90 OPT | 91 to 100 |
| 10 | 20 DOE + 100 OPT | none, kept for later analysis |

Cross-validation selects the L2 strength. It uses 5 folds over runs and one fold
partition. It takes the plain argmin of the mean held-out loss over
`[0, 1e-2, 1, 1e2]`. It applies no one-standard-error rule.

The driver never rescales an earlier row. A row keeps the estimate that the
vintage in force produced. `results.csv` records that vintage next to the
measured costs and the divisors. You can therefore recompute the whole history
under one fit without a new simulation:

    SSE_under_new_fit = SSE_measured / I_z(a_new, b_new)

## Records

`results/results.csv` and `results/init/results.csv` hold one row for each
evaluation. Each row carries the evaluation index, the timestamp, the phase, the
phi vintage, z, the measured costs before scaling, the two divisors, the scaled
costs, J, the solver time, the count of fmincon exit flags other than 1, whether
the divisor floor acted, five wall times, and the twelve theta components at
`%.17g`.

`results/registry/` holds four kinds of file.

- `manifest_<phase>.json`. The declared configuration, every optimizer flag, the
  seeds, the package versions, the git commit, and whether the working tree was
  dirty. On resume the driver compares the configuration against this file. Any
  difference stops the run.
- `working_tree_<phase>.diff`. The source diff, when the tree was dirty.
- `evaluations_<phase>.jsonl`. One line for each evaluation. It holds the
  acquisition value at the optimum and at the snapped candidate, the reference
  point, the runtime floor, the fitted GP hyperparameters, the acquisition seed,
  the vintage expected, the vintage applied, and the wall time of each step.
- `vintages/vintage_NN.json`. The shape parameters, the selected lambda, the
  whole cross-validation loss grid, the optimizer exit status, the evaluations
  that informed the fit, the iterations that the vintage governs, and the wall
  time of each stage of the fit.

`results/failures.csv` holds the evaluations that raised in MATLAB, with the
identifier and the message.

## Wall time

The pipeline measures and stores the time of each task.

MATLAB stores these times in `out.wall_s` and in the CSV: the whole
`simulate_nmpc` call, the case loop, the phi evaluation, the `build_nmpc` call,
and the `.mat` write. `runtime_s` remains the solver time summed over the
control steps.

Python stores these times in the ledger: the objective GP fit, the log-runtime
GP fit, the acquisition maximization, the whole proposal, the MATLAB round trip,
the ledger write, and the whole iteration. Each vintage record holds the load
time, the cross-validation time, the final fit time, the record write, and the
coefficient publish.

The end-of-run summary prints the totals and the share of wall time that the
driver used.

## Interruption

You can kill either process at any point. Restart the same two commands.

- The driver recovers a finished evaluation from `results.csv`, which MATLAB
  writes.
- The reader drops a ledger line that a kill truncated.
- A crash between the MATLAB row and the driver line leaves a gap. The driver
  folds the row back into the ledger and marks it `recovered`. It records the
  proposal metadata as lost. It does not invent it.
- The driver seeds each proposal from the iteration index. A resumed run
  proposes what an uninterrupted run would propose.
- The driver generates the whole Sobol design each time and indexes it by
  position. An interruption therefore does not shift the stream.
- The driver treats a `matlab.lock` from a killed process as abandoned after
  `lock_stale_s`. It does not wait on it forever.
- A MATLAB process that dies without writing ends in a timeout after
  `eval_timeout_s`.
- Five failures in a row stop the MATLAB serve loop. Past that point the fault
  is no longer specific to one candidate.

Two conditions stop a run instead of continuing.

- The configuration differs from the manifest that the run started under.
- A row comes back scaled by a vintage that the iteration did not require. This
  means that MATLAB read a stale coefficient file.

## To change a setting

Edit `run_config.py`. Do not edit the drivers. Every value in `run_config.py`
reaches the manifest.

The driver refuses a value that changed in the middle of a run. Move `results/`
aside to start a new run.

## Files

| Path | Role |
|---|---|
| `run_config.py` | the declared configuration, the cases, the budget, both z limits |
| `main.py` | the driver: design, proposals, refits, ledger |
| `provenance.py` | manifest, ledger, environment capture, reconciliation |
| `matlab_interface.py` | the request and result exchange |
| `J surrogate/runtime_surrogate/fit_beta_surrogate.py` | the fit of phi |
| `main_initialization.m`, `main_BO.m` | the MATLAB entry points |
| `dependencies/io/serve_requests.m` | the shared serve loop and failure handling |
| `dependencies/simulation/phi_eval.m` | phi(z) = I_z(a, b) |
| `dependencies/simulation/load_phi_coeffs.m` | coefficient load with validation |

The Chebyshev evaluator and its loader are removed. The analysis scripts under
`Result analysis/` still read the stored coefficient text file of the published
run.
