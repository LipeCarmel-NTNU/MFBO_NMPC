# Rerun procedure

Two processes exchange files. MATLAB runs the simulations. Python proposes the
candidates.

MATLAB is a server in both phases. It holds no budget and no design. It answers
requests until you stop it with Ctrl-C. Python owns the budget, and
`run_config.py` declares it in one place.

Start MATLAB first. The driver waits for the results file that MATLAB creates.

## Sequence

Two commands per case.

1. In MATLAB, run `main_initialization`.
2. In a shell, run `python run_pipeline.py --case case1`.

That is the whole run. The launcher runs the design phase, fits phi as vintage
0, waits ten seconds, and then runs the optimization phase.

MATLAB needs no second command. Every request carries a phase code.
`main_initialization` answers code 0. When the driver sends its first
optimization request, which carries code 1, `main_initialization` leaves its
loop without serving that request and calls `main_BO`. `main_BO` reads the same
request from the inbox and serves it, and the rest of the run.

Press Ctrl-C in MATLAB when the launcher reports that it is done.

For the second case, move `results/` aside and repeat with `--case case2`. The
two cases write the same paths.

To run one phase alone:

    python run_pipeline.py --case case1 --phase init
    python run_pipeline.py --case case1 --phase bo

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

The cadence is `refit_every` in `run_config.py`, which declares 5. Vintage v is
fitted on the DOE runs plus the first `v * refit_every` optimization runs, and it
governs the next `refit_every` iterations.

| Vintage | Fitted on | Governs optimization iterations |
|---|---|---|
| 0 | 20 DOE runs | 1 to 5 |
| 1 | 20 DOE + 5 OPT | 6 to 10 |
| ... | ... | ... |
| 19 | 20 DOE + 95 OPT | 96 to 100 |
| 20 | 20 DOE + 100 OPT | none, kept for later analysis |

The terminal vintage is fitted only when `refit_after_last` is true.

Cross-validation selects the L2 strength. It uses 5 folds over runs and one fold
partition. It takes the plain argmin of the mean held-out loss over
`[0, 1e-2, 1, 1e2]`. It applies no one-standard-error rule.

The driver never rewrites an earlier row. A row keeps the estimate that the
vintage in force at the time produced. `results.csv` records that vintage next to
the measured costs and the divisors. You can therefore recompute the whole
history under one fit without a new simulation:

    SSE_under_new_fit = SSE_measured / I_z(a_new, b_new)

That recomputation is what the driver itself does in memory before every
proposal. `history_tensors` in `pipeline/driver.py` re-measures every row under
the vintage in force, so the GP fits one function measured with one instrument
rather than a mixture left by successive vintages. Two consequences: the best
value so far is not monotone over a run, and replaying a trajectory needs the
vintage that each proposal saw, which is why the ledger carries
`phi_vintage_applied`.

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

## Watching a headless run

`run_supervised.py` runs MATLAB with no command window and sends what that
window would have shown to `<results>/logs/matlab_console.log`. That file is
appended for the whole campaign, so it reaches tens of megabytes.

Beside it, `<results>/logs/console_blocks/` holds the same output cut into
blocks of at most 1 MiB or 30 minutes, whichever comes first, with `index.csv`
naming the evaluation ids in each closed block. Read the recent output without
opening the whole log:

    python watch_matlab_log.py --blocks    # follow the block being written
    python watch_matlab_log.py --list      # the index, to find one evaluation
    python watch_matlab_log.py             # the raw log, as before

The blocks are a copy. The raw log is read and never written, so nothing the
mirror does can lose output. Set `MFBO_LOG_BLOCKS=0` to turn the mirror off, or
`MFBO_LOG_BLOCK_BYTES` and `MFBO_LOG_BLOCK_SECONDS` to change the limits. These
are environment variables and not `run_config.py` settings, because every value
in `run_config.py` reaches the manifest and a manifest difference refuses a
resume.

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
| `run_pipeline.py` | the launcher, and the only Python command you type |
| `run_config.py` | the declared configuration, the cases, the budget, both z limits |
| `main_initialization.m` | the MATLAB server for the design phase |
| `main_BO.m` | the MATLAB server for the optimization phase |
| `pipeline/driver.py` | the phase logic: design, proposals, refits, ledger |
| `pipeline/provenance.py` | manifest, ledger, environment capture, reconciliation |
| `pipeline/matlab_interface.py` | the request and result exchange |
| `pipeline/phi_surrogate.py` | the fit of phi |
| `dependencies/io/serve_requests.m` | the shared serve loop and failure handling |
| `dependencies/simulation/phi_eval.m` | phi(z) = I_z(a, b) |
| `dependencies/simulation/load_phi_coeffs.m` | coefficient load with validation |

The Chebyshev evaluator and its loader are removed. The analysis scripts under
`Result analysis/` still read the stored coefficient text file of the published
run.
