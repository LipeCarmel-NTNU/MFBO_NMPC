# MFBO_NMPC — `main_*` and `Result analysis/` Script Report

A high-level tour of every script in the repository root whose name starts with `main`, and every script in `Result analysis/`. For each main script: what it is trying to accomplish. For each analysis script: what it is testing, the images it generates, and what those images plus the numerical outputs let you conclude.

All paths are relative to the repository root.

## Root: `main*` scripts

### `main.py`

The Python side of the multi-fidelity Bayesian optimization loop. It drives the search: builds a scrambled Sobol DOE in optimization space, then iterates by fitting two GPs (one on the bi-objective vector $-[\text{SSE}, \text{SSdU}]$, one on $\log(\text{runtime})$), forming a `qLogNoisyExpectedHypervolumeImprovement` acquisition, and wrapping it in a custom `MFMOLogAcq` that adds a fidelity-bias term $w_z \log a_z(f)$ and a runtime-penalty term $-w_t \gamma \log \mathbb{E}[t]$. Communicates with MATLAB by writing $\theta$ to `inbox/theta.txt`, waiting for a new row in `results/results.csv`, and logging the same row plus acquisition diagnostics to `results/opt_trace.csv`. Integer dims $\theta_p, \theta_m$ are rounded after proposal, and three reserved input-weight dims are pinned to $-1000$ so $10^{r_u}$ underflows to zero (effectively disabling $R_u$).

### `main_BO.m`

The MATLAB side of the BO loop and the canonical evaluator. Initializes a "base" struct once (plant, setpoints, LQR data for terminal cost, fixed noise realization, Chebyshev surrogate coefficients) and runs in one of three modes: `external` polls `inbox/theta.txt` continuously (the live mode used in tandem with `main.py`); `doe` sweeps a predefined `ThetaDOE`; `single` runs one hardcoded $\theta$ for debugging. The shared evaluator `simulate_nmpc()` runs the NMPC closed loop on two initial conditions over $t_f = 10f$ h, then divides the truncated SSE/SSdU sums by `Cheb5(2f-1, c_*)` (floored at 0.01) to estimate the full-horizon totals. Each $\theta$ is written as one row in `results/results.csv` plus a per-run `.mat` payload.

### `main_NMPC_test_runs.m`

Replays a hand-picked list of Pareto-optimal controllers (two timestamp lists `t1`, `t2`) at full fidelity, no measurement noise. Forces `theta(1) = 1` so the simulation runs the full 10 h horizon and writes the resulting trajectories to `results/test_run/<run>_full_f1_no_noise/out_full_<ts>.mat`. This is what `analyze_test_run_metrics.m` later consumes for settling-time diagnostics.

### `main_benchmark_fix.m`

Runs a fixed reference NMPC controller ($m = 6$, $p = 61$, $Q = \text{diag}([10\;1\;1])$, $R_u = 0$, $R_{\Delta u} = \text{diag}([10\;10\;10])$, $P = 0$) at full fidelity in two scenarios: noiseless and same-noise. Output goes to `results/benchmark_fix/benchmark_full_f1_no_noise/` and `.../benchmark_full_f1_same_noise/`. This is the canonical "baseline NMPC" against which Pareto-selected controllers are compared. The "_fix" variant disables the input penalty by pushing $r^{(u)}$ to $-1000$.

### `main_benchmark_reference_controller.m`

Same idea as `main_benchmark_fix.m`, but keeps the default input penalty $R_u = \text{diag}([2\;2\;1])$ from `NMPC_terminal`. Outputs land under `results/benchmark_reference_controller/`. Pair this with `main_benchmark_fix.m` to study how $R_u$ alone changes the benchmark's behaviour.

### `main_final_fidelity.m`

Replays all non-DOE Pareto-frontier controllers (timestamps either loaded from `results/txt results/final_pareto_frontier_timestamps_only.txt` or recomputed from the run CSVs) at full fidelity ($f = 1$) with the **same** measurement-noise realization as the original BO runs. This is the "refined" evaluation set: it tests whether the BO-selected Pareto candidates retain their ranking when the truncation-and-extrapolation surrogate is replaced by a true full-horizon run. Output: `results/final_fidelity_same_noise/run{1,2}_full_f1_same_noise/`.

### `main_setpoint_schedule_controller_comparison.m`

Simulates each Pareto-selected controller plus the benchmark on a 30 h setpoint-schedule scenario where the biomass setpoint $X_{sp}$ steps through $7 \to 13 \to 16$ at $t = 10, 20$ h. Non-benchmark controllers rebuild the terminal cost $P$ at every setpoint change via the `TerminalCost/main_TerminalP.m` workflow; the benchmark keeps $P = 0$ throughout. This evaluates transferability of the BO-tuned Pareto set to a regime the BO never saw. Output: `results/setpoint_schedule_xsp_7_13_16/same_noise/`.

## `Result analysis/` scripts

### `main_pareto.m`

**What it tests.** The post-processing hub for the two BO runs (`run1`, `run2`). It loads `results.csv` from each run, enforces that DOE iterations $1..20$ are excluded from any Pareto/optimization conclusion, computes the Pareto mask in $(\text{SSE}, \text{SSdU})$ space, summarises runtime by phase, and produces both publication-ready figures and numerical tables. It also delegates the refined-frontier comparison and the noisy/noiseless metrics export.

**Images generated** (all under `results/graphical_results/`, also mirrored into `Overleaf---MFBO-NMPC/` by `run_analysis.m`):

- `pareto_samples_run1_run2.pdf` — combined Pareto-sample scatter, both runs overlaid in $(J_\text{TV}, J_\text{track})$. Shows a well-defined non-dominated front with a clear knee where reducing tracking error starts costing disproportionate control variation.

  A second observation is that the two runs specialize on different segments of the final front: **Case 1 (run1) populates the low-$J_\text{TV}$ tail** ($J_\text{TV}$ down to 0.026, but $J_\text{track}$ as high as 22,781), while **Case 2 (run2) populates the low-$J_\text{track}$ region** ($J_\text{track}$ as low as 12,221 but $J_\text{TV}$ in the 0.13–0.39 range). This is not coincidence: in Case 1's Pareto subset, 5/9 candidates have $N_p = 1$ alongside $N_c = 1$, whereas in Case 2 only 4/10 do — Case 2's Pareto controllers preferentially extend the prediction horizon while keeping $N_c = 1$. A pure $N_c = N_p = 1$ controller is essentially the LQR feedback law applied locally to the nonlinear model: the optimization horizon collapses to one stage cost plus the terminal $P$, so the result is dominated by the LQR solution that built $P$. This produces very smooth, conservative control (hence low $J_\text{TV}$) at the cost of slower convergence (hence higher $J_\text{track}$). When $N_c = 1$ but $N_p > 1$, a single control move is committed but evaluated against an $N_p$-step prediction of free-state evolution under that move, then closed off by the same terminal $P$; this makes the controller more anticipatory and reactive, lowering $J_\text{track}$ but encouraging larger move-to-move corrections. So the apparent "different degrees of freedom" between Case 1 and Case 2 is really the same $N_c = 1$ regime exploring two different sub-modes — LQR-dominated vs prediction-dominated — and BO's path dependence under a finite budget pushed each run into a different sub-mode.

  Numerically, the Pareto controllers also clarify how poorly the benchmark trades off the two objectives. The fixed-benchmark $f=1$ same-noise replays (10 h horizon, two initial conditions) give SSE $\approx 12{,}260$ and SSdU $\approx 3.5$ for the $R_u = 0$ variant and SSE $\approx 12{,}307$ with SSdU $\approx 3.0$ for the default-$R_u$ variant. The top Case 2 Pareto controllers reach essentially the **same tracking error** (SSE $\approx 12{,}221$–$12{,}299$) with **SSdU values of 0.16–0.39 — roughly 8–22$\times$ less control variation**. In other words, the benchmark is excessively aggressive on input movement without any tracking advantage: it sits strictly inside (dominated by) the BO-found Pareto set on this plant.
- `sse_vs_ssdu_side_by_side_z.pdf` — same objective scatter, side by side, points coloured by fidelity $z = f$. Low-$z$ points spread broadly over exploratory regions while high-$z$ points cluster near the front — the intended MFBO behaviour.
- `runtime_vs_iteration_side_by_side.pdf` — per-evaluation runtime against iteration, DOE region highlighted. The DOE block sits visibly higher than the optimization region, confirming that the runtime-aware acquisition shifts the search toward cheaper-but-still-informative points after the model is initialized.
- `runtime_cumulative_run1_run2.pdf` — cumulative wall-clock time for both runs. DOE consumes most of the budget early; the slope flattens noticeably during the BO phase.
- `refined_frontier_change.pdf` (via `run_refined_frontier_change`) — see entry for `plot_refined_frontier_change.m`.

**Numerical outputs** (under `results/numerical results/`):

- `resultssandbox_runtime_and_params.txt`: Case 1 DOE = 509.69 min (60.03%) + Opt = 339.36 min (39.97%), total 849.05 min; Case 2 DOE = 765.23 min (58.08%) + Opt = 552.21 min (41.92%), total 1317.44 min; combined DOE share = 58.85%. Mean per-eval time drops 31.87 → 4.41 min from DOE to BO. Mean optimization-phase fidelity is high (Case 1: 0.819, Case 2: 0.800, combined: 0.810). Optimization-phase $N_c = 1$ share is 98.0% / 80.2% (combined 89.1%); Pareto $m=1$ share is 9/9 in Case 1 and 10/10 in Case 2.
- `final_pareto_frontier_f1_noisy_noiseless_metrics.csv`: per-controller settling/IAE/objective values evaluated at $f=1$ both with and without noise.

**What you can conclude.** Runtime-aware MFBO works as intended: DOE dominates cost, and the BO phase preferentially samples cheaper candidates without collapsing to low fidelity. The recovered Pareto front is shaped by a clear knee, and the strong preference for $N_c = 1$ ($m = 1$) in Pareto solutions suggests that for this plant the dominant lever is not move-horizon freedom but the choice of $N_p$ and the weight balance.

### `resultssandbox.m`

**What it tests.** This is the predecessor of `main_pareto.m`. It performs the same per-run loading, DOE exclusion, Pareto-mask, runtime-summary, and side-by-side plotting workflow, but stops before the refined-frontier and noisy/noiseless metric exports.

**Images generated.** Same set as `main_pareto.m`'s first four PDFs: `pareto_samples_run1_run2.pdf`, `sse_vs_ssdu_side_by_side_z.pdf`, `runtime_vs_iteration_side_by_side.pdf`, `runtime_cumulative_run1_run2.pdf`. Numerical export covers the same `resultssandbox_runtime_and_params.txt` block as above. Conclusions are the same — this script remains useful as a lighter standalone entry point when you do not need the refinement step.

### `analyze_test_run_metrics.m`

**What it tests.** Settling-time, final-error, and IAE diagnostics computed from the `results/test_run/` MAT files produced by `main_NMPC_test_runs.m` (full horizon, no noise). Settling for state $i$ is the earliest time at which $|e_i| / |r_i| \le 5\%$ for all later times; NaN means final error stays above the tolerance. It also assembles the controller-by-case table that downstream Python plotting consumes.

**Images generated.** Two MATLAB figures shown interactively (`test_run_metric_summary` and `p_distribution_by_run_final_pareto`), but the publication PDFs come from the Python sibling script `recover_boxplot_data.py` (see below): `python_settling_time_boxplot.pdf`, `python_np_by_case_boxplot.pdf`, `python_nc_by_case_boxplot.pdf`. The settling-time boxplot makes the state-3 problem obvious: state 1 settles in most runs, state 2 in a clear majority, state 3 in only a minority. The $N_p$ and $N_c$ boxplots show that Case 2 reaches higher $N_p$ than Case 1 while both runs concentrate at $N_c = 1$.

The settling-time boxplot also carries a horizontal reference line at $t = 8.09603$ h — this is exactly $z \cdot T$, where $z = 0.81$ is the combined mean optimization-phase fidelity and $T = 10$ h is the full horizon. The fact that this line cuts through the middle of the settling-time distribution (state 1 mostly settles before it, state 2 straddles it, state 3 mostly does not) is informative: BO appears to have **converged on a fidelity that coincides with the point where most controllers reach steady state**. Past that point, additional simulation time mostly accumulates near-zero residual error and contributes little new information for ranking, so the runtime-aware acquisition correctly identifies it as a wasteful regime. In effect, BO learned an empirical truncation horizon that matches the system's natural settling timescale.

**Numerical outputs** (`results/numerical results/`):

- `analyze_test_run_metrics_summary.txt`: NaN-settling shares 12.50% / 25.00% / 62.50% for states 1/2/3; fraction settled before 8.09603 h is 87.50% / 72.50% / 25.00% for the three states and 25.00% for all-states-simultaneously; final-Pareto controller counts per case (3 in Case 1, 6 in Case 2, 9 total).
- `final_pareto_counts_by_case.txt`, `optimization_nc1_share.txt`: smaller standalone summaries.
- Boxplot CSVs (`boxplot_final_error_data.csv`, `boxplot_settling_time_data.csv`) consumed by Python.

**What you can conclude.** Minimizing $(J_\text{track}, J_\text{TV})$ buys aggregate performance but does **not** guarantee per-state settling — state 3 (substrate $S$) is the bottleneck and settles only 25% of the time at full fidelity. This is the key caveat that motivates considering settling as a constraint or additional objective in future iterations.

### `recover_boxplot_data.py`

**What it tests.** Reads the long-form CSVs that `analyze_test_run_metrics.m` exports and rebuilds the publication-quality boxplots in matplotlib using the Wong (2011) Nature Methods palette. Adds a 8.09603 h reference line (the mean BO simulation time used as a fair comparison threshold).

**Images generated** (saved to `results/graphical_results/`):

- `python_settling_time_boxplot.pdf` — settling time by state ($V, X, S$); shows the strong state-3 outlier pattern.
- `python_np_by_case_boxplot.pdf` — prediction horizon $N_p$ by case; Case 2 medians are visibly higher than Case 1.
- `python_nc_by_case_boxplot.pdf` — control horizon $N_c$ by case; concentrated at 1 in both runs.

**What you can conclude.** Visually re-confirms the numerical findings of `analyze_test_run_metrics.m` and provides the figures the paper cites in the Horizon Selection and Settling Behaviour subsections.

### `analyze_setpoint_schedule_metrics.m`

**What it tests.** Reads every `out_schedule_*.mat` produced by `main_setpoint_schedule_controller_comparison.m`, computes the per-case SSE/SSdU/IAE/settling totals over the full $7 \to 13 \to 16$ schedule, builds a 2-objective Pareto subset, and produces both schedule-trajectory plots and modified-vs-unmodified controller boxplots. It also identifies the top-3 (by relative SSE) non-benchmark controllers as the "coloured" controllers in subsequent figures.

**Images generated** (`results/graphics/`, with the two manuscript-named copies sitting directly under `results/`):

- `results/setpoint_tracking.pdf` (corresponds to `setpoint_schedule_all_cases_same_noise.pdf`) — every controller's state trajectories across both initial conditions and all three schedule segments. Pareto controllers and the selected controller are highlighted; non-Pareto controllers are shown as thin neutral lines for context.
- `results/setpoint_tracking_MV.pdf` (corresponds to `setpoint_schedule_selected_benchmark_fin_only_same_noise.pdf`) — $F_\text{in}$ time series for the selected controller and the benchmark, both cases. Lets you compare manipulated-variable activity directly between the two.
- `setpoint_schedule_black_controller_inputs_same_noise.pdf` + `..._mls_same_noise.pdf` — full input trace for the chosen black-marker controller, in L/h and mL/s units.
- `setpoint_schedule_selected_controller_inputs_same_noise.pdf` + `..._mls_same_noise.pdf` — same for the selected controller.
- `setpoint_schedule_modified_vs_unmodified_boxplots_same_noise.pdf` — side-by-side boxplots of objectives/settling for the modified vs unmodified controller variants.

**How the "selected" controller was chosen.** The script identifies candidates programmatically and then commits to one by hand:

1. Compute the 2-objective Pareto set on the schedule scenario (`get_schedule_pareto_2obj`) using the combined-case totals $(\text{SSE}_\text{total}, \text{SSdU}_\text{total})$.
2. Restrict to non-benchmark Pareto candidates, compute their $\text{SSE}_\text{total}$ relative to the benchmark, sort ascending, and keep the top three as `coloredControllerIds`.
3. Hard-code one of those candidates as `blackControllerId` — currently `ts_20260211_122653_modified`. The comments in the script document why this controller wins over the two runner-ups: against the benchmark it achieves **16.6× faster wall time** (6.02% of the benchmark's), a **4.49× smaller volume tracking error**, an **8.87% lower SSdU**, and **biomass IAE within 2%**, at the cost of only **11.54% worse substrate tracking**. The two commented-out alternatives (`ts_20260210_180826_modified` and `ts_20260210_151703_modified`) each beat the benchmark on subsets of those metrics but lose on others — typically by trading substrate quality for volume quality in a way the author judged less acceptable.

The discussion in `results/setpoint_tracking.pdf` and `results/setpoint_tracking_MV.pdf` should therefore be read as a deliberate "engineering pick" rather than a single-objective optimum: the script auto-narrows down to the top-3 by relative SSE, but the actual deployment recommendation is the candidate that simultaneously holds up on wall time, individual-state settling, and IAE balance. The MV plot is the supporting evidence — it visibly shows the selected controller producing smoother $F_\text{in}$ profiles than the benchmark while still tracking the schedule changes.

**Numerical outputs.** `setpoint_schedule_metrics_same_noise.csv` — one row per controller with `SSE_total`, `SSdU_total`, `wall_time_s`, per-case sub-totals, settling times per state per case, and IAE per state per case. The benchmark row shows much higher SSE and SSdU than the modified Pareto controllers, while a couple of modified controllers achieve SSE close to (or below) the benchmark with markedly lower SSdU.

**What you can conclude.** The Pareto controllers selected on the fixed setpoint generalize meaningfully to a multi-step schedule the BO never trained on: several "modified" variants beat the benchmark in both objectives, and modified-vs-unmodified boxplots show that the schedule-aware $P$ recomputation provides a tangible but bounded improvement. The selection process highlights that the BO Pareto set serves as a curated shortlist for a human decision-maker — once the optimization narrows the field to a handful of non-dominated candidates, picking the deployment controller is a small, transparent multi-criteria choice rather than another optimization problem. This is the strongest external-validity evidence in the project.

### `plot_refined_frontier_change.m`

**What it tests.** Direct before/after comparison of the optimization Pareto frontier (after DOE exclusion) against the refined frontier obtained by re-evaluating the same controllers at full fidelity with the same noise (output of `main_final_fidelity.m`). The fairness guards are explicit: only non-DOE original points are eligible, and matching is by `(run_key, timestamp)` rather than timestamp alone.

**Images generated.**

- `results/graphical_results/refined_frontier_change.pdf` — two-panel layout. Left: combined frontier from matched non-DOE points (context). Right: original points overlaid with refined evaluations, refined-Pareto points highlighted as purple diamonds, and "promoted" controllers (not Pareto originally but Pareto after refinement, with $z < 1$) marked as black-edged stars.

**Numerical outputs.**

- `results/txt results/refined_promoted_frontier_z_lt_1.txt` — list of controllers that became Pareto only after full-fidelity refinement.

**What you can conclude.** The shape of the original frontier is largely preserved after refinement, but a small number of controllers that were dominated under truncated evaluation become Pareto-optimal under full fidelity. That is exactly the kind of mismatch the multi-fidelity assumption forces you to confront: truncated evaluation is a useful ranker, but a handful of decisions near the frontier are sensitive to the extrapolation.

### `exploring_u.m`

**What it tests.** A purely post-hoc input-smoothing experiment for the selected controller. Inputs are partitioned at nodes every $\approx 0.1$ h. For each interval between consecutive nodes and each input channel, an LP enforces the cumulative-sum constraint $\sum_i u^\star_i = \sum_i u_i$ (so the integral over each block is preserved) while minimizing total absolute increments $\sum_i |u^\star_i - u^\star_{i-1}|$. The point is to ask: how much of the controller's TV is structurally necessary versus optimization-artefact noise?

**Images generated.**

- `results/graphics/exploring_u_selected_controller_same_noise.pdf` — original $U$ vs reconstructed $U^\star$ as staircase plots for all three inputs across both cases. The reconstructed series is visibly smoother but preserves the same coarse trajectory.

**Numerical outputs.** Printed reconstruction summary table with `TV_*`, `TVstar_*`, and max-absolute-difference per input per case, plus the largest node-residual error.

**What you can conclude.** Much of the controller's high-frequency input activity can be removed without changing the cumulative input applied at each 0.1-hour node. This suggests the BO-selected control trace contains structural slack — there is room for a tighter $R_{\Delta u}$ penalty (or a post-processing smoothing step) without sacrificing the tracking trajectory.

### `run_analysis.m`

**What it tests.** Nothing on its own — it is the orchestrator. It runs `analyze_test_run_metrics.m`, then `recover_boxplot_data.py`, then `main_pareto.m`, then `J surrogate/surrogate.m`, then `results/numerical results/compile_results_txt.m`, and finally copies every PDF in `results/graphical_results/` into `Overleaf---MFBO-NMPC/`. Use this as the one-button pipeline that regenerates the paper's figures and numerical tables after any data update.

## A note on the surrogate cost fit (`J surrogate/`)

Although these scripts live outside `Result analysis/`, `run_analysis.m` invokes them and they produce two of the figures cited above. They are worth a short explanation because they underpin every "estimated full-horizon" value in the BO loop.

For each stored evaluation, define the normalized cumulative cost trajectory $y(\tau) = C(\tau) / C(1)$, where $C(\tau)$ is the cumulative SSE (or SSdU) at relative time $\tau \in [0, 1]$ and $C(1)$ is the total accumulated at $\tau = 1$. This trajectory is by construction $0$ at $\tau = 0$ and $1$ at $\tau = 1$, and its shape encodes how cost accumulates over time. `J surrogate/surrogate.m` fits a single degree-5 Chebyshev polynomial $\phi(x)$ in $x = 2\tau - 1$ to the aggregated $y(\tau)$ values across all stored runs — one fit for SSE, one for SSdU. At runtime, given a truncated simulation that stopped at $\tau = f$, the extrapolated full-horizon cost is simply $\hat C(1) = C(f) / \phi(f)$. The aggregated training fits are very tight ($R^2 = 0.9998$ for SSE, $0.9919$ for SSdU), which is why the surrogate is an excellent **ranker** of candidates even at low fidelity. Coefficients are recorded in `results/numerical results/surrogate_cheb_coeffs.txt`.

`J surrogate/test_surrogate.m` validates extrapolation on held-out runs by reconstructing the predicted full-horizon cost from each truncated stored run and comparing it to the actually observed full-horizon total. It produces two PDFs:

- `surrogate_test_ratio_vs_time.pdf` — every run's $\widehat C(\tau)/C_\text{final}$ trajectory overlaid in faint colour, with the ensemble median in bold.
- `surrogate_test_ratio_vs_time_95band.pdf` — the pointwise 2.5th–97.5th percentile envelope across all trajectories, also showing the median.

The 95-band figure is the most informative diagnostic. Two observations stand out. First, the envelope is **narrow for $\tau$ ranging above $\approx 0.81$ but widens noticeably for $\tau < 0.81$**: extrapolation from a truncation point that lies within the active transient is much less reliable than from a truncation point near or past settling. The fact that BO's combined mean optimization-phase fidelity converged to $z \approx 0.81$ is therefore not an arbitrary number — it lines up almost exactly with the empirical reliability boundary of the surrogate. The runtime-aware acquisition implicitly discovered the lowest fidelity at which extrapolation remains trustworthy, which is also (per the previous section) approximately the system's natural settling time.

Second, the tracking-cost (SSE) panel of the band figure shows a **persistent modelling bias**: the median ratio sits noticeably below 1 across most of $\tau$, meaning the surrogate **systematically underestimates the true full-horizon tracking cost**. The control-variation (SSdU) panel is comparatively unbiased. The interpretation is that SSE accumulates most of its weight early in the transient — once states reach near-setpoint, the per-step contribution is small — and a single Chebyshev shape fit across heterogeneous controllers cannot fully capture controller-dependent transient profiles. In practice this bias is mostly harmless for ranking (it scales all candidates similarly), but it is the right place to look for residual MFBO error.

A natural extension of this workflow, **not implemented** in the current pipeline, would be to re-fit the surrogate coefficients periodically during the BO run as new evaluations accumulate. The surrogate currently uses one offline fit from a small set of runs; refitting online with each batch of new full-fidelity points would tighten both the bias and the low-$\tau$ envelope, and would likely allow BO to push the operating fidelity below $z \approx 0.81$ without losing ranking quality — a direct path to additional runtime savings.

A separate observation concerns the choice of model class itself. The surrogate $\phi(z)$ is what gets used at runtime to convert truncated cost into full-horizon cost, so its fit must be robust under online re-estimation with limited data. An unregularized degree-5 polynomial is a questionable choice for this role: six free parameters from a handful of (potentially noisy) cumulative-cost trajectories invites overfitting and Runge-style oscillations near the endpoints, which is exactly the regime where extrapolation matters most. A more defensible choice would be a low-parameter, regularized model — for example a degree-2 or degree-3 monotone spline, an isotonic-regression fit, or a Beta-CDF-style family with two or three shape parameters and an L2 penalty — and to **bake in the structural constraints that the surrogate must satisfy by construction**: $\phi(0) = 0$ (no cost has accumulated at zero time), $\phi(1) = 1$ (the normalization point), and $\phi'(z) \ge 0$ on $(0,1)$ (cumulative cost is non-decreasing in time). The current Chebyshev fit can violate the monotonicity property in principle, which is why the implementation defensively clips $\phi$ to $[0.01, 1]$ at evaluation time. A constrained fit removes the need for that clamp and tends to behave much better when refit on small online samples — both of which matter if the goal is to push the operating fidelity floor below the current $z \approx 0.81$.

## Conclusion

The repository is organized into a clean two-layer structure: a small set of `main_*` scripts at the root drives every kind of NMPC closed-loop simulation that the project needs (live MFBO loop, fixed Pareto replay at full fidelity with and without noise, fixed-benchmark controllers, and the setpoint-schedule transferability test), and the `Result analysis/` folder turns the resulting MAT/CSV artefacts into the figures and tables that the manuscript references.

Several findings hold up consistently across the analysis scripts. First, runtime-aware MFBO does what it promises: DOE consumes about 59% of the wall time, the optimization phase shifts to cheaper-but-still-high-fidelity points (mean $z \approx 0.81$), and the resulting Pareto front has a robust knee that overlaps cleanly between the two independent runs. The choice of $z \approx 0.81$ is not arbitrary — it coincides with the system's natural settling timescale (the $z \cdot T = 8.096$ h reference line on the settling-time boxplot cuts through the middle of the empirical settling distribution) and with the lower boundary at which the surrogate extrapolation remains tight in the 95% band figure. BO effectively discovered an empirical truncation horizon set by both physics and surrogate reliability.

Second, the Chebyshev fraction surrogate is an excellent ranker — $R^2 = 0.992$ for SSdU, $0.9998$ for SSE, with 98.76% of files within 1% endpoint error — but the validation figure also reveals a structural caveat: the tracking surrogate is biased low (it underestimates true full-horizon tracking cost) and its envelope widens significantly below $\tau \approx 0.81$. Periodically re-estimating the Chebyshev coefficients during the BO run (instead of using one offline fit) would likely correct the bias and push the usable fidelity floor lower, unlocking additional runtime savings. The refined-frontier comparison already shows that a small number of decisions near the Pareto boundary are sensitive to the extrapolation, and full-fidelity refinement does promote previously-dominated controllers onto the front.

Third, the two BO runs specialize on different segments of the Pareto front — Case 1 on low $J_\text{TV}$, Case 2 on low $J_\text{track}$ — and inspection of the horizon choices explains why: both runs concentrate at $N_c = 1$, but Case 1 frequently also picks $N_p = 1$ (collapsing the controller onto the LQR feedback law derived from the terminal cost $P$, which gives smooth conservative actions), whereas Case 2 retains $N_p > 1$ (committing one control move across an anticipatory prediction window, which gives sharper tracking at the cost of larger move-to-move corrections). The benchmark sits strictly inside this front: with SSE $\approx 12{,}260$–$12{,}307$ it matches the best Pareto controllers on tracking, but its SSdU of 3.0–3.5 is 8–22$\times$ larger than what the Pareto set achieves at the same tracking quality — the benchmark is excessively aggressive on inputs with no offsetting benefit.

Fourth, settling is the gap the analysis exposes: minimizing the bi-objective $(J_\text{track}, J_\text{TV})$ is not equivalent to settling all states, and state 3 (substrate) reaches the 5% band in only 25% of the full-horizon replays. The setpoint-schedule analysis is the most positive transferability evidence: several modified Pareto controllers outperform the benchmark on a 30 h schedule with multiple $X_{sp}$ steps that the BO never saw.

Operationally, the analysis scripts are well-instrumented for reuse — `run_analysis.m` regenerates everything in one shot, and the strict DOE-exclusion guards in `main_pareto.m`, `resultssandbox.m`, and `plot_refined_frontier_change.m` prevent the most common source of misreporting in Bayesian-optimization analyses.
