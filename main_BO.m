%% Phase 2: Bayesian optimisation runs
%
% Serves evaluation requests written to inbox/theta.txt by main.py and logs one
% row per evaluation. A run of length tf = 10*z accumulates only part of the
% full-horizon cost, and the fitted fraction f(z) = I_z(a, b) scales the
% measured partial costs up to full-horizon estimates.
%
% The surrogate is refitted by main.py during the run, every 10 optimisation
% iterations, and published by renaming a new coefficient file over
% results/surrogate/beta_coeffs.mat. This script reloads that file before every
% evaluation and stores the vintage it used in the results row, so an objective
% value can always be traced to the fit that produced it. Vintage 0 is the fit
% on the initialisation runs alone.
%
% Missing coefficients are an error rather than a fallback: a BO row scaled by
% no surrogate, or by one from an earlier run, is not comparable with the rest.
%
% Outputs:
%   results/results.csv           one summary row per evaluation
%   results/failures.csv          evaluations that raised, if any
%   results/out_<timestamp>.mat   full per-step trends
%   SIMULATIONS_LOG.txt           fmincon exit flags other than 1, and failures
%
% Reproducibility: rng(1) fixes the measurement-noise realisation, which is
% drawn once in nmpc_base and reused by every evaluation, and matches the seed
% used by main_initialization.m so that both phases see the same disturbance
% sequence.

clear all; close all; clc;

current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))

delete_if_exists('.lock')
delete_if_exists('matlab.lock')

rng(1)

USE_PARALLEL = true;
NumWorkers = 31;
configure_pool(USE_PARALLEL, NumWorkers);

%% Run configuration
cfg_run = struct();
cfg_run.theta_txt = fullfile("inbox", "theta.txt");
cfg_run.poll_s = 2.0;
cfg_run.out_dir = fullfile("results");
cfg_run.results_csv = fullfile("results", "results.csv");
cfg_run.failures_csv = fullfile("results", "failures.csv");
cfg_run.beta_coeffs = fullfile("results", "surrogate", "beta_coeffs.mat");
cfg_run.log_path = fullfile("SIMULATIONS_LOG.txt");
cfg_run.lock_path = "matlab.lock";
cfg_run.lock_stale_s = 6 * 3600;
cfg_run.max_consecutive_failures = 5;
cfg_run.sigma_y = [0.001 0.1 0.1];
cfg_run.NumWorkers = NumWorkers;
cfg_run.rng_seed = 1;

base = nmpc_base( ...
    sigma_y = cfg_run.sigma_y, ...
    beta_coeffs_path = cfg_run.beta_coeffs);

cfg_run.theta_len = 1 + 2 + base.nx + 2*base.nu;

fprintf("Loaded surrogate vintage %d (%s)\n", base.beta.vintage, base.beta.created_at);
fprintf("  SSE : a = %.6f, b = %.6f, lambda = %g\n", ...
    base.beta.SSE.a, base.beta.SSE.b, base.beta.SSE.lambda);
fprintf("  SSdU: a = %.6f, b = %.6f, lambda = %g\n", ...
    base.beta.SSdU.a, base.beta.SSdU.b, base.beta.SSdU.lambda);

ensure_dir(cfg_run.out_dir);
check_results_header(cfg_run.results_csv, cfg_run.theta_len);
init_results_csv(cfg_run.results_csv, cfg_run.theta_len);

n_done = count_results_rows(cfg_run.results_csv);
fprintf("Optimisation: %d rows already in %s.\n", n_done, cfg_run.results_csv);

%% External request loop
% main.py writes a request, waits for the matching eval_id to appear in the
% CSV, then writes the next one. The budget lives in main.py, so this loop runs
% until it is stopped.
serve_requests(cfg_run, @(req) run_and_log(cfg_run, base, req));


%% Local functions

function run_and_log(cfg_run, base, req)
%RUN_AND_LOG Evaluate one request against the current surrogate vintage.
% The coefficients are reloaded here rather than once at startup, because
% main.py refits during the run. Reloading before every evaluation, instead of
% only when a refit is expected, keeps this script indifferent to the refit
% cadence and to any interruption in either process.
    base.beta = load_beta_coeffs(cfg_run.beta_coeffs);

    ts = timestamp_compact();

    out = simulate_nmpc(base, req.theta, ...
        horizon = "fidelity", ...
        extrapolate = true, ...
        terminal_cost = "lqr", ...
        run_id = string(ts), ...
        log_path = cfg_run.log_path);

    % Saved in the default v7 format. The surrogate refit reads these files with
    % scipy.io.loadmat, which cannot open the HDF5-based v7.3 format, so the
    % format is a compatibility requirement rather than a preference.
    mat_path = fullfile(cfg_run.out_dir, "out_" + ts + ".mat");
    save(mat_path, "ts", "out", "cfg_run", "base");

    append_results_row(cfg_run.results_csv, req.eval_id, ts, "OPT", out, req.theta);

    fprintf("  OPT eval %d [vintage %d]: SSE=%.6g, SSdU=%.6g, z=%.4f, runtime=%.1fs\n", ...
        req.eval_id, out.beta_vintage, out.SSE, out.SSdU, out.cfg.f, out.runtime_s);
end
