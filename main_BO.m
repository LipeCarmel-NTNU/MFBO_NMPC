%% Phase 2: Bayesian optimization runs
%
% This script serves the evaluation requests that carry phase code 1. It logs
% one row for each evaluation.
%
% main_initialization calls this script when the design phase ends. You can also
% run it on its own to resume an interrupted optimization phase.
%
% A run of length tf = 10*z accumulates only part of the full-horizon cost. The
% fitted fraction phi(z) = I_z(a, b) scales the measured partial costs up to
% full-horizon estimates.
%
% The driver refits phi during the run, every refit_every iterations. It
% publishes each fit by renaming a new coefficient file over
% results/surrogate/phi_coeffs.mat. This script reloads that file before every
% evaluation. It stores the vintage that it used in the results row. You can
% therefore trace an objective value to the fit that produced it. Vintage 0 is
% the fit on the initialization runs alone.
%
% A missing coefficient file raises an error. The script applies no fallback. A
% row that no surrogate scaled, or that a surrogate from an earlier run scaled,
% does not compare with the rest.
%
% Outputs:
%   results/results.csv           one summary row for each evaluation
%   results/failures.csv          the evaluations that raised, if any
%   results/out_<timestamp>.mat   the full per-step trends
%   SIMULATIONS_LOG.txt           fmincon exit flags other than 1, and failures
%
% Reproducibility: rng(1) fixes the measurement-noise realization. nmpc_base
% draws it once and every evaluation reuses it. main_initialization.m uses the
% same seed, so both phases see the same disturbance sequence.

% The clear runs whether you start this script yourself or main_initialization
% calls it. Both scripts share the base workspace, so the clear also empties the
% caller. That is the intent: the handover then leaves the same state as closing
% MATLAB after the design phase and starting this script fresh. Nothing follows
% the call in main_initialization, so it needs none of its variables again.
clear all; close all; clc; %#ok<CLALL>

current_dir = fileparts(mfilename('fullpath'));

addpath(genpath(current_dir))

% The lock is cleared, but not inbox/theta.txt. main_initialization hands over
% with an unserved request already in the inbox, and this server has to read it.
delete_if_exists('.lock')
delete_if_exists('matlab.lock')

rng(1)

USE_PARALLEL = true;
NumWorkers = 31;
configure_pool(USE_PARALLEL, NumWorkers);

%% Run configuration
cfg_run = struct();
% Results tree, overridable with MFBO_RESULTS_DIR (default "results"); matches
% pipeline/matlab_interface.py RESULTS_DIR. The exchange stays at the root.
results_root = getenv("MFBO_RESULTS_DIR");
if isempty(results_root); results_root = "results"; end
cfg_run.theta_txt = fullfile("inbox", "theta.txt");
cfg_run.poll_s = 2.0;
cfg_run.out_dir = fullfile(results_root);
cfg_run.results_csv = fullfile(results_root, "results.csv");
cfg_run.failures_csv = fullfile(results_root, "failures.csv");
cfg_run.phi_coeffs = fullfile(results_root, "surrogate", "phi_coeffs.mat");
cfg_run.log_path = fullfile("SIMULATIONS_LOG.txt");
cfg_run.lock_path = "matlab.lock";
cfg_run.lock_stale_s = 6 * 3600;
cfg_run.max_consecutive_failures = 5;
cfg_run.sigma_y = [0.001 0.1 0.1];
cfg_run.NumWorkers = NumWorkers;
cfg_run.rng_seed = 1;
cfg_run.serves_phase = 1;      % 1 = optimization

base = nmpc_base( ...
    sigma_y = cfg_run.sigma_y, ...
    phi_coeffs_path = cfg_run.phi_coeffs);

cfg_run.theta_len = 1 + 2 + base.nx + 2*base.nu;

fprintf("Loaded phi vintage %d (%s)\n", base.phi.vintage, base.phi.created_at);
fprintf("  SSE : a = %.6f, b = %.6f, lambda = %g\n", ...
    base.phi.SSE.a, base.phi.SSE.b, base.phi.SSE.lambda);
fprintf("  SSdU: a = %.6f, b = %.6f, lambda = %g\n", ...
    base.phi.SSdU.a, base.phi.SSdU.b, base.phi.SSdU.lambda);

ensure_dir(cfg_run.out_dir);
check_results_header(cfg_run.results_csv, cfg_run.theta_len);
init_results_csv(cfg_run.results_csv, cfg_run.theta_len);

n_done = count_results_rows(cfg_run.results_csv);
fprintf("Optimisation: %d rows already in %s.\n", n_done, cfg_run.results_csv);

%% External request loop
% The driver writes a request, waits for the matching eval_id to appear in the
% CSV, then writes the next one. The budget lives in run_config.py, so this loop
% runs until you stop it.
serve_requests(cfg_run, @(req) run_and_log(cfg_run, base, req));


%% Local functions

function run_and_log(cfg_run, base, req)
%RUN_AND_LOG Evaluate one request against the current vintage of phi.
% This function reloads the coefficients on every call and not once at startup,
% because the driver refits phi during the run. A reload before every evaluation
% keeps this script independent of the refit cadence. It also survives an
% interruption of either process.
    t_load = tic;
    base.phi = load_phi_coeffs(cfg_run.phi_coeffs);
    wall_load_s = toc(t_load);

    ts = timestamp_compact();

    out = simulate_nmpc(base, req.theta, ...
        horizon = "fidelity", ...
        extrapolate = true, ...
        terminal_cost = "lqr", ...
        run_id = string(ts), ...
        log_path = cfg_run.log_path);

    out.wall_s.phi_load = wall_load_s;

    % The default v7 format. The refit of phi reads these files with
    % scipy.io.loadmat, which cannot open the HDF5-based v7.3 format. The format
    % is a compatibility requirement and not a preference.
    t_save = tic;
    mat_path = fullfile(cfg_run.out_dir, "out_" + ts + ".mat");
    save(mat_path, "ts", "out", "cfg_run", "base");
    wall_save_s = toc(t_save);

    append_results_row(cfg_run.results_csv, req.eval_id, ts, "OPT", out, ...
        req.theta, wall_save_s);

    fprintf(['  OPT eval %d [phi v%d]: SSE=%.6g, SSdU=%.6g, z=%.4f, ' ...
             'solver=%.1fs, wall=%.1fs (save %.2fs)\n'], ...
        req.eval_id, out.phi_vintage, out.SSE, out.SSdU, out.cfg.f, ...
        out.runtime_s, out.wall_s.total, wall_save_s);
end
