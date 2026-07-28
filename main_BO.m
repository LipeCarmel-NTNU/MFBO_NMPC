%% Phase 2: Bayesian optimisation runs
%
% Serves theta vectors written to inbox/theta.txt by main.py and logs one row
% per evaluation. Requires the fidelity surrogate fitted on the initialisation
% phase: a run length of tf = 10*f accumulates only part of the full-horizon
% cost, and the fitted fraction curves frac_SSE(f) and frac_SSdU(f) scale the
% measured partial costs up to full-horizon estimates.
%
% Run main_initialization.m and the surrogate fit first. Missing coefficients
% are an error rather than a fallback, so that every BO result is traceable to
% the initialisation set that produced its surrogate.
%
% Outputs:
%   results/results.csv           one summary row per theta (read by main.py)
%   results/out_<timestamp>.mat   full per-step trends
%   SIMULATIONS_LOG.txt           every fmincon exit flag other than 1
%
% Reproducibility: rng(1) fixes the noise realisation, which is drawn once in
% nmpc_base and reused by every theta evaluation.

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
cfg_run.cheb_coeffs = fullfile("results", "surrogate", "cheb_coeffs.mat");
cfg_run.log_path = fullfile("SIMULATIONS_LOG.txt");
cfg_run.sigma_y = [0.001 0.1 0.1];
cfg_run.NumWorkers = NumWorkers;

base = nmpc_base( ...
    sigma_y = cfg_run.sigma_y, ...
    cheb_coeffs_path = cfg_run.cheb_coeffs);

ensure_dir(cfg_run.out_dir);
ensure_parent_dir(cfg_run.results_csv);
theta_len = 1 + 2 + base.nx + 2*base.nu;
init_results_csv(cfg_run.results_csv, theta_len);

%% External theta loop
% main.py writes a theta, waits for the matching row to appear in the CSV, then
% writes the next one. The lock file signals that MATLAB is busy.
lock = 'matlab.lock';
last_signature = "";

while true
    [theta, signature, ok] = read_theta_from_txt(cfg_run.theta_txt);

    if ~ok
        warning("Failed to read theta. Time: %s", ...
            char(datetime('now', 'TimeZone', 'Europe/Oslo')));
        pause(cfg_run.poll_s);
        continue
    end

    if signature == last_signature
        disp('Stale theta signature')
        pause(cfg_run.poll_s);
        continue
    end
    last_signature = signature;

    fid = fopen(lock, 'w');
    if fid < 0
        warning('Failed to create lock (pwd = %s)', pwd);
    else
        fclose(fid);
    end

    try
        run_and_log(cfg_run, base, theta);
    catch ME
        delete_if_exists(lock);
        rethrow(ME);
    end

    delete_if_exists(lock);
end


%% Local functions

function run_and_log(cfg_run, base, theta)
%RUN_AND_LOG Evaluate one theta, append the CSV row and save the trends.
    ts = timestamp_compact();

    out = simulate_nmpc(base, theta, ...
        horizon = "fidelity", ...
        extrapolate = true, ...
        terminal_cost = "lqr", ...
        run_id = string(ts), ...
        log_path = cfg_run.log_path);

    J = out.SSE + 1e4 * out.SSdU;
    append_results_row(cfg_run.results_csv, ts, out.SSE, out.SSdU, J, out.runtime_s, theta);

    mat_path = fullfile(cfg_run.out_dir, "out_" + ts + ".mat");
    save(mat_path, "ts", "out", "cfg_run", "base");
end
