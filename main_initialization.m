%% Phase 1: initialisation runs for the fidelity surrogate
%
% Serves theta vectors written to inbox/theta.txt by main.py and stops once
% the initialisation budget is reached. The Sobol design lives in main.py
% (N_INIT points, sobol_unique_points, seed 1234), so this script does not
% generate a design of its own.
%
% Costs are reported at the fidelity that was actually simulated: the run
% length is tf = 10*f and no surrogate extrapolation is applied. The purpose
% of this phase is to produce the data the surrogate is fitted to, so applying
% a surrogate here would make the fit depend on its own output.
%
% Outputs:
%   results/init/results.csv          one summary row per theta
%   results/init/out_<timestamp>.mat  full per-step trends, including
%                                     case(k).partial_SSE and partial_SSdU
%
% Next step, once this script stops:
%   python "J surrogate/runtime_surrogate/fit_runtime_surrogate.py" results/init
% then run main_BO.m.

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
cfg_run.n_init = 20;                                            % must match N_INIT in main.py
cfg_run.theta_txt = fullfile("inbox", "theta.txt");
cfg_run.poll_s = 2.0;
cfg_run.out_dir = fullfile("results", "init");
cfg_run.results_csv = fullfile(cfg_run.out_dir, "results.csv");
cfg_run.log_path = fullfile("SIMULATIONS_LOG.txt");
cfg_run.sigma_y = [0.001 0.1 0.1];
cfg_run.NumWorkers = NumWorkers;

base = nmpc_base(sigma_y = cfg_run.sigma_y);

ensure_dir(cfg_run.out_dir);
theta_len = 1 + 2 + base.nx + 2*base.nu;
init_results_csv(cfg_run.results_csv, theta_len);

%% Initialisation loop
n_done = numel(load_results_csv_timestamps(cfg_run.results_csv));
fprintf("Initialisation: %d of %d points already complete.\n", n_done, cfg_run.n_init);

last_signature = "";
while n_done < cfg_run.n_init
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

    lock = 'matlab.lock';
    fid = fopen(lock, 'w');
    if fid < 0
        warning('Failed to create lock (pwd = %s)', pwd);
    else
        fclose(fid);
    end

    try
        ts = timestamp_compact();
        out = simulate_nmpc(base, theta, ...
            horizon = "fidelity", ...
            extrapolate = false, ...
            terminal_cost = "lqr", ...
            run_id = string(ts), ...
            log_path = cfg_run.log_path);

        append_results_row(cfg_run.results_csv, ts, out.SSE, out.SSdU, ...
            out.SSE + 1e4*out.SSdU, out.runtime_s, theta);
        save(fullfile(cfg_run.out_dir, "out_" + ts + ".mat"), "ts", "out", "theta", "cfg_run", "base");
    catch ME
        delete_if_exists(lock);
        rethrow(ME);
    end
    delete_if_exists(lock);

    n_done = n_done + 1;
    fprintf("Initialisation point %d of %d complete (SSE=%.6g, SSdU=%.6g, f=%.4f).\n", ...
        n_done, cfg_run.n_init, out.SSE, out.SSdU, out.cfg.f);
end

fprintf("\nInitialisation complete: %d points in %s\n", n_done, cfg_run.out_dir);
fprintf("Fit the fidelity surrogate before starting BO:\n");
fprintf('  python "J surrogate/runtime_surrogate/fit_runtime_surrogate.py" "%s"\n', cfg_run.out_dir);
fprintf("Then run main_BO.m.\n");
