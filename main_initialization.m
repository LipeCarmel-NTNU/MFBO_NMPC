%% Phase 1: initialization runs for the fidelity surrogate
%
% This script serves the evaluation requests that main.py writes to
% inbox/theta.txt. It stops when the initialization budget is complete.
%
% The Sobol design lives in main.py. This script generates no design of its own.
%
% The script reports the cost at the fidelity that it simulated. The run length
% is tf = 10*z and no surrogate scales the result. This phase produces the data
% that the fit of phi uses. A correction here would make the fit depend on its
% own output.
%
% Outputs:
%   results/init/results.csv          one summary row for each evaluation
%   results/init/failures.csv         the evaluations that raised, if any
%   results/init/out_<timestamp>.mat  the full per-step trends, which include
%                                     case(k).partial_SSE and partial_SSdU
%
% Next step: main.py fits phi as vintage 0. It then runs the optimization phase
% against main_BO.m.
%
% Reproducibility: rng(1) fixes the measurement-noise realization. nmpc_base
% draws it once and every evaluation of this run reuses it.

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
cfg_run.n_init = 20;                                    % must match N_INIT in main.py
cfg_run.theta_txt = fullfile("inbox", "theta.txt");
cfg_run.poll_s = 2.0;
cfg_run.out_dir = fullfile("results", "init");
cfg_run.results_csv = fullfile("results", "init", "results.csv");
cfg_run.failures_csv = fullfile("results", "init", "failures.csv");
cfg_run.log_path = fullfile("SIMULATIONS_LOG.txt");
cfg_run.lock_path = "matlab.lock";
cfg_run.lock_stale_s = 6 * 3600;
cfg_run.max_consecutive_failures = 5;
cfg_run.sigma_y = [0.001 0.1 0.1];
cfg_run.NumWorkers = NumWorkers;
cfg_run.rng_seed = 1;

base = nmpc_base(sigma_y = cfg_run.sigma_y);

cfg_run.theta_len = 1 + 2 + base.nx + 2*base.nu;

ensure_dir(cfg_run.out_dir);
check_results_header(cfg_run.results_csv, cfg_run.theta_len);
init_results_csv(cfg_run.results_csv, cfg_run.theta_len);

%% Initialisation loop
n_done = count_results_rows(cfg_run.results_csv);
fprintf("Initialisation: %d of %d points already complete.\n", n_done, cfg_run.n_init);

cfg_run.stop_after = max(0, cfg_run.n_init - n_done);

serve_requests(cfg_run, @(req) run_and_log(cfg_run, base, req));

fprintf("\nInitialisation complete: %d points in %s\n", ...
    count_results_rows(cfg_run.results_csv), cfg_run.out_dir);
fprintf("main.py fits vintage 0 from these runs and then starts the BO phase.\n");
fprintf("Run main_BO.m next.\n");


%% Local functions

function run_and_log(cfg_run, base, req)
%RUN_AND_LOG Evaluate one request, save the trends and append the CSV row.
% This function writes the .mat before the CSV row. The driver treats the CSV
% row as the record that an evaluation finished. The order therefore matters. A
% crash between the two writes leaves an unused .mat file. The reverse order
% would leave a results row whose trends file is missing.
    ts = timestamp_compact();

    out = simulate_nmpc(base, req.theta, ...
        horizon = "fidelity", ...
        extrapolate = false, ...
        terminal_cost = "lqr", ...
        run_id = string(ts), ...
        log_path = cfg_run.log_path);

    % The default v7 format. The fit of phi reads these files with
    % scipy.io.loadmat, which cannot open the HDF5-based v7.3 format. The format
    % is a compatibility requirement and not a preference.
    t_save = tic;
    mat_path = fullfile(cfg_run.out_dir, "out_" + ts + ".mat");
    save(mat_path, "ts", "out", "cfg_run", "base");
    wall_save_s = toc(t_save);

    append_results_row(cfg_run.results_csv, req.eval_id, ts, "DOE", out, ...
        req.theta, wall_save_s);

    fprintf(['  DOE eval %d: SSE=%.6g, SSdU=%.6g, z=%.4f, ' ...
             'solver=%.1fs, wall=%.1fs (save %.2fs)\n'], ...
        req.eval_id, out.SSE, out.SSdU, out.cfg.f, ...
        out.runtime_s, out.wall_s.total, wall_save_s);
end
