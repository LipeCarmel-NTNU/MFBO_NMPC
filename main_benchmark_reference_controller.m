%% Benchmark reference controller, with and without measurement noise
%
% Benchmark definition:
%   m = 6, p = 61, P = 0, Q = diag([10 1 1]), Ru = 0, Rdu = diag([10 10 10])
%
% Conditions match the full-fidelity tests: two initial-condition cases, full
% horizon tf = 10 h at Ts = 1/60 h, and the same model and solver stack as
% main_BO and main_NMPC_test_runs.
%
% Outputs are stored under results/benchmark_fix/<scenario>/, where <scenario>
% is benchmark_full_f1_no_noise or benchmark_full_f1_same_noise.

clear; close all; clc
rng(1)

current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir));

cfg = struct();
cfg.output_root = fullfile("results", "benchmark_fix");
cfg.log_path = fullfile("SIMULATIONS_LOG.txt");
cfg.tf_h = 10;
cfg.Ts = 1/60;
cfg.use_parallel = true;
cfg.num_workers = 2;

scenarios = [
    struct("name", "benchmark_full_f1_no_noise", "sigma_y", [0 0 0]);
    struct("name", "benchmark_full_f1_same_noise", "sigma_y", [0.001 0.1 0.1])
];

configure_pool(cfg.use_parallel, cfg.num_workers);

for s = 1:numel(scenarios)
    scenario = scenarios(s);
    fprintf("\n=== Running scenario: %s ===\n", scenario.name);

    out_dir = fullfile(cfg.output_root, scenario.name);
    ensure_dir(out_dir);
    existingMat = dir(fullfile(out_dir, "*.mat"));
    if ~isempty(existingMat)
        fprintf("Skipping scenario %s (existing MAT found in %s)\n", scenario.name, out_dir);
        continue
    end

    % The benchmark is defined with a zero terminal cost, so the LQR terminal
    % matrix is neither loaded nor built for this run.
    base = nmpc_base( ...
        sigma_y = scenario.sigma_y, ...
        Ts = cfg.Ts, ...
        tf = cfg.tf_h, ...
        load_lqr = false);

    theta = benchmark_theta();

    out = simulate_nmpc(base, theta, ...
        horizon = "full", ...
        extrapolate = false, ...
        terminal_cost = "zero", ...
        verbosity = "quiet", ...
        run_id = scenario.name, ...
        log_path = cfg.log_path);

    J = out.SSE + 1e4 * out.SSdU;

    ts = char(datetime("now", "TimeZone", "Europe/Oslo", "Format", "yyyyMMdd_HHmmss"));
    save(fullfile(out_dir, "out_benchmark.mat"), "ts", "out", "theta", "scenario", "cfg", "base");

    % One row, under the same schema as the optimisation runs. The run is at
    % full horizon with no surrogate, so beta_vintage is NaN and the two
    % fraction columns are 1: the row states that its costs were measured
    % rather than scaled.
    results_csv = fullfile(out_dir, "results_benchmark.csv");
    check_results_header(results_csv, numel(theta));
    init_results_csv(results_csv, numel(theta));
    append_results_row(results_csv, 1, ts, "BENCH", out, theta);

    write_metadata(out_dir, theta, scenario, out);
end


%% Local functions

function write_metadata(out_dir, theta, scenario, out)
%WRITE_METADATA Record the scenario definition next to its results.
    path = fullfile(out_dir, "benchmark_metadata.txt");
    fid = fopen(path, "w");
    if fid < 0
        warning("Could not open metadata file for writing: %s", path);
        return
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, "Scenario: %s\n", scenario.name);
    fprintf(fid, "sigma_y: [%.6g %.6g %.6g]\n", scenario.sigma_y(1), scenario.sigma_y(2), scenario.sigma_y(3));
    fprintf(fid, "Steps with fmincon exit flag other than 1: %d\n", out.n_flag_not_one);
    fprintf(fid, "Benchmark theta:\n");
    fprintf(fid, "%.17g\n", theta(:));
end
