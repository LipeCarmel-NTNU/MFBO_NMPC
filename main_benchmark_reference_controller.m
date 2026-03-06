%% Benchmark reference controller (full horizon) with and without measurement noise
%
% Benchmark definition:
%   m = 6, p = 61, P = 0, Q = diag([10 1 1]), Rdu = diag([10 10 10])
%   Ru uses the same default as NMPC_terminal: diag([2 2 1]).
%
% Conditions match the prior full-fidelity tests:
%   - Two initial-condition cases
%   - Full horizon tf = 10 h, Ts = 1/60 h
%   - Same model/solver stack as main_BO / main_NMPC_test_runs
%
% Outputs are stored under:
%   results/benchmark_reference_controller/<scenario>/
% where <scenario> is one of:
%   - benchmark_full_f1_no_noise
%   - benchmark_full_f1_same_noise

clear; close all; clc
rng(1)

current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir));

cfg = struct();
cfg.output_root = fullfile("results", "benchmark_reference_controller");
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

    base = init_base(scenario.sigma_y, cfg.Ts, cfg.tf_h);
    theta = benchmark_theta();

    out = simulate_benchmark(base, theta);
    SSE = out.SSE;
    SSdU = out.SSdU;
    runtime_s = out.runtime_s;
    J = SSE + 1e4 * SSdU;

    ts = char(datetime("now", "TimeZone", "Europe/Oslo", "Format", "yyyyMMdd_HHmmss"));
    save(fullfile(out_dir, "out_benchmark.mat"), "ts", "out", "theta", "scenario", "cfg", "base");
    write_results_csv(fullfile(out_dir, "results_benchmark.csv"), ts, SSE, SSdU, J, runtime_s, theta);
    write_metadata(out_dir, theta, scenario);
end


%% Local functions

function configure_pool(use_parallel, num_workers)
    p = gcp("nocreate");
    if use_parallel
        if isempty(p) || p.NumWorkers ~= num_workers
            if ~isempty(p)
                delete(p);
            end
            parpool("Processes", num_workers);
        end
    else
        if ~isempty(p)
            delete(p);
        end
    end
end

function theta = benchmark_theta()
% THETA layout:
% [f, theta_p, theta_m, log10(Qdiag_3), log10(Ru_diag_3), log10(Rdu_diag_3)]
    f = 1;
    m = 6;
    p = 61;
    theta_m = m - 1;
    theta_p = p - m;

    q_diag = [10 1 1];
    ru_diag = [2 2 1];
    rdu_diag = [10 10 10];

    theta = [f, theta_p, theta_m, log10(q_diag), log10(ru_diag), log10(rdu_diag)];
end

function base = init_base(sigma_y, Ts, tf_h)
    if nargin < 1 || isempty(sigma_y)
        sigma_y = [0.001 0.1 0.1];
    end
    if nargin < 2 || isempty(Ts)
        Ts = 1/60;
    end
    if nargin < 3 || isempty(tf_h)
        tf_h = 10;
    end

    get_par
    base.par = par;
    base.ode_opt = odeset("NonNegative", [2 3]);

    base.model = @(x, u) dilution_reduced(0, x, u(:)', base.par);
    base.plant = base.model;
    base.nx = 3;
    base.nu = 3;
    base.dt = Ts;
    base.tf = tf_h;
    base.N = ceil(base.tf / base.dt) + 1;
    base.T = (0:base.N-1).' * base.dt;
    base.tspan = [0 base.dt];
    base.sigma_y = sigma_y;
    base.noise = randn(base.N, base.nx) .* base.sigma_y;

    V_sp = 1;
    X_sp = 20;
    [xss, uss] = find_ss(V_sp, X_sp, base.par, base.model, base.ode_opt);
    base.xsp = xss;
    base.usp = uss;
end

function out = simulate_benchmark(base, theta)
    cfg = decode_theta(theta, base.nx, base.nu);
    pool_empty = isempty(gcp("nocreate"));

    NMPC = NMPC_terminal(base.model, base.nx, base.nu);
    NMPC.optimizer_options.MaxIterations = 100;
    NMPC.optimizer_options.UseParallel = ~pool_empty;
    NMPC.optimizer_options.Display = "final-detailed";
    NMPC.Ts = base.dt;

    NMPC.p = cfg.p;
    NMPC.m = cfg.m;
    NMPC.Q = cfg.Q;
    NMPC.Ru = cfg.Ru;
    NMPC.Rdu = cfg.Rdu;
    NMPC.constraints();
    NMPC.xsp = base.xsp;
    NMPC.usp = base.usp;
    NMPC.P = zeros(base.nx + base.nu); % benchmark requirement

    out = struct();
    out.theta = theta(:).';
    out.cfg = cfg;
    out.tf = base.tf;
    out.N = base.N;
    out.T = base.T;
    out.SSE = 0;
    out.SSdU = 0;
    out.runtime_s = 0;

    for case_id = 1:2
        if case_id == 1
            x0 = [1.0, 10, 0];
        else
            x0 = [1.1, 25, 5];
        end
        case_out = run_one_case(base, NMPC, x0);
        out.case(case_id) = case_out;
        out.SSE = out.SSE + case_out.SSE;
        out.SSdU = out.SSdU + case_out.SSdU;
        out.runtime_s = out.runtime_s + case_out.runtime_s;
    end
end

function case_out = run_one_case(base, NMPC, x0)
    N = base.N;
    noise = base.noise;

    uk = zeros(1, base.nu);
    Y = zeros(N, base.nx);
    Y_meas = zeros(N, base.nx);
    Ysp = repmat(NMPC.xsp(1:base.nx), N, 1);
    U = zeros(N, base.nu);
    RUNTIME = zeros(N, 1);

    xk = x0;
    timer = tic;
    for i = 1:N
        iter_timer = tic;
        Y(i, :) = xk;
        yk_meas = xk + noise(i, :);
        yk_meas(yk_meas < 0) = 0;
        Y_meas(i, :) = yk_meas;

        uk = NMPC.solve(yk_meas(:)', uk(:)');
        U(i, :) = uk;

        if N > 1 && i < N
            [~, y] = ode45(@(t, x) base.plant(x, uk), base.tspan, xk, base.ode_opt);
            xk = y(end, :);
        end
        RUNTIME(i) = toc(iter_timer);
    end
    runtime_s = toc(timer);

    E = Y - Ysp;
    E = E .* [10 1 1];
    SSE_vec = sum(E.^2, 2);
    dU = diff(U, 1, 1);
    SSdU_vec = sum(dU.^2, 2);
    SSE = sum(SSE_vec);
    SSdU = sum(SSdU_vec);

    case_out = struct();
    case_out.Y = Y;
    case_out.Y_meas = Y_meas;
    case_out.Ysp = Ysp;
    case_out.U = U;
    case_out.noise = noise;
    case_out.SSE = SSE;
    case_out.SSdU = SSdU;
    case_out.runtime_s = runtime_s;
    case_out.RUNTIME = RUNTIME;
    case_out.dt = base.dt;
    case_out.tf = base.tf;
    case_out.x0 = x0;
end

function cfg = decode_theta(theta, nx, nu)
    theta = theta(:).';
    expected_len = 1 + 2 + nx + nu + nu;
    if numel(theta) ~= expected_len
        error("theta must have length %d.", expected_len);
    end

    k = 1;
    cfg.f = theta(k); k = k + 1;
    theta_p = theta(k); k = k + 1;
    theta_m = theta(k); k = k + 1;

    cfg.m = theta_m + 1;
    cfg.p = theta_p + cfg.m;

    q_exp = theta(k:k+nx-1); k = k + nx;
    ru_exp = theta(k:k+nu-1); k = k + nu;
    rdu_exp = theta(k:k+nu-1);

    cfg.Q = diag(10.^q_exp);
    cfg.Ru = diag(10.^ru_exp);
    cfg.Rdu = diag(10.^rdu_exp);
end

function write_results_csv(path, ts, SSE, SSdU, J, runtime_s, theta)
    fid = fopen(path, "w");
    if fid < 0
        error("Could not open file for writing: %s", path);
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, "timestamp,SSE,SSdU,J,runtime_s");
    for k = 1:numel(theta)
        fprintf(fid, ",theta_%d", k);
    end
    fprintf(fid, "\n");

    fprintf(fid, "%s,%.17g,%.17g,%.17g,%.17g", ts, SSE, SSdU, J, runtime_s);
    for k = 1:numel(theta)
        fprintf(fid, ",%.17g", theta(k));
    end
    fprintf(fid, "\n");
end

function write_metadata(out_dir, theta, scenario)
    path = fullfile(out_dir, "benchmark_metadata.txt");
    fid = fopen(path, "w");
    if fid < 0
        warning("Could not open metadata file for writing: %s", path);
        return
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, "Scenario: %s\n", scenario.name);
    fprintf(fid, "sigma_y: [%.6g %.6g %.6g]\n", scenario.sigma_y(1), scenario.sigma_y(2), scenario.sigma_y(3));
    fprintf(fid, "Benchmark theta:\n");
    fprintf(fid, "%.17g\n", theta(:));
end

function ensure_dir(p)
    if ~isfolder(p)
        mkdir(p);
    end
end
