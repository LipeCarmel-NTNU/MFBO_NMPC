%% Setpoint-schedule comparison: benchmark vs timestamp-selected controllers
%
% Schedule on second-state setpoint Xsp:
%   0-10 h   : Xsp = 7
%   10-20 h  : Xsp = 13
%   20-30 h  : Xsp = 16
%
% Controllers simulated:
%   1) Benchmark reference controller
%      m = 6, p = 61, Q = diag([10 1 1]), Ru = 0, Rdu = diag([10 10 10]), P = 0
%   2) Controllers selected by timestamp file
%      (defaults to results/txt results/final_pareto_frontier_timestamps_only.txt)
%
% Important:
%   - Benchmark controller keeps P = 0 for all setpoint segments.
%   - Non-benchmark controllers recompute P whenever setpoint changes using
%     the terminal_P_xu_du workflow from TerminalCost/main_TerminalP.m.
%
% Output root:
%   results/setpoint_schedule_xsp_7_13_16/

clear; close all; clc
rng(1)

current_dir = fileparts(mfilename("fullpath"));
addpath(genpath(current_dir));
project_root = current_dir;

cfg = struct();
cfg.source_root = fullfile(project_root, "results");
cfg.output_root = fullfile(project_root, "results", "setpoint_schedule_xsp_7_13_16");
cfg.timestamp_file = fullfile(cfg.source_root, "txt results", "final_pareto_frontier_timestamps_only.txt");
cfg.sigma_y = [0.001 0.1 0.1];
cfg.use_parallel = true;
cfg.num_workers = 2;
cfg.Ts = 1/60;
cfg.tf_h = 30;

schedule = struct();
schedule.Vsp = 1;
schedule.segment_end_h = [10 20 30];
schedule.Xsp_values = [7 13 16];

% Same LQR tuning used in TerminalCost/main_TerminalP.m
lqr_tuning = [-1.9980 0.0003 1.4849 0.5267 -0.9742 0.0425 0.1074 -0.1175];

configure_pool(cfg.use_parallel, cfg.num_workers);
ensure_dir(cfg.output_root);

base = init_base(cfg.sigma_y, cfg.Ts, cfg.tf_h);
controllers = build_controller_set(cfg.source_root, cfg.timestamp_file);

if all(cfg.sigma_y == 0)
    error("This refactored workflow requires noise in both cases (cfg.sigma_y must be nonzero).");
end
scenario_name = "same_noise";
out_dir = fullfile(cfg.output_root, scenario_name);
ensure_dir(out_dir);

summary_csv = fullfile(out_dir, "results_schedule.csv");
init_summary_csv(summary_csv, 12);
write_selected_controller_list(out_dir, controllers, cfg.timestamp_file);

for i = 1:numel(controllers)
    ctrl = controllers(i);
    fprintf("\n=== Simulating controller: %s ===\n", ctrl.id);
    out_mat_path = fullfile(out_dir, "out_schedule_" + ctrl.id + ".mat");
    partial_mat_path = fullfile(out_dir, "out_schedule_" + ctrl.id + "_partial.mat");

    if isfile(out_mat_path)
        S = load(out_mat_path, "out");
        if isfield(S, "out") && is_complete_schedule_output(S.out, base)
            out = S.out;
            fprintf("Skipping simulation (existing result found): %s\n", out_mat_path);
        else
            warning("Existing final result is incomplete; recomputing/resuming: %s", out_mat_path);
            out = simulate_controller_schedule(base, schedule, ctrl, lqr_tuning, partial_mat_path);
            save(out_mat_path, "out", "ctrl", "schedule", "cfg", "base");
            if isfile(partial_mat_path), delete(partial_mat_path); end
        end
    else
        out = simulate_controller_schedule(base, schedule, ctrl, lqr_tuning, partial_mat_path);
        save(out_mat_path, "out", "ctrl", "schedule", "cfg", "base");
        if isfile(partial_mat_path), delete(partial_mat_path); end
    end

    SSE = out.SSE;
    SSdU = out.SSdU;
    runtime_s = out.runtime_s;
    J = SSE + 1e4 * SSdU;

    append_summary_row(summary_csv, ctrl, SSE, SSdU, J, runtime_s);
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

function base = init_base(sigma_y, Ts, tf_h)
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
end

function controllers = build_controller_set(source_root, timestamp_file)
    controllers = struct( ...
        "id", {}, ...
        "source", {}, ...
        "timestamp", {}, ...
        "theta", {}, ...
        "is_benchmark", {});

    controllers(end+1) = struct( ... %#ok<AGROW>
        "id", "benchmark_fix", ...
        "source", "benchmark_fix", ...
        "timestamp", "", ...
        "theta", benchmark_theta(), ...
        "is_benchmark", true);

    timestamps = load_frontier_timestamps(timestamp_file);
    for i = 1:numel(timestamps)
        ts = timestamps(i);
        [src_run, theta] = load_theta_by_timestamp(source_root, ts);
        controllers(end+1) = struct( ... %#ok<AGROW>
            "id", "ts_" + ts, ...
            "source", src_run, ...
            "timestamp", ts, ...
            "theta", theta, ...
            "is_benchmark", false);
    end

    modifiedTimestamps = [
        "20260210_171107";
        "20260201_192807";
        "20260202_000115";
        "20260210_154010";
        "20260201_232106";
        "20260201_192807";
        "20260211_134235";
        "20260211_122653";
    ];
    for i = 1:numel(modifiedTimestamps)
        ts = modifiedTimestamps(i);
        [src_run, theta] = load_theta_by_timestamp(source_root, ts);
        thetaMod = make_schedule_modified_theta(theta);
        controllers(end+1) = struct( ... %#ok<AGROW>
            "id", "ts_" + ts + "_modified", ...
            "source", src_run, ...
            "timestamp", ts, ...
            "theta", thetaMod, ...
            "is_benchmark", false);
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
    ru_exp = -1000 * ones(1, 3); % 10^-1000 underflows to 0 in double precision
    rdu_diag = [10 10 10];
    theta = [f, theta_p, theta_m, log10(q_diag), ru_exp, log10(rdu_diag)];
end

function timestamps = load_frontier_timestamps(path)
    if ~isfile(path)
        error("Timestamp file not found: %s", path);
    end
    lines = strip(string(splitlines(fileread(path))));
    is_ts = lines ~= "" & startsWith(lines, "20");
    timestamps = unique(lines(is_ts), "stable");
end

function [src_run, theta] = load_theta_by_timestamp(source_root, ts)
    mat_run1 = fullfile(source_root, "run1", "out_" + ts + ".mat");
    mat_run2 = fullfile(source_root, "run2", "out_" + ts + ".mat");
    if isfile(mat_run1)
        src_run = "run1";
        S = load(mat_run1, "out");
    elseif isfile(mat_run2)
        src_run = "run2";
        S = load(mat_run2, "out");
    else
        error("Timestamp %s not found in run1 or run2.", ts);
    end
    if ~isfield(S, "out") || ~isfield(S.out, "theta")
        error("Missing out.theta for timestamp %s.", ts);
    end
    theta = S.out.theta(:).';
    theta(1) = 1; % enforce full-horizon settings
end

function thetaMod = make_schedule_modified_theta(theta)
    % Flat gradient wasnt a problem NMPC.optimizer_options.FiniteDifferenceType="central";
% Set theta_7:theta_9 to theta_10:theta_12 - 1 for modified schedule runs.
thetaMod = theta(:).';
% if numel(thetaMod) < 12
%     error("Expected theta length >= 12 for modified schedule controller.");
% end
% thetaMod(7:9) = thetaMod(10:12) - 1;
% thetaMod(1) = 1; % preserve full-horizon setting explicitly
end

function out = simulate_controller_schedule(base, schedule, ctrl, lqr_tuning, partial_mat_path)
    cfg = decode_theta(ctrl.theta, base.nx, base.nu);
    pool_empty = isempty(gcp("nocreate"));

    NMPC = NMPC_terminal(base.model, base.nx, base.nu);
    NMPC.optimizer_options.MaxIterations = 100;
    NMPC.optimizer_options.UseParallel = ~pool_empty;
    NMPC.optimizer_options.Display = "final-detailed";
    NMPC.optimizer_options.FiniteDifferenceType="central";
    NMPC.Ts = base.dt;

    NMPC.p = cfg.p;
    NMPC.m = cfg.m;
    NMPC.Q = cfg.Q;
    NMPC.Ru = cfg.Ru;
    NMPC.Rdu = cfg.Rdu;
    NMPC.constraints();

    out = struct();
    out.controller = ctrl;
    out.cfg = cfg;
    out.schedule = schedule;
    out.T = base.T;
    out.SSE = 0;
    out.SSdU = 0;
    out.runtime_s = 0;
    out.case = struct([]);

    start_case = 1;
    resume_case_state = struct();
    partial = load_partial_state(partial_mat_path, ctrl.id);
    if partial.valid
        out = partial.out;
        start_case = partial.case_id;
        resume_case_state = partial.case_state;
        fprintf("Resuming from partial checkpoint (case %d): %s\n", start_case, partial_mat_path);
    end

    if start_case > 2
        out = recompute_out_totals(out);
        return
    end

    for case_id = start_case:2
        x0 = [1.0, 2.0, 2.0];

        if case_id == start_case && ~isempty(fieldnames(resume_case_state))
            case_resume = resume_case_state;
        else
            case_resume = struct();
        end

        case_out = run_one_case_schedule( ...
            base, NMPC, schedule, ctrl, lqr_tuning, x0, ...
            partial_mat_path, out, case_id, case_resume);

        if ~isfield(out, "case") || isempty(out.case)
            out.case = case_out;
        else
            out.case(case_id) = case_out;
        end

        out = recompute_out_totals(out);
        save_partial_state(partial_mat_path, ctrl.id, out, case_id + 1, struct());
    end

    out = recompute_out_totals(out);
end

function case_out = run_one_case_schedule(base, NMPC, schedule, ctrl, lqr_tuning, x0, partial_mat_path, out_prefix, case_id, case_resume)
    N = base.N;
    noise = base.noise;

    if nargin < 10
        case_resume = struct();
    end

    i_start = 1;
    uk = zeros(1, base.nu);
    Y = zeros(N, base.nx);
    Y_meas = zeros(N, base.nx);
    Ysp = zeros(N, base.nx);
    U = zeros(N, base.nu);
    RUNTIME = zeros(N, 1);
    Xsp2_trace = zeros(N, 1);
    xk = x0;
    currentXsp2 = NaN;

    if ~isempty(fieldnames(case_resume))
        i_start = max(1, min(N, double(case_resume.i_next)));
        if isfield(case_resume, "uk"), uk = double(case_resume.uk); end
        if isfield(case_resume, "xk"), xk = double(case_resume.xk); end
        if isfield(case_resume, "currentXsp2"), currentXsp2 = double(case_resume.currentXsp2); end
        if isfield(case_resume, "Y"), Y = double(case_resume.Y); end
        if isfield(case_resume, "Y_meas"), Y_meas = double(case_resume.Y_meas); end
        if isfield(case_resume, "Ysp"), Ysp = double(case_resume.Ysp); end
        if isfield(case_resume, "U"), U = double(case_resume.U); end
        if isfield(case_resume, "RUNTIME"), RUNTIME = double(case_resume.RUNTIME); end
        if isfield(case_resume, "Xsp2_trace"), Xsp2_trace = double(case_resume.Xsp2_trace); end
    end

    for i = i_start:N
        t_now = base.T(i);
        xsp2_target = scheduled_xsp2(schedule, t_now);

        if i == i_start || abs(xsp2_target - currentXsp2) > 1e-12
            [xss, uss, P] = update_setpoint_and_P(base, NMPC, ctrl, lqr_tuning, schedule.Vsp, xsp2_target);
            NMPC.xsp = xss;
            NMPC.usp = uss;
            NMPC.P = P;
            currentXsp2 = xsp2_target;
        end
        NMPC = apply_case_setpoint_policy(NMPC, case_id, xk, currentXsp2);

        iter_timer = tic;
        Y(i, :) = xk;
        Ysp(i, :) = NMPC.xsp(1:base.nx);
        Xsp2_trace(i) = currentXsp2;

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
disp(uk)
        if should_save_hourly_checkpoint(base.T(i), i, N)
            case_state = struct();
            case_state.i_next = min(i + 1, N + 1);
            case_state.uk = uk;
            case_state.xk = xk;
            case_state.currentXsp2 = currentXsp2;
            case_state.Y = Y;
            case_state.Y_meas = Y_meas;
            case_state.Ysp = Ysp;
            case_state.U = U;
            case_state.RUNTIME = RUNTIME;
            case_state.Xsp2_trace = Xsp2_trace;
            out_partial = out_prefix;
            if ~isfield(out_partial, "case") || isempty(out_partial.case)
                out_partial.case = finalize_case_out(base, noise, x0, Y, Y_meas, Ysp, Xsp2_trace, U, RUNTIME, i);
            else
                out_partial.case(case_id) = finalize_case_out(base, noise, x0, Y, Y_meas, Ysp, Xsp2_trace, U, RUNTIME, i);
            end
            out_partial = recompute_out_totals(out_partial);
            save_partial_state(partial_mat_path, ctrl.id, out_partial, case_id, case_state);
        end
    end

    case_out = finalize_case_out(base, noise, x0, Y, Y_meas, Ysp, Xsp2_trace, U, RUNTIME, N);
end

function NMPC = apply_case_setpoint_policy(NMPC, case_id, xk, xsp2_target)
% Case 1: keep steady-state setpoint from find_ss / terminal_P_xu_du.
% Case 2: overwrite xsp(3) with Ssp = min(3, 2*(Xsp - X)).
if case_id ~= 2
    return
end
X = xk(2);
Xsp = xsp2_target;
Ssp = min(3, 2 * (Xsp - X)); % can be negative by design
if numel(NMPC.xsp) >= 3
    NMPC.xsp(3) = Ssp;
end
end

function xsp2 = scheduled_xsp2(schedule, t_h)
    if t_h < schedule.segment_end_h(1)
        xsp2 = schedule.Xsp_values(1);
    elseif t_h < schedule.segment_end_h(2)
        xsp2 = schedule.Xsp_values(2);
    else
        xsp2 = schedule.Xsp_values(3);
    end
end

function [xss, uss, P] = update_setpoint_and_P(base, NMPC, ctrl, lqr_tuning, Vsp, Xsp)
    if ctrl.is_benchmark
        [xss, uss] = find_ss(Vsp, Xsp, base.par, base.model, base.ode_opt);
        P = zeros(base.nx + base.nu);
        return
    end

    [P, ~, ~, ~, xss, uss, ~] = terminal_P_xu_du_local( ...
        Vsp, Xsp, base.par, base.model, base.ode_opt, base.dt, ...
        lqr_tuning, NMPC.Q, NMPC.Ru, NMPC.Rdu);
end

function [P, K, Ai, Bi, xss, uss, LQR_data] = terminal_P_xu_du_local( ...
        Vsp, Xsp, par, model, ode_opt, Ts, log10w, Q_eval, R1_eval, R2_eval)
% Same workflow as TerminalCost/main_TerminalP.m
    [xss, uss] = find_ss(Vsp, Xsp, par, model, ode_opt);
    [A, B] = linearize(xss, uss, model);
    [Ai, Bi] = incremental(A, B, Ts);

    nx = size(A, 1);
    nu = size(B, 2);
    [K, ~, ~, ~] = build_LQR_full(log10w, Ai, Bi, nx, nu);
    Acl = Ai - Bi * K;

    Sx = [eye(nx), zeros(nx, nu)];
    Su = [zeros(nu, nx), eye(nu)];
    LQR_data = struct("Sx", Sx, "K", K, "Acl", Acl, "Su", Su);
    P = construct_P(LQR_data, Q_eval, R1_eval, R2_eval);
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

function init_summary_csv(path, theta_len)
    fid = fopen(path, "w");
    if fid < 0
        error("Could not open summary csv for writing: %s", path);
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, "controller_id,source,timestamp,is_benchmark,SSE,SSdU,J,runtime_s,m,p");
    for k = 1:theta_len
        fprintf(fid, ",theta_%d", k);
    end
    fprintf(fid, "\n");
end

function append_summary_row(path, ctrl, SSE, SSdU, J, runtime_s)
    cfg = decode_theta(ctrl.theta, 3, 3);
    fid = fopen(path, "a");
    if fid < 0
        error("Could not open summary csv for appending: %s", path);
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, "%s,%s,%s,%d,%.17g,%.17g,%.17g,%.17g,%d,%d", ...
        ctrl.id, ctrl.source, ctrl.timestamp, ctrl.is_benchmark, ...
        SSE, SSdU, J, runtime_s, cfg.m, cfg.p);
    for k = 1:numel(ctrl.theta)
        fprintf(fid, ",%.17g", ctrl.theta(k));
    end
    fprintf(fid, "\n");
end

function write_selected_controller_list(out_dir, controllers, timestamp_file)
    path = fullfile(out_dir, "selected_controllers.txt");
    fid = fopen(path, "w");
    if fid < 0
        warning("Could not open selected controller list for writing: %s", path);
        return
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, "Timestamp source: %s\n\n", timestamp_file);
    for i = 1:numel(controllers)
        c = controllers(i);
        fprintf(fid, "%s | source=%s | timestamp=%s | benchmark=%d\n", ...
            c.id, c.source, c.timestamp, c.is_benchmark);
    end
end

function tf = should_save_hourly_checkpoint(t_h, i, N)
    if i == N
        tf = true;
        return
    end
    tf = abs(t_h - round(t_h)) < 1e-12 && t_h > 0;
end

function case_out = finalize_case_out(base, noise, x0, Y, Y_meas, Ysp, Xsp2_trace, U, RUNTIME, i_last)
    Y_eval = Y(1:i_last, :);
    Ysp_eval = Ysp(1:i_last, :);
    U_eval = U(1:i_last, :);

    E = Y_eval - Ysp_eval;
    E = E .* [10 1 1];
    SSE_vec = sum(E.^2, 2);
    dU = diff(U_eval, 1, 1);
    SSdU_vec = sum(dU.^2, 2);
    SSE = sum(SSE_vec);
    SSdU = sum(SSdU_vec);
    runtime_s = sum(RUNTIME(1:i_last));

    case_out = struct();
    case_out.Y = Y;
    case_out.Y_meas = Y_meas;
    case_out.Ysp = Ysp;
    case_out.Xsp2_trace = Xsp2_trace;
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

function out = recompute_out_totals(out)
    out.SSE = 0;
    out.SSdU = 0;
    out.runtime_s = 0;
    if ~isfield(out, "case") || isempty(out.case)
        return
    end
    for k = 1:numel(out.case)
        c = out.case(k);
        if isfield(c, "SSE") && isfinite(c.SSE), out.SSE = out.SSE + c.SSE; end
        if isfield(c, "SSdU") && isfinite(c.SSdU), out.SSdU = out.SSdU + c.SSdU; end
        if isfield(c, "runtime_s") && isfinite(c.runtime_s), out.runtime_s = out.runtime_s + c.runtime_s; end
    end
end

function save_partial_state(path, ctrl_id, out, case_id, case_state)
    partial = struct();
    partial.ctrl_id = string(ctrl_id);
    partial.case_id = case_id;
    partial.out = out;
    partial.case_state = case_state;
    partial.saved_at = char(datetime("now"));
    save(path, "partial");
end

function partial = load_partial_state(path, ctrl_id)
    partial = struct("valid", false, "case_id", 1, "out", struct(), "case_state", struct());
    if ~isfile(path)
        return
    end
    S = load(path, "partial");
    if ~isfield(S, "partial")
        warning("Ignoring malformed partial file (missing struct): %s", path);
        return
    end
    p = S.partial;
    if ~isfield(p, "ctrl_id") || string(p.ctrl_id) ~= string(ctrl_id)
        warning("Ignoring partial file for different controller id: %s", path);
        return
    end
    if ~isfield(p, "case_id") || ~isfield(p, "out") || ~isfield(p, "case_state")
        warning("Ignoring malformed partial file (missing required fields): %s", path);
        return
    end
    partial.valid = true;
    partial.case_id = double(p.case_id);
    partial.out = p.out;
    partial.case_state = p.case_state;
end

function tf = is_complete_schedule_output(out, base)
tf = false;
if ~isstruct(out) || ~isfield(out, "SSE") || ~isfield(out, "SSdU") || ~isfield(out, "runtime_s")
    return
end
if ~isfield(out, "T") || numel(out.T) ~= base.N
    return
end
if ~isfield(out, "case") || numel(out.case) < 2
    return
end
for case_id = 1:2
    c = out.case(case_id);
    if ~isfield(c, "Y") || size(c.Y, 1) ~= base.N
        return
    end
    if ~isfield(c, "Ysp") || size(c.Ysp, 1) ~= base.N
        return
    end
    if ~isfield(c, "U") || size(c.U, 1) ~= base.N
        return
    end
    if ~isfield(c, "Xsp2_trace") || numel(c.Xsp2_trace) ~= base.N
        return
    end
end
tf = true;
end

function ensure_dir(p)
    if ~isfolder(p)
        mkdir(p);
    end
end
