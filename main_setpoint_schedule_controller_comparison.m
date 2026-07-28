%% Setpoint-schedule comparison: benchmark against timestamp-selected controllers
%
% Schedule on the second-state setpoint Xsp:
%   0-10 h  : Xsp = 7
%   10-20 h : Xsp = 13
%   20-30 h : Xsp = 16
%
% Controllers simulated:
%   1) Benchmark reference controller
%      m = 6, p = 61, Q = diag([10 1 1]), Ru = 0, Rdu = diag([10 10 10]), P = 0
%   2) Controllers selected by timestamp file
%
% The benchmark keeps P = 0 for every segment. The other controllers rebuild P
% whenever the setpoint changes, using terminal_P_xu_du.
%
% Output root: results/setpoint_schedule_xsp_7_13_16/

clear; close all; clc
rng(1)

current_dir = fileparts(mfilename("fullpath"));
addpath(genpath(current_dir));
project_root = current_dir;

cfg = struct();
cfg.source_root = fullfile(project_root, "results");
cfg.output_root = fullfile(project_root, "results", "setpoint_schedule_xsp_7_13_16");
cfg.timestamp_file = fullfile(cfg.source_root, "txt results", "final_pareto_frontier_timestamps_only.txt");
cfg.log_path = fullfile(project_root, "SIMULATIONS_LOG.txt");
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

if all(cfg.sigma_y == 0)
    error("This workflow requires noise in both cases (cfg.sigma_y must be nonzero).");
end

% Setpoints are recomputed for every schedule segment, so the base does not
% carry a fixed steady state.
base = nmpc_base( ...
    sigma_y = cfg.sigma_y, ...
    Ts = cfg.Ts, ...
    tf = cfg.tf_h, ...
    set_setpoint = false, ...
    load_lqr = false);

controllers = build_controller_set(cfg.source_root, cfg.timestamp_file);

scenario_name = "same_noise";
out_dir = fullfile(cfg.output_root, scenario_name);
ensure_dir(out_dir);

summary_csv = fullfile(out_dir, "results_schedule.csv");
init_summary_csv(summary_csv, 12);
write_selected_controller_list(out_dir, controllers, cfg.timestamp_file);

%% Both initial-condition cases start from the same state and differ in the
% sugar setpoint policy applied inside the schedule hook.
x0_cases = [1.0, 2.0, 2.0; 1.0, 2.0, 2.0];

for i = 1:numel(controllers)
    ctrl = controllers(i);
    fprintf("\n=== Simulating controller: %s ===\n", ctrl.id);
    out_mat_path = fullfile(out_dir, "out_schedule_" + ctrl.id + ".mat");
    partial_mat_path = fullfile(out_dir, "out_schedule_" + ctrl.id + "_partial.mat");

    skip = false;
    if isfile(out_mat_path)
        S = load(out_mat_path, "out");
        if isfield(S, "out") && is_complete_schedule_output(S.out, base)
            out = S.out;
            fprintf("Skipping simulation (existing result found): %s\n", out_mat_path);
            skip = true;
        else
            warning("Existing final result is incomplete; recomputing/resuming: %s", out_mat_path);
        end
    end

    if ~skip
        setpoint_fn = make_schedule_setpoint_fn(base, schedule, ctrl.is_benchmark, lqr_tuning);

        out = simulate_nmpc(base, ctrl.theta, ...
            horizon = "full", ...
            extrapolate = false, ...
            terminal_cost = "none", ...
            set_setpoint = false, ...
            x0 = x0_cases, ...
            setpoint_fn = setpoint_fn, ...
            verbosity = "control", ...
            run_id = ctrl.id, ...
            log_path = cfg.log_path, ...
            checkpoint_path = partial_mat_path);

        out.controller = ctrl;
        out.schedule = schedule;
        out = lift_schedule_trace(out);

        save(out_mat_path, "out", "ctrl", "schedule", "cfg", "base");
        if isfile(partial_mat_path), delete(partial_mat_path); end
    end

    J = out.SSE + 1e4 * out.SSdU;
    append_summary_row(summary_csv, ctrl, out.SSE, out.SSdU, J, out.runtime_s);
end


%% Local functions

function out = lift_schedule_trace(out)
%LIFT_SCHEDULE_TRACE Expose the scheduled Xsp as case(k).Xsp2_trace.
%   The hook accumulates the trace in sp_state; the analysis scripts read it
%   from Xsp2_trace.
    for k = 1:numel(out.case)
        if isfield(out.case(k), "sp_state") && isfield(out.case(k).sp_state, "Xsp2_trace")
            out.case(k).Xsp2_trace = out.case(k).sp_state.Xsp2_trace;
        end
    end
end

function controllers = build_controller_set(source_root, timestamp_file)
%BUILD_CONTROLLER_SET Benchmark plus every controller named in the list file.
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

    timestamps = load_timestamp_list(timestamp_file);
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

    % Controllers carrying the "_modified" suffix. Result analysis/
    % analyze_setpoint_schedule_metrics.m keeps the modified entry when both
    % exist and selects one of these ids for the reported figures, so the set
    % is produced here even though make_schedule_modified_theta currently
    % returns theta unchanged.
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
        controllers(end+1) = struct( ... %#ok<AGROW>
            "id", "ts_" + ts + "_modified", ...
            "source", src_run, ...
            "timestamp", ts, ...
            "theta", make_schedule_modified_theta(theta), ...
            "is_benchmark", false);
    end
end

function thetaMod = make_schedule_modified_theta(theta)
%MAKE_SCHEDULE_MODIFIED_THETA Theta for a "_modified" schedule controller.
%   The weight substitution is disabled, so the returned theta equals the
%   input and a "_modified" controller currently duplicates its unmodified
%   counterpart.
    thetaMod = theta(:).';
    % thetaMod(7:9) = thetaMod(10:12) - 1;
    % thetaMod(1) = 1;
end

function [src_run, theta] = load_theta_by_timestamp(source_root, ts)
%LOAD_THETA_BY_TIMESTAMP Find a stored theta in run1 or run2.
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
    theta(1) = 1;   % full-horizon settings
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

function tf = is_complete_schedule_output(out, base)
%IS_COMPLETE_SCHEDULE_OUTPUT True when every case ran the whole horizon.
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
        for fld = ["Y" "Ysp" "U"]
            if ~isfield(c, fld) || size(c.(fld), 1) ~= base.N
                return
            end
        end
        if ~isfield(c, "Xsp2_trace") || numel(c.Xsp2_trace) ~= base.N
            return
        end
    end
    tf = true;
end
