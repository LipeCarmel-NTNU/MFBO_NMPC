%% Fidelity-1 reruns of the Pareto frontier (full horizon, same noise as BO)
%
% Reruns every frontier controller with f forced to 1 under the same
% measurement noise as the BO phase, so the reported cost is measured over the
% full horizon instead of extrapolated from a truncated run.

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
cfg_run.source_root = "results";
cfg_run.output_root = fullfile("results", "final_fidelity_same_noise");
cfg_run.log_path = fullfile("SIMULATIONS_LOG.txt");
cfg_run.sigma_y = [0.001 0.1 0.1];          % same measurement noise as the BO phase
cfg_run.doe_last_iter = 20;                 % initialisation points precede the BO rows
cfg_run.NumWorkers = NumWorkers;

base = nmpc_base(sigma_y = cfg_run.sigma_y);

%% Frontier timestamps
% The stored list is cross-checked against the frontier recomputed from the
% run CSVs. On a mismatch the computed list wins, because the stored file can
% predate the most recent BO rows.
frontier_ts_file = fullfile(cfg_run.source_root, "txt results", "final_pareto_frontier_timestamps_only.txt");
all_ts_file = load_timestamp_list(frontier_ts_file);
all_ts_calc = load_combined_non_doe_pareto_timestamps(cfg_run.source_root, cfg_run.doe_last_iter);
[all_ts, tsSync] = resolve_frontier_timestamps(all_ts_file, all_ts_calc);
if ~tsSync.in_sync
    fprintf("Using computed combined non-DOE Pareto timestamps (count=%d) due to mismatch with file list.\n", numel(all_ts));
end
[t1, t2] = split_timestamps_by_run(cfg_run.source_root, all_ts);

% recover_csv_rows rebuilds a summary row when the output file survived an
% interruption but the CSV row did not.
rerun_stored_controllers(cfg_run, base, "run1", t1, ...
    out_suffix = "_full_f1_same_noise", recover_csv_rows = true, write_list = true);
rerun_stored_controllers(cfg_run, base, "run2", t2, ...
    out_suffix = "_full_f1_same_noise", recover_csv_rows = true, write_list = true);


%% Local functions

function [timestamps, info] = resolve_frontier_timestamps(tsFile, tsComputed)
%RESOLVE_FRONTIER_TIMESTAMPS Compare the stored list with the computed frontier.
%   Falls back to the computed list when the two disagree.
    tsFile = unique(string(tsFile), "stable");
    tsComputed = unique(string(tsComputed), "stable");

    missingInFile = setdiff(tsComputed, tsFile, "stable");
    extraInFile = setdiff(tsFile, tsComputed, "stable");

    info = struct();
    info.in_sync = isempty(missingInFile) && isempty(extraInFile);
    info.missing_in_file = missingInFile;
    info.extra_in_file = extraInFile;

    if info.in_sync
        timestamps = tsFile;
        fprintf("Frontier timestamp list is in sync with current combined non-DOE Pareto set (%d timestamps).\n", numel(timestamps));
        return
    end

    fprintf("Frontier timestamp list mismatch detected.\n");
    fprintf("  File list count: %d\n", numel(tsFile));
    fprintf("  Computed Pareto count: %d\n", numel(tsComputed));
    if ~isempty(missingInFile)
        fprintf("  Missing in file (present in computed Pareto):\n");
        for i = 1:numel(missingInFile)
            fprintf("    - %s\n", missingInFile(i));
        end
    end
    if ~isempty(extraInFile)
        fprintf("  Extra in file (not in computed Pareto):\n");
        for i = 1:numel(extraInFile)
            fprintf("    - %s\n", extraInFile(i));
        end
    end

    timestamps = tsComputed;
end

function timestamps = load_combined_non_doe_pareto_timestamps(source_root, doeLastIter)
%LOAD_COMBINED_NON_DOE_PARETO_TIMESTAMPS Combined Pareto set from the run CSVs.
    if nargin < 2
        doeLastIter = 20;
    end
    c1 = fullfile(source_root, "run1", "results.csv");
    c2 = fullfile(source_root, "run2", "results.csv");
    T1 = readtable(c1, "TextType", "string");
    T2 = readtable(c2, "TextType", "string");

    req = ["timestamp", "SSE", "SSdU"];
    for r = req
        if ~ismember(r, string(T1.Properties.VariableNames))
            error("Missing required column '%s' in %s", r, c1);
        end
        if ~ismember(r, string(T2.Properties.VariableNames))
            error("Missing required column '%s' in %s", r, c2);
        end
    end

    T1.iteration = (1:height(T1)).';
    T2.iteration = (1:height(T2)).';
    T1 = T1(T1.iteration > doeLastIter, :);
    T2 = T2(T2.iteration > doeLastIter, :);
    Tall = [T1(:, ["timestamp", "SSE", "SSdU"]); T2(:, ["timestamp", "SSE", "SSdU"])];

    isPareto = compute_pareto_mask(double(Tall.SSE), double(Tall.SSdU));
    Tp = Tall(isPareto, :);
    [~, ord] = sort(double(Tp.SSdU), "ascend");
    Tp = Tp(ord, :);
    timestamps = unique(string(Tp.timestamp), "stable");
end

function isPareto = compute_pareto_mask(J_track, J_TV)
%COMPUTE_PARETO_MASK Non-dominated mask for two objectives, lower is better.
    n = numel(J_track);
    isPareto = true(n, 1);
    for i = 1:n
        dominated = (J_track <= J_track(i)) & (J_TV <= J_TV(i)) & ...
                    ((J_track < J_track(i)) | (J_TV < J_TV(i)));
        dominated(i) = false;
        if any(dominated)
            isPareto(i) = false;
        end
    end
end

function [ts_run1, ts_run2] = split_timestamps_by_run(source_root, timestamps)
%SPLIT_TIMESTAMPS_BY_RUN Assign each timestamp to run1 or run2 by file presence.
    ts_run1 = strings(0,1);
    ts_run2 = strings(0,1);
    for i = 1:numel(timestamps)
        ts = timestamps(i);
        p1 = fullfile(source_root, "run1", "out_" + ts + ".mat");
        p2 = fullfile(source_root, "run2", "out_" + ts + ".mat");
        if isfile(p1)
            ts_run1(end+1,1) = ts; %#ok<AGROW>
        elseif isfile(p2)
            ts_run2(end+1,1) = ts; %#ok<AGROW>
        else
            warning("Timestamp not found in run1 or run2 folders: %s", ts);
        end
    end
end

