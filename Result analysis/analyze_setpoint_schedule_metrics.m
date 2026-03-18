%% Analyze setpoint-schedule controller comparison results
% Loads all out_schedule_*.mat files, plots all schedule cases, and computes
% objective + settling/IAE metrics for each controller.

clear; close all; clc

projectRoot = fileparts(fileparts(mfilename("fullpath")));
resultsRoot = fullfile(projectRoot, "results");

cfg = struct();
cfg.schedule_root = fullfile(resultsRoot, "setpoint_schedule_xsp_7_13_16");
cfg.scenario = "same_noise"; % "same_noise" or "no_noise"
cfg.settling_tol = 0.05;
cfg.graphics_dir = fullfile(resultsRoot, "graphics");
cfg.numerical_dir = fullfile(resultsRoot, "numerical results");

runDir = fullfile(cfg.schedule_root, cfg.scenario);
if ~isfolder(runDir)
    error("Setpoint schedule directory not found: %s", runDir);
end
if ~isfolder(cfg.graphics_dir), mkdir(cfg.graphics_dir); end
if ~isfolder(cfg.numerical_dir), mkdir(cfg.numerical_dir); end

files = dir(fullfile(runDir, "out_schedule_*.mat"));
if isempty(files)
    error("No schedule output files found in %s", runDir);
end

records = table();
plotData = struct("id", {}, "is_benchmark", {}, "T", {}, "cases", {});

for k = 1:numel(files)
    path = fullfile(files(k).folder, files(k).name);
    S = load(path, "out", "ctrl");
    if ~isfield(S, "out")
        warning("Missing 'out' in file, skipping: %s", path);
        continue
    end
    out = S.out;
    ctrl = struct();
    if isfield(S, "ctrl"), ctrl = S.ctrl; end

    [ctrlId, isBenchmark] = get_controller_meta(ctrl, files(k).name);
    if ~isfield(out, "case") || isempty(out.case)
        warning("Missing out.case in file, skipping: %s", path);
        continue
    end

    nCase = min(numel(out.case), 2);
    rec = table(string(ctrlId), logical(isBenchmark), ...
        double(read_num_field(out, "SSE", nan)), ...
        double(read_num_field(out, "SSdU", nan)), ...
        compute_schedule_wall_time(out), ...
        'VariableNames', {'controller_id','is_benchmark','SSE_total','SSdU_total','wall_time_s'});
    rec = append_controller_metadata(rec, ctrl);
    rec.J_total = rec.SSE_total + 1e4 * rec.SSdU_total;

    pdat = struct();
    pdat.id = string(ctrlId);
    pdat.is_benchmark = logical(isBenchmark);
    pdat.T = [];
    pdat.cases = struct([]);

    for c = 1:nCase
        caseData = out.case(c);
        M = summarize_case_metrics_schedule(caseData, cfg.settling_tol);

        rec.("SSE_c" + c) = M.SSE;
        rec.("SSdU_c" + c) = M.SSdU;
        rec.("J_c" + c) = M.SSE + 1e4 * M.SSdU;
        rec.("settle_case_h_c" + c) = M.settle_case_h;
        rec.("IAE_case_c" + c) = M.IAE_case;

        for s = 1:3
            rec.(sprintf("settle_x%d_h_c%d", s, c)) = M.settle_h(s);
            rec.(sprintf("IAE_x%d_c%d", s, c)) = M.IAE(s);
        end

        pdat.T = M.t;
        pdat.cases(c).Y = M.Y;
        pdat.cases(c).Ysp = M.Ysp;
    end

    records = [records; rec]; %#ok<AGROW>
    plotData(end+1) = pdat; %#ok<AGROW>
end

if isempty(records)
    error("No valid schedule outputs were parsed in %s", runDir);
end

records = sortrows(records, ["is_benchmark","J_total"], ["descend","ascend"]);

nonBenchmarkRecords = records(~records.is_benchmark, :);
topCount = min(5, height(nonBenchmarkRecords));
topByC1 = sortrows(nonBenchmarkRecords, "SSE_c1", "ascend");
topByC2 = sortrows(nonBenchmarkRecords, "SSE_c2", "ascend");
topIdsC1 = string(topByC1.controller_id(1:topCount));
topIdsC2 = string(topByC2.controller_id(1:topCount));
coloredControllerIds = intersect(topIdsC1, topIdsC2, "stable");
[peakX1Case2Id, peakX1Case2Val] = find_highest_case_peak_state(plotData, 2, 1, true);
[pareto3Mask, pareto3Ids] = get_schedule_pareto_3obj(records);

fprintf("Setpoint schedule metrics (%s):\n", cfg.scenario);
disp(records);
disp_schedule_pareto_table_3obj(records, pareto3Mask);
disp_case_peak_pick(records, 2, 1, "x1", peakX1Case2Id, peakX1Case2Val);
disp_selected_comparison_table(records, coloredControllerIds, peakX1Case2Id, pareto3Ids);

plot_setpoint_cases_all(plotData, cfg.scenario, cfg.graphics_dir, coloredControllerIds, peakX1Case2Id, pareto3Ids);

csvPath = fullfile(cfg.numerical_dir, "setpoint_schedule_metrics_" + cfg.scenario + ".csv");
writetable(records, csvPath);
fprintf("Saved: %s\n", csvPath);


function v = read_num_field(S, fieldName, defaultVal)
if isfield(S, fieldName)
    v = double(S.(fieldName));
else
    v = defaultVal;
end
end

function wallTime = compute_schedule_wall_time(out)
wallTime = nan;
if ~isstruct(out) || ~isfield(out, "case") || isempty(out.case)
    return
end
caseTimes = nan(numel(out.case), 1);
for k = 1:numel(out.case)
    if isfield(out.case(k), "runtime_s")
        caseTimes(k) = double(out.case(k).runtime_s);
    end
end
if any(isfinite(caseTimes))
    wallTime = sum(caseTimes, "omitnan");
end
end

function [ctrlId, isBenchmark] = get_controller_meta(ctrl, filename)
isBenchmark = false;
ctrlId = erase(string(filename), ".mat");
ctrlId = erase(ctrlId, "out_schedule_");
if ~isempty(fieldnames(ctrl))
    if isfield(ctrl, "id"), ctrlId = string(ctrl.id); end
    if isfield(ctrl, "is_benchmark"), isBenchmark = logical(ctrl.is_benchmark); end
end
if contains(lower(ctrlId), "benchmark")
    isBenchmark = true;
end
end

function rec = append_controller_metadata(rec, ctrl)
rec.source = "";
rec.timestamp = "";
rec.m = nan;
rec.p = nan;
for k = 1:12
    rec.(sprintf("theta_%d", k)) = nan;
end
if isempty(fieldnames(ctrl))
    return
end
if isfield(ctrl, "source"), rec.source = string(ctrl.source); end
if isfield(ctrl, "timestamp"), rec.timestamp = string(ctrl.timestamp); end
if isfield(ctrl, "theta")
    theta = double(ctrl.theta(:).');
    if numel(theta) >= 3
        rec.m = theta(3) + 1;
        rec.p = theta(2) + rec.m;
    end
    for k = 1:min(numel(theta), 12)
        rec.(sprintf("theta_%d", k)) = theta(k);
    end
end
end

function M = summarize_case_metrics_schedule(caseData, settlingTol)
M = struct();
M.Y = nan(0,3);
M.Ysp = nan(0,3);
M.t = nan(0,1);
M.SSE = nan;
M.SSdU = nan;
M.settle_h = nan(1,3);
M.settle_case_h = nan;
M.IAE = nan(1,3);
M.IAE_case = nan;

if ~isfield(caseData, "Y") || ~isfield(caseData, "Ysp")
    return
end
Y = double(caseData.Y);
Ysp = double(caseData.Ysp);
nState = min([size(Y,2), size(Ysp,2), 3]);
if nState < 1
    return
end
Y = Y(:, 1:nState);
Ysp = Ysp(:, 1:nState);

if isfield(caseData, "dt")
    dt = double(caseData.dt);
elseif isfield(caseData, "tf")
    dt = double(caseData.tf) / max(size(Y,1)-1, 1);
else
    dt = 1/60;
end
t = (0:size(Y,1)-1).' * dt;

if isfield(caseData, "SSE")
    M.SSE = double(caseData.SSE);
else
    E = (Y - Ysp) .* [10 1 1];
    M.SSE = sum(sum(E.^2, 2), "omitnan");
end
if isfield(caseData, "SSdU")
    M.SSdU = double(caseData.SSdU);
else
    if isfield(caseData, "U")
        dU = diff(double(caseData.U), 1, 1);
        M.SSdU = sum(sum(dU.^2, 2), "omitnan");
    else
        M.SSdU = nan;
    end
end

settle = nan(1,3);
iae = nan(1,3);
epsRef = 1e-9;
for s = 1:nState
    refVal = Ysp(end, s);
    relErr = abs(Y(:, s) - Ysp(:, s)) / max(abs(refVal), epsRef);
    if relErr(end) <= settlingTol
        idx = find_settling_index_schedule(relErr, settlingTol);
        if ~isempty(idx), settle(s) = t(idx); end
    end
    iae(s) = sum(abs(Y(:, s) - Ysp(:, s)), "omitnan") * dt;
end

M.Y = Y;
M.Ysp = Ysp;
M.t = t;
M.settle_h = settle;
M.IAE = iae;
if any(isfinite(settle))
    M.settle_case_h = max(settle(isfinite(settle)));
end
M.IAE_case = sum(iae, "omitnan");
end

function idx = find_settling_index_schedule(relErr, tol)
idx = [];
n = numel(relErr);
for k = 1:n
    if all(relErr(k:end) <= tol)
        idx = k;
        return
    end
end
end

function [controllerId, peakVal] = find_highest_case_peak_state(plotData, caseId, stateIdx, excludeBenchmark)
controllerId = strings(0,1);
peakVal = nan;
bestVal = -inf;
for i = 1:numel(plotData)
    if excludeBenchmark && plotData(i).is_benchmark
        continue
    end
    if numel(plotData(i).cases) < caseId
        continue
    end
    Y = plotData(i).cases(caseId).Y;
    if isempty(Y) || size(Y, 2) < stateIdx
        continue
    end
    peakNow = max(Y(:, stateIdx), [], "omitnan");
    if isfinite(peakNow) && peakNow > bestVal
        bestVal = peakNow;
        controllerId = string(plotData(i).id);
        peakVal = peakNow;
    end
end
end

function disp_case_peak_pick(records, caseId, stateIdx, stateLabel, controllerId, peakVal)
fprintf("Highest Case %d peak %s:\n", caseId, stateLabel);
if isempty(controllerId)
    fprintf("  none found\n");
    return
end
Tsel = records(string(records.controller_id) == string(controllerId), :);
if isempty(Tsel)
    fprintf("  controller %s not found in records table\n", controllerId);
    return
end
peakCol = sprintf("peak_%s_c%d", stateLabel, caseId);
Tsel.(peakCol) = repmat(peakVal, height(Tsel), 1);
Tsel = movevars(Tsel, peakCol, "After", "controller_id");
disp(Tsel);
end

function disp_selected_comparison_table(records, coloredControllerIds, blackControllerId, pareto3Ids)
benchmarkRows = records(records.is_benchmark, :);
benchmarkIds = intersect(string(benchmarkRows.controller_id), string(pareto3Ids), "stable");
coloredParetoIds = intersect(string(coloredControllerIds), string(pareto3Ids), "stable");
coloredIdsFinal = setdiff(coloredParetoIds, [string(blackControllerId); benchmarkIds], "stable");
keepIds = unique([benchmarkIds; string(blackControllerId); coloredIdsFinal], "stable");
Tcmp = records(ismember(string(records.controller_id), keepIds), :);
if isempty(Tcmp)
    fprintf("Benchmark/selected/colored comparison table: none\n");
    return
end
Tcmp = Tcmp(:, {'controller_id','SSE_total','SSdU_total','wall_time_s'});
Tcmp = sortrows(Tcmp, "SSE_total", "ascend");
fprintf("Benchmark, selected black controller, and colored controllers (sorted by SSE_total):\n");
disp(Tcmp);
end

function [isPareto3, paretoIds] = get_schedule_pareto_3obj(records)
x1IAETotal = double(records.IAE_x1_c1) + double(records.IAE_x1_c2);
objectives = [double(records.SSE_total), x1IAETotal, double(records.SSdU_total)];
isPareto3 = compute_pareto_mask_nd_local(objectives);
paretoIds = string(records.controller_id(isPareto3));
end

function disp_schedule_pareto_table_3obj(records, isPareto3)
x1IAETotal = double(records.IAE_x1_c1) + double(records.IAE_x1_c2);
Tpareto3 = records(isPareto3, :);
Tpareto3.x1_IAE_total = x1IAETotal(isPareto3);
Tpareto3 = movevars(Tpareto3, "x1_IAE_total", "After", "controller_id");
Tpareto3 = sortrows(Tpareto3, ["SSdU_total","SSE_total","x1_IAE_total"], ["ascend","ascend","ascend"]);
fprintf("Setpoint schedule 3-objective Pareto controllers (SSE_total, x1 IAE total, SSdU_total):\n");
disp(Tpareto3);
end

function isPareto = compute_pareto_mask_nd_local(objectives)
n = size(objectives, 1);
isPareto = false(n, 1);
for i = 1:n
    oi = objectives(i, :);
    if ~all(isfinite(oi))
        continue
    end
    dominated = false;
    for j = 1:n
        oj = objectives(j, :);
        if i == j || ~all(isfinite(oj))
            continue
        end
        if all(oj <= oi) && any(oj < oi)
            dominated = true;
            break
        end
    end
    isPareto(i) = ~dominated;
end
end

function plot_setpoint_cases_all(plotData, scenario, graphicsDir, coloredControllerIds, blackControllerId, pareto3Ids)
if isempty(plotData)
    return
end
nCase = 2;
nState = 3;
stateNames = ["x1","x2","x3"];
defaultColors = get(groot, "defaultAxesColorOrder");

fig = figure("Color", "w", "Name", "Setpoint schedule all cases");
tiledlayout(fig, nCase, nState, "TileSpacing", "compact", "Padding", "compact");
set(fig, "Position", [80 80 1400 760]);

for c = 1:nCase
    for s = 1:nState
        ax = nexttile; hold(ax, "on");
        setpointDrawn = false;
        for i = 1:numel(plotData)
            if numel(plotData(i).cases) < c
                continue
            end
            if ~ismember(string(plotData(i).id), pareto3Ids)
                continue
            end
            t = plotData(i).T;
            Y = plotData(i).cases(c).Y;
            Ysp = plotData(i).cases(c).Ysp;
            if isempty(t) || size(Y,2) < s || size(Ysp,2) < s
                continue
            end
            if string(plotData(i).id) == string(blackControllerId)
                plot(ax, t, Y(:, s), "-", "LineWidth", 1.6, "Color", [0 0 0]);
            elseif plotData(i).is_benchmark
                plot(ax, t, Y(:, s), "-", "LineWidth", 2.1, "Color", [1 0 0]);
            elseif ismember(string(plotData(i).id), coloredControllerIds)
                colorIdx = find(coloredControllerIds == string(plotData(i).id), 1, "first");
                colorVal = defaultColors(mod(colorIdx - 1, size(defaultColors, 1)) + 1, :);
                plot(ax, t, Y(:, s), "-", "LineWidth", 1.0, "Color", colorVal);
            else
                line(ax, t, Y(:, s), "LineStyle", "-", "LineWidth", 1.0, "Color", [0.65 0.65 0.65]);
            end
            if ~setpointDrawn
                plot(ax, t, Ysp(:, s), "--", "LineWidth", 1.2, "Color", [0 0 0]);
                setpointDrawn = true;
            end
        end
        grid(ax, "on");
        box(ax, "off");
        xlabel(ax, "Time [h]");
        ylabel(ax, stateNames(s));
        title(ax, sprintf("Case %d, %s", c, stateNames(s)));
    end
end

outStem = fullfile(graphicsDir, "setpoint_schedule_all_cases_" + scenario);
exportgraphics(fig, outStem + ".png", "Resolution", 300);
exportgraphics(fig, outStem + ".pdf", "ContentType", "vector");
fprintf("Saved: %s\n", outStem + ".png");
fprintf("Saved: %s\n", outStem + ".pdf");
end
