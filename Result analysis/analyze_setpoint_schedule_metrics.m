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
        'VariableNames', {'controller_id','is_benchmark','SSE_total','SSdU_total'});
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
[lowestId5h, lowestBiomass5h] = find_lowest_case1_biomass_at_hour(plotData, 5);
[lowestId20h, lowestBiomass20h] = find_lowest_case_biomass_at_hour(plotData, 2, 20);
specialGreyControllerIds = unique([lowestId5h; lowestId20h], "stable");

fprintf("Setpoint schedule metrics (%s):\n", cfg.scenario);
disp(records);
disp_case1_biomass_pick(records, "5 h", lowestId5h, lowestBiomass5h);
disp_case_biomass_pick(records, 2, "20 h", lowestId20h, lowestBiomass20h);

plot_setpoint_cases_all(plotData, cfg.scenario, cfg.graphics_dir, coloredControllerIds, specialGreyControllerIds);

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

function [controllerId, biomassVal] = find_lowest_case1_biomass_at_hour(plotData, targetHour)
 [controllerId, biomassVal] = find_lowest_case_biomass_at_hour(plotData, 1, targetHour);
end

function [controllerId, biomassVal] = find_lowest_case_biomass_at_hour(plotData, caseId, targetHour)
controllerId = strings(0,1);
biomassVal = nan;
bestVal = inf;
for i = 1:numel(plotData)
    if numel(plotData(i).cases) < caseId
        continue
    end
    t = plotData(i).T;
    Y = plotData(i).cases(caseId).Y;
    if isempty(t) || size(Y, 2) < 2
        continue
    end
    idx = find(abs(t - targetHour) < 1e-12, 1, "first");
    if isempty(idx)
        [~, idx] = min(abs(t - targetHour));
    end
    xVal = Y(idx, 2);
    if isfinite(xVal) && xVal < bestVal
        bestVal = xVal;
        controllerId = string(plotData(i).id);
        biomassVal = xVal;
    end
end
end

function disp_case1_biomass_pick(records, hourLabel, controllerId, biomassVal)
disp_case_biomass_pick(records, 1, hourLabel, controllerId, biomassVal);
end

function disp_case_biomass_pick(records, caseId, hourLabel, controllerId, biomassVal)
fprintf("Lowest Case %d biomass at %s:\n", caseId, hourLabel);
if isempty(controllerId)
    fprintf("  none found\n");
    return
end
Tsel = records(string(records.controller_id) == string(controllerId), :);
if isempty(Tsel)
    fprintf("  controller %s not found in records table\n", controllerId);
    return
end
Tsel.case1_biomass = repmat(biomassVal, height(Tsel), 1);
Tsel = movevars(Tsel, "case1_biomass", "After", "controller_id");
disp(Tsel);
end

function plot_setpoint_cases_all(plotData, scenario, graphicsDir, coloredControllerIds, specialGreyControllerIds)
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
            t = plotData(i).T;
            Y = plotData(i).cases(c).Y;
            Ysp = plotData(i).cases(c).Ysp;
            if isempty(t) || size(Y,2) < s || size(Ysp,2) < s
                continue
            end
            if ismember(string(plotData(i).id), specialGreyControllerIds)
                line(ax, t, Y(:, s), "LineStyle", "-", "LineWidth", 1.1, "Color", [0.55 0.55 0.55]);
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
