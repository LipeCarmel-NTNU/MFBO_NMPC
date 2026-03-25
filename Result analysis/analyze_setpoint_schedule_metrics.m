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
        pdat.cases(c).U = M.U;
    end

    records = [records; rec]; %#ok<AGROW>
    plotData(end+1) = pdat; %#ok<AGROW>
end

if isempty(records)
    error("No valid schedule outputs were parsed in %s", runDir);
end

% Drop unmodified duplicates: if both "foo" and "foo_modified" exist, keep
% only "foo_modified". Controllers with no modified counterpart are kept as-is.
allIds = string(records.controller_id);
modifiedIds = allIds(endsWith(allIds, "_modified"));
baseIdsWithModified = erase(modifiedIds, "_modified");
dropMask = ismember(allIds, baseIdsWithModified);
records(dropMask, :) = [];
plotData(dropMask) = [];

records = sortrows(records, ["is_benchmark","J_total"], ["descend","ascend"]);

[pareto2Mask, pareto2Ids] = get_schedule_pareto_2obj(records);
nonBenchmarkRecords = records(~records.is_benchmark, :);
benchmarkRows = records(records.is_benchmark, :);
coloredControllerIds = strings(0, 1);
if ~isempty(benchmarkRows)
    benchmarkSSE = double(benchmarkRows.SSE_total(1));
    paretoNonBenchmark = nonBenchmarkRecords(ismember(string(nonBenchmarkRecords.controller_id), pareto2Ids), :);
    if benchmarkSSE ~= 0 && ~isempty(paretoNonBenchmark)
        paretoNonBenchmark.rel_SSE = double(paretoNonBenchmark.SSE_total) ./ benchmarkSSE;
        paretoNonBenchmark = sortrows(paretoNonBenchmark, ["rel_SSE","SSE_total"], ["ascend","ascend"]);
        topCount = min(3, height(paretoNonBenchmark));
        coloredControllerIds = string(paretoNonBenchmark.controller_id(1:topCount));
    end
end


%% SELECTED

% blackControllerId = "benchmark_ref";
blackControllerId = "ts_20260211_122653_modified";
% Good alternative to benchmark
% as long as 11.54% worse S tracking is acceptable. 
% 6.0208% wall time (16.6 times faster)
% 4.4897 times smaller volume tracking error
% 8.8700% lower SSdU
% overall SSE  biomass IAE within 2%

% blackControllerId = "ts_20260210_180826_modified"; 
% Better than benchmark for most pratical cases
% 29.82% rel SSdU
% 3.8958% wall time
%  volume IAE 8% higher due to larger transient peak but no offset. 
% x2 2.2% worse (equivalent)
%  substrate 15.77% better 

% blackControllerId = "ts_20260210_151703_modified";
% 29.78% rel SSdU 
%  little compromise (4% SSE) 
% x1 to x3 tradeoff 35.64% lower x1 and 68.91% higher x3
% Higher substrate error could be problematic

%%



fprintf("Setpoint schedule metrics (%s):\n", cfg.scenario);
disp(sortrows(records, "SSE_total", "ascend"));
disp_schedule_pareto_table_2obj(records, pareto2Mask);
disp_schedule_pareto_relative_table(records, pareto2Mask);
disp_black_controller_pick(records, blackControllerId);
disp_selected_comparison_table(records, coloredControllerIds, blackControllerId, pareto2Ids);

plot_setpoint_cases_all(plotData, cfg.scenario, cfg.graphics_dir, coloredControllerIds, blackControllerId, pareto2Ids);
% plot_black_controller_inputs(plotData, cfg.scenario, cfg.graphics_dir, blackControllerId);
% plot_black_controller_inputs_mls(plotData, cfg.scenario, cfg.graphics_dir, blackControllerId);
plot_selected_benchmark_fin_only(plotData, cfg.scenario, cfg.graphics_dir, blackControllerId);
% plot_modified_unmodified_boxplots(records, cfg.scenario, cfg.graphics_dir);

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
M.U = nan(0,3);
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
if isfield(caseData, "U")
    M.U = double(caseData.U);
end
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

function disp_black_controller_pick(records, blackControllerId)
fprintf("Selected black controller parameters and metrics:\n");
controllerId = string(blackControllerId);
if isempty(controllerId)
    fprintf("  none found\n");
    return
end
Tsel = records(string(records.controller_id) == string(controllerId), :);
if isempty(Tsel)
    fprintf("  controller %s not found in records table\n", controllerId);
    return
end
disp(Tsel);
end

function disp_selected_comparison_table(records, coloredControllerIds, blackControllerId, pareto2Ids)
benchmarkRows = records(records.is_benchmark, :);
benchmarkIds = string(benchmarkRows.controller_id);
keepIds = unique([benchmarkIds; string(blackControllerId); string(coloredControllerIds(:))], "stable");
TcmpFull = records(ismember(string(records.controller_id), keepIds), :);
if isempty(TcmpFull)
    fprintf("Benchmark/selected/colored comparison table: none\n");
    return
end
Tcmp = TcmpFull;
Tcmp.x1_IAE_total = compute_state_iae_total(Tcmp, 1);
Tcmp.x2_IAE_total = compute_state_iae_total(Tcmp, 2);
Tcmp.x3_IAE_total = compute_state_iae_total(Tcmp, 3);
Tcmp.is_black = string(Tcmp.controller_id) == string(blackControllerId);
Tcmp.is_pareto2 = ismember(string(Tcmp.controller_id), string(pareto2Ids));
Tcmp = Tcmp(:, {'controller_id','SSE_total','SSdU_total','x1_IAE_total','x2_IAE_total','x3_IAE_total','wall_time_s','is_pareto2','is_black'});
Tcmp = sortrows(Tcmp, "SSE_total", "ascend");
fprintf("Benchmark, selected black controller, and colored controllers (sorted by SSE_total; is_pareto2 indicates SSE_total/SSdU_total Pareto membership):\n");
disp(Tcmp);

Ttune = sortrows(TcmpFull, "SSE_total", "ascend");
tuningCols = {'controller_id','m','p', ...
    'theta_1','theta_2','theta_3','theta_4','theta_5','theta_6', ...
    'theta_7','theta_8','theta_9','theta_10','theta_11','theta_12'};
presentCols = tuningCols(ismember(tuningCols, Ttune.Properties.VariableNames));
Ttune = Ttune(:, presentCols);
fprintf("Same controllers with tuning only (same SSE_total sorting):\n");
disp(Ttune);
end

function [isPareto2, paretoIds] = get_schedule_pareto_2obj(records)
objectives = [double(records.SSE_total), double(records.SSdU_total)];
isPareto2 = compute_pareto_mask_nd_local(objectives);
paretoIds = string(records.controller_id(isPareto2));
end

function disp_schedule_pareto_table_2obj(records, isPareto2)
Tpareto2 = records(isPareto2, :);
Tpareto2.x1_IAE_total = compute_state_iae_total(Tpareto2, 1);
Tpareto2.x2_IAE_total = compute_state_iae_total(Tpareto2, 2);
Tpareto2.x3_IAE_total = compute_state_iae_total(Tpareto2, 3);
Tpareto2 = movevars(Tpareto2, "x1_IAE_total", "After", "controller_id");
Tpareto2 = movevars(Tpareto2, "x2_IAE_total", "After", "x1_IAE_total");
Tpareto2 = movevars(Tpareto2, "x3_IAE_total", "After", "x2_IAE_total");
Tpareto2 = sortrows(Tpareto2, ["SSdU_total","SSE_total","x1_IAE_total","x2_IAE_total","x3_IAE_total"], ["ascend","ascend","ascend","ascend","ascend"]);
fprintf("Setpoint schedule 2-objective Pareto controllers (Pareto objectives: SSE_total, SSdU_total):\n");
disp(Tpareto2);
end

function disp_schedule_pareto_relative_table(records, isPareto2)
benchmarkRows = records(records.is_benchmark, :);
if isempty(benchmarkRows)
    fprintf("Setpoint schedule Pareto benchmark-relative table: benchmark not found\n");
    return
end

Trel = records(isPareto2, :);
if isempty(Trel)
    fprintf("Setpoint schedule Pareto benchmark-relative table: no Pareto controllers\n");
    return
end

benchmarkRow = benchmarkRows(1, :);
Trel.rel_SSE = double(Trel.SSE_total) ./ double(benchmarkRow.SSE_total(1));
Trel.rel_SSdU = double(Trel.SSdU_total) ./ double(benchmarkRow.SSdU_total(1));
Trel.rel_wall_time = double(Trel.wall_time_s) ./ double(benchmarkRow.wall_time_s(1));
Trel.rel_x1_IAE = compute_state_iae_total(Trel, 1) ./ compute_state_iae_total(benchmarkRow, 1);
Trel.rel_x2_IAE = compute_state_iae_total(Trel, 2) ./ compute_state_iae_total(benchmarkRow, 2);
Trel.rel_x3_IAE = compute_state_iae_total(Trel, 3) ./ compute_state_iae_total(benchmarkRow, 3);
Trel = Trel(:, {'controller_id','rel_SSE','rel_SSdU','rel_wall_time','rel_x1_IAE','rel_x2_IAE','rel_x3_IAE'});
Trel = sortrows(Trel, "rel_SSE", "ascend");
fprintf("Setpoint schedule 2-objective Pareto controllers relative to benchmark:\n");
disp(Trel);
end

function stateIAETotal = compute_state_iae_total(T, stateIdx)
stateIAETotal = double(T.(sprintf("IAE_x%d_c1", stateIdx))) + ...
    double(T.(sprintf("IAE_x%d_c2", stateIdx)));
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

function plot_setpoint_cases_all(plotData, scenario, graphicsDir, coloredControllerIds, blackControllerId, paretoIds)
if isempty(plotData)
    return
end
% coloredControllerIds is unused here; selected, benchmark, and Pareto
% membership already determine the plotted trajectories.
nCase = 2;
nState = 3;
fontSize = 14;
NATURE_COLOR = nature_methods_colors();
plotColors = nature_methods_colors(3);
stateLabels = ["$V$ [L]","$X$ [g/L]","$S$ [g/L]"];
casePanelLabels = ["a","b"];
thinLineColor = blend_color_with_white([0 0 0], 0.35);

fig = figure("Color", "w", "Name", "Setpoint schedule all cases");
tiledlayout(fig, nCase, nState, "TileSpacing", "compact", "Padding", "compact");
set(fig, "Position", [80 80 1400 760]);
stateAxes = gobjects(nCase, nState);

for c = 1:nCase
    for s = 1:nState
        ax = nexttile; hold(ax, "on");
        stateAxes(c, s) = ax;
        plot_state_panel_contents(ax, plotData, c, s, paretoIds, blackControllerId, plotColors, thinLineColor, NATURE_COLOR.Vermillion);
        if c < nCase
            xlabelText = "";
        else
            xlabelText = "Time [h]";
        end
        xlabel(ax, xlabelText, "Interpreter", "latex");
        ylabel(ax, stateLabels(s), "Interpreter", "latex");
        if s == 1
            title(ax, "$\mathbf{" + casePanelLabels(c) + "}$", "Interpreter", "latex");
            ax.TitleHorizontalAlignment = "left";
        else
            title(ax, "");
        end
        ax.TickLabelInterpreter = "latex";
        set(ax, "FontSize", fontSize);
        grid(ax, "off");
        box(ax, "off");
    end
end

drawnow;
for c = 1:nCase
    add_biomass_inset(fig, stateAxes(c, 2), plotData, c, paretoIds, blackControllerId, plotColors, thinLineColor, NATURE_COLOR.Vermillion);
end

outStem = fullfile(graphicsDir, "setpoint_schedule_all_cases_" + scenario);
exportgraphics(fig, outStem + ".png", "Resolution", 300);
exportgraphics(fig, outStem + ".pdf", "ContentType", "vector");
fprintf("Saved: %s\n", outStem + ".png");
fprintf("Saved: %s\n", outStem + ".pdf");
end

function plot_state_panel_contents(ax, plotData, caseIdx, stateIdx, paretoIds, selectedControllerId, plotColors, thinLineColor, setpointColor)
[tSetpoint, yspMedian] = compute_case_median_setpoint(plotData, caseIdx, stateIdx, paretoIds, selectedControllerId);

for pass = 1:3
    for i = 1:numel(plotData)
        if ~should_plot_schedule_controller(plotData(i), caseIdx, paretoIds, selectedControllerId)
            continue
        end
        t = plotData(i).T;
        Y = plotData(i).cases(caseIdx).Y;
        if isempty(t) || size(Y, 2) < stateIdx
            continue
        end

        controllerId = string(plotData(i).id);
        if controllerId == string(selectedControllerId)
            linePass = 3;
            lineColor = plotColors(1, :);
            lineWidth = 2.2;
            lineStyle = "-";
        elseif plotData(i).is_benchmark
            linePass = 2;
            lineColor = plotColors(2, :);
            lineWidth = 2.0;
            lineStyle = ":";
        else
            linePass = 1;
            lineColor = thinLineColor;
            lineWidth = 1.0;
            lineStyle = "-";
        end

        if pass == linePass
            plot(ax, t, Y(:, stateIdx), "LineStyle", lineStyle, "LineWidth", lineWidth, "Color", lineColor);
        end
    end
end

if ~isempty(tSetpoint) && ~isempty(yspMedian)
    plot(ax, tSetpoint, yspMedian, "--", "LineWidth", 1.2, "Color", setpointColor);
end
end

function add_biomass_inset(fig, mainAx, plotData, caseIdx, paretoIds, selectedControllerId, plotColors, thinLineColor, setpointColor)
switch caseIdx
    case 1
        xWindow = [17.75, 19.25];
        yWindow = [12.95, 13.05];
    case 2
        xWindow = [6, 10];
        yWindow = [6.95, 7.1];
    otherwise
        return
end

drawnow;
mainPos = mainAx.Position;
insetPos = [ ...
    mainPos(1) + 0.52 * mainPos(3), ...
    mainPos(2) + 0.06 * mainPos(4), ...
    0.43 * mainPos(3), ...
    0.40 * mainPos(4)];

insetAx = axes("Parent", fig, "Units", "normalized", "Position", insetPos, "Color", "w"); hold(insetAx, "on");
plot_state_panel_contents(insetAx, plotData, caseIdx, 2, paretoIds, selectedControllerId, plotColors, thinLineColor, setpointColor);
xlim(insetAx, xWindow);
ylim(insetAx, yWindow);
set(insetAx, "FontSize", 10);
insetAx.TickLabelInterpreter = "latex";
grid(insetAx, "off");
box(insetAx, "on");
end

function plot_selected_benchmark_fin_only(plotData, scenario, graphicsDir, selectedControllerId)
plotColors = nature_methods_colors(3);
selectedIdx = find(string({plotData.id}) == string(selectedControllerId), 1, "first");
if isempty(selectedIdx)
    warning("Selected controller not found for F_in plot: %s", selectedControllerId);
    return
end
benchmarkIdx = find([plotData.is_benchmark], 1, "first");
if isempty(benchmarkIdx)
    warning("Benchmark controller not found for F_in plot.");
    return
end

fontSize = 14;
casePanelLabels = ["a","b"];

fig = figure("Color", "w", "Name", "Selected and benchmark F_in");
tiledlayout(fig, 2, 1, "TileSpacing", "compact", "Padding", "compact");
set(fig, "Position", [160 120 760 760]);

for c = 1:2
    if numel(plotData(selectedIdx).cases) < c || numel(plotData(benchmarkIdx).cases) < c
        continue
    end

    tSelected = double(plotData(selectedIdx).T(:));
    USelected = double(plotData(selectedIdx).cases(c).U);
    tBenchmark = double(plotData(benchmarkIdx).T(:));
    UBenchmark = double(plotData(benchmarkIdx).cases(c).U);
    if isempty(tSelected) || isempty(USelected) || isempty(tBenchmark) || isempty(UBenchmark)
        continue
    end

    ax = nexttile; hold(ax, "on");
    plot(ax, tSelected, USelected(:, 1), "-", "LineWidth", 1.8, ...
        "Color", plotColors(1, :), "DisplayName", "Selected");

    markerIdx = compute_fin_peak_marker_indices(tBenchmark, UBenchmark(:, 1));
    benchmarkLine = plot(ax, tBenchmark, UBenchmark(:, 1), ":", ...
        "LineWidth", 1.8, "Color", plotColors(2, :), ...
        "Marker", ".", "MarkerSize", 14, "DisplayName", "Benchmark");
    if ~isempty(markerIdx)
        benchmarkLine.MarkerIndices = markerIdx;
    end

    if c < 2
        xlabelText = "";
    else
        xlabelText = "Time [h]";
    end
    xlabel(ax, xlabelText, "Interpreter", "latex");
    ylabel(ax, "$F_{in}$ [L/h]", "Interpreter", "latex");
    title(ax, sprintf('$\\mathbf{%s}$', char(casePanelLabels(c))), "Interpreter", "latex");
    ax.TitleHorizontalAlignment = "left";
    ax.TickLabelInterpreter = "latex";
    set(ax, "FontSize", fontSize);
    grid(ax, "off");
    box(ax, "off");
    if c == 1
        legend(ax, "Interpreter", "latex", "Location", "best", "Box", "off");
    end
end

outStem = fullfile(graphicsDir, "setpoint_schedule_selected_benchmark_fin_only_" + scenario);
exportgraphics(fig, outStem + ".png", "Resolution", 300);
exportgraphics(fig, outStem + ".pdf", "ContentType", "vector");
fprintf("Saved: %s\n", outStem + ".png");
fprintf("Saved: %s\n", outStem + ".pdf");
end

function markerIdx = compute_fin_peak_marker_indices(t, fin)
markerIdx = [];
if isempty(t) || isempty(fin)
    return
end

timeWindows = [9, 11; 19, 21];
for k = 1:size(timeWindows, 1)
    mask = t >= timeWindows(k, 1) & t <= timeWindows(k, 2);
    if ~any(mask)
        continue
    end
    localIdx = find(mask);
    [~, iMax] = max(fin(localIdx));
    markerIdx(end + 1) = localIdx(iMax); %#ok<AGROW>
end
markerIdx = unique(markerIdx);
end

function plot_black_controller_inputs(plotData, scenario, graphicsDir, blackControllerId)
plotColors = nature_methods_colors(3);
inputLabels = ["$F_{in}$ [L/h]","$F_m$ [L/h]","$F_{out}$ [L/h]"];
plot_black_controller_inputs_scaled(plotData, scenario, graphicsDir, blackControllerId, ...
    1, inputLabels, "setpoint_schedule_selected_controller_inputs_" + scenario, [0 0.4], ...
    "Selected and benchmark controller inputs", plotColors(1, :), plotColors(2, :));
end

function plot_black_controller_inputs_mls(plotData, scenario, graphicsDir, blackControllerId)
plotColors = nature_methods_colors(3);
inputLabels = ["$F_{in}$ [mL/s]","$F_m$ [mL/s]","$F_{out}$ [mL/s]"];
plot_black_controller_inputs_scaled(plotData, scenario, graphicsDir, blackControllerId, ...
    1000 / 3600, inputLabels, "setpoint_schedule_selected_controller_inputs_mls_" + scenario, ...
    [0 0.4] * (1000 / 3600), "Selected and benchmark controller inputs (mL/s)", ...
    plotColors(1, :), plotColors(2, :));
end

function plot_black_controller_inputs_scaled(plotData, scenario, graphicsDir, blackControllerId, unitScale, inputLabels, outBaseName, yLimits, figName, selectedColor, benchmarkColor)
selectedIdx = find(string({plotData.id}) == string(blackControllerId), 1, "first");
if isempty(selectedIdx)
    warning("Selected controller not found for input plot: %s", blackControllerId);
    return
end
benchmarkIdx = find([plotData.is_benchmark], 1, "first");

fontSize = 14;
casePanelLabels = ["a","b"];

fig = figure("Color", "w", "Name", figName);
tiledlayout(fig, 2, 3, "TileSpacing", "compact", "Padding", "compact");
set(fig, "Position", [120 120 1400 760]);

for c = 1:2
    if numel(plotData(selectedIdx).cases) < c
        continue
    end
    tSelected = plotData(selectedIdx).T;
    USelected = plotData(selectedIdx).cases(c).U;
    if isempty(tSelected) || isempty(USelected)
        continue
    end
    USelected = double(USelected) * unitScale;
    nu = min(size(USelected, 2), 3);

    hasBenchmark = ~isempty(benchmarkIdx) && numel(plotData(benchmarkIdx).cases) >= c;
    if hasBenchmark
        tBenchmark = plotData(benchmarkIdx).T;
        UBenchmark = plotData(benchmarkIdx).cases(c).U;
        hasBenchmark = ~isempty(tBenchmark) && ~isempty(UBenchmark);
        if hasBenchmark
            UBenchmark = double(UBenchmark) * unitScale;
        end
    end

    for u = 1:nu
        ax = nexttile((c - 1) * 3 + u); hold(ax, "on");
        if hasBenchmark && size(UBenchmark, 2) >= u
            plot(ax, tBenchmark, UBenchmark(:, u), ":", "LineWidth", 1.6, ...
                "Color", benchmarkColor, "DisplayName", "Benchmark");
        end
        plot(ax, tSelected, USelected(:, u), "-", "LineWidth", 1.6, ...
            "Color", selectedColor, "DisplayName", "Selected");
        ylim(ax, yLimits)
        if c < 2
            xlabelText = "";
        else
            xlabelText = "Time [h]";
        end
        xlabel(ax, xlabelText, "Interpreter", "latex");
        ylabel(ax, inputLabels(u), "Interpreter", "latex");
        if u == 1
            title(ax, "$\mathbf{" + casePanelLabels(c) + "}$", "Interpreter", "latex");
            ax.TitleHorizontalAlignment = "left";
        else
            title(ax, "");
        end
        if c == 1 && u == 1
            legend(ax, "Interpreter", "latex", "Location", "best", "Box", "off");
        end
        ax.TickLabelInterpreter = "latex";
        set(ax, "FontSize", fontSize);
        grid(ax, "off");
        box(ax, "off");
    end
end

outStem = fullfile(graphicsDir, outBaseName);
exportgraphics(fig, outStem + ".png", "Resolution", 300);
exportgraphics(fig, outStem + ".pdf", "ContentType", "vector");
fprintf("Saved: %s\n", outStem + ".png");
fprintf("Saved: %s\n", outStem + ".pdf");
end

function tf = should_plot_schedule_controller(controllerPlotData, caseIdx, paretoIds, selectedControllerId)
tf = false;
if numel(controllerPlotData.cases) < caseIdx
    return
end

controllerId = string(controllerPlotData.id);
tf = controllerPlotData.is_benchmark || ...
    controllerId == string(selectedControllerId) || ...
    ismember(controllerId, paretoIds);
end

function [tCommon, yspMedian] = compute_case_median_setpoint(plotData, caseIdx, stateIdx, paretoIds, selectedControllerId)
tCommon = [];
yspMedian = [];
yspMatrix = [];

for i = 1:numel(plotData)
    if ~should_plot_schedule_controller(plotData(i), caseIdx, paretoIds, selectedControllerId)
        continue
    end

    t = plotData(i).T;
    Ysp = plotData(i).cases(caseIdx).Ysp;
    if isempty(t) || size(Ysp, 2) < stateIdx
        continue
    end

    yspCol = double(Ysp(:, stateIdx));
    if isempty(tCommon)
        tCommon = double(t(:));
        yspMatrix = yspCol;
        continue
    end

    if numel(t) ~= numel(tCommon) || any(abs(double(t(:)) - tCommon) > 1e-12)
        continue
    end

    yspMatrix(:, end + 1) = yspCol; %#ok<AGROW>
end

if ~isempty(yspMatrix)
    yspMedian = median(yspMatrix, 2, "omitnan");
end
end

function blendedColor = blend_color_with_white(baseColor, alphaVal)
baseColor = double(baseColor(:)).';
alphaVal = double(alphaVal);
blendedColor = alphaVal * baseColor + (1 - alphaVal) * [1 1 1];
end

function plot_modified_unmodified_boxplots(records, scenario, graphicsDir)
if isempty(records)
    return
end

ids = string(records.controller_id);
isModified = endsWith(ids, "_modified");
isBenchmark = logical(records.is_benchmark);
unmodifiedIds = ids(~isModified & ~isBenchmark);
modifiedIds = ids(isModified);
baseIds = erase(modifiedIds, "_modified");
hasBase = ismember(baseIds, ids);

missingModifiedIds = setdiff(unmodifiedIds, baseIds, "stable");
if ~isempty(missingModifiedIds)
    warning("Unmodified controllers without modified counterpart are excluded from the paired boxplots: %s", ...
        strjoin(cellstr(missingModifiedIds), ", "));
end

modifiedIds = modifiedIds(hasBase);
baseIds = baseIds(hasBase);
if isempty(baseIds)
    warning("No modified/unmodified controller pairs found for boxplot comparison.");
    return
end

baseVals = nan(numel(baseIds), 3);
modVals = nan(numel(baseIds), 3);
metricVars = ["SSE_total", "SSdU_total", "wall_time_s"];
metricTitles = ["SSE", "SSdU", "Wall Time"];
metricYLabels = ["SSE", "SSdU", "Wall time [s]"];

for i = 1:numel(baseIds)
    baseRow = records(ids == baseIds(i), :);
    modRow = records(ids == modifiedIds(i), :);
    if isempty(baseRow) || isempty(modRow)
        continue
    end
    baseVals(i, :) = [double(baseRow.SSE_total(1)), double(baseRow.SSdU_total(1)), double(baseRow.wall_time_s(1))];
    modVals(i, :) = [double(modRow.SSE_total(1)), double(modRow.SSdU_total(1)), double(modRow.wall_time_s(1))];
end

fig = figure("Color", "w", "Name", "Modified vs unmodified controller boxplots");
tiledlayout(fig, 1, 3, "TileSpacing", "compact", "Padding", "compact");
set(fig, "Position", [140 140 1300 420]);

for m = 1:3
    ax = nexttile; hold(ax, "on");
    values = [baseVals(:, m); modVals(:, m)];
    groups = [repmat("Unmodified", size(baseVals, 1), 1); repmat("Modified", size(modVals, 1), 1)];
    validMask = isfinite(values);
    if ~any(validMask)
        title(ax, metricTitles(m));
        text(ax, 0.5, 0.5, "No valid paired data", "HorizontalAlignment", "center");
        axis(ax, "off");
        continue
    end
    boxplot(ax, values(validMask), cellstr(groups(validMask)));
    grid(ax, "on");
    box(ax, "off");
    title(ax, metricTitles(m));
    ylabel(ax, metricYLabels(m));
end

outStem = fullfile(graphicsDir, "setpoint_schedule_modified_vs_unmodified_boxplots_" + scenario);
exportgraphics(fig, outStem + ".png", "Resolution", 300);
exportgraphics(fig, outStem + ".pdf", "ContentType", "vector");
fprintf("Saved: %s\n", outStem + ".png");
fprintf("Saved: %s\n", outStem + ".pdf");
end
