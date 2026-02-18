% RESULTSSANDBOX
% Project post-processing script for MFBO-based NMPC tuning experiments.
%
% Purpose
% - Load optimisation results from two experiment runs (`run1`, `run2`).
% - Reconstruct derived tuning/diagnostic variables from `results.csv`.
% - Reproduce core analysis plots per run.
% - Compare Pareto frontiers across runs and compute a final frontier from
%   all points combined.
%
% Expected input data
% - `results/run1/results.csv`
% - `results/run2/results.csv`
% - Required columns:
%   `timestamp`, `SSE`, `SSdU`, `J`, `runtime_s`, `theta_1..theta_12`
%
% Main derived data
% - `T`: full processed table for a run, including:
%   `iteration`, `runtime_min`, `f`, `theta_p`, `theta_m`, `p`, `m`,
%   and mapped weights (`Q_x*`, `R_u_x*`, `R_du_x*`).
% - `Tp`: Pareto-optimal subset of `T` using minimisation of
%   `J_track = SSE` and `J_TV = SSdU`, sorted by `SSE` (descending).
% - `isPareto`: logical mask identifying Pareto points in `T`.
%


%% Dependencies
clear; close all; clc
addpath(genpath("dependencies"));

%% Paths
rootFolder = "results";
datasets = [
    struct("name","Case 1", "csvPath", fullfile(rootFolder, "run1", "results.csv"), "outDir", fullfile(rootFolder, "run1"));
    struct("name","Case 2", "csvPath", fullfile(rootFolder, "run2", "results.csv"), "outDir", fullfile(rootFolder, "run2"))
];

if ~isfolder(rootFolder)
    mkdir(rootFolder);
end

set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
fontSize = 14;

%% Colors
plotColors = good_colors(4);

figColors = figure("Color", "w");
axColors = axes(figColors); hold(axColors, "on");
xColor = 1:size(plotColors, 1);
for i = 1:size(plotColors, 1)
    scatter(axColors, xColor(i), 1, 1200, plotColors(i, :), "filled", ...
        "MarkerEdgeColor", "k", "LineWidth", 0.8);
end
xlim(axColors, [0.5, size(plotColors, 1) + 0.5]);
ylim(axColors, [0.7, 1.3]);
xticks(axColors, xColor);
yticks(axColors, []);
xlabel(axColors, "Color index");
set(axColors, "FontSize", fontSize);
grid(axColors, "off");
box(axColors, "off");

% Trimming
plotColors(4, :) = [];

%% 
allT = cell(numel(datasets), 1);
allTp = cell(numel(datasets), 1);
allPareto = cell(numel(datasets), 1);

for k = 1:numel(datasets)
    [T, Tp, isPareto] = load_results_table(datasets(k).csvPath);
    T = enrich_with_out_data(T, datasets(k).outDir);
    allT{k} = T;
    allTp{k} = Tp;
    allPareto{k} = isPareto;
    display_pareto_table(Tp);
    display_runtime_phase_summary(T, datasets(k).name, 20);
end

create_analysis_plots_side_by_side(allT, allPareto, datasets, rootFolder, fontSize, plotColors);
%% Combined Pareto
plot_combined_pareto_samples(allT{1}, allTp{1}, allT{2}, allTp{2}, fullfile(rootFolder, "pareto_samples_run1_run2.png"), fontSize, plotColors);
plot_cumulative_runtime_combined(allT, datasets, rootFolder, fontSize, plotColors);

% final_Tp.timestamp
% 
% ans = 
% 
%   11×1 string array
% 
%     "20260131_151035"
%     "20260201_111557"
%     "20260201_111928"
%     "20260201_192807"
%     "20260201_223337"
%     "20260201_232106"
%     "20260210_151703"
%     "20260210_171107"
%     "20260210_180826"
%     "20260211_122653"
%     "20260211_134235"


function [T, Tp, isPareto] = load_results_table(csvPath)
if nargin < 1 || strlength(string(csvPath)) == 0
    error("You must provide csvPath.");
end
if ~isfile(csvPath)
    error("CSV not found: %s", csvPath);
end

T = readtable(csvPath, "TextType", "string");

requiredCols = ["timestamp","SSE","SSdU","J","runtime_s"];
for c = requiredCols
    if ~ismember(c, string(T.Properties.VariableNames))
        error("Missing required column '%s' in %s", c, csvPath);
    end
end
for i = 1:12
    nm = "theta_" + i;
    if ~ismember(nm, string(T.Properties.VariableNames))
        error("Missing required column '%s' in %s", nm, csvPath);
    end
end

T.timestamp_dt = datetime(T.timestamp, "InputFormat", "yyyyMMdd_HHmmss");
T.iteration = (1:height(T)).';
T.runtime_min = T.runtime_s / 60;

T.f = double(T.theta_1);
T.theta_p = round(double(T.theta_2));
T.theta_m = round(double(T.theta_3));

T.q1_log10 = double(T.theta_4);
T.q2_log10 = double(T.theta_5);
T.q3_log10 = double(T.theta_6);

T.r_u1_log10 = double(T.theta_7);
T.r_u2_log10 = double(T.theta_8);
T.r_u3_log10 = double(T.theta_9);

T.r_du1_log10 = double(T.theta_10);
T.r_du2_log10 = double(T.theta_11);
T.r_du3_log10 = double(T.theta_12);

T.Q_x1 = 10.^T.q1_log10;
T.Q_x2 = 10.^T.q2_log10;
T.Q_x3 = 10.^T.q3_log10;

T.R_u_x1 = 10.^T.r_u1_log10;
T.R_u_x2 = 10.^T.r_u2_log10;
T.R_u_x3 = 10.^T.r_u3_log10;

T.R_du_x1 = 10.^T.r_du1_log10;
T.R_du_x2 = 10.^T.r_du2_log10;
T.R_du_x3 = 10.^T.r_du3_log10;

T.m = T.theta_m + 1;
T.p = T.theta_p + T.m;

isPareto = compute_pareto_mask(double(T.SSE), double(T.SSdU));

Tp = T(isPareto, :);
[~, ord_tune] = sort(double(Tp.SSE), "descend");
Tp = Tp(ord_tune, :);

Tp.runtime_per_f = NaN(height(Tp), 1);
if ismember("runtime_s", string(Tp.Properties.VariableNames)) && ismember("f", string(Tp.Properties.VariableNames))
    Tp.runtime_per_f = double(Tp.runtime_s) ./ max(double(Tp.f), eps);
elseif ismember("runtime_min", string(Tp.Properties.VariableNames)) && ismember("f", string(Tp.Properties.VariableNames))
    Tp.runtime_per_f = double(Tp.runtime_min) ./ max(double(Tp.f), eps);
end
end


function T = enrich_with_out_data(T, runFolder)
% Load simulation-level outputs (`out_*.mat`) linked by timestamp and
% append runtime metrics that are not present in CSV-only summaries.
n = height(T);
T.out_found = false(n,1);
T.out_runtime_s = NaN(n,1);
T.out_mean_case_runtime_s = NaN(n,1);
T.out_mean_step_runtime_ms = NaN(n,1);
T.out_total_steps = NaN(n,1);
T.out_case_count = NaN(n,1);

for i = 1:n
    ts = string(T.timestamp(i));
    outPath = fullfile(runFolder, "out_" + ts + ".mat");
    if ~isfile(outPath)
        continue;
    end

    S = load(outPath, "out");
    if ~isfield(S, "out")
        continue;
    end
    out = S.out;

    T.out_found(i) = true;
    if isfield(out, "runtime_s")
        T.out_runtime_s(i) = double(out.runtime_s);
    end

    if isfield(out, "case") && ~isempty(out.case)
        nCase = numel(out.case);
        caseRuntime = NaN(nCase,1);
        caseSteps = NaN(nCase,1);
        caseMeanStep = NaN(nCase,1);
        for c = 1:nCase
            ci = out.case(c);
            if isfield(ci, "runtime_s")
                caseRuntime(c) = double(ci.runtime_s);
            end
            if isfield(ci, "RUNTIME") && ~isempty(ci.RUNTIME)
                r = double(ci.RUNTIME(:));
                caseSteps(c) = numel(r);
                caseMeanStep(c) = mean(r, "omitnan");
            end
        end

        T.out_case_count(i) = nCase;
        T.out_mean_case_runtime_s(i) = mean(caseRuntime, "omitnan");
        T.out_total_steps(i) = sum(caseSteps, "omitnan");
        T.out_mean_step_runtime_ms(i) = 1000 * mean(caseMeanStep, "omitnan");
    end
end
end


function isPareto = compute_pareto_mask(J_track, J_TV)
n = numel(J_track);
isPareto = true(n,1);
for i = 1:n
    dominated = (J_track <= J_track(i)) & (J_TV <= J_TV(i)) & ...
                ((J_track < J_track(i)) | (J_TV < J_TV(i)));
    dominated(i) = false;
    if any(dominated)
        isPareto(i) = false;
    end
end
end


function create_analysis_plots(T, isPareto, outDir, fontSize, plotColors, runIdx)
if ~isfolder(outDir)
    mkdir(outDir);
end

outScatterPath = fullfile(outDir, "sse_vs_ssdu.png");
outRuntimePath = fullfile(outDir, "runtime_vs_iteration.png");
outCumRuntimePath = fullfile(outDir, "runtime_cumulative_vs_iteration.png");

J_track = double(T.SSE);
J_TV = double(T.SSdU);

fig1 = figure("Color", "w");
ax1 = axes(fig1); hold(ax1, "on");
scatter(ax1, J_TV, J_track, 80, double(T.J), "filled", "MarkerEdgeColor", "k", "LineWidth", 0.7);
scatter(ax1, J_TV(isPareto), J_track(isPareto), 150, "MarkerEdgeColor", "r", "MarkerFaceColor", "none", "LineWidth", 1.2);
set(ax1, "XScale", "log", "YScale", "log", "FontSize", fontSize);
xlabel(ax1, "$J_{\mathrm{TV}}$");
ylabel(ax1, "$J_{\mathrm{track}}$");
grid(ax1, "off");
box(ax1, "off");
cb = colorbar(ax1);
cb.Label.String = "$J$";
cb.Label.Interpreter = "latex";
cb.TickLabelInterpreter = "latex";
cb.FontSize = fontSize;
exportgraphics(fig1, outScatterPath, "Resolution", 300);

fig2 = figure("Color", "w");
ax2 = axes(fig2); hold(ax2, "on");
runtime_h = double(T.runtime_min) / 60;
yyaxis(ax2, "left");
plot(ax2, T.iteration, runtime_h, "-", "LineWidth", 2.0, "Color", plotColors(runIdx, :));
ax2.YColor = plotColors(runIdx, :);
xline(ax2, 20.5, "--", "LineWidth", 2.0);
xlabel(ax2, "$k$ (iteration)");
ylabel(ax2, "$t_{\mathrm{iter}}$ (h)");

yyaxis(ax2, "right");
plot(ax2, T.iteration, double(T.f), "o", "LineWidth", 2.0, "MarkerSize", 4, "Color", plotColors(3, :));
ax2.YColor = plotColors(3, :);
ylabel(ax2, "$z$");

xlim(ax2, [1, max(1, height(T))]);
set(ax2, "FontSize", fontSize);
grid(ax2, "off");
box(ax2, "off");
exportgraphics(fig2, outRuntimePath, "Resolution", 300);

fig3 = figure("Color", "w");
ax3 = axes(fig3); hold(ax3, "on");
plot(ax3, T.iteration, cumsum(runtime_h, "omitnan"), "-", "LineWidth", 2.0, "Color", plotColors(runIdx, :));
xline(ax3, 20.5, "--", "LineWidth", 2.0);
xlabel(ax3, "$k$ (iteration)");
ylabel(ax3, "$t_{\mathrm{run}}$ (h)");
xlim(ax3, [1, max(1, height(T))]);
set(ax3, "FontSize", fontSize);
grid(ax3, "off");
box(ax3, "off");
exportgraphics(fig3, outCumRuntimePath, "Resolution", 300);
end


function create_analysis_plots_side_by_side(allT, allPareto, datasets, outDir, fontSize, plotColors)
if numel(allT) < 2 || isempty(allT{1}) || isempty(allT{2})
    return
end
if ~isfolder(outDir)
    mkdir(outDir);
end

% SSE vs SSdU (side-by-side), color mapped by z; export multiple colormap variants.
cmapList = build_colormap_variants(256);
for ic = 1:numel(cmapList)
    fig1z = figure("Color", "w");
    tiledlayout(fig1z, 1, 2, "Padding", "compact", "TileSpacing", "compact");
    for k = 1:2
        T = allT{k};
        isPareto = allPareto{k};
        ax = nexttile; hold(ax, "on");
        scatter(ax, double(T.SSdU), double(T.SSE), 80, double(T.f), "filled", "MarkerEdgeColor", "k", "LineWidth", 0.7);
        hGuide = plot_pareto_continuum(ax, double(T.SSdU(isPareto)), double(T.SSE(isPareto)), ...
            plotColors(3, :), [1e-2, 1e2], [1e4, 1.3e5]);
        try
            hGuide.Color = [plotColors(3, :) 0.45];
        catch
            hGuide.Color = 0.55 * plotColors(3, :) + 0.45 * [1 1 1];
        end
        scatter(ax, double(T.SSdU(isPareto)), double(T.SSE(isPareto)), 170, plotColors(3,:), ...
            "o", "MarkerFaceColor", "none", "MarkerEdgeColor", plotColors(3,:), "LineWidth", 1.2);
        set(ax, "XScale", "log", "YScale", "log", "FontSize", fontSize);
        xlim(ax, [1e-2, 1e2]);
        ylim(ax, [1e4, 1.3e5]);
        colormap(ax, cmapList(ic).map);
        caxis(ax, [0, 1]);
        xlabel(ax, "$J_{\mathrm{TV}}$");
        ylabel(ax, "$J_{\mathrm{track}}$");
        title(ax, string(datasets(k).name), "Interpreter", "none");
        grid(ax, "off");
        box(ax, "off");
        cb = colorbar(ax);
        cb.Label.String = "$z$";
        cb.Label.Interpreter = "latex";
        cb.TickLabelInterpreter = "latex";
        cb.FontSize = fontSize;
    end
    sgtitle(fig1z, "Colormap: " + cmapList(ic).name, "Interpreter", "none");
    fileTag = regexprep(lower(char(cmapList(ic).name)), "[^a-z0-9]+", "_");
    exportgraphics(fig1z, fullfile(outDir, "sse_vs_ssdu_side_by_side_z_" + fileTag + ".png"), "Resolution", 300);
end

% Iteration runtime + z (side-by-side)
fig2 = figure("Color", "w");
tiledlayout(fig2, 1, 2, "Padding", "compact", "TileSpacing", "compact");
for k = 1:2
    T = allT{k};
    runtime_h = double(T.runtime_min) / 60;
    ax = nexttile; hold(ax, "on");
    yyaxis(ax, "left");
    plot(ax, T.iteration, runtime_h, "-", "LineWidth", 2.0, "Color", plotColors(k, :));
    ax.YColor = plotColors(k, :);
    ylim(ax, [0, 4]);
    xline(ax, 20.5, "--", "LineWidth", 2.0);
    ylabel(ax, "$t_{\mathrm{iter}}$ (h)");
    yyaxis(ax, "right");
    plot(ax, T.iteration, double(T.f), "o", "LineWidth", 2.0, "MarkerSize", 4, "Color", plotColors(3, :));
    ax.YColor = plotColors(3, :);
    ylim(ax, [0, 1]);
    ylabel(ax, "$z$");
    xlabel(ax, "$k$ (iteration)");
    xlim(ax, [1, max(1, height(T))]);
    title(ax, string(datasets(k).name), "Interpreter", "none");
    set(ax, "FontSize", fontSize);
    grid(ax, "off");
    box(ax, "off");
end
exportgraphics(fig2, fullfile(outDir, "runtime_vs_iteration_side_by_side.png"), "Resolution", 300);
end


function display_pareto_table(Tp)
tuningCols = ["timestamp","SSE","SSdU","J","p","m","f", ...
    "runtime_min","runtime_per_f", ...
    "Q_x1","Q_x2","Q_x3", ...
    "R_u_x1","R_u_x2","R_u_x3", ...
    "R_du_x1","R_du_x2","R_du_x3"];

tuningCols = tuningCols(ismember(tuningCols, string(Tp.Properties.VariableNames)));
ParetoTuningTable = Tp(:, tuningCols);

disp("Pareto frontier with tuning weights (Q, R, Rdu):");
disp(ParetoTuningTable);
end


function display_runtime_phase_summary(T, runName, optimizationStartIter)
if nargin < 3
    optimizationStartIter = 20; % last DOE iteration (1-based)
end

iter = double(T.iteration);
rtMin = double(T.runtime_min);

doeMask = iter <= optimizationStartIter;
optMask = iter > optimizationStartIter;

doeTime = sum(rtMin(doeMask), "omitnan");
optTime = sum(rtMin(optMask), "omitnan");
totalTime = sum(rtMin, "omitnan");

doeN = nnz(doeMask);
optN = nnz(optMask);
totalN = numel(iter);

if totalTime > 0
    doePct = 100 * doeTime / totalTime;
    optPct = 100 * optTime / totalTime;
else
    doePct = NaN;
    optPct = NaN;
end

phase = {'DOE'; 'Optimisation'; 'Total'};
iterations = [doeN; optN; totalN];
runtime_min = [doeTime; optTime; totalTime];
runtime_pct_total = [doePct; optPct; 100];

summaryTbl = table(phase, iterations, runtime_min, runtime_pct_total, ...
    'VariableNames', {'phase','iterations','runtime_min','runtime_pct_total'});

disp("Runtime phase summary - " + string(runName) + ":");
disp(summaryTbl);
end




function plot_combined_pareto_samples(T1, Tp1, T2, Tp2, outPath, fontSize, plotColors)
fig = figure("Color", "w");
ax = axes(fig); hold(ax, "on");

% All evaluated points as small grey dots
Tall = [T1; T2];
scatter(ax, double(Tall.SSdU), double(Tall.SSE), 18, [0.65 0.65 0.65], ...
    "filled", "MarkerEdgeColor", "none", "DisplayName", "All samples");

% One smooth curve through the final (combined) Pareto points.
finalMask = compute_pareto_mask(double(Tall.SSE), double(Tall.SSdU));
Tf = Tall(finalMask, :);
hCurve = plot_pareto_continuum(ax, double(Tf.SSdU), double(Tf.SSE), plotColors(3, :), [1e-2, 1e2], [1e4, 1.3e5]);
try
    hCurve.Color = [plotColors(3, :) 0.45];
catch
    hCurve.Color = 0.55 * plotColors(3, :) + 0.45 * [1 1 1];
end

% Pareto points per case (filled markers) drawn after curve to stay on top.
scatter(ax, double(Tp1.SSdU), double(Tp1.SSE), 40, plotColors(1,:), ...
    "o", "MarkerFaceColor", plotColors(1,:), "MarkerEdgeColor", plotColors(1,:), ...
    "LineWidth", 1.4);
scatter(ax, double(Tp2.SSdU), double(Tp2.SSE), 40, plotColors(2,:), ...
    "^", "MarkerFaceColor", plotColors(2,:), "MarkerEdgeColor", plotColors(2,:), ...
    "LineWidth", 1.4);

% Circle only final (combined) Pareto points with a visible gap.
scatter(ax, double(Tf.SSdU), double(Tf.SSE), 112, plotColors(3,:), ...
    "o", "MarkerFaceColor", "none", "MarkerEdgeColor", plotColors(3,:), ...
    "LineWidth", 1.1);

set(ax, "XScale", "log", "YScale", "log", "FontSize", fontSize);
xlim(ax, [1e-2, 2e0]);
ylim(ax, [1e4, 3.5e4]);
xlabel(ax, "$J_{\mathrm{TV}}$");
ylabel(ax, "$J_{\mathrm{track}}$");
grid(ax, "off");
box(ax, "off");

exportgraphics(fig, outPath, "Resolution", 300);
end


function h = plot_pareto_continuum(ax, x, y, curveColor, xBounds, yBounds)
[xSort, ord] = sort(x(:), "ascend");
ySort = y(ord);

if numel(xSort) < 3
    h = plot(ax, xSort, ySort, "-", "Color", curveColor, "LineWidth", 2.0);
    return;
end

% Interpolate in log-log domain for a smooth non-segmented visual guide.
lx = log10(xSort);
ly = log10(ySort);
lxqIn = linspace(min(lx), max(lx), 220);
lyqIn = pchip(lx, ly, lxqIn);
xq = (10.^lxqIn).';
yq = (10.^lyqIn).';

% Left extension: vertical at the leftmost Pareto point.
if nargin >= 6 && ~isempty(yBounds)
    yTop = max(yBounds);
    if yTop > yq(1)
        xq = [xq(1); xq];
        yq = [yTop; yq];
    end
end

% Right extension: horizontal from the rightmost Pareto point.
if nargin >= 5 && ~isempty(xBounds)
    xRight = max(xBounds);
    if xRight > xq(end)
        xq = [xq; xRight];
        yq = [yq; yq(end)];
    end
end

h = plot(ax, xq, yq, "-", "Color", curveColor, "LineWidth", 2.0);
end


function cmapList = build_colormap_variants(n)
% Curated colormaps: readable, publication-friendly, and colorblind-aware.
cmapList = struct("name", {}, "map", {});

% Three-color smooth map (dark blue -> teal -> yellow).
anchorsBlueTealYellow = [
    0.08 0.12 0.55
    0.13 0.55 0.57
    0.98 0.90 0.20
];
cmapList(end+1) = struct("name", "blue_teal_yellow", "map", interpolate_anchors(anchorsBlueTealYellow, n));

% Viridis-like (purple/blue -> green -> yellow), widely used in papers.
anchorsViridisLike = [
    0.27 0.00 0.33
    0.13 0.56 0.55
    0.99 0.91 0.14
];
cmapList(end+1) = struct("name", "viridis_like", "map", interpolate_anchors(anchorsViridisLike, n));
end


function cmap = interpolate_anchors(anchors, n)
tAnch = linspace(0, 1, size(anchors, 1));
t = linspace(0, 1, n);
cmap = zeros(n, 3);
for j = 1:3
    cmap(:, j) = interp1(tAnch, anchors(:, j), t, "pchip");
end
cmap = min(max(cmap, 0), 1);
end


function plot_available_paper_results(allT, allTp, datasets, outDir, optimizationStartIter, fontSize)
if nargin < 5
    optimizationStartIter = 20; % last DOE iteration (1-based)
end
if ~isfolder(outDir)
    mkdir(outDir);
end

nRuns = numel(datasets);
runNames = strings(nRuns,1);
for i = 1:nRuns
    runNames(i) = string(datasets(i).name);
end

%% Runtime phase comparison (DOE vs Optimisation)
% What: stacked runtime bars per run, split by DOE and optimisation phase.
% Why: shows where wall-clock time is spent and supports computational
% cost claims between the two tuning configurations.
doeHr = zeros(nRuns,1);
optHr = zeros(nRuns,1);
for i = 1:nRuns
    Ti = allT{i};
    iter = double(Ti.iteration);
    rt = double(Ti.runtime_min);
    doeHr(i) = sum(rt(iter <= optimizationStartIter), "omitnan") / 60;
    optHr(i) = sum(rt(iter > optimizationStartIter), "omitnan") / 60;
end

figRt = figure("Color","w");
axRt = axes(figRt); hold(axRt, "on");
bRt = bar(axRt, categorical(runNames), [doeHr, optHr], "stacked");
bRt(1).FaceColor = [0.35 0.35 0.35];
bRt(2).FaceColor = [0.10 0.10 0.10];
ylabel(axRt, "$t_{\mathrm{run}}$ (h)");
set(axRt, "FontSize", fontSize);
grid(axRt, "off");
box(axRt, "off");
exportgraphics(figRt, fullfile(outDir, "runtime_phase_comparison.png"), "Resolution", 300);

%% Runtime distribution across iterations
% What: per-run boxplot of per-iteration runtime.
% Why: highlights runtime variability/stability, not only average runtime.
rtAll = [];
grpAll = [];
for i = 1:nRuns
    % Standard outlier filtering: Tukey rule (outside [Q1-1.5*IQR, Q3+1.5*IQR]).
    rt = double(allT{i}.runtime_min) / 60;
    q1 = prctile(rt, 25);
    q3 = prctile(rt, 75);
    iqrVal = q3 - q1;
    lo = q1 - 1.5 * iqrVal;
    hi = q3 + 1.5 * iqrVal;
    keep = (rt >= lo) & (rt <= hi);

    rtAll = [rtAll; rt(keep)]; %#ok<AGROW>
    grpAll = [grpAll; repmat(runNames(i), nnz(keep), 1)]; %#ok<AGROW>
end

figBox = figure("Color","w");
axBox = axes(figBox); hold(axBox, "on");
boxplot(axBox, rtAll, grpAll);
ylabel(axBox, "$t_{\mathrm{iter}}$ (h)");
set(axBox, "FontSize", fontSize);
grid(axBox, "off");
box(axBox, "off");
exportgraphics(figBox, fullfile(outDir, "runtime_iteration_distribution.png"), "Resolution", 300);

%% Efficiency and objective-quality metrics
% What: bars for Pareto count, Pareto density, dominated fraction, and
% best objective values per run.
% Why: provides compact quantitative evidence of search efficiency and
% objective quality attained by each tuning strategy.
paretoCount = zeros(nRuns,1);
paretoDensity = zeros(nRuns,1);
dominatedFrac = zeros(nRuns,1);
bestSSE = zeros(nRuns,1);
bestSSdU = zeros(nRuns,1);

for i = 1:nRuns
    Ti = allT{i};
    Tpi = allTp{i};
    paretoCount(i) = height(Tpi);
    paretoDensity(i) = height(Tpi) / max(height(Ti), 1);
    dominatedFrac(i) = 1 - paretoDensity(i);
    bestSSE(i) = min(double(Ti.SSE));
    bestSSdU(i) = min(double(Ti.SSdU));
end

figEff = figure("Color","w");
tiledlayout(figEff, 2, 2, "Padding","compact", "TileSpacing","compact");

axE1 = nexttile; bar(axE1, categorical(runNames), paretoCount, "FaceColor", [0.15 0.15 0.15]);
ylabel(axE1, "$N_{\mathrm{Pareto}}$");
set(axE1, "FontSize", fontSize); grid(axE1, "off"); box(axE1, "off");

axE2 = nexttile; bar(axE2, categorical(runNames), paretoDensity, "FaceColor", [0.25 0.25 0.25]);
ylabel(axE2, "$\rho_{\mathrm{Pareto}}$");
set(axE2, "FontSize", fontSize); grid(axE2, "off"); box(axE2, "off");

axE3 = nexttile; bar(axE3, categorical(runNames), dominatedFrac, "FaceColor", [0.35 0.35 0.35]);
ylabel(axE3, "$\rho_{\mathrm{dom}}$");
set(axE3, "FontSize", fontSize); grid(axE3, "off"); box(axE3, "off");

axE4 = nexttile; hold(axE4, "on");
bar(axE4, categorical(runNames), [bestSSE, bestSSdU], "grouped");
ylabel(axE4, "$\min(J)$ components");
set(axE4, "FontSize", fontSize); grid(axE4, "off"); box(axE4, "off");

exportgraphics(figEff, fullfile(outDir, "run_efficiency_objective_metrics.png"), "Resolution", 300);

%% Pareto tuning-parameter distributions
% What: distributions of key tuning parameters on Pareto-optimal points.
% Why: reveals which parameter ranges are repeatedly selected as high-value
% solutions and helps interpret controller design tendencies.
keyVars = ["f","p","m","Q_x1","Q_x2","Q_x3","R_u_x1","R_u_x2","R_u_x3","R_du_x1","R_du_x2","R_du_x3"];
colors = lines(nRuns);

figPar = figure("Color","w");
tiledlayout(figPar, 3, 4, "Padding","compact", "TileSpacing","compact");
axFirst = [];
for v = 1:numel(keyVars)
    ax = nexttile; hold(ax, "on");
    if v == 1
        axFirst = ax;
    end
    varName = keyVars(v);
    for i = 1:nRuns
        if ismember(varName, string(allTp{i}.Properties.VariableNames))
            vals = double(allTp{i}.(varName));
            histogram(ax, vals, "DisplayStyle","stairs", "Normalization","probability", ...
                "NumBins", 12, "LineWidth",2.0, ...
                "EdgeColor", colors(i,:), "DisplayName", runNames(i));
        end
    end
    if startsWith(varName, "Q_") || startsWith(varName, "R_")
        set(ax, "XScale", "log");
    end
    xlabel(ax, "$" + strrep(varName, "_", "\_") + "$");
    set(ax, "FontSize", fontSize);
    grid(ax, "off");
    box(ax, "off");
end
exportgraphics(figPar, fullfile(outDir, "pareto_tuning_distributions.png"), "Resolution", 300);
end


function plot_out_runtime_vs_tuning(T, outDir, runName, fontSize)
hasOut = T.out_found & ~isnan(T.out_mean_step_runtime_ms);
if ~any(hasOut)
    return;
end

Tout = T(hasOut, :);

% What: mean NMPC step time per evaluated point versus key tuning
% parameters. Why: links computational burden directly to tuning decisions.
fig = figure("Color","w");
tiledlayout(fig, 2, 3, "Padding","compact", "TileSpacing","compact");

paramList = ["f","p","m","Q_x1","R_u_x1","R_du_x1"];
for i = 1:numel(paramList)
    ax = nexttile; hold(ax, "on");
    xv = double(Tout.(paramList(i)));
    yv = double(Tout.out_mean_step_runtime_ms);
    scatter(ax, xv, yv, 32, "filled", "MarkerFaceColor", [0.1 0.1 0.1], "MarkerFaceAlpha", 0.75);
    xlabel(ax, "$" + strrep(paramList(i), "_", "\_") + "$");
    ylabel(ax, "$\bar{t}_{\mathrm{step}}$ (ms)");
    if startsWith(paramList(i), "Q_") || startsWith(paramList(i), "R_")
        set(ax, "XScale", "log");
    end
    set(ax, "FontSize", fontSize);
    grid(ax, "off");
    box(ax, "off");
end

exportgraphics(fig, fullfile(outDir, "out_mean_step_runtime_vs_tuning_" + string(runName) + ".png"), "Resolution", 300);
end


function plot_out_runtime_vs_tuning_combined(allTout, datasets, outDir, fontSize)
Tall = table();
for i = 1:numel(allTout)
    Ti = allTout{i};
    Ti.run_id = repmat(string(datasets(i).name), height(Ti), 1);
    Tall = [Tall; Ti]; %#ok<AGROW>
end

hasOut = Tall.out_found & ~isnan(Tall.out_mean_step_runtime_ms);
if ~any(hasOut)
    return;
end
Tall = Tall(hasOut, :);

% What: aggregated simulation-step runtime against objective outcomes.
% Why: shows whether computationally expensive settings buy better trade-off.
fig = figure("Color","w");
tiledlayout(fig, 1, 2, "Padding","compact", "TileSpacing","compact");

ax1 = nexttile; hold(ax1, "on");
scatter(ax1, double(Tall.out_mean_step_runtime_ms), double(Tall.SSE), 26, "filled", ...
    "MarkerFaceColor", [0.2 0.2 0.2], "MarkerFaceAlpha", 0.7);
xlabel(ax1, "$\bar{t}_{\mathrm{step}}$ (ms)");
ylabel(ax1, "$J_{\mathrm{track}}$");
set(ax1, "YScale", "log", "FontSize", fontSize);
grid(ax1, "off"); box(ax1, "off");

ax2 = nexttile; hold(ax2, "on");
scatter(ax2, double(Tall.out_mean_step_runtime_ms), double(Tall.SSdU), 26, "filled", ...
    "MarkerFaceColor", [0.2 0.2 0.2], "MarkerFaceAlpha", 0.7);
xlabel(ax2, "$\bar{t}_{\mathrm{step}}$ (ms)");
ylabel(ax2, "$J_{\mathrm{TV}}$");
set(ax2, "YScale", "log", "FontSize", fontSize);
grid(ax2, "off"); box(ax2, "off");

exportgraphics(fig, fullfile(outDir, "out_mean_step_runtime_vs_objectives.png"), "Resolution", 300);
end


function plot_cumulative_runtime_combined(allT, datasets, outDir, fontSize, plotColors)
if numel(allT) < 2 || isempty(allT{1}) || isempty(allT{2})
    return
end

fig = figure("Color","w");
ax = axes(fig); hold(ax, "on");

runtime_h_1 = double(allT{1}.runtime_min) / 60;
runtime_h_2 = double(allT{2}.runtime_min) / 60;

plot(ax, double(allT{1}.iteration), cumsum(runtime_h_1, "omitnan"), "-", "LineWidth", 2.0, ...
    "Color", plotColors(1,:), ...
    "DisplayName", string(datasets(1).name));
plot(ax, double(allT{2}.iteration), cumsum(runtime_h_2, "omitnan"), "-.", "LineWidth", 2.0, ...
    "Color", plotColors(2,:), ...
    "DisplayName", string(datasets(2).name));
xline(ax, 20.5, "--", "LineWidth", 2.0);

xlabel(ax, "$k$ (iteration)");
ylabel(ax, "$t_{\mathrm{run}}$ (h)");
xlim(ax, [1, max(1, max(height(allT{1}), height(allT{2})))]);
set(ax, "FontSize", fontSize);
grid(ax, "off");
box(ax, "off");

exportgraphics(fig, fullfile(outDir, "runtime_cumulative_run1_run2.png"), "Resolution", 300);
end
