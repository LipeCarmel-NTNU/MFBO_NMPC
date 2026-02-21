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
scriptDir = fileparts(mfilename("fullpath"));
projectRoot = fileparts(scriptDir);
addpath(genpath(fullfile(projectRoot, "dependencies")));

%% Paths
rootFolder = fullfile(projectRoot, "results");
graphicsFolder = fullfile(rootFolder, "graphical_results");
numericalFolder = fullfile(rootFolder, "numerical results");
datasets = [
    struct("name","Case 1", "csvPath", fullfile(rootFolder, "run1", "results.csv"), "outDir", fullfile(rootFolder, "run1"));
    struct("name","Case 2", "csvPath", fullfile(rootFolder, "run2", "results.csv"), "outDir", fullfile(rootFolder, "run2"))
];

if ~isfolder(rootFolder)
    mkdir(rootFolder);
end
if ~isfolder(graphicsFolder)
    mkdir(graphicsFolder);
end
if ~isfolder(numericalFolder)
    mkdir(numericalFolder);
end

set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
fontSize = 18;

%% Colors
plotColors = nature_methods_colors(3); % Blue, Vermillion, Orange

% figColors = figure("Color", "w");
% axColors = axes(figColors); hold(axColors, "on");
% xColor = 1:size(plotColors, 1);
% for i = 1:size(plotColors, 1)
%     scatter(axColors, xColor(i), 1, 1200, plotColors(i, :), "filled", ...
%         "MarkerEdgeColor", "k", "LineWidth", 0.8);
% end
% xlim(axColors, [0.5, size(plotColors, 1) + 0.5]);
% ylim(axColors, [0.7, 1.3]);
% xticks(axColors, xColor);
% yticks(axColors, []);
% xlabel(axColors, "Color index");
% set(axColors, "FontSize", fontSize);
% grid(axColors, "off");
% box(axColors, "off");

% Palette entries: Case 1 (Blue), Case 2 (Vermillion), z/guideline (Orange).

%% 
allT = cell(numel(datasets), 1);
allTp = cell(numel(datasets), 1);
allPareto = cell(numel(datasets), 1);
runtimeSummaryTables = cell(numel(datasets), 1);
optStartIter = 20;

for k = 1:numel(datasets)
    [T, Tp, isPareto] = load_results_table(datasets(k).csvPath);
    T = enrich_with_out_data(T, datasets(k).outDir);
    allT{k} = T;
    allTp{k} = Tp;
    allPareto{k} = isPareto;
    display_pareto_table(Tp);
    runtimeSummaryTables{k} = display_runtime_phase_summary(T, datasets(k).name, optStartIter);
end

TallCombined = [allT{1}; allT{2}];
runtimeSummaryCombined = display_runtime_phase_summary(TallCombined, "Case 1 + Case 2 (combined)", optStartIter);
write_runtime_and_parameter_summary(allT, allTp, datasets, runtimeSummaryTables, runtimeSummaryCombined, numericalFolder, optStartIter);

create_analysis_plots_side_by_side(allT, allPareto, datasets, graphicsFolder, fontSize, plotColors);
%% Combined Pareto
plot_combined_pareto_samples(allT{1}, allTp{1}, allT{2}, allTp{2}, fullfile(graphicsFolder, "pareto_samples_run1_run2.png"), fontSize, plotColors);
plot_cumulative_runtime_combined(allT, datasets, graphicsFolder, fontSize, plotColors);


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


function create_analysis_plots_side_by_side(allT, allPareto, datasets, outDir, fontSize, plotColors)
if numel(allT) < 2 || isempty(allT{1}) || isempty(allT{2})
    return
end
if ~isfolder(outDir)
    mkdir(outDir);
end

% SSE vs SSdU (side-by-side), color mapped by z.
fig1z = figure("Color", "w");
tiledlayout(fig1z, 1, 2, "Padding", "compact", "TileSpacing", "compact");
viridisMap = make_viridis_like(256);
panelLabels = ["a", "b"];
for k = 1:2
    T = allT{k};
    isPareto = allPareto{k};
    ax = nexttile; hold(ax, "on");
    hGuide = plot_pareto_continuum(ax, double(T.SSdU(isPareto)), double(T.SSE(isPareto)), ...
        plotColors(3, :), [1e-2, 1e2], [1e4, 1.3e5]);
    try
        hGuide.Color = [plotColors(3, :) 0.45];
    catch
        hGuide.Color = 0.55 * plotColors(3, :) + 0.45 * [1 1 1];
    end
    scatter(ax, double(T.SSdU), double(T.SSE), 80, double(T.f), "filled", "MarkerEdgeColor", "k", "LineWidth", 0.7);
    scatter(ax, double(T.SSdU(isPareto)), double(T.SSE(isPareto)), 170, plotColors(3,:), ...
        "o", "MarkerFaceColor", "none", "MarkerEdgeColor", plotColors(3,:), "LineWidth", 1.2);
    set(ax, "XScale", "log", "YScale", "log", "FontSize", fontSize);
    xlim(ax, [1e-2, 1e2]);
    ylim(ax, [1e4, 1.3e5]);
    colormap(ax, viridisMap);
    caxis(ax, [0, 1]);
    xlabel(ax, "$J_{\mathrm{TV}}$");
    ylabel(ax, "$J_{\mathrm{track}}$");
    title(ax, "$\mathbf{" + panelLabels(k) + "}$", "Interpreter", "latex");
    ax.TitleHorizontalAlignment = "left";
    grid(ax, "off");
    box(ax, "off");
    cb = colorbar(ax);
    cb.Label.String = "$z$ (dimensionless)";
    cb.Label.Interpreter = "latex";
    cb.TickLabelInterpreter = "latex";
    cb.FontSize = fontSize;
end
save_plot_outputs(fig1z, fullfile(outDir, "sse_vs_ssdu_side_by_side_z.png"), fontSize, 1200, 460);

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
    ylabel(ax, "$z$ (dimensionless)");
    xlabel(ax, "$k$ (iteration)");
    xlim(ax, [1, max(1, height(T))]);
    title(ax, "$\mathbf{" + panelLabels(k) + "}$", "Interpreter", "latex");
    ax.TitleHorizontalAlignment = "left";
    set(ax, "FontSize", fontSize);
    grid(ax, "off");
    box(ax, "off");
    axes(ax);
    yyaxis(ax, "left");
    format_tick(0, 1);
    yyaxis(ax, "right");
    format_tick(0, 1);
end
save_plot_outputs(fig2, fullfile(outDir, "runtime_vs_iteration_side_by_side.png"), fontSize, 1200, 460);

% z, N_p, N_c by iteration (3x2)
fig3 = figure("Color", "w");
tiledlayout(fig3, 3, 2, "Padding", "compact", "TileSpacing", "compact");
panelLabels3x2 = ["a", "b", "c", "d", "e", "f"];

for k = 1:2
    T = allT{k};
    iter = double(T.iteration);

    % Row 1: z by k
    ax1 = nexttile((1 - 1) * 2 + k); hold(ax1, "on");
    plot(ax1, iter, double(T.f), "-", "LineWidth", 2.0, "Color", plotColors(3, :));
    plot(ax1, iter, double(T.f), "o", "MarkerSize", 3.5, "Color", plotColors(3, :));
    ylim(ax1, [0, 1]);
    xlim(ax1, [1, max(1, height(T))]);
    ylabel(ax1, "$z$");
    xlabel(ax1, "");
    title(ax1, "$\mathbf{" + panelLabels3x2((1 - 1) * 2 + k) + "}$", "Interpreter", "latex");
    ax1.TitleHorizontalAlignment = "left";
    set(ax1, "FontSize", fontSize);
    grid(ax1, "off"); box(ax1, "off");
    format_tick(0, 1);

    % Row 2: N_p by k
    ax2 = nexttile((2 - 1) * 2 + k); hold(ax2, "on");
    plot(ax2, iter, double(T.p), "-", "LineWidth", 2.0, "Color", plotColors(k, :));
    plot(ax2, iter, double(T.p), "o", "MarkerSize", 3.5, "Color", plotColors(k, :));
    xlim(ax2, [1, max(1, height(T))]);
    yminP = min(double(T.p), [], "omitnan");
    ymaxP = max(double(T.p), [], "omitnan");
    if isfinite(yminP) && isfinite(ymaxP)
        ylim(ax2, [max(0, floor(yminP) - 1), ceil(ymaxP) + 1]);
    end
    ylabel(ax2, "$N_p$");
    xlabel(ax2, "");
    title(ax2, "$\mathbf{" + panelLabels3x2((2 - 1) * 2 + k) + "}$", "Interpreter", "latex");
    ax2.TitleHorizontalAlignment = "left";
    set(ax2, "FontSize", fontSize);
    grid(ax2, "off"); box(ax2, "off");
    format_tick(0, 0);

    % Row 3: N_c by k (stored as m)
    ax3 = nexttile((3 - 1) * 2 + k); hold(ax3, "on");
    plot(ax3, iter, double(T.m), "-", "LineWidth", 2.0, "Color", plotColors(k, :));
    plot(ax3, iter, double(T.m), "o", "MarkerSize", 3.5, "Color", plotColors(k, :));
    xlim(ax3, [1, max(1, height(T))]);
    yminM = min(double(T.m), [], "omitnan");
    ymaxM = max(double(T.m), [], "omitnan");
    if isfinite(yminM) && isfinite(ymaxM)
        ylim(ax3, [max(0, floor(yminM) - 1), ceil(ymaxM) + 1]);
    end
    ylabel(ax3, "$N_c$");
    xlabel(ax3, "$k$ (iteration)");
    title(ax3, "$\mathbf{" + panelLabels3x2((3 - 1) * 2 + k) + "}$", "Interpreter", "latex");
    ax3.TitleHorizontalAlignment = "left";
    set(ax3, "FontSize", fontSize);
    grid(ax3, "off"); box(ax3, "off");
    format_tick(0, 0);
end

save_plot_outputs(fig3, fullfile(outDir, "z_np_nc_vs_iteration_3x2.png"), fontSize, 1200, 980);

% N_c by N_p density maps (side-by-side)
fig4 = figure("Color", "w");
tiledlayout(fig4, 1, 2, "Padding", "compact", "TileSpacing", "compact");
panelLabelsNcNp = ["a", "b"];

for k = 1:2
    T = allT{k};
    p = double(T.p);
    m = double(T.m); % N_c
    valid = isfinite(p) & isfinite(m);
    p = p(valid);
    m = m(valid);

    ax = nexttile; hold(ax, "on");
    if isempty(p)
        text(ax, 0.5, 0.5, "No valid $(N_p, N_c)$ data", ...
            "Units", "normalized", "HorizontalAlignment", "center", "Interpreter", "latex");
        axis(ax, "off");
        continue
    end

    pMin = floor(min(p)); pMax = ceil(max(p));
    mMin = floor(min(m)); mMax = ceil(max(m));
    xEdges = (pMin - 0.5):(pMax + 0.5);
    yEdges = (mMin - 0.5):(mMax + 0.5);
    counts = histcounts2(p, m, xEdges, yEdges);

    xCenters = xEdges(1:end-1) + 0.5;
    yCenters = yEdges(1:end-1) + 0.5;
    imagesc(ax, xCenters, yCenters, counts.');
    set(ax, "YDir", "normal");
    colormap(ax, parula(256));
    cb = colorbar(ax);
    cb.Label.String = "Count";
    cb.Label.Interpreter = "latex";
    cb.TickLabelInterpreter = "latex";
    cb.FontSize = fontSize;

    % Overlay points to preserve discrete-location visibility.
    scatter(ax, p, m, 26, "k", "filled", "MarkerFaceAlpha", 0.35, "MarkerEdgeColor", "none");

    xlabel(ax, "$N_p$");
    ylabel(ax, "$N_c$");
    title(ax, "$\mathbf{" + panelLabelsNcNp(k) + "}$", "Interpreter", "latex");
    ax.TitleHorizontalAlignment = "left";
    set(ax, "FontSize", fontSize);
    xlim(ax, [pMin - 0.5, pMax + 0.5]);
    ylim(ax, [mMin - 0.5, mMax + 0.5]);
    xticks(ax, pMin:pMax);
    yticks(ax, mMin:mMax);
    grid(ax, "off");
    box(ax, "off");
end

save_plot_outputs(fig4, fullfile(outDir, "nc_by_np_density_side_by_side.png"), fontSize, 1100, 470);
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


function summaryTbl = display_runtime_phase_summary(T, runName, optimizationStartIter)
arguments
    T table
    runName
    optimizationStartIter (1,1) double = 20 % last DOE iteration (1-based)
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

% All evaluated points as small black dots
Tall = [T1; T2];
scatter(ax, double(Tall.SSdU), double(Tall.SSE), 18, [0 0 0], ...
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
axes(ax);
format_tick(1, 1);

save_plot_outputs(fig, outPath, fontSize, 920, 520);
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


function cmap = make_viridis_like(n)
anchorsViridisLike = [
    0.27 0.00 0.33
    0.13 0.56 0.55
    0.99 0.91 0.14
];
cmap = interpolate_anchors(anchorsViridisLike, n);
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


function write_runtime_and_parameter_summary(allT, allTp, datasets, runSummaryTables, combinedSummaryTable, outDir, optimizationStartIter)
outPath = fullfile(outDir, "resultssandbox_runtime_and_params.txt");
fid = fopen(outPath, "w");
if fid == -1
    warning("Unable to write numerical summary: %s", outPath);
    return
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

for k = 1:numel(datasets)
    fprintf(fid, "Runtime phase summary - %s:\n", char(string(datasets(k).name)));
    write_runtime_table(fid, runSummaryTables{k});
    fprintf(fid, "\n");
end

fprintf(fid, "Runtime phase summary - Case 1 + Case 2 (combined):\n");
write_runtime_table(fid, combinedSummaryTable);
fprintf(fid, "\n");

meanZCase = NaN(numel(allT), 1);
for k = 1:numel(allT)
    meanZCase(k) = compute_mean_z_optimization(allT{k}, optimizationStartIter);
end
Tall = vertcat(allT{:});
meanZCombined = compute_mean_z_optimization(Tall, optimizationStartIter);

fprintf("Mean Fidelity z during optimization iterations:\n");
for k = 1:numel(datasets)
    fprintf("  %s: %.6g\n", string(datasets(k).name), meanZCase(k));
end
fprintf("  Case 1 + Case 2 (combined): %.6g\n", meanZCombined);

fprintf(fid, "Mean Fidelity z during optimization iterations:\n");
for k = 1:numel(datasets)
    fprintf(fid, "  %s: %.6g\n", char(string(datasets(k).name)), meanZCase(k));
end
fprintf(fid, "  Case 1 + Case 2 (combined): %.6g\n\n", meanZCombined);

for k = 1:numel(allTp)
    T = allTp{k};
    validM = isfinite(T.m);
    validP = isfinite(T.p);
    nRows = height(T);
    countM1 = nnz(T.m(validM) == 1);
    countP1 = nnz(T.p(validP) == 1);

    fprintf(fid, "%s Pareto counts:\n", char(string(datasets(k).name)));
    fprintf(fid, "  m = 1: %d/%d\n", countM1, nRows);
    fprintf(fid, "  p = 1: %d/%d\n", countP1, nRows);

    ssdU = double(T.SSdU);
    minSSdU = min(ssdU, [], "omitnan");
    topMask = abs(ssdU - minSSdU) <= 10 * eps(max(1, abs(minSSdU)));
    topCount = nnz(topMask);
    topPM1 = nnz(topMask & T.p == 1 & T.m == 1);
    fprintf(fid, "  Top (lowest SSdU) rows: %d\n", topCount);
    fprintf(fid, "  Top (lowest SSdU) with p = 1 and m = 1: %d\n\n", topPM1);
end
end


function meanZ = compute_mean_z_optimization(T, optimizationStartIter)
if isempty(T) || ~ismember("f", string(T.Properties.VariableNames))
    meanZ = NaN;
    return
end
mask = double(T.iteration) > optimizationStartIter;
meanZ = mean(double(T.f(mask)), "omitnan");
end


function write_runtime_table(fid, summaryTbl)
for i = 1:height(summaryTbl)
    fprintf(fid, "  %s | iterations=%d | runtime_min=%.6g | runtime_pct_total=%.6g\n", ...
        char(string(summaryTbl.phase{i})), summaryTbl.iterations(i), ...
        summaryTbl.runtime_min(i), summaryTbl.runtime_pct_total(i));
end
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

save_plot_outputs(fig, fullfile(outDir, "runtime_cumulative_run1_run2.png"), fontSize, 920, 520);
end


function save_plot_outputs(figHandle, pngPath, fontSize, figWidthPx, figHeightPx)
arguments
    figHandle
    pngPath
    fontSize (1,1) double = 14
    figWidthPx (1,1) double = 900
    figHeightPx (1,1) double = 500
end

% Apply consistent figure size and font before exporting.
figure(figHandle);
set_fig_size(figWidthPx, figHeightPx);
set_font_size(fontSize);
exportgraphics(figHandle, pngPath, "Resolution", 300);

% Save a vector PDF with the same base filename.
[folderPath, fileStem] = fileparts(pngPath);
pdfPath = fullfile(folderPath, strcat(fileStem, ".pdf"));
save_figure(pdfPath, NaN, false);
end
