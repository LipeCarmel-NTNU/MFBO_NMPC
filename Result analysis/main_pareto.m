function main_pareto()
% MAIN_PARETO
% Post-processing hub for MFBO-NMPC optimization outcomes.
%
% Intent
% - Provide one coherent analysis view of both optimization runs.
% - Separate exploration setup time (DOE) from true optimization behavior.
% - Summarize the tradeoff between tracking quality and control variation.
% - Produce publication-ready plots + numeric tables from one entry point.
%
% Data contract
% - Inputs: results/run1/results.csv and results/run2/results.csv
% - Required optimization columns:
%   timestamp, SSE, SSdU, J, runtime_s, theta_1..theta_12
%
% Core policy
% - DOE exclusion is mandatory for optimization conclusions:
%   only iterations k > 20 are eligible for Pareto/optimal-controller
%   analysis and Pareto scatter views.
% - Built-in guards raise errors if DOE points leak into Pareto outputs.
%
% Main outputs
% - Figures in results/graphical_results (PNG + PDF).
% - Numerical summaries in results/numerical results.
% - Console diagnostics for runtime split and Pareto composition.
%


%% Dependencies
close all; clc
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
fontSize = 20;

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
allTopt = cell(numel(datasets), 1);
allTp = cell(numel(datasets), 1);
allPareto = cell(numel(datasets), 1);
runtimeSummaryTables = cell(numel(datasets), 1);
optStartIter = 20;

for k = 1:numel(datasets)
    [T, Topt, Tp, isPareto] = load_results_table(datasets(k).csvPath, optStartIter);
    validate_no_doe_in_pareto(Tp, optStartIter, datasets(k).name);
    T = enrich_with_out_data(T, datasets(k).outDir);
    allT{k} = T;
    allTopt{k} = Topt;
    allTp{k} = Tp;
    allPareto{k} = isPareto;
    display_pareto_table(Tp);
    runtimeSummaryTables{k} = display_runtime_phase_summary(T, datasets(k).name, optStartIter);
end

TallCombined = [allT{1}; allT{2}];
runtimeSummaryCombined = display_runtime_phase_summary(TallCombined, "Case 1 + Case 2 (combined)", optStartIter);
write_runtime_and_parameter_summary(allT, allTp, datasets, runtimeSummaryTables, runtimeSummaryCombined, numericalFolder, optStartIter);

create_analysis_plots_side_by_side(allT, allTopt, allPareto, datasets, graphicsFolder, fontSize, plotColors);
%% Combined Pareto
plot_combined_pareto_samples(allTopt{1}, allTp{1}, allTopt{2}, allTp{2}, fullfile(graphicsFolder, "pareto_samples_run1_run2.png"), fontSize, plotColors);
plot_cumulative_runtime_combined(allT, datasets, graphicsFolder, fontSize, plotColors);

%% Refined comparison (former plot_refined_frontier_change.m logic).
run_refined_frontier_change(projectRoot, graphicsFolder, optStartIter);
report_final_frontier_f1_metrics(projectRoot, numericalFolder, optStartIter);
end

%% FUNCTIONS
function [T, Topt, Tp, isPareto] = load_results_table(csvPath, optimizationStartIter)
%LOAD_RESULTS_TABLE Build one run table and derive optimization-only subset.
if nargin < 1 || strlength(string(csvPath)) == 0
    error("You must provide csvPath.");
end
if nargin < 2
    optimizationStartIter = 20;
end
if optimizationStartIter ~= 20
    error("DOE must be exactly 20 iterations. Received optimizationStartIter=%d.", optimizationStartIter);
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

% Exclude DOE points (iterations 1..optimizationStartIter) from all
% Pareto/optimal-controller analyses.
optMask = double(T.iteration) > optimizationStartIter;
Topt = T(optMask, :);

doeCount = nnz(double(T.iteration) <= optimizationStartIter);
if doeCount ~= 20
    error("DOE must be exactly 20 rows in %s. Found %d.", csvPath, doeCount);
end
if height(Topt) ~= 101
    error("Optimization phase must be exactly 101 rows in %s. Found %d.", csvPath, height(Topt));
end

isPareto = compute_pareto_mask(double(Topt.SSE), double(Topt.SSdU));

Tp = Topt(isPareto, :);
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
%ENRICH_WITH_OUT_DATA Add simulation-level runtime diagnostics when available.
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
%COMPUTE_PARETO_MASK Mark non-dominated tradeoff points (lower is better).
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


function create_analysis_plots_side_by_side(allT, allTopt, allPareto, datasets, outDir, fontSize, plotColors)
%CREATE_ANALYSIS_PLOTS_SIDE_BY_SIDE Generate primary manuscript comparison figures.
if numel(allT) < 2 || isempty(allT{1}) || isempty(allT{2})
    return
end
if ~isfolder(outDir)
    mkdir(outDir);
end

% SSE vs SSdU (side-by-side), color mapped by z.
fig1z = figure("Color", "w", "Name", "Pareto SSE vs SSdU by Case");
tiledlayout(fig1z, 1, 2, "Padding", "compact", "TileSpacing", "compact");
seqMap = load_navia_colormap(256);
panelLabels = ["a", "b"];
for k = 1:2
    T = allTopt{k};
    isPareto = allPareto{k};
    ax = nexttile; hold(ax, "on");
    if isempty(T)
        text(ax, 0.5, 0.5, "No optimization points after DOE filtering", ...
            "Units", "normalized", "HorizontalAlignment", "center", "Interpreter", "latex");
    else
        hGuide = plot_pareto_continuum(ax, double(T.SSdU(isPareto)), double(T.SSE(isPareto)), ...
            plotColors(3, :), [1e-2, 1e2], [1e4, 1.3e5]);
        hGuide.Color = plotColors(3, :);
        scatter(ax, double(T.SSdU), double(T.SSE), 80, double(T.f), "filled", "MarkerEdgeColor", "k", "LineWidth", 0.7);
        scatter(ax, double(T.SSdU(isPareto)), double(T.SSE(isPareto)), 170, plotColors(3,:), ...
            "o", "MarkerFaceColor", "none", "MarkerEdgeColor", plotColors(3,:), "LineWidth", 1.2);
    end
    set(ax, "XScale", "log", "YScale", "log", "FontSize", fontSize);
    xlim(ax, [1e-2, 1e2]);
    ylim(ax, [1e4, 1.3e5]);
    colormap(ax, seqMap);
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
fig2 = figure("Color", "w", "Name", "Iteration Runtime and Fidelity by Case");
tiledlayout(fig2, 1, 2, "Padding", "compact", "TileSpacing", "compact");
for k = 1:2
    T = allT{k};
    runtime_h = double(T.runtime_min) / 60;
    ax = nexttile; hold(ax, "on");
    yyaxis(ax, "left");
    plot(ax, T.iteration, runtime_h, "-", "LineWidth", 2.0, "Color", plotColors(k, :));
    ax.YColor = "k";
    ylim(ax, [0, 4]);
    xline(ax, 20, "--", "LineWidth", 2.0, "Color", "k", "Alpha", 1);
    ylabel(ax, "$t_{\mathrm{iter}}$ (h)");
    yyaxis(ax, "right");
    plot(ax, T.iteration, double(T.f), "o", "LineWidth", 2.0, "MarkerSize", 4, "Color", plotColors(3, :));
    ax.YColor = "k";
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
    yticks(ax, 0:1:4);
    format_tick(0, 1);
    yyaxis(ax, "right");
    yticks(ax, 0:0.2:1);
    format_tick(0, 1);
end
save_plot_outputs(fig2, fullfile(outDir, "runtime_vs_iteration_side_by_side.png"), fontSize, 1200, 460);
end


function display_pareto_table(Tp)
%DISPLAY_PARETO_TABLE Print optimization-eligible Pareto controller settings.
tuningCols = ["timestamp","SSE","SSdU","J","p","m","f", ...
    "runtime_min","runtime_per_f", ...
    "Q_x1","Q_x2","Q_x3", ...
    "R_u_x1","R_u_x2","R_u_x3", ...
    "R_du_x1","R_du_x2","R_du_x3"];

tuningCols = tuningCols(ismember(tuningCols, string(Tp.Properties.VariableNames)));
ParetoTuningTable = Tp(:, tuningCols);

disp("Pareto frontier with tuning weights (Q, R, Rdu) [optimization phase only, DOE excluded]:");
disp(ParetoTuningTable);
end


function summaryTbl = display_runtime_phase_summary(T, runName, optimizationStartIter)
%DISPLAY_RUNTIME_PHASE_SUMMARY Explain runtime share spent in DOE vs optimization.
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
%PLOT_COMBINED_PARETO_SAMPLES Show global optimization tradeoff across both runs.
fig = figure("Color", "w", "Name", "Combined Pareto Samples");
ax = axes(fig); hold(ax, "on");
plot_combined_pareto_base(ax, T1, Tp1, T2, Tp2, plotColors);

apply_combined_axes_style(ax, fontSize);

save_plot_outputs(fig, outPath, fontSize, 920, 520);
end


function h = plot_pareto_continuum(ax, x, y, curveColor, xBounds, yBounds)
%PLOT_PARETO_CONTINUUM Draw a smooth monotone visual guide of frontier trend.
[xSort, ord] = sort(x(:), "ascend");
ySort = y(ord);

if numel(xSort) < 3
    h = plot(ax, xSort, ySort, "-o", ...
        "Color", curveColor, "LineWidth", 2.0, ...
        "MarkerSize", 7, "MarkerFaceColor", "w", "MarkerEdgeColor", curveColor);
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
plot(ax, xSort, ySort, "o", ...
    "Color", curveColor, ...
    "MarkerSize", 7, ...
    "MarkerFaceColor", "w", ...
    "MarkerEdgeColor", curveColor, ...
    "LineWidth", 1.4);
end


function h = plot_pareto_continuum_line_only(ax, x, y, curveColor, xBounds, yBounds)
%PLOT_PARETO_CONTINUUM_LINE_ONLY Smooth pchip Pareto guide without point markers.
[xSort, ord] = sort(x(:), "ascend");
ySort = y(ord);

if numel(xSort) < 3
    h = plot(ax, xSort, ySort, "-", "Color", curveColor, "LineWidth", 2.0);
    return;
end

lx = log10(xSort);
ly = log10(ySort);
lxqIn = linspace(min(lx), max(lx), 220);
lyqIn = pchip(lx, ly, lxqIn);
xq = (10.^lxqIn).';
yq = (10.^lyqIn).';

if nargin >= 6 && ~isempty(yBounds)
    yTop = max(yBounds);
    if yTop > yq(1)
        xq = [xq(1); xq];
        yq = [yTop; yq];
    end
end

if nargin >= 5 && ~isempty(xBounds)
    xRight = max(xBounds);
    if xRight > xq(end)
        xq = [xq; xRight];
        yq = [yq; yq(end)];
    end
end

h = plot(ax, xq, yq, "-", "Color", curveColor, "LineWidth", 2.0);
end


function [finalMask, Tf] = plot_combined_pareto_base(ax, T1, Tp1, T2, Tp2, plotColors)
%PLOT_COMBINED_PARETO_BASE Shared base drawing used by Figure 3 and Figure 5a/5b.
Tall = [T1; T2];
if ismember("iteration", string(Tall.Properties.VariableNames)) && any(double(Tall.iteration) <= 20)
    error("Combined Pareto plot received DOE rows (iteration <= 20), which is not allowed.");
end
scatter(ax, double(Tall.SSdU), double(Tall.SSE), 18, ...
    "filled", "MarkerFaceColor", "k", "MarkerEdgeColor", "none", "DisplayName", "Optimization samples");

finalMask = compute_pareto_mask(double(Tall.SSE), double(Tall.SSdU));
Tf = Tall(finalMask, :);
plot_pareto_continuum(ax, double(Tf.SSdU), double(Tf.SSE), plotColors(3, :), [1e-2, 1e2], [1e4, 1.3e5]);

scatter(ax, double(Tp1.SSdU), double(Tp1.SSE), 80, plotColors(1,:), ...
    "o", "MarkerFaceColor", plotColors(1,:), "MarkerEdgeColor", plotColors(1,:), ...
    "LineWidth", 1.4);
scatter(ax, double(Tp2.SSdU), double(Tp2.SSE), 90, plotColors(2,:), ...
    "^", "MarkerFaceColor", plotColors(2,:), "MarkerEdgeColor", plotColors(2,:), ...
    "LineWidth", 1.4);

scatter(ax, double(Tf.SSdU), double(Tf.SSE), 300, plotColors(3,:), ...
    "o", "MarkerFaceColor", "none", "MarkerEdgeColor", plotColors(3,:), ...
    "LineWidth", 2);
end


function plot_combined_samples_no_guide(ax, T1, Tp1, T2, Tp2, plotColors)
%PLOT_COMBINED_SAMPLES_NO_GUIDE Combined samples + case markers, without guideline/pareto ring.
Tall = [T1; T2];
if ismember("iteration", string(Tall.Properties.VariableNames)) && any(double(Tall.iteration) <= 20)
    error("Combined Pareto plot received DOE rows (iteration <= 20), which is not allowed.");
end
scatter(ax, double(Tall.SSdU), double(Tall.SSE), 18, ...
    "filled", "MarkerFaceColor", "k", "MarkerEdgeColor", "none");

scatter(ax, double(Tp1.SSdU), double(Tp1.SSE), 80, plotColors(1,:), ...
    "o", "MarkerFaceColor", plotColors(1,:), "MarkerEdgeColor", plotColors(1,:), ...
    "LineWidth", 1.4);
scatter(ax, double(Tp2.SSdU), double(Tp2.SSE), 90, plotColors(2,:), ...
    "^", "MarkerFaceColor", plotColors(2,:), "MarkerEdgeColor", plotColors(2,:), ...
    "LineWidth", 1.4);
end


function apply_combined_axes_style(ax, fontSize)
%APPLY_COMBINED_AXES_STYLE Match Figure 3 axes/tick styling.
set(ax, "XScale", "log", "YScale", "log", "FontSize", fontSize);
xlim(ax, [1e-2, 2e0]);
ylim(ax, [1e4, 3.5e4]);
xlabel(ax, "$J_{\mathrm{TV}}$");
ylabel(ax, "$J_{\mathrm{track}}$");
grid(ax, "off");
box(ax, "off");
axes(ax);
format_tick(1, 1);
end


function [isParetoCombined, Tf] = plot_original_frontier_from_non_doe(ax, T, plotColors, doeIterationsPerRun)
%PLOT_ORIGINAL_FRONTIER_FROM_NON_DOE Plot combined original frontier with strict DOE exclusion.
if nargin < 4
    doeIterationsPerRun = 20;
end
if ismember("iteration", string(T.Properties.VariableNames)) && any(double(T.iteration) <= doeIterationsPerRun)
    error("Pareto pool contains DOE rows (iteration <= %d), which is not allowed.", doeIterationsPerRun);
end

scatter(ax, double(T.SSdU), double(T.SSE), 18, ...
    "filled", "MarkerFaceColor", "k", "MarkerEdgeColor", "none");

isParetoCombined = compute_pareto_mask(double(T.SSE), double(T.SSdU));
Tf = T(isParetoCombined, :);
plot_pareto_polyline_with_markers_refined(ax, double(Tf.SSdU), double(Tf.SSE), plotColors(3, :), 2.0, 8);

if ismember("run_key", string(T.Properties.VariableNames))
    isRun1Front = isParetoCombined & (string(T.run_key) == "run1");
    isRun2Front = isParetoCombined & (string(T.run_key) == "run2");
    scatter(ax, double(T.SSdU(isRun1Front)), double(T.SSE(isRun1Front)), 80, plotColors(1,:), ...
        "o", "MarkerFaceColor", plotColors(1,:), "MarkerEdgeColor", plotColors(1,:), "LineWidth", 1.4);
    scatter(ax, double(T.SSdU(isRun2Front)), double(T.SSE(isRun2Front)), 90, plotColors(2,:), ...
        "^", "MarkerFaceColor", plotColors(2,:), "MarkerEdgeColor", plotColors(2,:), "LineWidth", 1.4);
end
end


function cmap = load_navia_colormap(n)
%LOAD_NAVIA_COLORMAP Load the sequential map used for fidelity-color plots.
matPath = which("navia.mat");
if strlength(string(matPath)) == 0
    error("Unable to locate navia.mat on MATLAB path.");
end
S = load(matPath, "navia");
if ~isfield(S, "navia")
    error("File %s does not contain variable 'navia'.", matPath);
end
cmap = double(S.navia);
if size(cmap, 2) ~= 3
    error("Variable 'navia' in %s must have 3 columns (RGB).", matPath);
end
if size(cmap, 1) ~= n
    x = linspace(0, 1, size(cmap, 1));
    xq = linspace(0, 1, n);
    cmap = interp1(x, cmap, xq, "linear");
end
cmap = min(max(cmap, 0), 1);
end


function write_runtime_and_parameter_summary(allT, allTp, datasets, runSummaryTables, combinedSummaryTable, outDir, optimizationStartIter)
%WRITE_RUNTIME_AND_PARAMETER_SUMMARY Persist manuscript-facing runtime/tuning summaries.
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

% Share of optimization-phase points where N_c (=m) is 1.
nc1PctCase = NaN(numel(allT), 1);
nc1CountCase = zeros(numel(allT), 1);
optRowsCase = zeros(numel(allT), 1);
for k = 1:numel(allT)
    [nc1PctCase(k), nc1CountCase(k), optRowsCase(k)] = ...
        compute_m1_share_optimization(allT{k}, optimizationStartIter);
end
[nc1PctCombined, nc1CountCombined, optRowsCombined] = ...
    compute_m1_share_optimization(Tall, optimizationStartIter);

fprintf("Percentage of optimization points with N_c = 1:\n");
for k = 1:numel(datasets)
    fprintf("  %s: %d/%d (%.6g%%)\n", string(datasets(k).name), ...
        nc1CountCase(k), optRowsCase(k), nc1PctCase(k));
end
fprintf("  Case 1 + Case 2 (combined): %d/%d (%.6g%%)\n", ...
    nc1CountCombined, optRowsCombined, nc1PctCombined);

fprintf(fid, "Percentage of optimization points with N_c = 1:\n");
for k = 1:numel(datasets)
    fprintf(fid, "  %s: %d/%d (%.6g%%)\n", char(string(datasets(k).name)), ...
        nc1CountCase(k), optRowsCase(k), nc1PctCase(k));
end
fprintf(fid, "  Case 1 + Case 2 (combined): %d/%d (%.6g%%)\n\n", ...
    nc1CountCombined, optRowsCombined, nc1PctCombined);

% Export a compact standalone txt for manuscript bookkeeping.
nc1Path = fullfile(outDir, "optimization_nc1_share.txt");
fidNc1 = fopen(nc1Path, "w");
if fidNc1 ~= -1
    cleanupNc1 = onCleanup(@() fclose(fidNc1)); %#ok<NASGU>
    fprintf(fidNc1, "Percentage of optimization points with N_c = 1:\n");
    for k = 1:numel(datasets)
        fprintf(fidNc1, "  %s: %d/%d (%.6g%%)\n", char(string(datasets(k).name)), ...
            nc1CountCase(k), optRowsCase(k), nc1PctCase(k));
    end
    fprintf(fidNc1, "  Case 1 + Case 2 (combined): %d/%d (%.6g%%)\n", ...
        nc1CountCombined, optRowsCombined, nc1PctCombined);
else
    warning("Unable to write N_c=1 optimization summary: %s", nc1Path);
end

for k = 1:numel(allTp)
    T = allTp{k};
    validM = isfinite(T.m);
    validP = isfinite(T.p);
    nRows = height(T);
    countM1 = nnz(T.m(validM) == 1);
    countP1 = nnz(T.p(validP) == 1);

    fprintf(fid, "%s Pareto counts (optimization phase only, DOE excluded):\n", char(string(datasets(k).name)));
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


function [pctM1, countM1, nOpt] = compute_m1_share_optimization(T, optimizationStartIter)
%COMPUTE_M1_SHARE_OPTIMIZATION Quantify how often N_c=1 appears in optimization.
if isempty(T) || ...
        ~ismember("m", string(T.Properties.VariableNames)) || ...
        ~ismember("iteration", string(T.Properties.VariableNames))
    pctM1 = NaN;
    countM1 = 0;
    nOpt = 0;
    return
end

iter = double(T.iteration);
mVals = double(T.m);
optMask = iter > optimizationStartIter;
nOpt = nnz(optMask);
if nOpt == 0
    pctM1 = NaN;
    countM1 = 0;
    return
end

countM1 = nnz(optMask & (mVals == 1));
pctM1 = 100 * countM1 / nOpt;
end


function meanZ = compute_mean_z_optimization(T, optimizationStartIter)
%COMPUTE_MEAN_Z_OPTIMIZATION Track average fidelity level during optimization.
if isempty(T) || ~ismember("f", string(T.Properties.VariableNames))
    meanZ = NaN;
    return
end
mask = double(T.iteration) > optimizationStartIter;
meanZ = mean(double(T.f(mask)), "omitnan");
end


function write_runtime_table(fid, summaryTbl)
%WRITE_RUNTIME_TABLE Write one runtime summary table to disk.
for i = 1:height(summaryTbl)
    fprintf(fid, "  %s | iterations=%d | runtime_min=%.6g | runtime_pct_total=%.6g\n", ...
        char(string(summaryTbl.phase{i})), summaryTbl.iterations(i), ...
        summaryTbl.runtime_min(i), summaryTbl.runtime_pct_total(i));
end
end


function validate_no_doe_in_pareto(Tp, optimizationStartIter, runName)
%VALIDATE_NO_DOE_IN_PARETO Safety guard: prevent DOE leakage into Pareto outputs.
if isempty(Tp) || ~ismember("iteration", string(Tp.Properties.VariableNames))
    return
end
if any(double(Tp.iteration) <= optimizationStartIter)
    error("DOE point detected in Pareto table for %s.", string(runName));
end
fprintf("Verified Pareto DOE exclusion for %s: %d rows, all with iteration > %d.\n", ...
    string(runName), height(Tp), optimizationStartIter);
end


function plot_cumulative_runtime_combined(allT, datasets, outDir, fontSize, plotColors)
%PLOT_CUMULATIVE_RUNTIME_COMBINED Compare wall-clock accumulation across runs.
if numel(allT) < 2 || isempty(allT{1}) || isempty(allT{2})
    return
end

fig = figure("Color","w", "Name", "Cumulative Runtime by Case");
ax = axes(fig); hold(ax, "on");

runtime_h_1 = double(allT{1}.runtime_min) / 60;
runtime_h_2 = double(allT{2}.runtime_min) / 60;

plot(ax, double(allT{1}.iteration), cumsum(runtime_h_1, "omitnan"), "-", "LineWidth", 2.0, ...
    "Color", plotColors(1,:), ...
    "DisplayName", string(datasets(1).name));
plot(ax, double(allT{2}.iteration), cumsum(runtime_h_2, "omitnan"), "-.", "LineWidth", 2.0, ...
    "Color", plotColors(2,:), ...
    "DisplayName", string(datasets(2).name));
xline(ax, 20, "--", "LineWidth", 2.0, "Color", "k", "Alpha", 1);

xlabel(ax, "$k$ (iteration)");
ylabel(ax, "$t_{\mathrm{run}}$ (h)");
xlim(ax, [1, max(1, max(height(allT{1}), height(allT{2})))]);
set(ax, "FontSize", fontSize);
grid(ax, "off");
box(ax, "off");

save_plot_outputs(fig, fullfile(outDir, "runtime_cumulative_run1_run2.png"), fontSize, 920, 520);
end


function save_plot_outputs(figHandle, pngPath, fontSize, figWidthPx, figHeightPx)
%SAVE_PLOT_OUTPUTS Export standardized figure files for reporting pipelines.
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


function run_refined_frontier_change(projectRoot, graphicsFolder, doeIterationsPerRun)
%RUN_REFINED_FRONTIER_CHANGE Compare original (non-DOE) and refined Pareto points.
    originalFiles = [
        struct("runKey","run1", "runLabel","Case 1", "path", fullfile(projectRoot, "results", "run1", "results.csv"));
        struct("runKey","run2", "runLabel","Case 2", "path", fullfile(projectRoot, "results", "run2", "results.csv"))
    ];
    refinedFiles = [
        struct("runKey","run1", "runLabel","Case 1 refined", "path", fullfile(projectRoot, "results", "final_fidelity_same_noise", "run1_full_f1_same_noise", "results_full.csv"));
        struct("runKey","run2", "runLabel","Case 2 refined", "path", fullfile(projectRoot, "results", "final_fidelity_same_noise", "run2_full_f1_same_noise", "results_full.csv"))
    ];

    print_runtime_cfg_refined(originalFiles, refinedFiles, doeIterationsPerRun);

    T_orig_full = load_results_for_refined(originalFiles, 0);
    T_orig = T_orig_full(T_orig_full.iteration > doeIterationsPerRun, :);
    T_ref = load_results_for_refined(refinedFiles, 0);
    if isempty(T_orig)
        error("No original rows remain after DOE filtering.");
    end
    if isempty(T_ref)
        error("No rows loaded from refined result files.");
    end

    isParetoOrig = compute_pareto_mask(double(T_orig.SSE), double(T_orig.SSdU));
    [T_ref, refFilterInfo] = keep_refined_from_original_pareto_refined(T_ref, T_orig, isParetoOrig, T_orig_full);
    print_refined_filter_summary_refined(refFilterInfo);
    if isempty(T_ref)
        error("No refined rows remain after filtering to original non-DOE Pareto points.");
    end
    origParetoKeys = unique(string(T_orig.match_key(isParetoOrig)), "stable");
    refKeys = unique(string(T_ref.match_key), "stable");
    missingRefined = setdiff(origParetoKeys, refKeys, "stable");
    if ~isempty(missingRefined)
        error("Missing z=1 refined results for %d Pareto controller(s): %s", ...
            numel(missingRefined), strjoin(missingRefined, ", "));
    end

    T_jtv = compute_jtv_change_table_refined(T_orig, T_ref);
    fprintf("\n=== J_TV change (original -> refined), sorted by |delta %%| descending ===\n");
    disp(T_jtv);
    T_jtrack = compute_jtrack_change_table_refined(T_orig, T_ref);
    fprintf("\n=== J_track change (original -> refined), sorted by |delta %%| descending ===\n");
    disp(T_jtrack);
    print_top1_cfg_refined(projectRoot, T_jtv);

    isParetoRef = compute_pareto_mask(double(T_ref.SSE), double(T_ref.SSdU));
    [commonKeys, idxOrig, idxRef] = intersect(string(T_orig.match_key), string(T_ref.match_key), "stable");
    [matchedRunKey, matchedTs] = split_match_key_refined(commonKeys);

    NATURE_COLOR = nature_methods_colors();
    plotColors = nature_methods_colors(3); % Blue, BluishGreen, ReddishPurple
    colRefAll = [0.60 0.82 0.98];
    colPromoted = plotColors(2, :);
    % Match dependencies/plot_utils/good_colors.m: C.orange = [242, 133, 34].
    benchmarkColor = [242, 133, 34] / 255;
    [benchJtrack, benchJTV, benchFound] = load_noisy_benchmark_point(projectRoot);

    % Promotion bookkeeping.
    origParetoMatched = isParetoOrig(idxOrig);
    refParetoMatched = isParetoRef(idxRef);
    zOrigMatched = double(T_orig.z_eval(idxOrig));
    promotedMask = ~origParetoMatched & refParetoMatched;
    promotedZlt1Mask = promotedMask & isfinite(zOrigMatched) & (zOrigMatched < 1 - 1e-12);
    promotedIdxOrig = idxOrig(promotedZlt1Mask);
    promotedIdxRef = idxRef(promotedZlt1Mask);
    promotedRunKey = matchedRunKey(promotedZlt1Mask);
    promotedTs = matchedTs(promotedZlt1Mask);
    matchedOrigParetoCount = nnz(origParetoMatched);
    matchedRefParetoCount = nnz(refParetoMatched);

    T_orig_run1 = T_orig(string(T_orig.run_key) == "run1", :);
    T_orig_run2 = T_orig(string(T_orig.run_key) == "run2", :);
    Tp_orig_run1 = T_orig_run1(compute_pareto_mask(double(T_orig_run1.SSE), double(T_orig_run1.SSdU)), :);
    Tp_orig_run2 = T_orig_run2(compute_pareto_mask(double(T_orig_run2.SSE), double(T_orig_run2.SSdU)), :);

    fig = figure("Color", "w", "Toolbar", "none", "Name", "Pareto Frontier Change");
    tiledlayout(fig, 1, 2, "Padding", "compact", "TileSpacing", "compact");
    set(fig, "Position", [80 80 1400 560]);

    axL = nexttile; hold(axL, "on");
    [isParetoLeft, ~] = plot_combined_pareto_base(axL, T_orig_run1, Tp_orig_run1, T_orig_run2, Tp_orig_run2, plotColors);
    if benchFound
        scatter(axL, benchJTV, benchJtrack, 130, "s", ...
            "MarkerFaceColor", benchmarkColor, "MarkerEdgeColor", "none", "LineWidth", 1.0);
    end

    axR = nexttile; hold(axR, "on");
    plot_combined_samples_no_guide(axR, T_orig_run1, Tp_orig_run1, T_orig_run2, Tp_orig_run2, plotColors);
    scatter(axR, double(T_ref.SSdU), double(T_ref.SSE), 42, ...
        "filled", "MarkerFaceColor", colRefAll, "MarkerEdgeColor", "none");
    TrefPareto = T_ref(isParetoRef, :);
    scatter(axR, double(TrefPareto.SSdU), double(TrefPareto.SSE), 80, ...
        "d", "MarkerFaceColor", NATURE_COLOR.ReddishPurple, ...
        "MarkerEdgeColor", NATURE_COLOR.ReddishPurple, "LineWidth", 1.0);
    if ~isempty(promotedIdxRef)
        scatter(axR, double(T_ref.SSdU(promotedIdxRef)), double(T_ref.SSE(promotedIdxRef)), 140, ...
            "p", "MarkerFaceColor", colPromoted, "MarkerEdgeColor", "k", "LineWidth", 0.7);
    end
    if benchFound
        scatter(axR, benchJTV, benchJtrack, 130, "s", ...
            "MarkerFaceColor", benchmarkColor, "MarkerEdgeColor", "none", "LineWidth", 1.0);
    end

    if benchFound
        [nOrigSuperior, superiorMaskOrig] = count_benchmark_strict_superior_refined( ...
            double(T_orig.SSE(isParetoOrig)), double(T_orig.SSdU(isParetoOrig)), benchJtrack, benchJTV);
        [nRefSuperior, superiorMaskRef] = count_benchmark_strict_superior_refined( ...
            double(TrefPareto.SSE), double(TrefPareto.SSdU), benchJtrack, benchJTV);
        fprintf("Pareto controllers strictly better than noisy benchmark (J_track and J_TV both lower):\n");
        fprintf("  Original Pareto (non-DOE): %d/%d\n", nOrigSuperior, nnz(isParetoOrig));
        print_benchmark_ratio_stats_refined("Original Pareto (non-DOE)", ...
            double(T_orig.SSE(isParetoOrig)), double(T_orig.SSdU(isParetoOrig)), superiorMaskOrig, benchJtrack, benchJTV);
        fprintf("  Refined Pareto: %d/%d\n", nRefSuperior, height(TrefPareto));
        print_benchmark_ratio_stats_refined("Refined Pareto", ...
            double(TrefPareto.SSE), double(TrefPareto.SSdU), superiorMaskRef, benchJtrack, benchJTV);
        TorigBench = build_benchmark_comparison_table_refined( ...
            string(T_orig.run_key(isParetoOrig)), string(T_orig.timestamp(isParetoOrig)), ...
            double(T_orig.SSE(isParetoOrig)), double(T_orig.SSdU(isParetoOrig)), benchJtrack, benchJTV);
        fprintf("  Original Pareto benchmark comparison table:\n");
        disp(TorigBench);
        TrefBench = build_benchmark_comparison_table_refined( ...
            string(TrefPareto.run_key), string(TrefPareto.timestamp), ...
            double(TrefPareto.SSE), double(TrefPareto.SSdU), benchJtrack, benchJTV);
        fprintf("  Refined Pareto benchmark comparison table:\n");
        disp(TrefBench);
    end

    apply_combined_axes_style(axL, 20);
    % Keep previous panel-b custom limits for quick swap:
    % [xLimR, yLimR] = refined_pareto_axis_limits_10pct(TrefPareto);
    % xLimR(1) = 2e-2;
    xLimR = [1e-2, 2e0];
    yLimR = [1e4, 3.5e4];
    if benchFound
        xLimR = expand_log_limits_to_include_point_refined(xLimR, benchJTV, 1.10);
        yLimR = expand_log_limits_to_include_point_refined(yLimR, benchJtrack, 1.10);
    end
    plot_pareto_continuum_line_only(axR, double(TrefPareto.SSdU), double(TrefPareto.SSE), ...
        NATURE_COLOR.ReddishPurple, xLimR, yLimR);
    apply_combined_axes_style(axR, 20);
    xlim(axL, xLimR);
    ylim(axL, yLimR);
    xlim(axR, xLimR);
    ylim(axR, yLimR);
    axes(axR);
    format_tick(1, 1);

    title(axL, "$\mathbf{a}$", "Interpreter", "latex");
    title(axR, "$\mathbf{b}$", "Interpreter", "latex");
    axL.TitleHorizontalAlignment = "left";
    axR.TitleHorizontalAlignment = "left";

    outStem = fullfile(graphicsFolder, "refined_frontier_change");
    exportgraphics(fig, outStem + ".png", "Resolution", 300);
    exportgraphics(fig, outStem + ".pdf", "ContentType", "vector");

    outTxtDir = fullfile(projectRoot, "results", "txt results");
    if ~isfolder(outTxtDir)
        mkdir(outTxtDir);
    end
    reportPath = fullfile(outTxtDir, "refined_promoted_frontier_z_lt_1.txt");
    write_promoted_report_refined(reportPath, promotedRunKey, promotedTs, T_orig, T_ref, promotedIdxOrig, promotedIdxRef);

    fprintf("Saved: %s\n", outStem + ".png");
    fprintf("Saved: %s\n", outStem + ".pdf");
    fprintf("Saved: %s\n", reportPath);
    fprintf("J_TV rows compared: %d\n", height(T_jtv));
    fprintf("Matched points (run_key + timestamp): %d\n", numel(idxOrig));
    fprintf("Left-panel frontier pool size (all original non-DOE points): %d\n", height(T_orig));
    fprintf("Left-panel frontier points (combined Pareto): %d\n", nnz(isParetoLeft));
    fprintf("Original Pareto points (all original rows after DOE filter): %d\n", nnz(isParetoOrig));
    fprintf("Refined Pareto points (all refined rows): %d\n", nnz(isParetoRef));
    fprintf("Original Pareto points in matched set: %d\n", matchedOrigParetoCount);
    fprintf("Refined Pareto points in matched set: %d\n", matchedRefParetoCount);
    fprintf("Promoted to Pareto with z < 1: %d\n", numel(promotedIdxRef));
end


function limOut = expand_log_limits_to_include_point_refined(limIn, v, padFactor)
%EXPAND_LOG_LIMITS_TO_INCLUDE_POINT_REFINED Expand positive log-axis limits to include point v.
limOut = limIn;
if ~(isfinite(v) && v > 0)
    return
end
if nargin < 3 || ~isfinite(padFactor) || padFactor <= 1
    padFactor = 1.10;
end
if v < limOut(1)
    limOut(1) = v / padFactor;
end
if v > limOut(2)
    limOut(2) = v * padFactor;
end
end


function [nSuperior, superiorMask] = count_benchmark_strict_superior_refined(Jtrack, JTV, benchJtrack, benchJTV)
%COUNT_BENCHMARK_STRICT_SUPERIOR_REFINED Strictly better than benchmark in both objectives.
superiorMask = isfinite(Jtrack) & isfinite(JTV) & ...
    (Jtrack < benchJtrack) & (JTV < benchJTV);
nSuperior = nnz(superiorMask);
end


function print_benchmark_ratio_stats_refined(label, Jtrack, JTV, superiorMask, benchJtrack, benchJTV)
%PRINT_BENCHMARK_RATIO_STATS_REFINED Print benchmark-relative ratios for superior controllers.
if nargin < 4 || isempty(superiorMask) || ~any(superiorMask)
    fprintf("    %s ratios: n/a (no strictly superior controllers)\n", label);
    return
end
JtrackSup = double(Jtrack(superiorMask));
JTVSup = double(JTV(superiorMask));
trackPct = 100 * mean(JtrackSup / benchJtrack, "omitnan");
tvPct = 100 * mean(JTVSup / benchJTV, "omitnan");
trackTimes = mean(benchJtrack ./ JtrackSup, "omitnan");
tvTimes = mean(benchJTV ./ JTVSup, "omitnan");
fprintf("    %s mean (J_better/J_bench): J_track=%.1f%%, J_TV=%.1f%%\n", label, trackPct, tvPct);
fprintf("    %s mean benchmark higher factor: J_track=%.1fx, J_TV=%.1fx\n", label, trackTimes, tvTimes);
end


function Tcmp = build_benchmark_comparison_table_refined(runKey, timestamp, Jtrack, JTV, benchJtrack, benchJTV)
%BUILD_BENCHMARK_COMPARISON_TABLE_REFINED Per-controller benchmark-relative metrics.
runKey = string(runKey(:));
timestamp = string(timestamp(:));
Jtrack = double(Jtrack(:));
JTV = double(JTV(:));
n = numel(Jtrack);

ratioTrackPct = nan(n,1);
ratioTVPct = nan(n,1);
timesTrack = nan(n,1);
timesTV = nan(n,1);
isSuperior = false(n,1);

valid = isfinite(Jtrack) & isfinite(JTV) & (Jtrack > 0) & (JTV > 0);
ratioTrackPct(valid) = 100 * (Jtrack(valid) / benchJtrack);
ratioTVPct(valid) = 100 * (JTV(valid) / benchJTV);
timesTrack(valid) = benchJtrack ./ Jtrack(valid);
timesTV(valid) = benchJTV ./ JTV(valid);
isSuperior(valid) = (Jtrack(valid) < benchJtrack) & (JTV(valid) < benchJTV);

Tcmp = table(ratioTrackPct, ratioTVPct, timesTrack, timesTV, isSuperior, ...
    'VariableNames', {'J_track_ratio_pct','J_TV_ratio_pct', ...
    'J_track_bench_higher_x','J_TV_bench_higher_x','strictly_better'});
Tcmp = Tcmp(Tcmp.strictly_better, :);
Tcmp = Tcmp(:, {'J_track_ratio_pct','J_TV_ratio_pct','J_track_bench_higher_x','J_TV_bench_higher_x'});
if ~isempty(Tcmp)
    Tcmp = sortrows(Tcmp, 'J_TV_ratio_pct', 'ascend');
end
end


function [xLim, yLim] = refined_pareto_axis_limits_10pct(TrefPareto)
%REFINED_PARETO_AXIS_LIMITS_10PCT Axis limits set to +/-10% of refined Pareto cost range.
xVals = double(TrefPareto.SSdU);
yVals = double(TrefPareto.SSE);
xVals = xVals(isfinite(xVals) & (xVals > 0));
yVals = yVals(isfinite(yVals) & (yVals > 0));
if isempty(xVals) || isempty(yVals)
    xLim = [1e-2, 2e0];
    yLim = [1e4, 3.5e4];
    return
end
xLim = [0.9 * min(xVals), 1.1 * max(xVals)];
yLim = [0.9 * min(yVals), 1.1 * max(yVals)];
end


function [Jtrack, JTV, found] = load_noisy_benchmark_point(projectRoot)
%LOAD_NOISY_BENCHMARK_POINT Load noisy benchmark aggregate objectives.
benchPath = fullfile(projectRoot, "results", "benchmark_reference_controller", ...
    "benchmark_full_f1_same_noise_fix", "out_benchmark.mat");
Jtrack = nan;
JTV = nan;
found = false;
if ~isfile(benchPath)
    warning("Noisy benchmark file not found for refined frontier overlay: %s", benchPath);
    return
end
S = load(benchPath, "out");
if ~isfield(S, "out") || ~isfield(S.out, "SSE") || ~isfield(S.out, "SSdU")
    warning("Noisy benchmark file missing required out.SSE/out.SSdU fields: %s", benchPath);
    return
end
Jtrack = double(S.out.SSE);
JTV = double(S.out.SSdU);
found = isfinite(Jtrack) && isfinite(JTV) && Jtrack > 0 && JTV > 0;
if ~found
    warning("Noisy benchmark values are invalid for refined frontier overlay: %s", benchPath);
end
end


function print_runtime_cfg_refined(originalFiles, refinedFiles, doeIterationsPerRun)
%PRINT_RUNTIME_CFG_REFINED Print reproducibility settings for refined comparison.
    fprintf("=== refined frontier settings ===\n");
    fprintf("MATLAB version: %s\n", version);
    fprintf("Current folder: %s\n", pwd);
    fprintf("DOE exclusion in original runs: drop first %d iterations per run\n", doeIterationsPerRun);
    p = gcp("nocreate");
    if isempty(p)
        fprintf("Parallel pool: not running\n");
        workers = 0;
    else
        workers = p.NumWorkers;
        fprintf("Parallel pool: running (%d workers)\n", workers);
    end
    fprintf("Configured worker count (observed): %d\n", workers);
    fprintf("Original CSV sources:\n");
    for k = 1:numel(originalFiles)
        fprintf("  - %s\n", originalFiles(k).path);
    end
    fprintf("Refined CSV sources:\n");
    for k = 1:numel(refinedFiles)
        fprintf("  - %s\n", refinedFiles(k).path);
    end
end


function T = load_results_for_refined(fileDefs, dropFirstNRows)
%LOAD_RESULTS_FOR_REFINED Load and normalize CSV rows for refined comparison.
    T = table();
    for k = 1:numel(fileDefs)
        p = fileDefs(k).path;
        if ~isfile(p)
            warning("Missing file: %s", p);
            continue
        end
        Tk = readtable(p, "TextType", "string");
        required = ["timestamp", "SSE", "SSdU"];
        for c = required
            if ~ismember(c, string(Tk.Properties.VariableNames))
                error("Missing column '%s' in %s", c, p);
            end
        end
        Tk.run_label = repmat(string(fileDefs(k).runLabel), height(Tk), 1);
        Tk.run_key = repmat(string(fileDefs(k).runKey), height(Tk), 1);
        Tk.iteration = (1:height(Tk)).';
        if dropFirstNRows > 0
            Tk = Tk(Tk.iteration > dropFirstNRows, :);
        end
        if ismember("theta_1", string(Tk.Properties.VariableNames))
            Tk.z_eval = double(Tk.theta_1);
        else
            Tk.z_eval = nan(height(Tk), 1);
        end
        Tk.match_key = Tk.run_key + "|" + string(Tk.timestamp);
        Tk = Tk(:, ["run_label", "run_key", "match_key", "iteration", "timestamp", "SSE", "SSdU", "z_eval"]);
        T = [T; Tk]; %#ok<AGROW>
    end
end


function [TrefKeep, info] = keep_refined_from_original_pareto_refined(T_ref, T_orig, isParetoOrig, T_orig_full)
%KEEP_REFINED_FROM_ORIGINAL_PARETO_REFINED Keep only refined points linked to original Pareto.
    refKeys = string(T_ref.match_key);
    origKeys = string(T_orig.match_key);
    origFullKeys = string(T_orig_full.match_key);
    origParetoKeys = string(T_orig.match_key(isParetoOrig));

    inOrigFull = ismember(refKeys, origFullKeys);
    inOrig = ismember(refKeys, origKeys);
    inOrigPareto = ismember(refKeys, origParetoKeys);

    info = struct();
    info.n_ref_input = height(T_ref);
    info.n_drop_doe = nnz(inOrigFull & ~inOrig);
    info.n_drop_not_in_orig = nnz(~inOrigFull);
    info.n_drop_in_orig_not_pareto = nnz(inOrig & ~inOrigPareto);
    info.dropped_doe_keys = refKeys(inOrigFull & ~inOrig);
    TrefKeep = T_ref(inOrigPareto, :);
    if ~isempty(TrefKeep)
        [~, uniqIdx] = unique(string(TrefKeep.match_key), "stable");
        info.n_drop_duplicate_refined = height(TrefKeep) - numel(uniqIdx);
        TrefKeep = TrefKeep(uniqIdx, :);
    else
        info.n_drop_duplicate_refined = 0;
    end
    info.n_keep_final = height(TrefKeep);
end


function print_refined_filter_summary_refined(info)
%PRINT_REFINED_FILTER_SUMMARY_REFINED Print keep/drop diagnostics.
    fprintf("\n=== Refined sample eligibility check ===\n");
    fprintf("Input refined rows: %d\n", info.n_ref_input);
    fprintf("Dropped (mapped to original DOE rows): %d\n", info.n_drop_doe);
    fprintf("Dropped (not found in original run datasets): %d\n", info.n_drop_not_in_orig);
    fprintf("Dropped (original but not Pareto): %d\n", info.n_drop_in_orig_not_pareto);
    fprintf("Dropped duplicate refined rows: %d\n", info.n_drop_duplicate_refined);
    fprintf("Kept refined rows (original non-DOE Pareto-linked): %d\n", info.n_keep_final);
    if info.n_drop_doe > 0
        fprintf("Dropped DOE-linked keys:\n");
        for i = 1:numel(info.dropped_doe_keys)
            fprintf("  - %s\n", string(info.dropped_doe_keys(i)));
        end
    end
end


function T = compute_jtv_change_table_refined(T_orig, T_ref)
%COMPUTE_JTV_CHANGE_TABLE_REFINED Build delta-JTV table for matched points.
    [commonKeys, idxOrig, idxRef] = intersect(string(T_orig.match_key), string(T_ref.match_key), "stable");
    if isempty(commonKeys)
        error("No matching (run_key, timestamp) points between original and refined tables.");
    end
    [runKey, ts] = split_match_key_refined(commonKeys);
    T = table(runKey, ts, double(T_orig.SSdU(idxOrig)), double(T_ref.SSdU(idxRef)), ...
        'VariableNames', {'run_key', 'timestamp', 'JTV_original', 'JTV_refined'});
    T.delta_JTV = T.JTV_refined - T.JTV_original;
    T.delta_pct = 100 * T.delta_JTV ./ max(abs(T.JTV_original), eps);
    T.abs_delta_JTV = abs(T.delta_JTV);
    T.abs_delta_pct = abs(T.delta_pct);
    T = sortrows(T, "abs_delta_pct", "descend");
end


function T = compute_jtrack_change_table_refined(T_orig, T_ref)
%COMPUTE_JTRACK_CHANGE_TABLE_REFINED Build delta-J_track table for matched points.
    [commonKeys, idxOrig, idxRef] = intersect(string(T_orig.match_key), string(T_ref.match_key), "stable");
    if isempty(commonKeys)
        error("No matching (run_key, timestamp) points between original and refined tables.");
    end
    [runKey, ts] = split_match_key_refined(commonKeys);
    T = table(runKey, ts, double(T_orig.SSE(idxOrig)), double(T_ref.SSE(idxRef)), ...
        'VariableNames', {'run_key', 'timestamp', 'Jtrack_original', 'Jtrack_refined'});
    T.delta_Jtrack = T.Jtrack_refined - T.Jtrack_original;
    T.delta_pct = 100 * T.delta_Jtrack ./ max(abs(T.Jtrack_original), eps);
    T.abs_delta_Jtrack = abs(T.delta_Jtrack);
    T.abs_delta_pct = abs(T.delta_pct);
    T = sortrows(T, "abs_delta_pct", "descend");
end


function print_top1_cfg_refined(projectRoot, T_jtv)
%PRINT_TOP1_CFG_REFINED Print original/refined configs for the largest JTV change.
    if isempty(T_jtv)
        return
    end
    runKey = string(T_jtv.run_key(1));
    ts = string(T_jtv.timestamp(1));
    origMat = find_mat_for_timestamp_refined(projectRoot, runKey, ts, true);
    refMat = find_mat_for_timestamp_refined(projectRoot, runKey, ts, false);
    if strlength(origMat) == 0 || strlength(refMat) == 0
        fprintf("\nTop-1 cfg print skipped (missing MAT files) for %s | %s.\n", runKey, ts);
        return
    end
    Sorig = load_cfg_snapshot_refined(origMat);
    Sref = load_cfg_snapshot_refined(refMat);
    fprintf("\n=== Top-1 point by |delta J_TV|: %s | %s ===\n", runKey, ts);
    fprintf("\n--- Original cfg (%s) ---\n", origMat);
    print_cfg_struct_refined(Sorig);
    fprintf("\n--- Refined cfg (%s) ---\n", refMat);
    print_cfg_struct_refined(Sref);
end


function matPath = find_mat_for_timestamp_refined(projectRoot, runKey, ts, isOriginal)
%FIND_MAT_FOR_TIMESTAMP_REFINED Resolve original/refined MAT path for one controller.
    if isOriginal
        cand = fullfile(projectRoot, "results", runKey, "out_" + ts + ".mat");
    else
        cand = fullfile(projectRoot, "results", "final_fidelity_same_noise", runKey + "_full_f1_same_noise", "out_full_" + ts + ".mat");
    end
    if isfile(cand)
        matPath = string(cand);
    else
        matPath = "";
    end
end


function S = load_cfg_snapshot_refined(matPath)
%LOAD_CFG_SNAPSHOT_REFINED Load only cfg-relevant MAT fields.
    S = struct();
    vars = string(who("-file", matPath));
    if any(vars == "out")
        tmp = load(matPath, "out");
        S.out = tmp.out;
    end
    if any(vars == "cfg_run")
        tmp = load(matPath, "cfg_run");
        S.cfg_run = tmp.cfg_run;
    end
end


function print_cfg_struct_refined(S)
%PRINT_CFG_STRUCT_REFINED Print compact configuration information.
    if isfield(S, "out") && isfield(S.out, "cfg")
        cfg = S.out.cfg;
        if isfield(cfg, "f"), fprintf("f = %.12g\n", cfg.f); end
        if isfield(cfg, "m"), fprintf("m = %d\n", cfg.m); end
        if isfield(cfg, "p"), fprintf("p = %d\n", cfg.p); end
        if isfield(cfg, "Q"), fprintf("Q =\n"); disp(cfg.Q); end
        if isfield(cfg, "Ru"), fprintf("Ru =\n"); disp(cfg.Ru); end
        if isfield(cfg, "Rdu"), fprintf("Rdu =\n"); disp(cfg.Rdu); end
    end
    if isfield(S, "out") && isfield(S.out, "theta")
        fprintf("theta =\n");
        disp(S.out.theta);
    end
    if isfield(S, "cfg_run")
        cfgRun = S.cfg_run;
        if isfield(cfgRun, "sigma_y"), fprintf("sigma_y = "); disp(cfgRun.sigma_y); end
        if isfield(cfgRun, "mode"), fprintf("cfg_run.mode = %s\n", string(cfgRun.mode)); end
        if isfield(cfgRun, "source_root"), fprintf("cfg_run.source_root = %s\n", string(cfgRun.source_root)); end
        if isfield(cfgRun, "output_root"), fprintf("cfg_run.output_root = %s\n", string(cfgRun.output_root)); end
        if isfield(cfgRun, "results_csv"), fprintf("cfg_run.results_csv = %s\n", string(cfgRun.results_csv)); end
        if isfield(cfgRun, "out_dir"), fprintf("cfg_run.out_dir = %s\n", string(cfgRun.out_dir)); end
        if isfield(cfgRun, "NumWorkers"), fprintf("cfg_run.NumWorkers = %d\n", double(cfgRun.NumWorkers)); end
    end
end


function [runKey, ts] = split_match_key_refined(matchKey)
%SPLIT_MATCH_KEY_REFINED Split "runKey|timestamp" identifiers.
    n = numel(matchKey);
    runKey = strings(n, 1);
    ts = strings(n, 1);
    for i = 1:n
        parts = split(string(matchKey(i)), "|");
        runKey(i) = parts(1);
        ts(i) = parts(2);
    end
end


function h = plot_pareto_polyline_with_markers_refined(ax, x, y, colorVal, lineWidth, markerSize)
%PLOT_PARETO_POLYLINE_WITH_MARKERS_REFINED Draw frontier as one '-o' object.
    if isempty(x) || isempty(y)
        h = gobjects(0);
        return
    end
    [xSort, ord] = sort(x(:), "ascend");
    ySort = y(ord);
    h = plot(ax, xSort, ySort, "-o", ...
        "Color", colorVal, ...
        "LineWidth", lineWidth, ...
        "MarkerSize", markerSize, ...
        "MarkerFaceColor", "w", ...
        "MarkerEdgeColor", colorVal);
end


function write_promoted_report_refined(reportPath, promotedRunKey, promotedTs, T_orig, T_ref, promotedIdxOrig, promotedIdxRef)
%WRITE_PROMOTED_REPORT_REFINED Save promoted-point table for traceability.
    fid = fopen(reportPath, "w");
    if fid < 0
        warning("Could not write promoted-point report: %s", reportPath);
        return
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, "Promoted to refined Pareto frontier with original z < 1\n");
    fprintf(fid, "Original dataset excludes DOE rows (iterations 1..20 in run1/run2).\n");
    fprintf(fid, "Definition: original non-Pareto -> refined Pareto AND original theta_1 < 1\n");
    fprintf(fid, "Count: %d\n\n", numel(promotedTs));
    fprintf(fid, "run_key,timestamp,orig_SSE,orig_SSdU,orig_z,refined_SSE,refined_SSdU\n");
    for i = 1:numel(promotedTs)
        io = promotedIdxOrig(i);
        ir = promotedIdxRef(i);
        fprintf(fid, "%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g\n", ...
            promotedRunKey(i), promotedTs(i), ...
            double(T_orig.SSE(io)), double(T_orig.SSdU(io)), double(T_orig.z_eval(io)), ...
            double(T_ref.SSE(ir)), double(T_ref.SSdU(ir)));
    end
end


function report_final_frontier_f1_metrics(projectRoot, numericalFolder, optimizationStartIter)
%REPORT_FINAL_FRONTIER_F1_METRICS Build final Pareto f=1 table with noisy/noiseless metrics.
rootFolder = fullfile(projectRoot, "results");
runDefs = [
    struct("run_key","run1","csvPath", fullfile(rootFolder, "run1", "results.csv"));
    struct("run_key","run2","csvPath", fullfile(rootFolder, "run2", "results.csv"))
];

Tall = table();
for k = 1:numel(runDefs)
    [~, Topt, ~, ~] = load_results_table(runDefs(k).csvPath, optimizationStartIter);
    Topt.run_key = repmat(string(runDefs(k).run_key), height(Topt), 1);
    Tall = [Tall; Topt]; %#ok<AGROW>
end
if isempty(Tall)
    fprintf("Final Pareto f=1 metrics table: no optimization rows available.\n");
    return
end

isFront = compute_pareto_mask(double(Tall.SSE), double(Tall.SSdU));
Tf = Tall(isFront, :);
if isempty(Tf)
    fprintf("Final Pareto f=1 metrics table: no Pareto rows available.\n");
    return
end

n = height(Tf);
noisy_Jtrack = nan(n,1); noisy_JTV = nan(n,1);
noiseless_Jtrack = nan(n,1); noiseless_JTV = nan(n,1);

noisy_settle_x1_h_c1 = nan(n,1); noisy_settle_x2_h_c1 = nan(n,1); noisy_settle_x3_h_c1 = nan(n,1);
noisy_settle_x1_h_c2 = nan(n,1); noisy_settle_x2_h_c2 = nan(n,1); noisy_settle_x3_h_c2 = nan(n,1);
noisy_settle_case_h_c1 = nan(n,1); noisy_settle_case_h_c2 = nan(n,1);
noisy_IAE_x1_c1 = nan(n,1); noisy_IAE_x2_c1 = nan(n,1); noisy_IAE_x3_c1 = nan(n,1);
noisy_IAE_x1_c2 = nan(n,1); noisy_IAE_x2_c2 = nan(n,1); noisy_IAE_x3_c2 = nan(n,1);
noisy_IAE_case_c1 = nan(n,1); noisy_IAE_case_c2 = nan(n,1);

noiseless_settle_x1_h_c1 = nan(n,1); noiseless_settle_x2_h_c1 = nan(n,1); noiseless_settle_x3_h_c1 = nan(n,1);
noiseless_settle_x1_h_c2 = nan(n,1); noiseless_settle_x2_h_c2 = nan(n,1); noiseless_settle_x3_h_c2 = nan(n,1);
noiseless_settle_case_h_c1 = nan(n,1); noiseless_settle_case_h_c2 = nan(n,1);
noiseless_IAE_x1_c1 = nan(n,1); noiseless_IAE_x2_c1 = nan(n,1); noiseless_IAE_x3_c1 = nan(n,1);
noiseless_IAE_x1_c2 = nan(n,1); noiseless_IAE_x2_c2 = nan(n,1); noiseless_IAE_x3_c2 = nan(n,1);
noiseless_IAE_case_c1 = nan(n,1); noiseless_IAE_case_c2 = nan(n,1);

for i = 1:n
    run_key = string(Tf.run_key(i));
    ts = string(Tf.timestamp(i));
    noisyOutPath = fullfile(rootFolder, "final_fidelity_same_noise", run_key + "_full_f1_same_noise", "out_full_" + ts + ".mat");
    noiselessOutPath = fullfile(rootFolder, "test_run", run_key + "_full_f1_no_noise", "out_full_" + ts + ".mat");
    if ~isfile(noisyOutPath)
        error("Missing z=1 refined MAT for Pareto controller: %s | %s (%s)", run_key, ts, noisyOutPath);
    end
    if ~isfile(noiselessOutPath)
        warning("Missing no-noise MAT for Pareto controller: %s | %s (%s)", run_key, ts, noiselessOutPath);
    end

    noisyM = compute_out_metrics_by_path(noisyOutPath, 0.05);
    noiselessM = compute_out_metrics_by_path(noiselessOutPath, 0.05);

    noisy_Jtrack(i) = noisyM.Jtrack;
    noisy_JTV(i) = noisyM.JTV;
    noiseless_Jtrack(i) = noiselessM.Jtrack;
    noiseless_JTV(i) = noiselessM.JTV;

    noisy_settle_x1_h_c1(i) = noisyM.settle_h(1,1); noisy_settle_x2_h_c1(i) = noisyM.settle_h(1,2); noisy_settle_x3_h_c1(i) = noisyM.settle_h(1,3);
    noisy_settle_x1_h_c2(i) = noisyM.settle_h(2,1); noisy_settle_x2_h_c2(i) = noisyM.settle_h(2,2); noisy_settle_x3_h_c2(i) = noisyM.settle_h(2,3);
    noisy_settle_case_h_c1(i) = noisyM.settle_case_h(1); noisy_settle_case_h_c2(i) = noisyM.settle_case_h(2);
    noisy_IAE_x1_c1(i) = noisyM.IAE(1,1); noisy_IAE_x2_c1(i) = noisyM.IAE(1,2); noisy_IAE_x3_c1(i) = noisyM.IAE(1,3);
    noisy_IAE_x1_c2(i) = noisyM.IAE(2,1); noisy_IAE_x2_c2(i) = noisyM.IAE(2,2); noisy_IAE_x3_c2(i) = noisyM.IAE(2,3);
    noisy_IAE_case_c1(i) = noisyM.IAE_case(1); noisy_IAE_case_c2(i) = noisyM.IAE_case(2);

    noiseless_settle_x1_h_c1(i) = noiselessM.settle_h(1,1); noiseless_settle_x2_h_c1(i) = noiselessM.settle_h(1,2); noiseless_settle_x3_h_c1(i) = noiselessM.settle_h(1,3);
    noiseless_settle_x1_h_c2(i) = noiselessM.settle_h(2,1); noiseless_settle_x2_h_c2(i) = noiselessM.settle_h(2,2); noiseless_settle_x3_h_c2(i) = noiselessM.settle_h(2,3);
    noiseless_settle_case_h_c1(i) = noiselessM.settle_case_h(1); noiseless_settle_case_h_c2(i) = noiselessM.settle_case_h(2);
    noiseless_IAE_x1_c1(i) = noiselessM.IAE(1,1); noiseless_IAE_x2_c1(i) = noiselessM.IAE(1,2); noiseless_IAE_x3_c1(i) = noiselessM.IAE(1,3);
    noiseless_IAE_x1_c2(i) = noiselessM.IAE(2,1); noiseless_IAE_x2_c2(i) = noiselessM.IAE(2,2); noiseless_IAE_x3_c2(i) = noiselessM.IAE(2,3);
    noiseless_IAE_case_c1(i) = noiselessM.IAE_case(1); noiseless_IAE_case_c2(i) = noiselessM.IAE_case(2);
end

Treport = table( ...
    string(Tf.run_key), string(Tf.timestamp), double(Tf.p), double(Tf.m), ...
    double(Tf.Q_x1), double(Tf.Q_x2), double(Tf.Q_x3), ...
    double(Tf.R_u_x1), double(Tf.R_u_x2), double(Tf.R_u_x3), ...
    double(Tf.R_du_x1), double(Tf.R_du_x2), double(Tf.R_du_x3), ...
    noisy_Jtrack, noisy_JTV, ...
    noisy_settle_x1_h_c1, noisy_settle_x2_h_c1, noisy_settle_x3_h_c1, ...
    noisy_settle_x1_h_c2, noisy_settle_x2_h_c2, noisy_settle_x3_h_c2, ...
    noisy_settle_case_h_c1, noisy_settle_case_h_c2, ...
    noisy_IAE_x1_c1, noisy_IAE_x2_c1, noisy_IAE_x3_c1, ...
    noisy_IAE_x1_c2, noisy_IAE_x2_c2, noisy_IAE_x3_c2, ...
    noisy_IAE_case_c1, noisy_IAE_case_c2, ...
    noiseless_Jtrack, noiseless_JTV, ...
    noiseless_settle_x1_h_c1, noiseless_settle_x2_h_c1, noiseless_settle_x3_h_c1, ...
    noiseless_settle_x1_h_c2, noiseless_settle_x2_h_c2, noiseless_settle_x3_h_c2, ...
    noiseless_settle_case_h_c1, noiseless_settle_case_h_c2, ...
    noiseless_IAE_x1_c1, noiseless_IAE_x2_c1, noiseless_IAE_x3_c1, ...
    noiseless_IAE_x1_c2, noiseless_IAE_x2_c2, noiseless_IAE_x3_c2, ...
    noiseless_IAE_case_c1, noiseless_IAE_case_c2, ...
    'VariableNames', { ...
    'run_key','timestamp','p','m', ...
    'Q_x1','Q_x2','Q_x3','R_u_x1','R_u_x2','R_u_x3','R_du_x1','R_du_x2','R_du_x3', ...
    'noisy_Jtrack','noisy_JTV', ...
    'noisy_settle_x1_h_c1','noisy_settle_x2_h_c1','noisy_settle_x3_h_c1', ...
    'noisy_settle_x1_h_c2','noisy_settle_x2_h_c2','noisy_settle_x3_h_c2', ...
    'noisy_settle_case_h_c1','noisy_settle_case_h_c2', ...
    'noisy_IAE_x1_c1','noisy_IAE_x2_c1','noisy_IAE_x3_c1', ...
    'noisy_IAE_x1_c2','noisy_IAE_x2_c2','noisy_IAE_x3_c2', ...
    'noisy_IAE_case_c1','noisy_IAE_case_c2', ...
    'noiseless_Jtrack','noiseless_JTV', ...
    'noiseless_settle_x1_h_c1','noiseless_settle_x2_h_c1','noiseless_settle_x3_h_c1', ...
    'noiseless_settle_x1_h_c2','noiseless_settle_x2_h_c2','noiseless_settle_x3_h_c2', ...
    'noiseless_settle_case_h_c1','noiseless_settle_case_h_c2', ...
    'noiseless_IAE_x1_c1','noiseless_IAE_x2_c1','noiseless_IAE_x3_c1', ...
    'noiseless_IAE_x1_c2','noiseless_IAE_x2_c2','noiseless_IAE_x3_c2', ...
    'noiseless_IAE_case_c1','noiseless_IAE_case_c2'});

Treport = sortrows(Treport, 'noisy_JTV', 'descend');
fprintf("\nFinal Pareto frontier f=1 controller table (sorted by decreasing noisy J_TV):\n");
disp(Treport);

% Controllers simultaneously in top-3 (lowest) for noisy J_track and all individual noisy settling/IAE metrics.
noisyJtrackSel = double(Treport.noisy_Jtrack);
noisyJtrackSel(~isfinite(noisyJtrackSel)) = inf;

settleColNames = { ...
    'noisy_settle_x1_h_c1','noisy_settle_x2_h_c1','noisy_settle_x3_h_c1', ...
    'noisy_settle_x1_h_c2','noisy_settle_x2_h_c2','noisy_settle_x3_h_c2'};
iaeColNames = { ...
    'noisy_IAE_x1_c1','noisy_IAE_x2_c1','noisy_IAE_x3_c1', ...
    'noisy_IAE_x1_c2','noisy_IAE_x2_c2','noisy_IAE_x3_c2'};

nTop = min(3, height(Treport));
[~, ordJtrackSel] = sort(noisyJtrackSel, "ascend");
topJtrackMask = false(height(Treport),1);
topJtrackMask(ordJtrackSel(1:nTop)) = true;
topSettleAllMask = true(height(Treport),1);
topIAEAllMask = true(height(Treport),1);

for c = 1:numel(settleColNames)
    v = double(Treport.(settleColNames{c}));
    v(~isfinite(v)) = inf;
    [~, ord] = sort(v, "ascend");
    m = false(height(Treport),1);
    m(ord(1:nTop)) = true;
    topSettleAllMask = topSettleAllMask & m;
end
for c = 1:numel(iaeColNames)
    v = double(Treport.(iaeColNames{c}));
    v(~isfinite(v)) = inf;
    [~, ord] = sort(v, "ascend");
    m = false(height(Treport),1);
    m(ord(1:nTop)) = true;
    topIAEAllMask = topIAEAllMask & m;
end
topAllMask = topJtrackMask & topSettleAllMask & topIAEAllMask;

showCols = [{'run_key','timestamp','noisy_Jtrack'}, settleColNames, iaeColNames];
TtopAll = Treport(topAllMask, showCols);
if ~isempty(TtopAll)
    TtopAll = sortrows(TtopAll, 'noisy_Jtrack', 'ascend');
end
fprintf("Controllers simultaneously top-3 lowest in noisy J_track and all individual noisy settling/IAE metrics (NaN treated as Inf):\n");
disp(TtopAll);

[benchJtrack, benchJTV, benchFound] = load_noisy_benchmark_point(projectRoot);
Tf1Bench = table();
if benchFound
    [nSuperiorF1, isSuperiorF1] = count_benchmark_strict_superior_refined( ...
        Treport.noisy_Jtrack, Treport.noisy_JTV, benchJtrack, benchJTV);
    fprintf("Final Pareto f=1 controllers strictly better than noisy benchmark (J_track and J_TV both lower): %d/%d\n", ...
        nSuperiorF1, height(Treport));
    print_benchmark_ratio_stats_refined("Final Pareto f=1", ...
        Treport.noisy_Jtrack, Treport.noisy_JTV, isSuperiorF1, benchJtrack, benchJTV);
    if nSuperiorF1 > 0
        fprintf("Timestamps strictly better than benchmark:\n");
        disp(Treport.timestamp(isSuperiorF1));
    end
    Tf1Bench = build_benchmark_comparison_table_refined( ...
        string(Treport.run_key), string(Treport.timestamp), ...
        double(Treport.noisy_Jtrack), double(Treport.noisy_JTV), benchJtrack, benchJTV);
    fprintf("Final Pareto f=1 benchmark comparison table:\n");
    disp(Tf1Bench);
else
    warning("Benchmark superiority verification skipped: noisy benchmark not available.");
end

if ~isfolder(numericalFolder)
    mkdir(numericalFolder);
end
writetable(Treport, fullfile(numericalFolder, "final_pareto_frontier_f1_noisy_noiseless_metrics.csv"));

outTxtDir = fullfile(projectRoot, "results", "txt results");
if ~isfolder(outTxtDir)
    mkdir(outTxtDir);
end
if benchFound
    benchTxtPath = fullfile(outTxtDir, "final_pareto_f1_benchmark_comparison_table.txt");
    writetable(Tf1Bench, benchTxtPath, "FileType", "text", "Delimiter", "\t");
    fprintf("Saved: %s\n", benchTxtPath);
end
end


function M = compute_out_metrics_by_path(matPath, settlingTol)
%COMPUTE_OUT_METRICS_BY_PATH Load one out_full MAT and compute objective/settling/IAE metrics.
M = struct();
M.Jtrack = nan;
M.JTV = nan;
M.settle_h = nan(2,3);
M.settle_case_h = nan(2,1);
M.IAE = nan(2,3);
M.IAE_case = nan(2,1);
if ~isfile(matPath)
    return
end
S = load(matPath, "out");
if ~isfield(S, "out")
    return
end
out = S.out;
if isfield(out, "SSE"), M.Jtrack = double(out.SSE); end
if isfield(out, "SSdU"), M.JTV = double(out.SSdU); end
if ~isfield(out, "case") || isempty(out.case)
    return
end
nCase = min(numel(out.case), 2);
for c = 1:nCase
    caseData = out.case(c);
    [settle_h, iae] = summarize_case_metrics_simple(caseData, settlingTol);
    nState = min(numel(settle_h), 3);
    M.settle_h(c,1:nState) = settle_h(1:nState);
    M.IAE(c,1:nState) = iae(1:nState);
    if any(isfinite(settle_h))
        M.settle_case_h(c) = max(settle_h(isfinite(settle_h)));
    else
        M.settle_case_h(c) = nan;
    end
    M.IAE_case(c) = sum(iae, "omitnan");
end
end


function [settlingTimes_h, IAEByState] = summarize_case_metrics_simple(caseStruct, settlingTol)
%SUMMARIZE_CASE_METRICS_SIMPLE Compute settling times and IAE from one case struct.
if ~isfield(caseStruct, "Y") || ~isfield(caseStruct, "Ysp")
    settlingTimes_h = nan(1,3);
    IAEByState = nan(1,3);
    return
end
Y = double(caseStruct.Y);
Ysp = double(caseStruct.Ysp);
nState = min(size(Y,2), size(Ysp,2));
Y = Y(:,1:nState);
Ysp = Ysp(:,1:nState);
if isfield(caseStruct, "dt")
    dt = double(caseStruct.dt);
elseif isfield(caseStruct, "tf")
    dt = double(caseStruct.tf) / max(size(Y,1)-1, 1);
else
    dt = 1/60;
end
t = (0:size(Y,1)-1).' * dt;
settlingTimes_h = nan(1,nState);
IAEByState = sum(abs(Y - Ysp), 1) * dt;
epsRef = 1e-9;
for i = 1:nState
    refVal = Ysp(end, i);
    relErr = abs(Y(:,i) - Ysp(:,i)) / max(abs(refVal), epsRef);
    if relErr(end) > settlingTol
        settlingTimes_h(i) = nan;
        continue
    end
    settleIdx = find_settling_index(relErr, settlingTol);
    if isempty(settleIdx)
        settlingTimes_h(i) = nan;
    else
        settlingTimes_h(i) = t(settleIdx);
    end
end
if nState < 3
    settlingTimes_h(1, end+1:3) = nan;
    IAEByState(1, end+1:3) = nan;
end
end


function idx = find_settling_index(relErr, tol)
%FIND_SETTLING_INDEX First index where relative error stays below tol thereafter.
n = numel(relErr);
idx = [];
for k = 1:n
    if all(relErr(k:end) <= tol)
        idx = k;
        return
    end
end
end
