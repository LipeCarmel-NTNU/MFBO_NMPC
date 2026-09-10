function main_pareto()
% MAIN_PARETO
% Post-processing hub for MFBO-NMPC optimization outcomes (preprocessed pipeline).
%
% Data contract
% - Input: Result analysis/storage/preprocessed.mat, built by
%   initial_preprocessing.m (run parse_registry.py first). Variables used:
%     evals : table, one row per BO evaluation (analysis object)
%     doe   : table, one row per DOE initialisation evaluation (same schema)
%     meta  : provenance struct
% - Cases are whatever categories `evals.case` carries (driven by cases.txt),
%   e.g. results_case1, results_case2, and later a baseline campaign.
% - Costs: SSE (J_track) and SSdU (J_TV) are the extrapolated full-horizon
%   costs the BO loop optimizes. SSE_measured/SSdU_measured exist in the
%   schema but are not used here.
%
% Core policy
% - DOE exclusion is structural: Pareto/optimal-controller analyses use only
%   `evals` (the BO phase). `doe` enters runtime accounting and the iteration
%   timelines only.
%
% Guarded sections (skipped with a warning when their data is absent)
% - Refined frontier change: z=1 re-evaluation CSVs expected at
%   results/final_fidelity_same_noise/<case>_full_f1_same_noise/results_full.csv
% - Noisy benchmark overlay: results/benchmark_reference_controller/
%   benchmark_full_f1_same_noise_fix/out_benchmark.mat
% - Final-frontier f=1 metrics: out_full_<ts>.mat under the refined folder and
%   results/test_run/<case>_full_f1_no_noise (no-noise counterpart).
%
% Main outputs
% - Figures in results/graphical_results (PNG 300 dpi + vector PDF).
% - Numerical summaries in results/numerical results
%   (guarded sections also write to results/txt results).
% - Console diagnostics for runtime split and Pareto composition.
%
% Color scheme: see COLOR_SCHEME.md (Wong categorical + navia sequential).

%% Dependencies
close all; clc
here      = fileparts(mfilename('fullpath'));
repo_root = fileparts(here);
addpath(genpath(fullfile(repo_root, 'dependencies')));

%% Paths
resultsRoot  = fullfile(repo_root, 'results');
graphicsDir  = fullfile(resultsRoot, 'graphical_results');
numericalDir = fullfile(resultsRoot, 'numerical results');
if ~isfolder(graphicsDir);  mkdir(graphicsDir);  end
if ~isfolder(numericalDir); mkdir(numericalDir); end

%% Load preprocessed data
storePath = fullfile(here, 'storage', 'preprocessed.mat');
if ~isfile(storePath)
    error('main_pareto:noData', ...
        'Missing %s. Run parse_registry.py, then initial_preprocessing.m.', storePath);
end
S     = load(storePath, 'evals', 'doe', 'meta');
evals = S.evals;
doe   = S.doe;
meta  = S.meta; %#ok<NASGU>
if isempty(evals)
    error('main_pareto:emptyEvals', 'evals table is empty in %s.', storePath);
end

%% Style (COLOR_SCHEME.md)
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
fontSize = 20;

plotColors  = nature_methods_colors(3);   % Blue (Case 1), BluishGreen (Case 2), ReddishPurple (accent)
accentColor = plotColors(3, :);
caseMarkers = ["o", "^", "d"];
caseLines   = ["-", "-.", ":"];
seqMap      = load_navia_colormap(256);   % sequential map for fidelity z

%% Per-case split (BO-phase evals; DOE kept separately for runtime/timelines)
evals.case = removecats(categorical(evals.case));
doe.case   = removecats(categorical(doe.case));
caseNames  = categories(evals.case);
nCases     = numel(caseNames);
if nCases > numel(caseMarkers)
    error('main_pareto:tooManyCases', ...
        'Only %d case styles defined; found %d cases.', numel(caseMarkers), nCases);
end
caseLabels = arrayfun(@pretty_case, string(caseNames));

E        = cell(nCases, 1);   % BO evaluations per case
D        = cell(nCases, 1);   % DOE evaluations per case
Tp       = cell(nCases, 1);   % per-case Pareto subset (sorted by SSE desc)
isPareto = cell(nCases, 1);

for k = 1:nCases
    E{k} = sortrows(evals(evals.case == caseNames{k}, :), 'iter');
    D{k} = sortrows(doe(doe.case == caseNames{k},     :), 'iter');
    if height(D{k}) ~= 20
        warning('main_pareto:doeCount', ...
            '%s: expected 20 DOE rows, found %d.', caseLabels(k), height(D{k}));
    end
    if height(E{k}) ~= 100
        warning('main_pareto:boCount', ...
            '%s: expected 100 BO rows, found %d.', caseLabels(k), height(E{k}));
    end

    isPareto{k} = compute_pareto_mask(double(E{k}.SSE), double(E{k}.SSdU));
    Tp{k}       = E{k}(isPareto{k}, :);
    [~, ord]    = sort(double(Tp{k}.SSE), 'descend');
    Tp{k}       = Tp{k}(ord, :);

    fprintf('%s: %d BO evaluations, %d on the case frontier (DOE excluded structurally).\n', ...
        caseLabels(k), height(E{k}), height(Tp{k}));
    display_pareto_table(Tp{k}, caseLabels(k));
end

%% Runtime phase summaries (DOE vs BO, from t_total)
runtimeSummaryTables = cell(nCases, 1);
for k = 1:nCases
    runtimeSummaryTables{k} = display_runtime_phase_summary(D{k}, E{k}, caseLabels(k));
end
runtimeSummaryCombined = display_runtime_phase_summary(doe, evals, "All cases (combined)");

write_runtime_and_parameter_summary(E, Tp, caseLabels, ...
    runtimeSummaryTables, runtimeSummaryCombined, numericalDir);

%% Shared axis limits from the data (log-padded, shared across panels)
xLimAll = padded_log_limits(double(evals.SSdU), 1.25);
yLimAll = padded_log_limits(double(evals.SSE),  1.25);
zLo     = min(0.5, min(double(evals.z)));       % colorbar floor for fidelity

%% Figure 1: SSE vs SSdU per case, color mapped by fidelity z
fig1 = figure('Color', 'w', 'Name', 'Pareto SSE vs SSdU by Case');
tiledlayout(fig1, 1, nCases, 'Padding', 'compact', 'TileSpacing', 'compact');
for k = 1:nCases
    T  = E{k};
    ax = nexttile; hold(ax, 'on');
    plot_pareto_continuum(ax, double(T.SSdU(isPareto{k})), double(T.SSE(isPareto{k})), ...
        accentColor, xLimAll, yLimAll);
    scatter(ax, double(T.SSdU), double(T.SSE), 80, double(T.z), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
    scatter(ax, double(T.SSdU(isPareto{k})), double(T.SSE(isPareto{k})), 170, accentColor, ...
        'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', accentColor, 'LineWidth', 1.2);

    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', fontSize);
    xlim(ax, xLimAll); ylim(ax, yLimAll);
    colormap(ax, seqMap);
    caxis(ax, [zLo, 1]);
    xlabel(ax, '$J_{\mathrm{TV}}$');
    ylabel(ax, '$J_{\mathrm{track}}$');
    title(ax, "$\mathbf{" + char('a' + k - 1) + "}$", 'Interpreter', 'latex');
    ax.TitleHorizontalAlignment = 'left';
    grid(ax, 'off'); box(ax, 'off');
    cb = colorbar(ax);
    cb.Label.String = '$z$ (dimensionless)';
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = fontSize;
end
save_plot_outputs(fig1, fullfile(graphicsDir, 'sse_vs_ssdu_side_by_side_z.png'), ...
    fontSize, 600 * nCases, 460);

%% Figure 2: iteration runtime + fidelity per case (DOE + BO timeline)
tMaxH = 0;
for k = 1:nCases
    tMaxH = max(tMaxH, max([double(D{k}.t_total); double(E{k}.t_total)]) / 3600);
end
fig2 = figure('Color', 'w', 'Name', 'Iteration Runtime and Fidelity by Case');
tiledlayout(fig2, 1, nCases, 'Padding', 'compact', 'TileSpacing', 'compact');
for k = 1:nCases
    nDoe    = height(D{k});
    iterAll = [double(D{k}.iter); nDoe + double(E{k}.iter)];
    tAllH   = [double(D{k}.t_total); double(E{k}.t_total)] / 3600;
    zAll    = [double(D{k}.z); double(E{k}.z)];

    ax = nexttile; hold(ax, 'on');
    yyaxis(ax, 'left');
    plot(ax, iterAll, tAllH, '-', 'LineWidth', 2.0, 'Color', plotColors(k, :));
    ax.YColor = 'k';
    ylim(ax, [0, 1.05 * max(tMaxH, eps)]);
    xline(ax, nDoe, '--', 'LineWidth', 2.0, 'Color', 'k', 'Alpha', 1);
    ylabel(ax, '$t_{\mathrm{iter}}$ (h)');
    yyaxis(ax, 'right');
    plot(ax, iterAll, zAll, 'o', 'LineWidth', 2.0, 'MarkerSize', 4, 'Color', accentColor);
    ax.YColor = 'k';
    ylim(ax, [0, 1]);
    ylabel(ax, '$z$ (dimensionless)');
    xlabel(ax, '$k$ (iteration)');
    xlim(ax, [1, max(1, max(iterAll))]);
    title(ax, "$\mathbf{" + char('a' + k - 1) + "}$", 'Interpreter', 'latex');
    ax.TitleHorizontalAlignment = 'left';
    set(ax, 'FontSize', fontSize);
    grid(ax, 'off'); box(ax, 'off');
    axes(ax); %#ok<LAXES>
    yyaxis(ax, 'left');  format_tick(0, 1);
    yyaxis(ax, 'right'); yticks(ax, 0:0.2:1); format_tick(0, 1);
end
save_plot_outputs(fig2, fullfile(graphicsDir, 'runtime_vs_iteration_side_by_side.png'), ...
    fontSize, 600 * nCases, 460);

%% Figure 3: combined Pareto samples across all cases
fig3 = figure('Color', 'w', 'Name', 'Combined Pareto Samples');
ax = axes(fig3); hold(ax, 'on');
plot_combined_pareto_base(ax, E, Tp, plotColors, caseMarkers, accentColor, xLimAll, yLimAll);
apply_combined_axes_style(ax, fontSize, xLimAll, yLimAll);
save_plot_outputs(fig3, fullfile(graphicsDir, 'pareto_samples_combined.png'), fontSize, 920, 520);

%% Figure 4: cumulative runtime per case
fig4 = figure('Color', 'w', 'Name', 'Cumulative Runtime by Case');
ax = axes(fig4); hold(ax, 'on');
xMax = 1;
for k = 1:nCases
    nDoe    = height(D{k});
    iterAll = [double(D{k}.iter); nDoe + double(E{k}.iter)];
    cumH    = cumsum([double(D{k}.t_total); double(E{k}.t_total)], 'omitnan') / 3600;
    plot(ax, iterAll, cumH, caseLines(k), 'LineWidth', 2.0, ...
        'Color', plotColors(k, :), 'DisplayName', caseLabels(k));
    xMax = max(xMax, max(iterAll));
end
xline(ax, height(D{1}), '--', 'LineWidth', 2.0, 'Color', 'k', 'Alpha', 1, ...
    'HandleVisibility', 'off');
xlabel(ax, '$k$ (iteration)');
ylabel(ax, '$t_{\mathrm{run}}$ (h)');
xlim(ax, [1, xMax]);
legend(ax, 'Location', 'northwest');
set(ax, 'FontSize', fontSize);
grid(ax, 'off'); box(ax, 'off');
save_plot_outputs(fig4, fullfile(graphicsDir, 'runtime_cumulative.png'), fontSize, 920, 520);

%% Guarded: refined frontier comparison and final f=1 metrics
run_refined_frontier_change(E, caseNames, caseLabels, repo_root, graphicsDir, ...
    plotColors, caseMarkers, accentColor, fontSize);
report_final_frontier_f1_metrics(E, caseNames, repo_root, numericalDir);
end

%% ===================== CORE FUNCTIONS =====================

function label = pretty_case(name)
%PRETTY_CASE Turn a folder name like results_case1 into "Case 1".
tok = regexp(name, 'case(\d+)\s*$', 'tokens', 'once');
if isempty(tok)
    label = string(strrep(strrep(name, 'results_', ''), '_', ' '));
else
    label = "Case " + string(tok{1});
end
end


function ts = ts_string(dt)
%TS_STRING Filename timestamp string yyyyMMdd_HHmmss from a datetime column.
ts = string(datetime(dt, 'Format', 'yyyyMMdd_HHmmss'));
end


function isPareto = compute_pareto_mask(J_track, J_TV)
%COMPUTE_PARETO_MASK Mark non-dominated tradeoff points (lower is better).
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


function lim = padded_log_limits(vals, padFactor)
%PADDED_LOG_LIMITS Data-driven limits for a positive log axis.
if nargin < 2 || ~isfinite(padFactor) || padFactor <= 1
    padFactor = 1.25;
end
vals = vals(isfinite(vals) & (vals > 0));
if isempty(vals)
    lim = [1e-2, 1e2];
    return
end
lim = [min(vals) / padFactor, max(vals) * padFactor];
end


function display_pareto_table(Tp, caseLabel)
%DISPLAY_PARETO_TABLE Print BO-phase Pareto controller settings for one case.
tuningCols = ["id", "iter", "SSE", "SSdU", "z", "Np", "Nc", "t_total", ...
    "Q1", "Q2", "Q3", "Ru1", "Ru2", "Ru3", "Rdu1", "Rdu2", "Rdu3"];
tuningCols = tuningCols(ismember(tuningCols, string(Tp.Properties.VariableNames)));
disp(caseLabel + " Pareto frontier with tuning weights (BO phase only, DOE excluded):");
disp(Tp(:, tuningCols));
end


function summaryTbl = display_runtime_phase_summary(Ddoe, Ebo, runName)
%DISPLAY_RUNTIME_PHASE_SUMMARY Runtime share spent in DOE vs BO (t_total).
doeMin   = sum(double(Ddoe.t_total), 'omitnan') / 60;
boMin    = sum(double(Ebo.t_total),  'omitnan') / 60;
totalMin = doeMin + boMin;

if totalMin > 0
    doePct = 100 * doeMin / totalMin;
    boPct  = 100 * boMin  / totalMin;
else
    doePct = NaN;
    boPct  = NaN;
end

phase             = {'DOE'; 'Optimisation'; 'Total'};
iterations        = [height(Ddoe); height(Ebo); height(Ddoe) + height(Ebo)];
runtime_min       = [doeMin; boMin; totalMin];
runtime_pct_total = [doePct; boPct; 100];

summaryTbl = table(phase, iterations, runtime_min, runtime_pct_total, ...
    'VariableNames', {'phase', 'iterations', 'runtime_min', 'runtime_pct_total'});

disp("Runtime phase summary - " + string(runName) + ":");
disp(summaryTbl);
end


function write_runtime_and_parameter_summary(E, Tp, caseLabels, ...
    runSummaryTables, combinedSummaryTable, outDir)
%WRITE_RUNTIME_AND_PARAMETER_SUMMARY Persist manuscript-facing runtime/tuning summaries.
nCases  = numel(E);
outPath = fullfile(outDir, 'runtime_and_params.txt');
fid = fopen(outPath, 'w');
if fid == -1
    warning('Unable to write numerical summary: %s', outPath);
    return
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

for k = 1:nCases
    fprintf(fid, 'Runtime phase summary - %s:\n', char(caseLabels(k)));
    write_runtime_table(fid, runSummaryTables{k});
    fprintf(fid, '\n');
end
fprintf(fid, 'Runtime phase summary - All cases (combined):\n');
write_runtime_table(fid, combinedSummaryTable);
fprintf(fid, '\n');

% Mean fidelity z during the BO phase.
Eall  = vertcat(E{:});
meanZ = cellfun(@(T) mean(double(T.z), 'omitnan'), E);
fprintf('Mean fidelity z during BO evaluations:\n');
fprintf(fid, 'Mean fidelity z during BO evaluations:\n');
for k = 1:nCases
    fprintf('  %s: %.6g\n', caseLabels(k), meanZ(k));
    fprintf(fid, '  %s: %.6g\n', char(caseLabels(k)), meanZ(k));
end
meanZAll = mean(double(Eall.z), 'omitnan');
fprintf('  All cases (combined): %.6g\n', meanZAll);
fprintf(fid, '  All cases (combined): %.6g\n\n', meanZAll);

% Share of BO points with N_c = 1.
fprintf('Percentage of BO points with N_c = 1:\n');
fprintf(fid, 'Percentage of BO points with N_c = 1:\n');
nc1Lines = strings(nCases + 1, 1);
for k = 1:nCases
    n1 = nnz(double(E{k}.Nc) == 1);
    nT = height(E{k});
    nc1Lines(k) = sprintf('  %s: %d/%d (%.6g%%)', caseLabels(k), n1, nT, 100 * n1 / max(nT, 1));
    fprintf('%s\n', nc1Lines(k));
    fprintf(fid, '%s\n', nc1Lines(k));
end
n1 = nnz(double(Eall.Nc) == 1);
nT = height(Eall);
nc1Lines(end) = sprintf('  All cases (combined): %d/%d (%.6g%%)', n1, nT, 100 * n1 / max(nT, 1));
fprintf('%s\n', nc1Lines(end));
fprintf(fid, '%s\n\n', nc1Lines(end));

% Compact standalone txt for manuscript bookkeeping.
nc1Path = fullfile(outDir, 'optimization_nc1_share.txt');
fidNc1 = fopen(nc1Path, 'w');
if fidNc1 ~= -1
    cleanupNc1 = onCleanup(@() fclose(fidNc1)); %#ok<NASGU>
    fprintf(fidNc1, 'Percentage of BO points with N_c = 1:\n');
    fprintf(fidNc1, '%s\n', nc1Lines);
else
    warning('Unable to write N_c=1 share summary: %s', nc1Path);
end

% Per-case Pareto horizon composition.
for k = 1:nCases
    T = Tp{k};
    nRows   = height(T);
    countNc1 = nnz(double(T.Nc) == 1);
    countNp1 = nnz(double(T.Np) == 1);
    fprintf(fid, '%s Pareto counts (BO phase only, DOE excluded):\n', char(caseLabels(k)));
    fprintf(fid, '  N_c = 1: %d/%d\n', countNc1, nRows);
    fprintf(fid, '  N_p = 1: %d/%d\n', countNp1, nRows);

    ssdU    = double(T.SSdU);
    minSSdU = min(ssdU, [], 'omitnan');
    topMask = abs(ssdU - minSSdU) <= 10 * eps(max(1, abs(minSSdU)));
    fprintf(fid, '  Top (lowest SSdU) rows: %d\n', nnz(topMask));
    fprintf(fid, '  Top (lowest SSdU) with N_p = 1 and N_c = 1: %d\n\n', ...
        nnz(topMask & double(T.Np) == 1 & double(T.Nc) == 1));
end
fprintf('Saved: %s\n', outPath);
end


function write_runtime_table(fid, summaryTbl)
%WRITE_RUNTIME_TABLE Write one runtime summary table to disk.
for i = 1:height(summaryTbl)
    fprintf(fid, '  %s | iterations=%d | runtime_min=%.6g | runtime_pct_total=%.6g\n', ...
        char(string(summaryTbl.phase{i})), summaryTbl.iterations(i), ...
        summaryTbl.runtime_min(i), summaryTbl.runtime_pct_total(i));
end
end

%% ===================== PLOTTING HELPERS =====================

function h = plot_pareto_continuum(ax, x, y, curveColor, xBounds, yBounds)
%PLOT_PARETO_CONTINUUM Smooth monotone visual guide of the frontier trend.
[xSort, ord] = sort(x(:), 'ascend');
ySort = y(ord);

if numel(xSort) < 3
    h = plot(ax, xSort, ySort, '-o', ...
        'Color', curveColor, 'LineWidth', 2.0, ...
        'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', curveColor);
    return
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

h = plot(ax, xq, yq, '-', 'Color', curveColor, 'LineWidth', 2.0);
plot(ax, xSort, ySort, 'o', ...
    'Color', curveColor, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', curveColor, 'LineWidth', 1.4);
end


function h = plot_pareto_continuum_line_only(ax, x, y, curveColor, xBounds, yBounds)
%PLOT_PARETO_CONTINUUM_LINE_ONLY Smooth pchip Pareto guide without point markers.
[xSort, ord] = sort(x(:), 'ascend');
ySort = y(ord);

if numel(xSort) < 3
    h = plot(ax, xSort, ySort, '-', 'Color', curveColor, 'LineWidth', 2.0);
    return
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

h = plot(ax, xq, yq, '-', 'Color', curveColor, 'LineWidth', 2.0);
end


function [finalMask, Tall] = plot_combined_pareto_base(ax, E, Tp, plotColors, caseMarkers, accentColor, xBounds, yBounds)
%PLOT_COMBINED_PARETO_BASE All-case sample cloud + per-case frontiers + combined frontier.
Tall = vertcat(E{:});   % BO evaluations only; DOE excluded structurally
scatter(ax, double(Tall.SSdU), double(Tall.SSE), 18, ...
    'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'Optimization samples');

finalMask = compute_pareto_mask(double(Tall.SSE), double(Tall.SSdU));
Tf = Tall(finalMask, :);
plot_pareto_continuum(ax, double(Tf.SSdU), double(Tf.SSE), accentColor, xBounds, yBounds);

for k = 1:numel(Tp)
    scatter(ax, double(Tp{k}.SSdU), double(Tp{k}.SSE), 80 + 10 * (k - 1), plotColors(k, :), ...
        caseMarkers(k), 'MarkerFaceColor', plotColors(k, :), ...
        'MarkerEdgeColor', plotColors(k, :), 'LineWidth', 1.4);
end

scatter(ax, double(Tf.SSdU), double(Tf.SSE), 300, accentColor, ...
    'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', accentColor, 'LineWidth', 2);
end


function plot_combined_samples_no_guide(ax, E, Tp, plotColors, caseMarkers)
%PLOT_COMBINED_SAMPLES_NO_GUIDE Sample cloud + per-case frontier markers, no guideline/ring.
Tall = vertcat(E{:});
scatter(ax, double(Tall.SSdU), double(Tall.SSE), 18, ...
    'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
for k = 1:numel(Tp)
    scatter(ax, double(Tp{k}.SSdU), double(Tp{k}.SSE), 80 + 10 * (k - 1), plotColors(k, :), ...
        caseMarkers(k), 'MarkerFaceColor', plotColors(k, :), ...
        'MarkerEdgeColor', plotColors(k, :), 'LineWidth', 1.4);
end
end


function apply_combined_axes_style(ax, fontSize, xLim, yLim)
%APPLY_COMBINED_AXES_STYLE Shared log-log axes/tick styling for combined views.
set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', fontSize);
xlim(ax, xLim);
ylim(ax, yLim);
xlabel(ax, '$J_{\mathrm{TV}}$');
ylabel(ax, '$J_{\mathrm{track}}$');
grid(ax, 'off');
box(ax, 'off');
axes(ax);
format_tick(1, 1);
end


function cmap = load_navia_colormap(n)
%LOAD_NAVIA_COLORMAP Load the sequential map used for fidelity-color plots.
matPath = which('navia.mat');
if strlength(string(matPath)) == 0
    error('Unable to locate navia.mat on MATLAB path.');
end
S = load(matPath, 'navia');
if ~isfield(S, 'navia')
    error('File %s does not contain variable ''navia''.', matPath);
end
cmap = double(S.navia);
if size(cmap, 2) ~= 3
    error('Variable ''navia'' in %s must have 3 columns (RGB).', matPath);
end
if size(cmap, 1) ~= n
    x  = linspace(0, 1, size(cmap, 1));
    xq = linspace(0, 1, n);
    cmap = interp1(x, cmap, xq, 'linear');
end
cmap = min(max(cmap, 0), 1);
end


function save_plot_outputs(figHandle, pngPath, fontSize, figWidthPx, figHeightPx)
%SAVE_PLOT_OUTPUTS Export standardized figure files (PNG 300 dpi + vector PDF).
arguments
    figHandle
    pngPath
    fontSize (1,1) double = 14
    figWidthPx (1,1) double = 900
    figHeightPx (1,1) double = 500
end
figure(figHandle);
set_fig_size(figWidthPx, figHeightPx);
set_font_size(fontSize);
exportgraphics(figHandle, pngPath, 'Resolution', 300);
[folderPath, fileStem] = fileparts(pngPath);
pdfPath = fullfile(folderPath, strcat(fileStem, '.pdf'));
save_figure(pdfPath, NaN, false);
end

%% ===================== GUARDED: REFINED FRONTIER =====================

function run_refined_frontier_change(E, caseNames, caseLabels, projectRoot, graphicsFolder, ...
    plotColors, caseMarkers, accentColor, fontSize)
%RUN_REFINED_FRONTIER_CHANGE Compare original (BO) and refined (z=1 re-eval) Pareto points.
% Guarded: skips with a warning when no refined result CSVs exist.
nCases = numel(E);
resultsRoot = fullfile(projectRoot, 'results');

refinedFiles = strings(nCases, 1);
for k = 1:nCases
    refinedFiles(k) = fullfile(resultsRoot, 'final_fidelity_same_noise', ...
        string(caseNames{k}) + "_full_f1_same_noise", 'results_full.csv');
end
present = arrayfun(@isfile, refinedFiles);
if ~any(present)
    warning('main_pareto:noRefined', ...
        ['Refined frontier section skipped: no z=1 re-evaluation results found.\n' ...
         'Expected e.g. %s'], refinedFiles(1));
    return
end
if ~all(present)
    warning('main_pareto:partialRefined', ...
        'Refined results missing for %d case(s); comparing the available ones only.', ...
        nnz(~present));
end

% Original pool: all BO evaluations, keyed case|timestamp.
T_orig = table();
for k = 1:nCases
    Tk = E{k};
    Tk.run_key   = repmat(string(caseNames{k}), height(Tk), 1);
    Tk.run_label = repmat(caseLabels(k), height(Tk), 1);
    Tk.ts        = ts_string(Tk.timestamp);
    Tk.match_key = Tk.run_key + "|" + Tk.ts;
    Tk.z_eval    = double(Tk.z);
    T_orig = [T_orig; Tk(:, ["run_key", "run_label", "match_key", "ts", ...
        "iter", "SSE", "SSdU", "z_eval"])]; %#ok<AGROW>
end

% Refined pool from the re-evaluation CSVs.
T_ref = table();
for k = 1:nCases
    if ~present(k); continue; end
    Tk = readtable(refinedFiles(k), 'TextType', 'string');
    required = ["timestamp", "SSE", "SSdU"];
    for c = required
        if ~ismember(c, string(Tk.Properties.VariableNames))
            error('Missing column ''%s'' in %s', c, refinedFiles(k));
        end
    end
    Tk.run_key   = repmat(string(caseNames{k}), height(Tk), 1);
    Tk.run_label = repmat(caseLabels(k), height(Tk), 1);
    Tk.ts        = string(Tk.timestamp);
    Tk.match_key = Tk.run_key + "|" + Tk.ts;
    T_ref = [T_ref; Tk(:, ["run_key", "run_label", "match_key", "ts", "SSE", "SSdU"])]; %#ok<AGROW>
end
if isempty(T_ref)
    warning('main_pareto:emptyRefined', 'Refined result files were found but contained no rows.');
    return
end

isParetoOrig = compute_pareto_mask(double(T_orig.SSE), double(T_orig.SSdU));
[T_ref, refFilterInfo] = keep_refined_from_original_pareto(T_ref, T_orig, isParetoOrig);
print_refined_filter_summary(refFilterInfo);
if isempty(T_ref)
    warning('main_pareto:noMatchedRefined', ...
        'No refined rows match original Pareto points; refined section skipped.');
    return
end
origParetoKeys = unique(string(T_orig.match_key(isParetoOrig)), 'stable');
missingRefined = setdiff(origParetoKeys, unique(string(T_ref.match_key), 'stable'), 'stable');
if ~isempty(missingRefined)
    warning('main_pareto:missingRefined', ...
        'Missing z=1 refined results for %d Pareto controller(s): %s', ...
        numel(missingRefined), strjoin(missingRefined, ', '));
end

T_jtv = compute_change_table(T_orig, T_ref, 'SSdU', 'JTV');
fprintf('\n=== J_TV change (original -> refined), sorted by |delta %%| descending ===\n');
disp(T_jtv);
T_jtrack = compute_change_table(T_orig, T_ref, 'SSE', 'Jtrack');
fprintf('\n=== J_track change (original -> refined), sorted by |delta %%| descending ===\n');
disp(T_jtrack);
print_top1_cfg(projectRoot, T_jtv);

% Promotion bookkeeping.
isParetoRef = compute_pareto_mask(double(T_ref.SSE), double(T_ref.SSdU));
[commonKeys, idxOrig, idxRef] = intersect(string(T_orig.match_key), string(T_ref.match_key), 'stable');
[matchedRunKey, matchedTs] = split_match_key(commonKeys);

origParetoMatched = isParetoOrig(idxOrig);
refParetoMatched  = isParetoRef(idxRef);
zOrigMatched      = double(T_orig.z_eval(idxOrig));
promotedMask      = ~origParetoMatched & refParetoMatched;
promotedZlt1Mask  = promotedMask & isfinite(zOrigMatched) & (zOrigMatched < 1 - 1e-12);
promotedIdxOrig   = idxOrig(promotedZlt1Mask);
promotedIdxRef    = idxRef(promotedZlt1Mask);
promotedRunKey    = matchedRunKey(promotedZlt1Mask);
promotedTs        = matchedTs(promotedZlt1Mask);

% Per-case Pareto subsets of the original pool (for the base plot).
Ecase = cell(nCases, 1);
Tpcase = cell(nCases, 1);
for k = 1:nCases
    Tk = T_orig(string(T_orig.run_key) == string(caseNames{k}), :);
    Ecase{k} = Tk;
    Tpcase{k} = Tk(compute_pareto_mask(double(Tk.SSE), double(Tk.SSdU)), :);
end

NATURE_COLOR = nature_methods_colors();
colRefAll   = [0.60, 0.82, 0.98];          % re-evaluated reference cloud
colPromoted = plotColors(2, :);            % "promoted" points
benchmarkColor = [242, 133, 34] / 255;     % good_colors.m C.orange
[benchJtrack, benchJTV, benchFound] = load_noisy_benchmark_point(projectRoot);

% Axis limits from the pooled original + refined data.
xLimR = padded_log_limits([double(T_orig.SSdU); double(T_ref.SSdU)], 1.25);
yLimR = padded_log_limits([double(T_orig.SSE);  double(T_ref.SSE)],  1.25);
if benchFound
    xLimR = expand_log_limits_to_include_point(xLimR, benchJTV, 1.10);
    yLimR = expand_log_limits_to_include_point(yLimR, benchJtrack, 1.10);
end

fig = figure('Color', 'w', 'Toolbar', 'none', 'Name', 'Pareto Frontier Change');
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
set(fig, 'Position', [80, 80, 1400, 560]);

axL = nexttile; hold(axL, 'on');
[isParetoLeft, ~] = plot_combined_pareto_base(axL, Ecase, Tpcase, plotColors, caseMarkers, ...
    accentColor, xLimR, yLimR);
if benchFound
    scatter(axL, benchJTV, benchJtrack, 130, 's', ...
        'MarkerFaceColor', benchmarkColor, 'MarkerEdgeColor', 'none', 'LineWidth', 1.0);
end

axR = nexttile; hold(axR, 'on');
plot_combined_samples_no_guide(axR, Ecase, Tpcase, plotColors, caseMarkers);
scatter(axR, double(T_ref.SSdU), double(T_ref.SSE), 42, ...
    'filled', 'MarkerFaceColor', colRefAll, 'MarkerEdgeColor', 'none');
TrefPareto = T_ref(isParetoRef, :);
scatter(axR, double(TrefPareto.SSdU), double(TrefPareto.SSE), 80, ...
    'd', 'MarkerFaceColor', NATURE_COLOR.ReddishPurple, ...
    'MarkerEdgeColor', NATURE_COLOR.ReddishPurple, 'LineWidth', 1.0);
if ~isempty(promotedIdxRef)
    scatter(axR, double(T_ref.SSdU(promotedIdxRef)), double(T_ref.SSE(promotedIdxRef)), 140, ...
        'p', 'MarkerFaceColor', colPromoted, 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
end
if benchFound
    scatter(axR, benchJTV, benchJtrack, 130, 's', ...
        'MarkerFaceColor', benchmarkColor, 'MarkerEdgeColor', 'none', 'LineWidth', 1.0);
end

if benchFound
    [nOrigSuperior, superiorMaskOrig] = count_benchmark_strict_superior( ...
        double(T_orig.SSE(isParetoOrig)), double(T_orig.SSdU(isParetoOrig)), benchJtrack, benchJTV);
    [nRefSuperior, superiorMaskRef] = count_benchmark_strict_superior( ...
        double(TrefPareto.SSE), double(TrefPareto.SSdU), benchJtrack, benchJTV);
    fprintf('Pareto controllers strictly better than noisy benchmark (J_track and J_TV both lower):\n');
    fprintf('  Original Pareto (BO phase): %d/%d\n', nOrigSuperior, nnz(isParetoOrig));
    print_benchmark_ratio_stats('Original Pareto (BO phase)', ...
        double(T_orig.SSE(isParetoOrig)), double(T_orig.SSdU(isParetoOrig)), ...
        superiorMaskOrig, benchJtrack, benchJTV);
    fprintf('  Refined Pareto: %d/%d\n', nRefSuperior, height(TrefPareto));
    print_benchmark_ratio_stats('Refined Pareto', ...
        double(TrefPareto.SSE), double(TrefPareto.SSdU), superiorMaskRef, benchJtrack, benchJTV);
    TorigBench = build_benchmark_comparison_table( ...
        double(T_orig.SSE(isParetoOrig)), double(T_orig.SSdU(isParetoOrig)), benchJtrack, benchJTV);
    fprintf('  Original Pareto benchmark comparison table:\n');
    disp(TorigBench);
    TrefBench = build_benchmark_comparison_table( ...
        double(TrefPareto.SSE), double(TrefPareto.SSdU), benchJtrack, benchJTV);
    fprintf('  Refined Pareto benchmark comparison table:\n');
    disp(TrefBench);
end

plot_pareto_continuum_line_only(axR, double(TrefPareto.SSdU), double(TrefPareto.SSE), ...
    NATURE_COLOR.ReddishPurple, xLimR, yLimR);
apply_combined_axes_style(axL, fontSize, xLimR, yLimR);
apply_combined_axes_style(axR, fontSize, xLimR, yLimR);

title(axL, '$\mathbf{a}$', 'Interpreter', 'latex');
title(axR, '$\mathbf{b}$', 'Interpreter', 'latex');
axL.TitleHorizontalAlignment = 'left';
axR.TitleHorizontalAlignment = 'left';

outStem = fullfile(graphicsFolder, 'refined_frontier_change');
exportgraphics(fig, outStem + ".png", 'Resolution', 300);
exportgraphics(fig, outStem + ".pdf", 'ContentType', 'vector');

outTxtDir = fullfile(projectRoot, 'results', 'txt results');
if ~isfolder(outTxtDir); mkdir(outTxtDir); end
reportPath = fullfile(outTxtDir, 'refined_promoted_frontier_z_lt_1.txt');
write_promoted_report(reportPath, promotedRunKey, promotedTs, T_orig, T_ref, ...
    promotedIdxOrig, promotedIdxRef);

fprintf('Saved: %s\n', outStem + ".png");
fprintf('Saved: %s\n', outStem + ".pdf");
fprintf('Saved: %s\n', reportPath);
fprintf('J_TV rows compared: %d\n', height(T_jtv));
fprintf('Matched points (case + timestamp): %d\n', numel(idxOrig));
fprintf('Left-panel frontier pool size (all BO points): %d\n', height(T_orig));
fprintf('Left-panel frontier points (combined Pareto): %d\n', nnz(isParetoLeft));
fprintf('Original Pareto points (all BO rows): %d\n', nnz(isParetoOrig));
fprintf('Refined Pareto points (all refined rows): %d\n', nnz(isParetoRef));
fprintf('Promoted to Pareto with z < 1: %d\n', numel(promotedIdxRef));
end


function [TrefKeep, info] = keep_refined_from_original_pareto(T_ref, T_orig, isParetoOrig)
%KEEP_REFINED_FROM_ORIGINAL_PARETO Keep only refined points linked to original Pareto.
refKeys        = string(T_ref.match_key);
origKeys       = string(T_orig.match_key);
origParetoKeys = string(T_orig.match_key(isParetoOrig));

inOrig       = ismember(refKeys, origKeys);
inOrigPareto = ismember(refKeys, origParetoKeys);

info = struct();
info.n_ref_input            = height(T_ref);
info.n_drop_not_in_orig     = nnz(~inOrig);      % includes DOE-linked rows (not in evals)
info.n_drop_in_orig_not_pareto = nnz(inOrig & ~inOrigPareto);
TrefKeep = T_ref(inOrigPareto, :);
if ~isempty(TrefKeep)
    [~, uniqIdx] = unique(string(TrefKeep.match_key), 'stable');
    info.n_drop_duplicate_refined = height(TrefKeep) - numel(uniqIdx);
    TrefKeep = TrefKeep(uniqIdx, :);
else
    info.n_drop_duplicate_refined = 0;
end
info.n_keep_final = height(TrefKeep);
end


function print_refined_filter_summary(info)
%PRINT_REFINED_FILTER_SUMMARY Print keep/drop diagnostics.
fprintf('\n=== Refined sample eligibility check ===\n');
fprintf('Input refined rows: %d\n', info.n_ref_input);
fprintf('Dropped (not a BO evaluation: DOE-linked or unknown): %d\n', info.n_drop_not_in_orig);
fprintf('Dropped (BO evaluation but not Pareto): %d\n', info.n_drop_in_orig_not_pareto);
fprintf('Dropped duplicate refined rows: %d\n', info.n_drop_duplicate_refined);
fprintf('Kept refined rows (BO Pareto-linked): %d\n', info.n_keep_final);
end


function T = compute_change_table(T_orig, T_ref, col, label)
%COMPUTE_CHANGE_TABLE Build a delta table (original -> refined) for one cost column.
[commonKeys, idxOrig, idxRef] = intersect(string(T_orig.match_key), string(T_ref.match_key), 'stable');
if isempty(commonKeys)
    error('No matching (case, timestamp) points between original and refined tables.');
end
[runKey, ts] = split_match_key(commonKeys);
T = table(runKey, ts, double(T_orig.(col)(idxOrig)), double(T_ref.(col)(idxRef)), ...
    'VariableNames', ["run_key", "timestamp", label + "_original", label + "_refined"]);
T.("delta_" + label) = T.(label + "_refined") - T.(label + "_original");
T.delta_pct = 100 * T.("delta_" + label) ./ max(abs(T.(label + "_original")), eps);
T.("abs_delta_" + label) = abs(T.("delta_" + label));
T.abs_delta_pct = abs(T.delta_pct);
T = sortrows(T, 'abs_delta_pct', 'descend');
end


function print_top1_cfg(projectRoot, T_jtv)
%PRINT_TOP1_CFG Print original/refined configs for the largest JTV change.
if isempty(T_jtv)
    return
end
runKey = string(T_jtv.run_key(1));
ts     = string(T_jtv.timestamp(1));
origMat = find_mat_for_timestamp(projectRoot, runKey, ts, true);
refMat  = find_mat_for_timestamp(projectRoot, runKey, ts, false);
if strlength(origMat) == 0 || strlength(refMat) == 0
    fprintf('\nTop-1 cfg print skipped (missing MAT files) for %s | %s.\n', runKey, ts);
    return
end
Sorig = load_cfg_snapshot(origMat);
Sref  = load_cfg_snapshot(refMat);
fprintf('\n=== Top-1 point by |delta J_TV|: %s | %s ===\n', runKey, ts);
fprintf('\n--- Original cfg (%s) ---\n', origMat);
print_cfg_struct(Sorig);
fprintf('\n--- Refined cfg (%s) ---\n', refMat);
print_cfg_struct(Sref);
end


function matPath = find_mat_for_timestamp(projectRoot, runKey, ts, isOriginal)
%FIND_MAT_FOR_TIMESTAMP Resolve original/refined MAT path for one controller.
if isOriginal
    cand = fullfile(projectRoot, 'results', runKey, "out_" + ts + ".mat");
else
    cand = fullfile(projectRoot, 'results', 'final_fidelity_same_noise', ...
        runKey + "_full_f1_same_noise", "out_full_" + ts + ".mat");
end
if isfile(cand)
    matPath = string(cand);
else
    matPath = "";
end
end


function S = load_cfg_snapshot(matPath)
%LOAD_CFG_SNAPSHOT Load only cfg-relevant MAT fields.
S = struct();
vars = string(who('-file', matPath));
if any(vars == "out")
    tmp = load(matPath, 'out');
    S.out = tmp.out;
end
if any(vars == "cfg_run")
    tmp = load(matPath, 'cfg_run');
    S.cfg_run = tmp.cfg_run;
end
end


function print_cfg_struct(S)
%PRINT_CFG_STRUCT Print compact configuration information.
if isfield(S, 'out') && isfield(S.out, 'cfg')
    cfg = S.out.cfg;
    if isfield(cfg, 'f'), fprintf('f = %.12g\n', cfg.f); end
    if isfield(cfg, 'm'), fprintf('m = %d\n', cfg.m); end
    if isfield(cfg, 'p'), fprintf('p = %d\n', cfg.p); end
    if isfield(cfg, 'Q'),   fprintf('Q =\n');   disp(cfg.Q);   end
    if isfield(cfg, 'Ru'),  fprintf('Ru =\n');  disp(cfg.Ru);  end
    if isfield(cfg, 'Rdu'), fprintf('Rdu =\n'); disp(cfg.Rdu); end
end
if isfield(S, 'out') && isfield(S.out, 'theta')
    fprintf('theta =\n');
    disp(S.out.theta);
end
if isfield(S, 'cfg_run')
    cfgRun = S.cfg_run;
    if isfield(cfgRun, 'sigma_y'), fprintf('sigma_y = '); disp(cfgRun.sigma_y); end
    if isfield(cfgRun, 'mode'), fprintf('cfg_run.mode = %s\n', string(cfgRun.mode)); end
    if isfield(cfgRun, 'source_root'), fprintf('cfg_run.source_root = %s\n', string(cfgRun.source_root)); end
    if isfield(cfgRun, 'output_root'), fprintf('cfg_run.output_root = %s\n', string(cfgRun.output_root)); end
    if isfield(cfgRun, 'results_csv'), fprintf('cfg_run.results_csv = %s\n', string(cfgRun.results_csv)); end
    if isfield(cfgRun, 'out_dir'), fprintf('cfg_run.out_dir = %s\n', string(cfgRun.out_dir)); end
    if isfield(cfgRun, 'NumWorkers'), fprintf('cfg_run.NumWorkers = %d\n', double(cfgRun.NumWorkers)); end
end
end


function [runKey, ts] = split_match_key(matchKey)
%SPLIT_MATCH_KEY Split "case|timestamp" identifiers.
n = numel(matchKey);
runKey = strings(n, 1);
ts = strings(n, 1);
for i = 1:n
    parts = split(string(matchKey(i)), '|');
    runKey(i) = parts(1);
    ts(i) = parts(2);
end
end


function write_promoted_report(reportPath, promotedRunKey, promotedTs, T_orig, T_ref, ...
    promotedIdxOrig, promotedIdxRef)
%WRITE_PROMOTED_REPORT Save promoted-point table for traceability.
fid = fopen(reportPath, 'w');
if fid < 0
    warning('Could not write promoted-point report: %s', reportPath);
    return
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'Promoted to refined Pareto frontier with original z < 1\n');
fprintf(fid, 'Original pool is the BO phase only (DOE excluded structurally).\n');
fprintf(fid, 'Definition: original non-Pareto -> refined Pareto AND original z < 1\n');
fprintf(fid, 'Count: %d\n\n', numel(promotedTs));
fprintf(fid, 'case,timestamp,orig_SSE,orig_SSdU,orig_z,refined_SSE,refined_SSdU\n');
for i = 1:numel(promotedTs)
    io = promotedIdxOrig(i);
    ir = promotedIdxRef(i);
    fprintf(fid, '%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g\n', ...
        promotedRunKey(i), promotedTs(i), ...
        double(T_orig.SSE(io)), double(T_orig.SSdU(io)), double(T_orig.z_eval(io)), ...
        double(T_ref.SSE(ir)), double(T_ref.SSdU(ir)));
end
end


function limOut = expand_log_limits_to_include_point(limIn, v, padFactor)
%EXPAND_LOG_LIMITS_TO_INCLUDE_POINT Expand positive log-axis limits to include point v.
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

%% ===================== GUARDED: BENCHMARK =====================

function [Jtrack, JTV, found] = load_noisy_benchmark_point(projectRoot)
%LOAD_NOISY_BENCHMARK_POINT Load noisy benchmark aggregate objectives (guarded).
benchPath = fullfile(projectRoot, 'results', 'benchmark_reference_controller', ...
    'benchmark_full_f1_same_noise_fix', 'out_benchmark.mat');
Jtrack = nan;
JTV = nan;
found = false;
if ~isfile(benchPath)
    warning('Noisy benchmark file not found; benchmark overlays skipped: %s', benchPath);
    return
end
S = load(benchPath, 'out');
if ~isfield(S, 'out') || ~isfield(S.out, 'SSE') || ~isfield(S.out, 'SSdU')
    warning('Noisy benchmark file missing required out.SSE/out.SSdU fields: %s', benchPath);
    return
end
Jtrack = double(S.out.SSE);
JTV = double(S.out.SSdU);
found = isfinite(Jtrack) && isfinite(JTV) && Jtrack > 0 && JTV > 0;
if ~found
    warning('Noisy benchmark values are invalid; benchmark overlays skipped: %s', benchPath);
end
end


function [nSuperior, superiorMask] = count_benchmark_strict_superior(Jtrack, JTV, benchJtrack, benchJTV)
%COUNT_BENCHMARK_STRICT_SUPERIOR Strictly better than benchmark in both objectives.
superiorMask = isfinite(Jtrack) & isfinite(JTV) & ...
    (Jtrack < benchJtrack) & (JTV < benchJTV);
nSuperior = nnz(superiorMask);
end


function print_benchmark_ratio_stats(label, Jtrack, JTV, superiorMask, benchJtrack, benchJTV)
%PRINT_BENCHMARK_RATIO_STATS Print benchmark-relative ratios for superior controllers.
if nargin < 4 || isempty(superiorMask) || ~any(superiorMask)
    fprintf('    %s ratios: n/a (no strictly superior controllers)\n', label);
    return
end
JtrackSup = double(Jtrack(superiorMask));
JTVSup = double(JTV(superiorMask));
trackPct = 100 * mean(JtrackSup / benchJtrack, 'omitnan');
tvPct = 100 * mean(JTVSup / benchJTV, 'omitnan');
trackTimes = mean(benchJtrack ./ JtrackSup, 'omitnan');
tvTimes = mean(benchJTV ./ JTVSup, 'omitnan');
fprintf('    %s mean (J_better/J_bench): J_track=%.1f%%, J_TV=%.1f%%\n', label, trackPct, tvPct);
fprintf('    %s mean benchmark higher factor: J_track=%.1fx, J_TV=%.1fx\n', label, trackTimes, tvTimes);
end


function Tcmp = build_benchmark_comparison_table(Jtrack, JTV, benchJtrack, benchJTV)
%BUILD_BENCHMARK_COMPARISON_TABLE Per-controller benchmark-relative metrics.
Jtrack = double(Jtrack(:));
JTV = double(JTV(:));
n = numel(Jtrack);

ratioTrackPct = nan(n, 1);
ratioTVPct = nan(n, 1);
timesTrack = nan(n, 1);
timesTV = nan(n, 1);
isSuperior = false(n, 1);

valid = isfinite(Jtrack) & isfinite(JTV) & (Jtrack > 0) & (JTV > 0);
ratioTrackPct(valid) = 100 * (Jtrack(valid) / benchJtrack);
ratioTVPct(valid) = 100 * (JTV(valid) / benchJTV);
timesTrack(valid) = benchJtrack ./ Jtrack(valid);
timesTV(valid) = benchJTV ./ JTV(valid);
isSuperior(valid) = (Jtrack(valid) < benchJtrack) & (JTV(valid) < benchJTV);

Tcmp = table(ratioTrackPct, ratioTVPct, timesTrack, timesTV, isSuperior, ...
    'VariableNames', {'J_track_ratio_pct', 'J_TV_ratio_pct', ...
    'J_track_bench_higher_x', 'J_TV_bench_higher_x', 'strictly_better'});
Tcmp = Tcmp(Tcmp.strictly_better, :);
Tcmp = Tcmp(:, {'J_track_ratio_pct', 'J_TV_ratio_pct', 'J_track_bench_higher_x', 'J_TV_bench_higher_x'});
if ~isempty(Tcmp)
    Tcmp = sortrows(Tcmp, 'J_TV_ratio_pct', 'ascend');
end
end

%% ===================== GUARDED: FINAL f=1 METRICS =====================

function report_final_frontier_f1_metrics(E, caseNames, projectRoot, numericalFolder)
%REPORT_FINAL_FRONTIER_F1_METRICS Final Pareto f=1 table with noisy/noiseless metrics.
% Guarded: skips with a warning when the z=1 re-evaluation MAT folders are absent.
rootFolder = fullfile(projectRoot, 'results');
refinedRoot = fullfile(rootFolder, 'final_fidelity_same_noise');
if ~isfolder(refinedRoot)
    warning('main_pareto:noF1Metrics', ...
        'Final f=1 metrics section skipped: %s not found.', refinedRoot);
    return
end

nCases = numel(E);
Tall = table();
for k = 1:nCases
    Tk = E{k};
    Tk.run_key = repmat(string(caseNames{k}), height(Tk), 1);
    Tk.ts = ts_string(Tk.timestamp);
    Tall = [Tall; Tk]; %#ok<AGROW>
end
if isempty(Tall)
    fprintf('Final Pareto f=1 metrics table: no BO rows available.\n');
    return
end

isFront = compute_pareto_mask(double(Tall.SSE), double(Tall.SSdU));
Tf = Tall(isFront, :);
if isempty(Tf)
    fprintf('Final Pareto f=1 metrics table: no Pareto rows available.\n');
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

nMissingNoisy = 0;
for i = 1:n
    run_key = string(Tf.run_key(i));
    ts = string(Tf.ts(i));
    noisyOutPath = fullfile(refinedRoot, run_key + "_full_f1_same_noise", "out_full_" + ts + ".mat");
    noiselessOutPath = fullfile(rootFolder, 'test_run', run_key + "_full_f1_no_noise", "out_full_" + ts + ".mat");
    if ~isfile(noisyOutPath)
        nMissingNoisy = nMissingNoisy + 1;
        warning('Missing z=1 refined MAT for Pareto controller: %s | %s (%s)', run_key, ts, noisyOutPath);
    end
    if ~isfile(noiselessOutPath)
        warning('Missing no-noise MAT for Pareto controller: %s | %s (%s)', run_key, ts, noiselessOutPath);
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

if nMissingNoisy == n
    warning('main_pareto:noF1Data', ...
        'Final f=1 metrics section skipped: no z=1 re-evaluation MATs found under %s.', refinedRoot);
    return
end

Treport = table( ...
    string(Tf.run_key), string(Tf.ts), double(Tf.Np), double(Tf.Nc), ...
    double(Tf.Q1), double(Tf.Q2), double(Tf.Q3), ...
    double(Tf.Ru1), double(Tf.Ru2), double(Tf.Ru3), ...
    double(Tf.Rdu1), double(Tf.Rdu2), double(Tf.Rdu3), ...
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
    'run_key','timestamp','Np','Nc', ...
    'Q1','Q2','Q3','Ru1','Ru2','Ru3','Rdu1','Rdu2','Rdu3', ...
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
fprintf('\nFinal Pareto frontier f=1 controller table (sorted by decreasing noisy J_TV):\n');
disp(Treport);

% Controllers simultaneously in top-3 (lowest) for noisy J_track and all
% individual noisy settling/IAE metrics.
noisyJtrackSel = double(Treport.noisy_Jtrack);
noisyJtrackSel(~isfinite(noisyJtrackSel)) = inf;

settleColNames = { ...
    'noisy_settle_x1_h_c1','noisy_settle_x2_h_c1','noisy_settle_x3_h_c1', ...
    'noisy_settle_x1_h_c2','noisy_settle_x2_h_c2','noisy_settle_x3_h_c2'};
iaeColNames = { ...
    'noisy_IAE_x1_c1','noisy_IAE_x2_c1','noisy_IAE_x3_c1', ...
    'noisy_IAE_x1_c2','noisy_IAE_x2_c2','noisy_IAE_x3_c2'};

nTop = min(3, height(Treport));
[~, ordJtrackSel] = sort(noisyJtrackSel, 'ascend');
topJtrackMask = false(height(Treport), 1);
topJtrackMask(ordJtrackSel(1:nTop)) = true;
topSettleAllMask = true(height(Treport), 1);
topIAEAllMask = true(height(Treport), 1);

for c = 1:numel(settleColNames)
    v = double(Treport.(settleColNames{c}));
    v(~isfinite(v)) = inf;
    [~, ord] = sort(v, 'ascend');
    m = false(height(Treport), 1);
    m(ord(1:nTop)) = true;
    topSettleAllMask = topSettleAllMask & m;
end
for c = 1:numel(iaeColNames)
    v = double(Treport.(iaeColNames{c}));
    v(~isfinite(v)) = inf;
    [~, ord] = sort(v, 'ascend');
    m = false(height(Treport), 1);
    m(ord(1:nTop)) = true;
    topIAEAllMask = topIAEAllMask & m;
end
topAllMask = topJtrackMask & topSettleAllMask & topIAEAllMask;

showCols = [{'run_key','timestamp','noisy_Jtrack'}, settleColNames, iaeColNames];
TtopAll = Treport(topAllMask, showCols);
if ~isempty(TtopAll)
    TtopAll = sortrows(TtopAll, 'noisy_Jtrack', 'ascend');
end
fprintf('Controllers simultaneously top-3 lowest in noisy J_track and all individual noisy settling/IAE metrics (NaN treated as Inf):\n');
disp(TtopAll);

[benchJtrack, benchJTV, benchFound] = load_noisy_benchmark_point(projectRoot);
Tf1Bench = table();
if benchFound
    [nSuperiorF1, isSuperiorF1] = count_benchmark_strict_superior( ...
        Treport.noisy_Jtrack, Treport.noisy_JTV, benchJtrack, benchJTV);
    fprintf('Final Pareto f=1 controllers strictly better than noisy benchmark (J_track and J_TV both lower): %d/%d\n', ...
        nSuperiorF1, height(Treport));
    print_benchmark_ratio_stats('Final Pareto f=1', ...
        Treport.noisy_Jtrack, Treport.noisy_JTV, isSuperiorF1, benchJtrack, benchJTV);
    if nSuperiorF1 > 0
        fprintf('Timestamps strictly better than benchmark:\n');
        disp(Treport.timestamp(isSuperiorF1));
    end
    Tf1Bench = build_benchmark_comparison_table( ...
        double(Treport.noisy_Jtrack), double(Treport.noisy_JTV), benchJtrack, benchJTV);
    fprintf('Final Pareto f=1 benchmark comparison table:\n');
    disp(Tf1Bench);
else
    warning('Benchmark superiority verification skipped: noisy benchmark not available.');
end

if ~isfolder(numericalFolder)
    mkdir(numericalFolder);
end
writetable(Treport, fullfile(numericalFolder, 'final_pareto_frontier_f1_noisy_noiseless_metrics.csv'));

outTxtDir = fullfile(projectRoot, 'results', 'txt results');
if ~isfolder(outTxtDir)
    mkdir(outTxtDir);
end
if benchFound
    benchTxtPath = fullfile(outTxtDir, 'final_pareto_f1_benchmark_comparison_table.txt');
    writetable(Tf1Bench, benchTxtPath, 'FileType', 'text', 'Delimiter', '\t');
    fprintf('Saved: %s\n', benchTxtPath);
end
end


function M = compute_out_metrics_by_path(matPath, settlingTol)
%COMPUTE_OUT_METRICS_BY_PATH Load one out_full MAT and compute objective/settling/IAE metrics.
M = struct();
M.Jtrack = nan;
M.JTV = nan;
M.settle_h = nan(2, 3);
M.settle_case_h = nan(2, 1);
M.IAE = nan(2, 3);
M.IAE_case = nan(2, 1);
if ~isfile(matPath)
    return
end
S = load(matPath, 'out');
if ~isfield(S, 'out')
    return
end
out = S.out;
if isfield(out, 'SSE'), M.Jtrack = double(out.SSE); end
if isfield(out, 'SSdU'), M.JTV = double(out.SSdU); end
if ~isfield(out, 'case') || isempty(out.case)
    return
end
nCase = min(numel(out.case), 2);
for c = 1:nCase
    caseData = out.case(c);
    [settle_h, iae] = summarize_case_metrics_simple(caseData, settlingTol);
    nState = min(numel(settle_h), 3);
    M.settle_h(c, 1:nState) = settle_h(1:nState);
    M.IAE(c, 1:nState) = iae(1:nState);
    if any(isfinite(settle_h))
        M.settle_case_h(c) = max(settle_h(isfinite(settle_h)));
    else
        M.settle_case_h(c) = nan;
    end
    M.IAE_case(c) = sum(iae, 'omitnan');
end
end


function [settlingTimes_h, IAEByState] = summarize_case_metrics_simple(caseStruct, settlingTol)
%SUMMARIZE_CASE_METRICS_SIMPLE Compute settling times and IAE from one case struct.
if ~isfield(caseStruct, 'Y') || ~isfield(caseStruct, 'Ysp')
    settlingTimes_h = nan(1, 3);
    IAEByState = nan(1, 3);
    return
end
Y = double(caseStruct.Y);
Ysp = double(caseStruct.Ysp);
nState = min(size(Y, 2), size(Ysp, 2));
Y = Y(:, 1:nState);
Ysp = Ysp(:, 1:nState);
if isfield(caseStruct, 'dt')
    dt = double(caseStruct.dt);
elseif isfield(caseStruct, 'tf')
    dt = double(caseStruct.tf) / max(size(Y, 1) - 1, 1);
else
    dt = 1/60;
end
t = (0:size(Y, 1) - 1).' * dt;
settlingTimes_h = nan(1, nState);
IAEByState = sum(abs(Y - Ysp), 1) * dt;
epsRef = 1e-9;
for i = 1:nState
    refVal = Ysp(end, i);
    relErr = abs(Y(:, i) - Ysp(:, i)) / max(abs(refVal), epsRef);
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
