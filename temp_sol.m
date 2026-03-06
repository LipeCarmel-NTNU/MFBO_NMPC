function main_pareto()
% MAIN_PARETO  Post-processing hub for MFBO-NMPC optimization outcomes.
%
% Structure
%   1. Setup      – paths, global style, color palette
%   2. Data       – load + enrich each run table
%   3. Compute    – Pareto masks, runtime splits, fidelity stats
%   4. Report     – console summaries + numerical output files
%   5. Plot       – all figures, saved as PNG + PDF
%
% Data contract
%   results/run1/results.csv  and  results/run2/results.csv
%   Required columns: timestamp, SSE, SSdU, J, runtime_s, theta_1..theta_12
%
% Pareto policy
%   Only iterations k > DOE_END (20) are eligible for Pareto analysis.
%   A guard raises an error if any DOE point leaks into Pareto outputs.

close all; clc

% ── 1. SETUP ─────────────────────────────────────────────────────────────────

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
addpath(genpath(fullfile(projectRoot, 'dependencies')));

paths = build_paths(projectRoot);
ensure_dirs({paths.graphics, paths.numerical, paths.root});

apply_global_latex_style();
FONT_SIZE = 20;
COLORS    = nature_methods_colors(3);   % Blue | Vermillion | Orange
DOE_END   = 20;                         % last DOE iteration (inclusive)

datasets = [ ...
    struct('name', 'Case 1', 'csv', fullfile(paths.root, 'run1', 'results.csv'), ...
           'outDir', fullfile(paths.root, 'run1')); ...
    struct('name', 'Case 2', 'csv', fullfile(paths.root, 'run2', 'results.csv'), ...
           'outDir', fullfile(paths.root, 'run2')) ];

% ── 2. DATA ──────────────────────────────────────────────────────────────────

nRuns = numel(datasets);
T_all    = cell(nRuns, 1);   % full tables (DOE + optimization)
T_opt    = cell(nRuns, 1);   % optimization-phase only  (k > DOE_END)
T_pareto = cell(nRuns, 1);   % non-dominated subset of T_opt
isPareto = cell(nRuns, 1);   % logical index into T_opt

for k = 1:nRuns
    [T_all{k}, T_opt{k}, T_pareto{k}, isPareto{k}] = ...
        load_run(datasets(k).csv, DOE_END);
    T_all{k} = enrich_with_sim_data(T_all{k}, datasets(k).outDir);
    validate_no_doe_in_pareto(T_pareto{k}, DOE_END, datasets(k).name);
end

% ── 3. COMPUTE ───────────────────────────────────────────────────────────────

T_combined = vertcat(T_all{:});

runtime_tables  = cellfun(@(T, d) runtime_phase_summary(T, d.name, DOE_END), ...
                           T_all, num2cell(datasets), 'UniformOutput', false);
runtime_combined = runtime_phase_summary(T_combined, 'Case 1 + Case 2', DOE_END);

% ── 4. REPORT ────────────────────────────────────────────────────────────────

for k = 1:nRuns
    disp_pareto_table(T_pareto{k});
    disp(runtime_tables{k});
end
disp(runtime_combined);

write_numerical_summary(T_all, T_pareto, datasets, ...
    runtime_tables, runtime_combined, paths.numerical, DOE_END);

% ── 5. PLOT ──────────────────────────────────────────────────────────────────

plot_pareto_scatter_panels(T_opt, isPareto, datasets, ...
    paths.graphics, FONT_SIZE, COLORS);

plot_runtime_panels(T_all, datasets, paths.graphics, FONT_SIZE, COLORS);

plot_combined_pareto(T_opt{1}, T_pareto{1}, T_opt{2}, T_pareto{2}, ...
    fullfile(paths.graphics, 'pareto_samples_run1_run2.png'), FONT_SIZE, COLORS);

plot_cumulative_runtime(T_all, datasets, paths.graphics, FONT_SIZE, COLORS);

plot_refined_frontier(projectRoot, paths.graphics, DOE_END);
end


% ═══════════════════════════════════════════════════════════════════════════════
%  DATA
% ═══════════════════════════════════════════════════════════════════════════════

function [T, T_opt, T_pareto, isPareto] = load_run(csvPath, doeEnd)
% Load one results CSV, decode theta columns, and split DOE from optimization.

assert(isfile(csvPath), 'CSV not found: %s', csvPath);
T = readtable(csvPath, 'TextType', 'string');

required = {'timestamp','SSE','SSdU','J','runtime_s'};
for c = required
    assert(ismember(c{1}, T.Properties.VariableNames), ...
        'Missing column ''%s'' in %s', c{1}, csvPath);
end
for i = 1:12
    assert(ismember(sprintf('theta_%d',i), T.Properties.VariableNames), ...
        'Missing column ''theta_%d'' in %s', i, csvPath);
end

% ── index & timing ──
T.iteration  = (1:height(T)).';
T.runtime_min = T.runtime_s / 60;

% ── decode theta parameters ──
T.f       = double(T.theta_1);
T.theta_p = round(double(T.theta_2));
T.theta_m = round(double(T.theta_3));

T.q1_log10   = double(T.theta_4);   T.q2_log10   = double(T.theta_5);   T.q3_log10   = double(T.theta_6);
T.r_u1_log10 = double(T.theta_7);   T.r_u2_log10 = double(T.theta_8);   T.r_u3_log10 = double(T.theta_9);
T.r_du1_log10= double(T.theta_10);  T.r_du2_log10= double(T.theta_11);  T.r_du3_log10= double(T.theta_12);

T.Q_x1   = 10.^T.q1_log10;    T.Q_x2   = 10.^T.q2_log10;    T.Q_x3   = 10.^T.q3_log10;
T.R_u_x1 = 10.^T.r_u1_log10;  T.R_u_x2 = 10.^T.r_u2_log10;  T.R_u_x3 = 10.^T.r_u3_log10;
T.R_du_x1= 10.^T.r_du1_log10; T.R_du_x2= 10.^T.r_du2_log10; T.R_du_x3= 10.^T.r_du3_log10;

T.m = T.theta_m + 1;
T.p = T.theta_p + T.m;

% ── split at DOE boundary ──
T_opt  = T(T.iteration > doeEnd, :);
isPareto = pareto_mask(T_opt.SSE, T_opt.SSdU);
T_pareto = T_opt(isPareto, :);
T_pareto = sortrows(T_pareto, 'SSE', 'descend');

T_pareto.runtime_per_f = T_pareto.runtime_s ./ max(T_pareto.f, eps);
end


function T = enrich_with_sim_data(T, runDir)
% Append per-iteration simulation diagnostics from out_<timestamp>.mat files.
n = height(T);
T.out_found              = false(n,1);
T.out_runtime_s          = NaN(n,1);
T.out_mean_case_runtime_s= NaN(n,1);
T.out_mean_step_runtime_ms = NaN(n,1);
T.out_total_steps        = NaN(n,1);
T.out_case_count         = NaN(n,1);

for i = 1:n
    matFile = fullfile(runDir, sprintf('out_%s.mat', T.timestamp(i)));
    if ~isfile(matFile), continue; end

    S = load(matFile, 'out');
    if ~isfield(S, 'out'), continue; end
    out = S.out;

    T.out_found(i) = true;
    if isfield(out, 'runtime_s'),  T.out_runtime_s(i) = out.runtime_s; end
    if ~isfield(out, 'case') || isempty(out.case), continue; end

    nC = numel(out.case);
    caseRT   = NaN(nC, 1);
    stepMean = NaN(nC, 1);
    nSteps   = NaN(nC, 1);
    for c = 1:nC
        ci = out.case(c);
        if isfield(ci, 'runtime_s'),  caseRT(c)   = ci.runtime_s;      end
        if isfield(ci, 'RUNTIME') && ~isempty(ci.RUNTIME)
            r = double(ci.RUNTIME(:));
            nSteps(c)   = numel(r);
            stepMean(c) = mean(r, 'omitnan');
        end
    end
    T.out_case_count(i)           = nC;
    T.out_mean_case_runtime_s(i)  = mean(caseRT, 'omitnan');
    T.out_total_steps(i)          = sum(nSteps,  'omitnan');
    T.out_mean_step_runtime_ms(i) = 1000 * mean(stepMean, 'omitnan');
end
end


% ═══════════════════════════════════════════════════════════════════════════════
%  COMPUTATIONS
% ═══════════════════════════════════════════════════════════════════════════════

function mask = pareto_mask(SSE, SSdU)
% Return logical vector: true where point (SSE_i, SSdU_i) is non-dominated.
SSE  = double(SSE(:));
SSdU = double(SSdU(:));
n    = numel(SSE);
mask = true(n, 1);
for i = 1:n
    dominated = SSE <= SSE(i) & SSdU <= SSdU(i) & ...
               (SSE < SSE(i)  | SSdU < SSdU(i));
    dominated(i) = false;
    if any(dominated), mask(i) = false; end
end
end


function tbl = runtime_phase_summary(T, label, doeEnd)
% Summarize wall-clock time spent in DOE vs. optimization phases.
doeMask = T.iteration <= doeEnd;
optMask = ~doeMask;
rt = double(T.runtime_min);

phase            = {'DOE'; 'Optimisation'; 'Total'};
iterations       = [nnz(doeMask); nnz(optMask); height(T)];
runtime_min      = [sum(rt(doeMask),'omitnan'); sum(rt(optMask),'omitnan'); sum(rt,'omitnan')];
runtime_pct      = [NaN; NaN; 100];
if runtime_min(3) > 0
    runtime_pct(1:2) = 100 * runtime_min(1:2) / runtime_min(3);
end

tbl = table(phase, iterations, runtime_min, runtime_pct, ...
    'VariableNames', {'phase','iterations','runtime_min','runtime_pct_total'});
fprintf('\nRuntime summary – %s:\n', label);
disp(tbl);
end


function z = mean_z_opt(T, doeEnd)
% Average fidelity z over optimization-phase iterations only.
mask = double(T.iteration) > doeEnd;
z    = mean(double(T.f(mask)), 'omitnan');
end


function [pct, count, nOpt] = nc1_share_opt(T, doeEnd)
% Fraction of optimization iterations where N_c (= m) equals 1.
mask = double(T.iteration) > doeEnd;
nOpt  = nnz(mask);
count = nnz(mask & (double(T.m) == 1));
pct   = 100 * count / max(nOpt, 1);
if nOpt == 0, pct = NaN; end
end


function T_delta = jtv_change_table(T_orig, T_ref)
% Compute J_TV delta between matched original and refined evaluation points.
[keys, iO, iR] = intersect(string(T_orig.match_key), string(T_ref.match_key), 'stable');
assert(~isempty(keys), 'No (run_key, timestamp) matches between original and refined.');

[runKey, ts] = split_match_key(keys);
T_delta = table(runKey, ts, ...
    double(T_orig.SSdU(iO)), double(T_ref.SSdU(iR)), ...
    'VariableNames', {'run_key','timestamp','JTV_original','JTV_refined'});
T_delta.delta_JTV     = T_delta.JTV_refined - T_delta.JTV_original;
T_delta.delta_pct     = 100 * T_delta.delta_JTV ./ T_delta.JTV_original;
T_delta.abs_delta_JTV = abs(T_delta.delta_JTV);
T_delta = sortrows(T_delta, 'abs_delta_JTV', 'descend');
end


% ═══════════════════════════════════════════════════════════════════════════════
%  REPORTING
% ═══════════════════════════════════════════════════════════════════════════════

function disp_pareto_table(Tp)
% Print Pareto controller settings (optimization phase only).
cols = {'timestamp','SSE','SSdU','J','p','m','f','runtime_min','runtime_per_f', ...
        'Q_x1','Q_x2','Q_x3','R_u_x1','R_u_x2','R_u_x3','R_du_x1','R_du_x2','R_du_x3'};
cols = cols(ismember(cols, Tp.Properties.VariableNames));
disp('Pareto frontier (optimization phase, DOE excluded):');
disp(Tp(:, cols));
end


function validate_no_doe_in_pareto(Tp, doeEnd, label)
% Guard: raise an error if any DOE point leaked into the Pareto table.
if isempty(Tp) || ~ismember('iteration', Tp.Properties.VariableNames), return; end
assert(~any(double(Tp.iteration) <= doeEnd), ...
    'DOE point detected in Pareto table for %s.', label);
fprintf('DOE exclusion verified for %s: %d Pareto rows, all k > %d.\n', ...
    label, height(Tp), doeEnd);
end


function write_numerical_summary(T_all, T_pareto, datasets, ...
    runtime_tables, runtime_combined, outDir, doeEnd)
% Write runtime splits, mean fidelity, and N_c statistics to a text file.
outFile = fullfile(outDir, 'runtime_and_params.txt');
fid = fopen(outFile, 'w');
assert(fid ~= -1, 'Cannot open: %s', outFile);
cleanup = onCleanup(@() fclose(fid));

nRuns = numel(datasets);

% ── runtime splits ──
for k = 1:nRuns
    write_runtime_block(fid, runtime_tables{k}, datasets(k).name);
end
write_runtime_block(fid, runtime_combined, 'Case 1 + Case 2 (combined)');

% ── mean fidelity z ──
fprintf(fid, '\nMean fidelity z (optimization phase):\n');
for k = 1:nRuns
    z = mean_z_opt(T_all{k}, doeEnd);
    fprintf(fid, '  %s: %.6g\n', datasets(k).name, z);
end
z_comb = mean_z_opt(vertcat(T_all{:}), doeEnd);
fprintf(fid, '  Combined: %.6g\n', z_comb);

% ── N_c = 1 share ──
fprintf(fid, '\nOptimization points with N_c = 1:\n');
for k = 1:nRuns
    [pct, cnt, n] = nc1_share_opt(T_all{k}, doeEnd);
    fprintf(fid, '  %s: %d/%d (%.6g%%)\n', datasets(k).name, cnt, n, pct);
end
[pct, cnt, n] = nc1_share_opt(vertcat(T_all{:}), doeEnd);
fprintf(fid, '  Combined: %d/%d (%.6g%%)\n', cnt, n, pct);

% ── Pareto horizon counts ──
fprintf(fid, '\nPareto point composition (optimization phase, DOE excluded):\n');
for k = 1:nRuns
    Tp = T_pareto{k};
    fprintf(fid, '  %s: m=1 -> %d/%d | p=1 -> %d/%d\n', datasets(k).name, ...
        nnz(Tp.m == 1), height(Tp), nnz(Tp.p == 1), height(Tp));
end
end


function write_runtime_block(fid, tbl, label)
fprintf(fid, '\nRuntime summary – %s:\n', label);
for i = 1:height(tbl)
    fprintf(fid, '  %-14s | iters=%-4d | min=%.4g | pct=%.2f%%\n', ...
        tbl.phase{i}, tbl.iterations(i), tbl.runtime_min(i), tbl.runtime_pct_total(i));
end
end


% ═══════════════════════════════════════════════════════════════════════════════
%  PLOTS
% ═══════════════════════════════════════════════════════════════════════════════

function plot_pareto_scatter_panels(T_opt, isPareto, datasets, outDir, fontSize, clr)
% Side-by-side SSE vs SSdU scatter, coloured by fidelity z.
fig = figure('Color', 'w', 'Name', 'Pareto – SSE vs SSdU by Case');
tl  = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
cmap = load_navia_colormap(256);
labels = {'a', 'b'};
XLim   = [1e-2, 1e2];
YLim   = [1e4, 1.3e5];

for k = 1:2
    T  = T_opt{k};
    ip = isPareto{k};
    ax = nexttile(tl); hold(ax, 'on');

    if isempty(T)
        text(ax, 0.5, 0.5, 'No optimization points', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center');
    else
        draw_pareto_guide(ax, T.SSdU(ip), T.SSE(ip), clr(3,:), XLim, YLim);
        scatter(ax, double(T.SSdU), double(T.SSE), 80, double(T.f), ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
        scatter(ax, double(T.SSdU(ip)), double(T.SSE(ip)), 170, clr(3,:), ...
            'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr(3,:), 'LineWidth', 1.2);
    end

    style_log_axes(ax, XLim, YLim, fontSize, labels{k});
    colormap(ax, cmap); caxis(ax, [0, 1]);
    xlabel(ax, '$J_{\mathrm{TV}}$');
    ylabel(ax, '$J_{\mathrm{track}}$');
    cb = colorbar(ax);
    cb.Label.String = '$z$ (dimensionless)';
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = fontSize;
end

save_figure(fig, fullfile(outDir, 'sse_vs_ssdu_side_by_side_z.png'), fontSize, 1200, 460);
end


function plot_runtime_panels(T_all, datasets, outDir, fontSize, clr)
% Side-by-side per-iteration runtime and fidelity z.
fig = figure('Color', 'w', 'Name', 'Iteration Runtime and Fidelity');
tl  = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
labels = {'a', 'b'};

for k = 1:2
    T  = T_all{k};
    ax = nexttile(tl); hold(ax, 'on');

    yyaxis(ax, 'left');
    plot(ax, T.iteration, double(T.runtime_min)/60, '-', ...
        'LineWidth', 2, 'Color', clr(k,:));
    ylim(ax, [0, 4]); yticks(ax, 0:1:4);
    ylabel(ax, '$t_{\mathrm{iter}}$ (h)');
    ax.YColor = 'k';

    xline(ax, 20, '--k', 'LineWidth', 2, 'Alpha', 1);

    yyaxis(ax, 'right');
    plot(ax, T.iteration, double(T.f), 'o', ...
        'LineWidth', 2, 'MarkerSize', 4, 'Color', clr(3,:));
    ylim(ax, [0, 1]); yticks(ax, 0:0.2:1);
    ylabel(ax, '$z$ (dimensionless)');
    ax.YColor = 'k';

    xlabel(ax, '$k$ (iteration)');
    xlim(ax, [1, height(T)]);
    set(ax, 'FontSize', fontSize);
    title_label(ax, labels{k});
    grid(ax, 'off'); box(ax, 'off');
end

save_figure(fig, fullfile(outDir, 'runtime_vs_iteration_side_by_side.png'), fontSize, 1200, 460);
end


function plot_combined_pareto(T_opt1, Tp1, T_opt2, Tp2, outPath, fontSize, clr)
% Combined scatter of all optimization samples with the global Pareto frontier.
fig = figure('Color', 'w', 'Name', 'Combined Pareto Samples');
ax  = axes(fig); hold(ax, 'on');
XLim = [1e-2, 2e0];
YLim = [1e4,  3.5e4];

Tall = [T_opt1; T_opt2];
assert(~any(double(Tall.iteration) <= 20), ...
    'Combined Pareto plot received DOE rows (iteration <= 20).');

% all optimization samples (grey background)
scatter(ax, double(Tall.SSdU), double(Tall.SSE), 18, ...
    'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

% global Pareto guide curve
globalPareto = pareto_mask(Tall.SSE, Tall.SSdU);
Tf = Tall(globalPareto, :);
draw_pareto_guide(ax, Tf.SSdU, Tf.SSE, clr(3,:), XLim, YLim);

% per-case Pareto markers
scatter(ax, double(Tp1.SSdU), double(Tp1.SSE), 80, clr(1,:), ...
    'o', 'MarkerFaceColor', clr(1,:), 'MarkerEdgeColor', clr(1,:), 'LineWidth', 1.4);
scatter(ax, double(Tp2.SSdU), double(Tp2.SSE), 90, clr(2,:), ...
    '^', 'MarkerFaceColor', clr(2,:), 'MarkerEdgeColor', clr(2,:), 'LineWidth', 1.4);

% global Pareto ring overlay
scatter(ax, double(Tf.SSdU), double(Tf.SSE), 300, clr(3,:), ...
    'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', clr(3,:), 'LineWidth', 2);

style_log_axes(ax, XLim, YLim, fontSize, '');
xlabel(ax, '$J_{\mathrm{TV}}$');
ylabel(ax, '$J_{\mathrm{track}}$');
save_figure(fig, outPath, fontSize, 920, 520);
end


function plot_cumulative_runtime(T_all, datasets, outDir, fontSize, clr)
% Cumulative wall-clock time for both runs on one axis.
fig = figure('Color', 'w', 'Name', 'Cumulative Runtime');
ax  = axes(fig); hold(ax, 'on');

lineStyles = {'-', '-.'};
for k = 1:2
    T = T_all{k};
    plot(ax, double(T.iteration), cumsum(double(T.runtime_min)/60, 'omitnan'), ...
        lineStyles{k}, 'LineWidth', 2, 'Color', clr(k,:), ...
        'DisplayName', datasets(k).name);
end
xline(ax, 20, '--k', 'LineWidth', 2, 'Alpha', 1);

xlabel(ax, '$k$ (iteration)');
ylabel(ax, '$t_{\mathrm{run}}$ (h)');
xlim(ax, [1, max(height(T_all{1}), height(T_all{2}))]);
set(ax, 'FontSize', fontSize);
legend(ax, 'Location', 'northwest', 'Interpreter', 'latex');
grid(ax, 'off'); box(ax, 'off');

save_figure(fig, fullfile(outDir, 'runtime_cumulative_run1_run2.png'), fontSize, 920, 520);
end


function plot_refined_frontier(projectRoot, outDir, doeEnd)
% Compare original non-DOE Pareto points vs. refined (full-fidelity) evaluations.
origFiles = [ ...
    struct('runKey','run1','runLabel','Case 1', ...
           'path', fullfile(projectRoot,'results','run1','results.csv')); ...
    struct('runKey','run2','runLabel','Case 2', ...
           'path', fullfile(projectRoot,'results','run2','results.csv')) ];

refFiles = [ ...
    struct('runKey','run1','runLabel','Case 1 refined', ...
           'path', fullfile(projectRoot,'results','final_fidelity_same_noise','run1_full_f1_same_noise','results_full.csv')); ...
    struct('runKey','run2','runLabel','Case 2 refined', ...
           'path', fullfile(projectRoot,'results','final_fidelity_same_noise','run2_full_f1_same_noise','results_full.csv')) ];

T_orig_full = load_run_for_refined(origFiles, 0);
T_orig      = T_orig_full(T_orig_full.iteration > doeEnd, :);
T_ref       = load_run_for_refined(refFiles, 0);

assert(~isempty(T_orig), 'No original rows after DOE filtering.');
assert(~isempty(T_ref),  'No refined rows loaded.');

isParetoO = pareto_mask(T_orig.SSE, T_orig.SSdU);
[T_ref, info] = filter_refined_to_original_pareto(T_ref, T_orig, isParetoO, T_orig_full);
print_refined_filter_summary(info);
assert(~isempty(T_ref), 'No refined rows after filtering to original Pareto.');

T_delta = jtv_change_table(T_orig, T_ref);
fprintf('\n=== J_TV change (original -> refined), sorted by |delta| ===\n');
disp(T_delta);

isParetoR = pareto_mask(T_ref.SSE, T_ref.SSdU);

[commonKeys, iO, iR] = intersect(string(T_orig.match_key), string(T_ref.match_key), 'stable');
origParetoBit = isParetoO(iO);
refParetoBit  = isParetoR(iR);
zOrig = double(T_orig.z_eval(iO));

promotedMask = ~origParetoBit & refParetoBit & isfinite(zOrig) & (zOrig < 1 - 1e-12);
iOProm = iO(promotedMask);
iRProm = iR(promotedMask);
[pRunKey, pTs] = split_match_key(commonKeys(promotedMask));

% ── matched-pool frontier (left panel) ──
T1 = T_orig(T_orig.run_key == 'run1', :);
T2 = T_orig(T_orig.run_key == 'run2', :);
poolMask1 = ismember(string(T1.match_key), string(commonKeys));
poolMask2 = ismember(string(T2.match_key), string(commonKeys));
Tpool = [T1(poolMask1,:); T2(poolMask2,:)];
leftFrontierMask = pareto_mask(Tpool.SSE, Tpool.SSdU);
Tf = Tpool(leftFrontierMask, :);
doeKeys = string(T_orig_full.match_key(T_orig_full.iteration <= doeEnd));
assert(~any(ismember(string(Tf.match_key), doeKeys)), ...
    'Left-panel frontier unexpectedly contains DOE-derived points.');

clr = nature_methods_colors(3);
PURPLE = [0.60 0.82 0.98];

fig = figure('Color', 'w', 'Position', [80 80 1400 560]);
tl  = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
XLim = [1e-2, 2e0];
YLim = [1e4,  3.5e4];

% ── Left panel: matched-pool frontier ──
axL = nexttile(tl); hold(axL, 'on');
scatter(axL, double([T1.SSdU;T2.SSdU]), double([T1.SSEs;T2.SSE]), 18, ...
    'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
draw_frontier_polyline(axL, Tf.SSdU, Tf.SSE, clr(3,:), 2, 8);
scatter(axL, double(T1.SSdU(poolMask1)), double(T1.SSE(poolMask1)), 80, clr(1,:), ...
    'o', 'MarkerFaceColor', clr(1,:), 'MarkerEdgeColor', clr(1,:), 'LineWidth', 1.4);
scatter(axL, double(T2.SSdU(poolMask2)), double(T2.SSE(poolMask2)), 90, clr(2,:), ...
    '^', 'MarkerFaceColor', clr(2,:), 'MarkerEdgeColor', clr(2,:), 'LineWidth', 1.4);

% ── Right panel: original vs refined ──
axR = nexttile(tl); hold(axR, 'on');
scatter(axR, double(T1.SSdU), double(T1.SSE), 16, ...
    'o', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
scatter(axR, double(T2.SSdU), double(T2.SSE), 20, ...
    '^', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
scatter(axR, double(T_ref.SSdU), double(T_ref.SSE), 42, ...
    'filled', 'MarkerFaceColor', PURPLE, 'MarkerEdgeColor', 'none');
TrefP = T_ref(isParetoR, :);
scatter(axR, double(TrefP.SSdU), double(TrefP.SSE), 80, ...
    'd', 'MarkerFaceColor', clr(2,:), 'MarkerEdgeColor', clr(2,:), 'LineWidth', 1);
if ~isempty(iRProm)
    scatter(axR, double(T_ref.SSdU(iRProm)), double(T_ref.SSE(iRProm)), 140, ...
        'p', 'MarkerFaceColor', clr(2,:), 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
end

for ax = [axL, axR]
    style_log_axes(ax, XLim, YLim, 16, '');
    xlabel(ax, '$J_{\mathrm{TV}}$', 'Interpreter', 'latex');
    ylabel(ax, '$J_{\mathrm{track}}$', 'Interpreter', 'latex');
end
title_label(axL, 'a'); title_label(axR, 'b');

stem = fullfile(outDir, 'refined_frontier_change');
exportgraphics(fig, [stem, '.png'], 'Resolution', 300);
exportgraphics(fig, [stem, '.pdf'], 'ContentType', 'vector');

% ── promoted-point report ──
rptDir = fullfile(projectRoot, 'results', 'txt results');
if ~isfolder(rptDir), mkdir(rptDir); end
write_promoted_report(fullfile(rptDir, 'refined_promoted_z_lt_1.txt'), ...
    pRunKey, pTs, T_orig, T_ref, iOProm, iRProm);

fprintf('Refined frontier figure saved to: %s\n', stem);
fprintf('Promoted points (z < 1 -> refined Pareto): %d\n', numel(iRProm));
end


% ═══════════════════════════════════════════════════════════════════════════════
%  REFINED ANALYSIS HELPERS
% ═══════════════════════════════════════════════════════════════════════════════

function T = load_run_for_refined(fileDefs, dropFirst)
% Load CSV files for refined comparison; add run_key, match_key, z_eval columns.
T = table();
for k = 1:numel(fileDefs)
    p = fileDefs(k).path;
    if ~isfile(p), warning('Missing: %s', p); continue; end
    Tk = readtable(p, 'TextType', 'string');
    for c = {'timestamp','SSE','SSdU'}
        assert(ismember(c{1}, Tk.Properties.VariableNames), ...
            'Missing ''%s'' in %s', c{1}, p);
    end
    Tk.run_label = repmat(string(fileDefs(k).runLabel), height(Tk), 1);
    Tk.run_key   = repmat(string(fileDefs(k).runKey),   height(Tk), 1);
    Tk.iteration = (1:height(Tk)).';
    if dropFirst > 0, Tk = Tk(Tk.iteration > dropFirst, :); end
    Tk.z_eval    = double(getfield_or_nan(Tk, 'theta_1'));
    Tk.match_key = Tk.run_key + "|" + string(Tk.timestamp);
    Tk = Tk(:, {'run_label','run_key','match_key','iteration','timestamp','SSE','SSdU','z_eval'});
    T  = [T; Tk]; %#ok<AGROW>
end
end


function v = getfield_or_nan(T, col)
if ismember(col, T.Properties.VariableNames)
    v = double(T.(col));
else
    v = NaN(height(T), 1);
end
end


function [TrefKeep, info] = filter_refined_to_original_pareto(T_ref, T_orig, isParetoO, T_orig_full)
% Retain only refined rows whose match_key maps to an original Pareto point.
refKeys       = string(T_ref.match_key);
origKeys      = string(T_orig.match_key);
origFullKeys  = string(T_orig_full.match_key);
origParetoKeys= string(T_orig.match_key(isParetoO));

inFull      = ismember(refKeys, origFullKeys);
inOrig      = ismember(refKeys, origKeys);
inPareto    = ismember(refKeys, origParetoKeys);

info = struct( ...
    'n_input',           height(T_ref), ...
    'n_drop_doe',        nnz(inFull & ~inOrig), ...
    'n_drop_not_orig',   nnz(~inFull), ...
    'n_drop_not_pareto', nnz(inOrig & ~inPareto), ...
    'dropped_doe_keys',  {refKeys(inFull & ~inOrig)});

TrefKeep = T_ref(inPareto, :);
[~, ui]  = unique(string(TrefKeep.match_key), 'stable');
info.n_drop_duplicate = height(TrefKeep) - numel(ui);
TrefKeep = TrefKeep(ui, :);
info.n_keep = height(TrefKeep);
end


function print_refined_filter_summary(info)
fprintf('\n=== Refined eligibility filter ===\n');
fprintf('  Input rows:              %d\n', info.n_input);
fprintf('  Dropped (DOE-linked):    %d\n', info.n_drop_doe);
fprintf('  Dropped (not in orig):   %d\n', info.n_drop_not_orig);
fprintf('  Dropped (not Pareto):    %d\n', info.n_drop_not_pareto);
fprintf('  Dropped (duplicates):    %d\n', info.n_drop_duplicate);
fprintf('  Kept:                    %d\n', info.n_keep);
if info.n_drop_doe > 0
    fprintf('  DOE-linked keys dropped:\n');
    for i = 1:numel(info.dropped_doe_keys)
        fprintf('    %s\n', info.dropped_doe_keys(i));
    end
end
end


function write_promoted_report(path, runKey, ts, T_orig, T_ref, iO, iR)
% Write CSV-style report of points promoted to refined Pareto with z < 1.
fid = fopen(path, 'w');
if fid < 0, warning('Cannot write: %s', path); return; end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Promoted to refined Pareto (original z < 1): %d points\n', numel(ts));
fprintf(fid, 'DOE excluded from original (iterations 1..20).\n\n');
fprintf(fid, 'run_key,timestamp,orig_SSE,orig_SSdU,orig_z,refined_SSE,refined_SSdU\n');
for i = 1:numel(ts)
    fprintf(fid, '%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g\n', ...
        runKey(i), ts(i), ...
        double(T_orig.SSE(iO(i))), double(T_orig.SSdU(iO(i))), double(T_orig.z_eval(iO(i))), ...
        double(T_ref.SSE(iR(i))),  double(T_ref.SSdU(iR(i))));
end
end


% ═══════════════════════════════════════════════════════════════════════════════
%  PLOT UTILITIES
% ═══════════════════════════════════════════════════════════════════════════════

function draw_pareto_guide(ax, x, y, color, XLim, YLim)
% Smooth pchip guide curve through Pareto frontier in log-log space.
x = double(x(:)); y = double(y(:));
if numel(x) < 2, return; end

[xs, ord] = sort(x, 'ascend'); ys = y(ord);
if numel(xs) >= 3
    lx = log10(xs); ly = log10(ys);
    lxq = linspace(min(lx), max(lx), 220);
    lyq = pchip(lx, ly, lxq);
    xs = 10.^lxq'; ys = 10.^lyq';
end

% vertical extension to top of Y range
if ys(1) < YLim(2), xs = [xs(1); xs]; ys = [YLim(2); ys]; end
% horizontal extension to right of X range
if xs(end) < XLim(2), xs = [xs; XLim(2)]; ys = [ys; ys(end)]; end

plot(ax, xs, ys, '-', 'Color', color, 'LineWidth', 2);
end


function draw_frontier_polyline(ax, x, y, color, lw, ms)
% Simple sorted polyline with circle markers on Pareto points.
[xs, ord] = sort(double(x(:)), 'ascend'); ys = double(y(ord));
plot(ax, xs, ys, '-o', 'Color', color, 'LineWidth', lw, ...
    'MarkerSize', ms, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', color);
end


function style_log_axes(ax, XLim, YLim, fontSize, panelLabel)
% Apply consistent log-scale style and optional bold panel label.
set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', fontSize);
xlim(ax, XLim); ylim(ax, YLim);
grid(ax, 'off'); box(ax, 'off');
if ~isempty(panelLabel)
    title_label(ax, panelLabel);
end
end


function title_label(ax, label)
% Left-aligned bold panel label (e.g. "a", "b").
title(ax, sprintf('$\\mathbf{%s}$', label), 'Interpreter', 'latex');
ax.TitleHorizontalAlignment = 'left';
end


function save_figure(fig, pngPath, fontSize, widthPx, heightPx)
% Set figure size, font, then export PNG and PDF.
figure(fig);
set_fig_size(widthPx, heightPx);
set_font_size(fontSize);
exportgraphics(fig, pngPath, 'Resolution', 300);
[d, stem] = fileparts(pngPath);
exportgraphics(fig, fullfile(d, [stem, '.pdf']), 'ContentType', 'vector');
end


% ═══════════════════════════════════════════════════════════════════════════════
%  GENERAL UTILITIES
% ═══════════════════════════════════════════════════════════════════════════════

function paths = build_paths(projectRoot)
paths.root      = fullfile(projectRoot, 'results');
paths.graphics  = fullfile(paths.root, 'graphical_results');
paths.numerical = fullfile(paths.root, 'numerical results');
end


function ensure_dirs(dirList)
for i = 1:numel(dirList)
    if ~isfolder(dirList{i}), mkdir(dirList{i}); end
end
end


function apply_global_latex_style()
set(groot, 'defaultTextInterpreter',          'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');
end


function cmap = load_navia_colormap(n)
% Load the sequential Navia colormap; interpolate to n entries.
matPath = which('navia.mat');
assert(~isempty(matPath), 'Cannot locate navia.mat on MATLAB path.');
S = load(matPath, 'navia');
assert(isfield(S,'navia') && size(S.navia,2)==3, 'Invalid navia.mat format.');
cmap = S.navia;
if size(cmap,1) ~= n
    x = linspace(0,1,size(cmap,1));
    cmap = interp1(x, double(cmap), linspace(0,1,n), 'linear');
end
cmap = min(max(cmap, 0), 1);
end


function [runKey, ts] = split_match_key(keys)
% Split "runKey|timestamp" strings into two arrays.
n = numel(keys);
runKey = strings(n,1); ts = strings(n,1);
for i = 1:n
    p = split(string(keys(i)), '|');
    runKey(i) = p(1); ts(i) = p(2);
end
end