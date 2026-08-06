function runtime_consistency()
% RUNTIME_CONSISTENCY
% Check whether the runtime used in the analysis (t_total) is trustworthy by
% comparing it against the accumulated NMPC solve time (t_nmpc) for every
% evaluation (DOE + BO). The two should be close: times in hours, gap in
% minutes. A large gap means t_total contains something other than solving.
%
% Timing semantics (from dependencies/simulation/*.m)
% - t_nmpc  = out.runtime_s: sum over control steps of the per-step timer in
%   nmpc_run_case.m (fmincon solve + one-step ode45 plant integration + the
%   capped flag logging), summed over both setpoint scenarios.
% - t_total = out.wall_s.total: toc over the whole simulate_nmpc call. On top
%   of the timed steps it contains phi evaluation, build_nmpc, setpoint
%   updates, checkpoint writes and case finalisation.
% - gap = t_total - t_nmpc, expected small and positive.
%     gap >> 0 : wall-clock overhead that is NOT solver difficulty
%                (e.g. checkpoint/log writes on a OneDrive-synced folder).
%     gap < 0  : structurally impossible in one call -> the evaluation was
%                resumed from a checkpoint. On resume, RUNTIME (and thus
%                t_nmpc) accumulates across calls, but wall_s.total is
%                re-tic'ed, so t_total undercounts the true wall time.
%
% Interpretation for the crash/runtime question (main_pareto.m):
% - If evals with runtime jumps have a small gap: the jump is genuine solve
%   time -> slow = difficult to optimize = likely to crash is plausible.
% - If the jump lives in the gap: the clock, not the solver -> t_total is
%   the wrong runtime measure for those evaluations.
%
% Input : Result analysis/storage/preprocessed.mat (evals, doe, meta)
% Output: results/graphical_results/runtime_consistency.png/.pdf
%         results/numerical results/runtime_consistency.txt
%         console diagnostics (per-case summary, correlations, offenders)
%
% Color scheme: see COLOR_SCHEME.md (Wong categorical, accent ReddishPurple).

%% Dependencies and paths
close all; clc
here      = fileparts(mfilename('fullpath'));
repo_root = fileparts(here);
addpath(genpath(fullfile(repo_root, 'dependencies')));

graphicsDir  = fullfile(repo_root, 'results', 'graphical_results');
numericalDir = fullfile(repo_root, 'results', 'numerical results');
if ~isfolder(graphicsDir);  mkdir(graphicsDir);  end
if ~isfolder(numericalDir); mkdir(numericalDir); end

%% Load preprocessed data
storePath = fullfile(here, 'storage', 'preprocessed.mat');
if ~isfile(storePath)
    error('runtime_consistency:noData', ...
        'Missing %s. Run parse_registry.py, then initial_preprocessing.m.', storePath);
end
S     = load(storePath, 'evals', 'doe');
evals = S.evals;
doe   = S.doe;

%% Style (COLOR_SCHEME.md)
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
fontSize = 20;

plotColors  = nature_methods_colors(3);   % Blue, BluishGreen, ReddishPurple
accentColor = plotColors(3, :);
caseMarkers = ["o", "^", "d"];

%% Combine DOE + BO into one execution timeline per case
evals.case = removecats(categorical(evals.case));
doe.case   = removecats(categorical(doe.case));
caseNames  = categories(evals.case);
nCases     = numel(caseNames);
if nCases > numel(caseMarkers)
    error('runtime_consistency:tooManyCases', ...
        'Only %d case styles defined; found %d cases.', numel(caseMarkers), nCases);
end
caseLabels = arrayfun(@pretty_case, string(caseNames));

A = cell(nCases, 1);          % per-case timeline table (DOE first, then BO)
doeCount = zeros(nCases, 1);
for k = 1:nCases
    Dk = sortrows(doe(doe.case == caseNames{k},     :), 'iter');
    Ek = sortrows(evals(evals.case == caseNames{k}, :), 'iter');
    doeCount(k) = height(Dk);

    T = [Dk; Ek];
    T.phase_label = [repmat("DOE", height(Dk), 1); repmat("BO", height(Ek), 1)];
    T.gidx    = (1:height(T))';                              % execution order
    T.totalH  = double(T.t_total) / 3600;                    % hours
    T.nmpcH   = double(T.t_nmpc)  / 3600;                    % hours
    T.gapMin  = (double(T.t_total) - double(T.t_nmpc)) / 60; % minutes
    T.nBad    = double(T.n_flag0) + double(T.n_flag_neg);    % failed solves
    A{k} = T;
end

%% Figure: solve time vs wall time, gap timeline, gap vs crashes
fig = figure('Color', 'w');
tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- (a) t_nmpc vs t_total, y = x reference -----------------------------
ax1 = nexttile; hold(ax1, 'on'); grid(ax1, 'off'); box(ax1, 'off');
lim = [0, 1.05 * max(cellfun(@(T) max([T.totalH; T.nmpcH]), A))];
plot(ax1, lim, lim, '-', 'Color', accentColor, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
for k = 1:nCases
    scatter(ax1, A{k}.nmpcH, A{k}.totalH, 46, plotColors(k, :), ...
        caseMarkers(k), 'LineWidth', 1.1, 'DisplayName', caseLabels(k));
end
xlim(ax1, lim); ylim(ax1, lim); axis(ax1, 'square');
xlabel(ax1, 'NMPC solve time $t_{\mathrm{nmpc}}$ (h)');
ylabel(ax1, 'Wall time $t_{\mathrm{total}}$ (h)');
title(ax1, '$\mathbf{a}$');
ax1.TitleHorizontalAlignment = 'left';
legend(ax1, 'Location', 'northwest', 'Box', 'off');

% --- (b) gap vs execution order ------------------------------------------
ax2 = nexttile; hold(ax2, 'on'); grid(ax2, 'off'); box(ax2, 'off');
for k = 1:nCases
    plot(ax2, A{k}.gidx, A{k}.gapMin, '-', 'Color', plotColors(k, :), ...
        'LineWidth', 1.1, 'Marker', caseMarkers(k), 'MarkerSize', 4.5, ...
        'HandleVisibility', 'off');
end
for c = unique(doeCount(doeCount > 0))'
    xline(ax2, c + 0.5, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
end
yline(ax2, 0, 'k-', 'LineWidth', 0.75, 'HandleVisibility', 'off');
xlabel(ax2, 'Evaluation (execution order, DOE then BO)');
ylabel(ax2, '$t_{\mathrm{total}} - t_{\mathrm{nmpc}}$ (min)');
title(ax2, '$\mathbf{b}$');
ax2.TitleHorizontalAlignment = 'left';

% --- (c) gap vs number of failed solves ----------------------------------
ax3 = nexttile; hold(ax3, 'on'); grid(ax3, 'off'); box(ax3, 'off');
for k = 1:nCases
    scatter(ax3, A{k}.nBad, A{k}.gapMin, 46, plotColors(k, :), ...
        caseMarkers(k), 'LineWidth', 1.1, 'HandleVisibility', 'off');
end
yline(ax3, 0, 'k-', 'LineWidth', 0.75, 'HandleVisibility', 'off');
xlabel(ax3, 'Failed solves per eval. (flag $\le 0$)');
ylabel(ax3, '$t_{\mathrm{total}} - t_{\mathrm{nmpc}}$ (min)');
title(ax3, '$\mathbf{c}$');
ax3.TitleHorizontalAlignment = 'left';

save_plot_outputs(fig, fullfile(graphicsDir, 'runtime_consistency.png'), ...
    fontSize, 1400, 460);

%% Numerical summary (console + txt)
txtPath = fullfile(numericalDir, 'runtime_consistency.txt');
fid = fopen(txtPath, 'w');
if fid == -1
    warning('runtime_consistency:txt', 'Could not open %s; console only.', txtPath);
    fid = 1;
end
emit(fid, 'Runtime consistency: t_total (wall_s.total) vs t_nmpc (runtime_s)');
emit(fid, 'gap = t_total - t_nmpc. Expected: small positive (minutes).');
emit(fid, 'gap >> 0 -> non-solver wall overhead; gap < 0 -> checkpoint resume artifact.');
emit(fid, '');

for k = 1:nCases
    T = A{k};
    emit(fid, '=== %s (%d evaluations: %d DOE + %d BO) ===', ...
        caseLabels(k), height(T), doeCount(k), height(T) - doeCount(k));
    emit(fid, '  t_total : median %.2f h | max %.2f h | sum %.1f h', ...
        median(T.totalH, 'omitnan'), max(T.totalH), sum(T.totalH, 'omitnan'));
    emit(fid, '  t_nmpc  : median %.2f h | max %.2f h | sum %.1f h', ...
        median(T.nmpcH, 'omitnan'), max(T.nmpcH), sum(T.nmpcH, 'omitnan'));
    emit(fid, '  gap     : median %.2f min | p95 %.2f min | max %.2f min | min %.2f min', ...
        median(T.gapMin, 'omitnan'), prctile_safe(T.gapMin, 95), ...
        max(T.gapMin), min(T.gapMin));
    emit(fid, '  gap/t_total: median %.2f %% | max %.2f %%', ...
        100 * median(T.gapMin * 60 ./ double(T.t_total), 'omitnan'), ...
        100 * max(T.gapMin * 60 ./ double(T.t_total)));

    % Correlations: is the wall time jump explained by solving or by the gap?
    rt = rank_corr(T.nBad, T.totalH);
    rn = rank_corr(T.nBad, T.nmpcH);
    rg = rank_corr(T.nBad, T.gapMin);
    emit(fid, '  Spearman rho vs failed solves: t_total %.3f | t_nmpc %.3f | gap %.3f', ...
        rt, rn, rg);

    % Offenders: gap above max(10 min, 10 %% of t_total), or negative gap.
    thr = max(10, 0.10 * double(T.t_total) / 60);
    bad = T.gapMin > thr | T.gapMin < -1;
    if any(bad)
        emit(fid, '  %d evaluation(s) with inconsistent timing:', nnz(bad));
        B = sortrows(T(bad, :), 'gapMin', 'descend');
        for i = 1:height(B)
            emit(fid, ['    id=%d | %s iter %3d | z=%.3f | t_total %6.2f h | ' ...
                't_nmpc %6.2f h | gap %8.2f min | failed solves %d%s'], ...
                B.id(i), string(B.phase_label(i)), B.iter(i), B.z(i), ...
                B.totalH(i), B.nmpcH(i), B.gapMin(i), B.nBad(i), ...
                resume_note(B.gapMin(i)));
        end
    else
        emit(fid, '  No inconsistent evaluations: t_total is a faithful runtime.');
    end
    emit(fid, '');
end

if fid ~= 1
    fclose(fid);
    fprintf('Wrote %s\n', txtPath);
end
end


%% ===================== helpers =====================

function emit(fid, fmt, varargin)
%EMIT Print one line to the file and mirror it on the console.
line = sprintf(fmt, varargin{:});
fprintf(fid, '%s\n', line);
if fid ~= 1
    fprintf('%s\n', line);
end
end


function r = rank_corr(x, y)
%RANK_CORR Spearman rank correlation without toolbox dependencies.
x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
if numel(x) < 3 || all(x == x(1)) || all(y == y(1))
    r = NaN;
    return
end
r = corrcoef(tied_rank(x), tied_rank(y));
r = r(1, 2);
end


function rk = tied_rank(v)
%TIED_RANK Average ranks with ties.
[~, ord] = sort(v);
rk = zeros(size(v));
rk(ord) = 1:numel(v);
u = unique(v);
for i = 1:numel(u)
    m = v == u(i);
    rk(m) = mean(rk(m));
end
end


function p = prctile_safe(v, q)
%PRCTILE_SAFE Percentile without the Statistics Toolbox.
v = sort(v(isfinite(v)));
if isempty(v)
    p = NaN;
    return
end
idx = max(1, min(numel(v), round(q / 100 * numel(v))));
p = v(idx);
end


function s = resume_note(gapMin)
%RESUME_NOTE Tag negative gaps as checkpoint-resume artifacts.
if gapMin < -1
    s = '  <- resumed run: t_total undercounts';
else
    s = '';
end
end


function label = pretty_case(name)
%PRETTY_CASE Turn a folder name like results_case1 into "Case 1".
tok = regexp(name, 'case(\d+)\s*$', 'tokens', 'once');
if isempty(tok)
    label = string(strrep(strrep(name, 'results_', ''), '_', ' '));
else
    label = "Case " + string(tok{1});
end
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
