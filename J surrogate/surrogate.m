% fit_surrogates_all.m
% Fit on all surrogate_data_<i>.mat files:
%   SSdU: Cheb5 in x on grid f = (1:N)/N  (length N)
%   SSE : Cheb5 in x on grid f = (1:N)/N  (length N)
%
% Figure-generation method summary
% - Build normalized cumulative trajectories for each surrogate dataset.
% - Fit global Chebyshev polynomials by least squares (shared across files).
% - Plot per-fidelity envelopes (min/max + mean) and fitted curves.
% - Export side-by-side and stacked diagnostics to results/graphical_results.

clear; close all; clc
rng(1)
scriptDir = fileparts(mfilename("fullpath"));
projectRoot = fileparts(scriptDir);
addpath(genpath(fullfile(projectRoot, "dependencies")));

set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
fontSize = 18;
NATURE_COLOR = nature_methods_colors();
plotColors = [NATURE_COLOR.Blue; [0 0 0]];

files = dir(fullfile(scriptDir, "surrogate_data_*.mat"));
if isempty(files)
    error('No files matching surrogate_data_*.mat in the current folder.');
end

% -----------------------------
% Aggregate all data
% -----------------------------
x_all     = [];
SSdU_all  = [];
SSE_all   = [];
file_id   = [];
f_cases    = {};
SSdU_cases = {};
SSE_cases  = {};

for kf = 1:numel(files)
    S = load(fullfile(files(kf).folder, files(kf).name));
    if ~isfield(S, "out")
        error("File %s does not contain variable 'out'.", files(kf).name);
    end
    out = S.out;

    % Grid
    N  = length(out.case(1).SSE);
    f  = (1:N).'/N;
    x  = 2*f - 1;

    % Signals
    SSdU = [0; cumsum(out.case(1).SSdU)] + [0; cumsum(out.case(2).SSdU)];
    SSdU = SSdU / SSdU(end);

    SSE  = cumsum(out.case(1).SSE) + cumsum(out.case(2).SSE);
    SSE  = SSE / SSE(end);

    % Sanity checks (your intended construction)
    if numel(SSdU) ~= N
        error("File %s: expected SSdU length N, got %d (N=%d).", files(kf).name, numel(SSdU), N);
    end
    if numel(SSE) ~= N
        error("File %s: expected SSE length N, got %d (N=%d).", files(kf).name, numel(SSE), N);
    end

    % Append
    x_all    = [x_all;    x];
    SSdU_all = [SSdU_all; SSdU];
    SSE_all  = [SSE_all;  SSE];
    file_id  = [file_id;  repmat(kf, N, 1)];

    f_cases{end+1} = f;
    SSdU_cases{end+1} = SSdU;
    SSE_cases{end+1} = SSE;
end

% -----------------------------
% Optimisation options
% -----------------------------
opts = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'Display','off', ...
    'MaxIterations', 2000, ...
    'OptimalityTolerance', 1e-10, ...
    'StepTolerance', 1e-12);

% -----------------------------
% Fit SSdU and SSE with Cheb5
% -----------------------------
c5_0   = zeros(6,1);
c_SSdU = fminunc(@(c) obj_Cheb5(c, x_all, SSdU_all), c5_0, opts);
c_SSE  = fminunc(@(c) obj_Cheb5(c, x_all, SSE_all ), c5_0, opts);

% Goodness-of-fit on the aggregated training set
SSdU_hat_all = clamp01(Cheb5(x_all, c_SSdU));
SSE_hat_all  = clamp01(Cheb5(x_all, c_SSE));
R2_SSdU = compute_r2(SSdU_all, SSdU_hat_all);
R2_SSE  = compute_r2(SSE_all,  SSE_hat_all);

% Check rng(1) results. If non destructive changes are made to the code, then
% the results should be the same within rouding error
c_SSdU_check = [6.442657e-01   4.682368e-01  -1.455242e-01   2.457724e-02   2.198219e-02  -1.598959e-02]';
c_SSE_check = [7.827330e-01   3.771709e-01  -2.433549e-01   1.091471e-01  -2.846259e-02   1.358572e-03]';
exp_SSdU = floor(log10(abs(c_SSdU_check)));
exp_SSE  = floor(log10(abs(c_SSE_check)));

% Absolute rounding tolerance per entry (within significant digits)
significant_digits = 6;
tol_SSdU = 0.5 * 10.^(exp_SSdU - significant_digits);
tol_SSE  = 0.5 * 10.^(exp_SSE  - significant_digits);

if all(abs(c_SSdU - c_SSdU_check) <= tol_SSdU) && ...
   all(abs(c_SSE  - c_SSE_check ) <= tol_SSE)

    fprintf("Chebyshev coefficient check passed (within printed rounding tolerance).\n");

else
    error("Chebyshev coefficient check failed: deviation exceeds printed rounding tolerance.");
end

% -----------------------------
% Print coefficients
% -----------------------------
fprintf('Cheb5 coeffs SSdU (c0..c5):\n');
fprintf('% .6e  ', c_SSdU); fprintf('\n\n');

fprintf('Cheb5 coeffs SSE  (c0..c5):\n');
fprintf('% .6e  ', c_SSE); fprintf('\n\n');

fprintf('\nR^2 on aggregated fit data:\n');
fprintf('R^2(SSdU) = %.6f\n', R2_SSdU);
fprintf('R^2(SSE)  = %.6f\n', R2_SSE);

% Export numerical summary for SSdU/SSE surrogate coefficients.
numDir = fullfile(projectRoot, "results", "numerical results");
if ~isfolder(numDir)
    mkdir(numDir);
end
numTxtPath = fullfile(numDir, "surrogate_cheb_coeffs.txt");
fid = fopen(numTxtPath, "w");
if fid == -1
    error("Unable to open file for writing: %s", numTxtPath);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, "=== Surrogate Coefficients (Cheb5) ===\n");
fprintf(fid, "SSdU c0..c5:\n");
fprintf(fid, "% .10e\n", c_SSdU);
fprintf(fid, "\nSSE c0..c5:\n");
fprintf(fid, "% .10e\n", c_SSE);
fprintf(fid, "\nR^2:\n");
fprintf(fid, "R^2(SSdU) = %.10f\n", R2_SSdU);
fprintf(fid, "R^2(SSE)  = %.10f\n", R2_SSE);

%% -----------------------------
% Quick diagnostics plots (optional)
% -----------------------------
[f_ref, SSdU_mat, SSE_mat] = build_case_matrices(f_cases, SSdU_cases, SSE_cases);
x_ref = 2*f_ref - 1;
SSdU_hat_ref = clamp01(Cheb5(x_ref, c_SSdU));
SSE_hat_ref  = clamp01(Cheb5(x_ref, c_SSE));

SSdU_min = min(SSdU_mat, [], 2);
SSdU_max = max(SSdU_mat, [], 2);
SSdU_mean = mean(SSdU_mat, 2);

SSE_min = min(SSE_mat, [], 2);
SSE_max = max(SSE_mat, [], 2);
SSE_mean = mean(SSE_mat, 2);

fig_side = figure;
subplot(1,2,1);
plot_case_envelope_with_fit(f_ref, SSdU_min, SSdU_max, SSdU_mean, SSdU_hat_ref, ...
    'Fidelity $z$ (dimensionless)', '$\frac{J_{\mathrm{TV}}(z)}{\left.J_{\mathrm{TV}}\right|_{z=1}}$', ...
    'a', plotColors, fontSize, true);

subplot(1,2,2);
plot_case_envelope_with_fit(f_ref, SSE_min, SSE_max, SSE_mean, SSE_hat_ref, ...
    'Fidelity $z$ (dimensionless)', '$\frac{J_{\mathrm{track}}(z)}{\left.J_{\mathrm{track}}\right|_{z=1}}$', ...
    'b', plotColors, fontSize, true);
set_fig_size(1300, 520);

fig_stack = figure;
subplot(2,1,1);
plot_case_envelope_with_fit(f_ref, SSdU_min, SSdU_max, SSdU_mean, SSdU_hat_ref, ...
    'Fidelity $z$ (dimensionless)', '$\frac{J_{\mathrm{TV}}(z)}{\left.J_{\mathrm{TV}}\right|_{z=1}}$', ...
    'a', plotColors, fontSize, false);

subplot(2,1,2);
plot_case_envelope_with_fit(f_ref, SSE_min, SSE_max, SSE_mean, SSE_hat_ref, ...
    'Fidelity $z$ (dimensionless)', '$\frac{J_{\mathrm{track}}(z)}{\left.J_{\mathrm{track}}\right|_{z=1}}$', ...
    'b', plotColors, fontSize, true);
set_fig_size(920, 900);

plotDir = fullfile(projectRoot, "results", "graphical_results");
if ~isfolder(plotDir)
    mkdir(plotDir);
end
save_figure_pair(fig_side, fullfile(plotDir, "surrogate_diagnostics_side_by_side"));
save_figure_pair(fig_stack, fullfile(plotDir, "surrogate_diagnostics_stacked"));

% =============================
% Objectives
% =============================
function J = obj_Cheb5(c, x, y)
%OBJ_CHEB5 Least-squares objective for Chebyshev coefficient fitting.
    r = y - Cheb5(x, c);
    J = r.'*r;
    if ~isfinite(J); J = realmax; end
end

% =============================
% Models
% =============================
function y = Cheb5(x, c)
%CHEB5 Evaluate degree-5 Chebyshev polynomial using c0..c5 coefficients.
    c = c(:);
    if numel(c) ~= 6
        error('Cheb5 expects 6 coefficients (c0..c5).');
    end
    T0 = ones(size(x));
    T1 = x;
    T2 = 2*x.^2 - 1;
    T3 = 4*x.^3 - 3*x;
    T4 = 8*x.^4 - 8*x.^2 + 1;
    T5 = 16*x.^5 - 20*x.^3 + 5*x;

    y = c(1)*T0 + c(2)*T1 + c(3)*T2 + c(4)*T3 + c(5)*T4 + c(6)*T5;
end

function y = clamp01(y)
%CLAMP01 Clip values to [0, 1].
    y = min(max(y, 0), 1);
end

function r2 = compute_r2(y, y_hat)
%COMPUTE_R2 Compute coefficient of determination for fit diagnostics.
    y = y(:);
    y_hat = y_hat(:);
    sse_res = sum((y - y_hat).^2);
    sse_tot = sum((y - mean(y)).^2);
    if sse_tot <= eps
        r2 = NaN;
    else
        r2 = 1 - sse_res / sse_tot;
    end
end

function [f_ref, SSdU_mat, SSE_mat] = build_case_matrices(f_cases, SSdU_cases, SSE_cases)
%BUILD_CASE_MATRICES Align case trajectories onto one reference fidelity grid.
    K = numel(f_cases);
    if K == 0
        error('No cases available for diagnostics plotting.');
    end

    f_ref = f_cases{1}(:);
    Nref = numel(f_ref);
    SSdU_mat = nan(Nref, K);
    SSE_mat = nan(Nref, K);

    for k = 1:K
        fk = f_cases{k}(:);
        SSdUk = SSdU_cases{k}(:);
        SSEk = SSE_cases{k}(:);

        if numel(fk) == Nref && max(abs(fk - f_ref)) < 1e-12
            SSdU_mat(:, k) = SSdUk;
            SSE_mat(:, k) = SSEk;
        else
            SSdU_mat(:, k) = interp1(fk, SSdUk, f_ref, 'linear', 'extrap');
            SSE_mat(:, k) = interp1(fk, SSEk, f_ref, 'linear', 'extrap');
        end
    end
end

function plot_case_envelope_with_fit(f, y_min, y_max, y_mean, y_hat, xlab, ylab, titleLabel, plotColors, fontSize, showXLabel)
%PLOT_CASE_ENVELOPE_WITH_FIT Plot min-max envelope, mean trend, and fitted curve.
    if nargin < 11
        showXLabel = true;
    end

    if size(plotColors, 1) < 1 || size(plotColors, 2) ~= 3
        error("plotColors must be an N-by-3 RGB matrix with at least one row.");
    end
    fitColorIdx = min(3, size(plotColors, 1));

    fill([f; flipud(f)], [y_min; flipud(y_max)], plotColors(1,:), ...
        'FaceAlpha', 0.20, 'EdgeColor', 'none'); hold on
    plot(f, y_mean, '-', 'LineWidth', 2.2, 'Color', plotColors(1,:));
    plot(f, y_hat, '--', 'LineWidth', 2.2, 'Color', plotColors(fitColorIdx,:));

    xlim([0, 1]);
    ylim([0, 1.005]);
    if showXLabel
        xlabel(xlab);
    else
        xlabel("");
    end
    ylabel(ylab);
    t = title(sprintf('\\textbf{%s}', titleLabel), 'Interpreter', 'latex');
    t.Units = 'normalized';
    t.Position(1) = 0;
    t.HorizontalAlignment = 'left';
    grid off; box off
    set_font_size(fontSize);
    format_tick(1, 1);
end

function save_figure_pair(figHandle, fileStem)
%SAVE_FIGURE_PAIR Save one figure to PNG and PDF as a single combined layout.
    exportgraphics(figHandle, fileStem + ".png", "Resolution", 300);
    exportgraphics(figHandle, fileStem + ".pdf", "ContentType", "vector");
end
