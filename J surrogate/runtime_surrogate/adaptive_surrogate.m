%% Runtime adaptive Chebyshev surrogate fitting with ridge regression
%
% Pipeline
%   1. Load all out_*.mat simulation pairs (each contains case 1 and 2).
%   2. For every pair, build cumulative SSdU and SSE in fidelity z, then
%      linearly interpolate onto the fixed grid z = 0:0.05:1, restricted to
%      the z range actually covered by the simulation.
%   3. Initialise the Chebyshev fit on the first 20 simulation pairs, warm
%      started from c = [0, 1, 0, 0, 0] (c0 derived, c1 = 1, c2..c4 = 0).
%   4. Refit after every new pair is added, warm starting from the previous
%      coefficients.
%   5. Every refit selects lambda via K = 10 fold CV, with folds grouped by
%      simulation pair so the held-out pairs are fully held out.
%   6. Diagnostics tracked per step:
%        lambda*       chosen regularisation strength
%        coefficients  c0..c4 (c0 derived from f(0) = 0 constraint)
%        seen R^2      coefficient of determination on training pool
%        seen loss     residual SSE on training pool
%        unseen R^2    coefficient of determination on remaining pairs
%        unseen loss   residual SSE on remaining pairs
%
% Model
%   y(z) = sum_{i=0..4} c_i * T_i(2z - 1),  x = 2z - 1 in [-1, 1].
%   T_k built via recurrence: T_0 = 1, T_1 = x, T_{k+1} = 2x*T_k - T_{k-1}.
%
% Endpoint constraints (passed as nonlinear equalities via fmincon nonlcon):
%   f(z=0) = 0  ->  ceq(1) = X(z=0) * c     = 0
%   f(z=1) = 1  ->  ceq(2) = X(z=1) * c - 1 = 0
%
% Ridge fit
%   c* = argmin_c  || y - X * c ||^2 + lambda * || c ||^2
%        subject to  Aeq * c = beq  (the two endpoint constraints above)
%   Solved as a general NLP with fmincon (sqp), analytic gradient.

clear; clc; close all;

%% Configuration
fullHorizonHours   = 10.0;                  % T in z = t / T
zStep              = 0.05;                  % uniform z resampling step
nInit              = 5;                    % initial number of pairs
K_cv               = 5;                    % K-fold CV folds
chebOrder          = 20;                    % nth order Chebyshev
lambdaGrid         = logspace(-3, 3, 7);    % CV grid for lambda
c0WarmStart        = [0.5; 0.5; zeros(chebOrder - 1, 1)];  % satisfies f(0)=0 and f(1)=1 for any order
rngSeed            = 1;                     % reproducibility of CV folds

%% Paths
scriptDir  = fileparts(mfilename("fullpath"));
surrogDir  = fileparts(scriptDir);
projectRoot = fileparts(surrogDir);
addpath(genpath(fullfile(projectRoot, "dependencies")));

set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");

fontSize = 14;
hasNatureColors = exist('nature_methods_colors', 'file') == 2;
if hasNatureColors
    NC = nature_methods_colors();
    colSSdU = NC.Blue;
    colSSE  = NC.BluishGreen;
    colAccent = NC.ReddishPurple;
else
    colSSdU = [0.00 0.45 0.74];
    colSSE  = [0.00 0.62 0.45];
    colAccent = [0.80 0.40 0.69];
end

%% Load data
runDirs = {
    fullfile(projectRoot, "results", "run1");
    fullfile(projectRoot, "results", "run2");
};

zGridFull = (0:zStep:1).';

pairs = load_all_pairs(runDirs, fullHorizonHours, zGridFull);
nPairs = numel(pairs);
if nPairs < nInit + 1
    error("Need at least %d pairs (got %d) to run init + one update.", ...
        nInit + 1, nPairs);
end
fprintf("Loaded %d simulation pairs (initial %d, then %d adaptive updates).\n", ...
    nPairs, nInit, nPairs - nInit);

%% Pre-compute Cheb features for the full z grid (for plotting)
Xgrid = cheb_features(zGridFull, chebOrder);

%% Progressive fitting
% Step k = 1 uses pairs 1..nInit (initial fit).
% Step k > 1 adds pair nInit + k - 1.
nSteps = nPairs - nInit + 1;

history = init_history(nSteps, chebOrder, numel(lambdaGrid));

c_prev_dU = c0WarmStart;
c_prev_E  = c0WarmStart;

rng(rngSeed);

t0 = tic;
for step = 1:nSteps

    nSeen     = nInit + step - 1;
    seenIdx   = 1:nSeen;
    unseenIdx = (nSeen + 1):nPairs;

    % SSdU
    [c_dU, lam_dU, cv_dU] = fit_with_cv( ...
        pairs(seenIdx), 'ratioSSdU', c_prev_dU, chebOrder, lambdaGrid, K_cv);

    history.step(step)        = step;
    history.nSeen(step)       = nSeen;
    history.lambda_dU(step)   = lam_dU;
    history.coef_dU(:, step)  = c_dU;
    history.cvLoss_dU(:, step) = cv_dU(:);

    [history.seenR2_dU(step), history.seenLoss_dU(step)] = eval_on_pairs( ...
        pairs(seenIdx), 'ratioSSdU', c_dU);

    if ~isempty(unseenIdx)
        [history.unseenR2_dU(step), history.unseenLoss_dU(step)] = ...
            eval_on_pairs(pairs(unseenIdx), 'ratioSSdU', c_dU);
    else
        history.unseenR2_dU(step)   = NaN;
        history.unseenLoss_dU(step) = NaN;
    end
    c_prev_dU = c_dU;

    % SSE
    [c_E, lam_E, cv_E] = fit_with_cv( ...
        pairs(seenIdx), 'ratioSSE', c_prev_E, chebOrder, lambdaGrid, K_cv);

    history.lambda_E(step)   = lam_E;
    history.coef_E(:, step)  = c_E;
    history.cvLoss_E(:, step) = cv_E(:);

    [history.seenR2_E(step), history.seenLoss_E(step)] = eval_on_pairs( ...
        pairs(seenIdx), 'ratioSSE', c_E);

    if ~isempty(unseenIdx)
        [history.unseenR2_E(step), history.unseenLoss_E(step)] = ...
            eval_on_pairs(pairs(unseenIdx), 'ratioSSE', c_E);
    else
        history.unseenR2_E(step)   = NaN;
        history.unseenLoss_E(step) = NaN;
    end
    c_prev_E = c_E;

    if mod(step, 1) == 0 || step == 1 || step == nSteps
        fprintf("Step %3d/%d | seen=%3d unseen=%3d | lambda_dU=%.2e lambda_E=%.2e | R2u_dU=%.3f R2u_E=%.3f\n", ...
            step, nSteps, nSeen, nPairs - nSeen, ...
            lam_dU, lam_E, history.unseenR2_dU(step), history.unseenR2_E(step));
    end
end
fprintf("Total adaptive runtime: %.2f s\n", toc(t0));

%% Diagnostics plots
plot_progression(history, "ratioSSdU", "$J_{TV}$", colSSdU, fontSize, ...
    "Adaptive Chebyshev surrogate -- SSdU");
plot_progression(history, "ratioSSE",  "$J_{track}$", colSSE,  fontSize, ...
    "Adaptive Chebyshev surrogate -- SSE");
plot_r2_comparison(history, colSSdU, colSSE, fontSize);
%%
plot_final_fit(zGridFull, Xgrid, pairs, history, "ratioSSdU", "$J_{TV}$", ...
    colSSdU, colAccent, fontSize, ...
    "Final adaptive Chebyshev fit -- SSdU");
plot_final_fit(zGridFull, Xgrid, pairs, history, "ratioSSE",  "$J_{track}$", ...
    colSSE,  colAccent, fontSize, ...
    "Final adaptive Chebyshev fit -- SSE");

%% Save numeric summary
outDir = fullfile(projectRoot, "results", "numerical results");
if ~isfolder(outDir); mkdir(outDir); end

txtPath = fullfile(outDir, "surrogate_adaptive_runtime_summary.txt");
fid = fopen(txtPath, "w");
if fid == -1
    warning("Could not open %s for writing.", txtPath);
else
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, "=== Adaptive Chebyshev Surrogate Runtime Summary ===\n");
    fprintf(fid, "Pairs total: %d\nInit pairs: %d\nK-fold: %d\n", ...
        nPairs, nInit, K_cv);
    fprintf(fid, "Lambda grid: %d points in [%.1e, %.1e]\n\n", ...
        numel(lambdaGrid), lambdaGrid(1), lambdaGrid(end));

    fprintf(fid, "Final SSdU coefficients (c0..c4):\n");
    fprintf(fid, "% .6e\n", history.coef_dU(:, end));
    fprintf(fid, "Final SSE coefficients (c0..c4):\n");
    fprintf(fid, "% .6e\n", history.coef_E(:, end));

    fprintf(fid, "\nstep,nSeen,lambda_dU,lambda_E,seenR2_dU,seenR2_E,unseenR2_dU,unseenR2_E\n");
    for s = 1:nSteps
        fprintf(fid, "%d,%d,%.6e,%.6e,%.6f,%.6f,%.6f,%.6f\n", ...
            history.step(s), history.nSeen(s), ...
            history.lambda_dU(s), history.lambda_E(s), ...
            history.seenR2_dU(s), history.seenR2_E(s), ...
            history.unseenR2_dU(s), history.unseenR2_E(s));
    end
    fprintf("Saved summary to %s\n", txtPath);
end

%% Local functions

function pairs = load_all_pairs(runDirs, fullHorizonHours, zGridFull)
%LOAD_ALL_PAIRS Walk run directories, build one entry per simulation pair.
    pairs = struct('z', {}, 'ratioSSdU', {}, 'ratioSSE', {}, ...
                   'zEnd', {}, 'fileName', {}, 'runLabel', {});

    for r = 1:numel(runDirs)
        runPath = runDirs{r};
        if ~isfolder(runPath)
            warning("Run directory missing: %s", runPath);
            continue
        end
        files = dir(fullfile(runPath, "out_*.mat"));
        [~, idx] = sort({files.name});
        files = files(idx);
        runLabel = sprintf("run%d", r);

        for k = 1:numel(files)
            fpath = fullfile(files(k).folder, files(k).name);
            S = load(fpath, "out");
            if ~isfield(S, "out"); continue; end
            out = S.out;
            if ~isfield(out, "case") || numel(out.case) < 2; continue; end

            [z, rdU, rE, zEnd] = interpolate_pair(out, fullHorizonHours, zGridFull);
            if isempty(z); continue; end

            pairs(end+1) = struct( ... %#ok<AGROW>
                'z',         z, ...
                'ratioSSdU', rdU, ...
                'ratioSSE',  rE, ...
                'zEnd',      zEnd, ...
                'fileName',  files(k).name, ...
                'runLabel',  runLabel);
        end
    end
end

function [zOut, rdU, rE, zEnd] = interpolate_pair(out, T_full, zGridFull)
%INTERPOLATE_PAIR Reproduce the cumulative ratio construction of
% test_surrogate.m, then linearly interpolate onto the uniform z grid,
% anchored at z = 0 with ratio = 0.
    zOut = []; rdU = []; rE = []; zEnd = NaN;

    if ~isfield(out, "T"); return; end
    p1_dU = as_col(out.case(1).partial_SSdU);
    p2_dU = as_col(out.case(2).partial_SSdU);
    p1_e  = as_col(out.case(1).partial_SSE);
    p2_e  = as_col(out.case(2).partial_SSE);
    t     = as_col(out.T);

    NdU = min([numel(p1_dU), numel(p2_dU), max(numel(t) - 1, 0)]);
    NSE = min([numel(p1_e),  numel(p2_e),  numel(t)]);
    if NdU < 3 || NSE < 3; return; end

    p1_dU = p1_dU(1:NdU); p2_dU = p2_dU(1:NdU);
    p1_e  = p1_e(1:NSE);  p2_e  = p2_e(1:NSE);
    tSSdU = t(2:(NdU + 1));
    tSSE  = t(1:NSE);

    cumSSdU = cumsum(p1_dU) + cumsum(p2_dU);
    cumSSE  = cumsum(p1_e)  + cumsum(p2_e);

    final_dU = safe_scalar(out, "SSdU", NaN);
    final_E  = safe_scalar(out, "SSE",  NaN);
    if ~isfinite(final_dU) || final_dU <= 0; return; end
    if ~isfinite(final_E)  || final_E  <= 0; return; end

    ratioRaw_dU = cumSSdU / final_dU;
    ratioRaw_E  = cumSSE  / final_E;

    fSSdU = tSSdU / T_full;
    fSSE  = tSSE  / T_full;

    zEnd_dU = max(fSSdU);
    zEnd_E  = max(fSSE);
    zEnd    = min(zEnd_dU, zEnd_E);

    if ~isfinite(zEnd) || zEnd <= zGridFull(2)
        return
    end

    mask = zGridFull <= zEnd + 1e-12;
    zOut = zGridFull(mask);
    if numel(zOut) < 3
        zOut = []; return
    end

    % Anchor at z = 0 with ratio = 0 (no cost accumulated yet).
    rdU = interp1(fSSdU, ratioRaw_dU, zOut, 'linear', 'extrap');
    rE  = interp1(fSSE,  ratioRaw_E,  zOut, 'linear', 'extrap');

    keep = isfinite(rdU) & isfinite(rE);
    zOut = zOut(keep);
    rdU  = rdU(keep);
    rE   = rE(keep);

    if numel(zOut) < 3
        zOut = []; rdU = []; rE = []; return
    end
end

function v = safe_scalar(S, fieldName, fallback)
%SAFE_SCALAR Tolerant scalar extraction.
    if isfield(S, fieldName)
        v = double(S.(fieldName));
        if numel(v) ~= 1 || ~isfinite(v)
            v = fallback;
        end
    else
        v = fallback;
    end
end

function x = as_col(x)
%AS_COL Force column-vector double.
    x = double(x(:));
end

function X = cheb_features(z, order)
%CHEB_FEATURES Build the Chebyshev design matrix on the [0,1] -> [-1,1] map.
%
% Evaluates T_0(x) .. T_order(x) with x = 2*z - 1.
% Columns built via the three-term recurrence:
%   T_0 = 1,  T_1 = x,  T_{k+1} = 2*x*T_k - T_{k-1}.
    z = z(:);
    x = 2 * z - 1;
    X = zeros(numel(z), order + 1);
    X(:, 1) = 1;
    if order < 1; return; end
    X(:, 2) = x;
    for k = 2:order
        X(:, k+1) = 2 * x .* X(:, k) - X(:, k-1);
    end
end

function [X, y] = stack_pairs(pairs, yField, chebOrder)
%STACK_PAIRS Concatenate (X, y) across the given simulation pairs.
    if isempty(pairs); X = []; y = []; return; end
    zCells = cell(numel(pairs), 1);
    yCells = cell(numel(pairs), 1);
    for i = 1:numel(pairs)
        zCells{i} = pairs(i).z;
        yCells{i} = pairs(i).(yField);
    end
    z = vertcat(zCells{:});
    y = vertcat(yCells{:});
    X = cheb_features(z, chebOrder);
end

function c = fit_ridge_cheb(X, y, lambda, cWarm, chebOrder)
%FIT_RIDGE_CHEB Ridge regression Chebyshev fit with two endpoint constraints.
%
% Constraints enforced via fmincon nonlcon (nonlinear equalities):
%   f(z=0) = 0  ->  ceq(1) = X(z=0) * c       = 0
%   f(z=1) = 1  ->  ceq(2) = X(z=1) * c - 1   = 0
%
% Objective (general NLP, solved with fmincon sqp):
%   min_c  || y - X * c ||^2 + lambda * || c ||^2
%   Analytic gradient supplied.

    X0 = cheb_features(0, chebOrder);   % 1 x (order+1), evaluated at z = 0
    X1 = cheb_features(1, chebOrder);   % 1 x (order+1), evaluated at z = 1

    opts = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'off', ...
        'SpecifyObjectiveGradient', true, ...
        'OptimalityTolerance', 1e-8, ...
        'StepTolerance', 1e-8, ...
        'MaxIterations', 400, ...
        'MaxFunctionEvaluations', 4000);

    c = fmincon(@(c) ridge_obj(c, X, y, lambda), cWarm, ...
        [], [], [], [], [], [], @(c) endpoint_con(c, X0, X1), opts);
end

function [cineq, ceq] = endpoint_con(c, X0, X1)
%ENDPOINT_CON Nonlinear equality constraints: f(z=0) = 0, f(z=1) = 1.
    c     = c(:);
    cineq = [];
    ceq   = [X0 * c; X1 * c - 1];
end

function [val, grad] = ridge_obj(c, X, y, lambda)
%RIDGE_OBJ L2 loss + L2 regularisation with analytic gradient.
    c = c(:);
    r = y - X * c;
    val = r' * r + lambda * (c' * c);
    if nargout > 1
        grad = -2 * (X' * r) + 2 * lambda * c;
    end
end

function [c, lambdaBest, cvLossPerLambda] = fit_with_cv( ...
    pairsSeen, yField, cWarm, chebOrder, lambdaGrid, K)
%FIT_WITH_CV Choose lambda by K-fold CV (grouped by simulation pair),
% then refit on the full seen set with the chosen lambda.
%
% Folds are assigned by simulation pair so every point from a held-out
% pair is excluded together. Within each fold the lambda grid is swept
% from largest to smallest so each fit warm-starts from the previous
% (more regularised) solution along the regularisation path.
%
% CV is skipped when:
%   numel(lambdaGrid) == 1  ->  that lambda is used directly.
%   K < 2                   ->  lambda = 0 is used.
    nP   = numel(pairsSeen);
    nLam = numel(lambdaGrid);

    [Xall, yall] = stack_pairs(pairsSeen, yField, chebOrder);

    if nLam == 1
        lambdaBest      = lambdaGrid(1);
        cvLossPerLambda = NaN;
        c = fit_ridge_cheb(Xall, yall, lambdaBest, cWarm, chebOrder);
        return
    end

    if K < 2
        lambdaBest      = 0;
        cvLossPerLambda = NaN(nLam, 1);
        c = fit_ridge_cheb(Xall, yall, lambdaBest, cWarm, chebOrder);
        return
    end

    Keff = min(K, nP);

    % Shuffle pair indices and assign fold labels.
    perm     = randperm(nP);
    foldSize = floor(nP / Keff);
    leftover = nP - foldSize * Keff;
    foldAssign = zeros(nP, 1);
    cursor = 1;
    for f = 1:Keff
        extra = double(f <= leftover);
        nThis = foldSize + extra;
        foldAssign(perm(cursor:cursor + nThis - 1)) = f;
        cursor = cursor + nThis;
    end

    % Sweep lambdas largest-to-smallest (regularisation path ordering).
    [lamSorted, sortIdx] = sort(lambdaGrid, 'descend');

    cvLossPerLambda = zeros(nLam, 1);

    for f = 1:Keff
        valIdx = find(foldAssign == f);
        trnIdx = find(foldAssign ~= f);
        if isempty(valIdx) || isempty(trnIdx); continue; end

        [Xtr, ytr] = stack_pairs(pairsSeen(trnIdx), yField, chebOrder);
        [Xva, yva] = stack_pairs(pairsSeen(valIdx), yField, chebOrder);
        if isempty(ytr) || isempty(yva); continue; end

        % Fit along the lambda path for this fold.
        cPath  = cWarm;
        cByLam = zeros(numel(cWarm), nLam);
        for li = 1:nLam
            cPath = fit_ridge_cheb(Xtr, ytr, lamSorted(li), cPath, chebOrder);
            cByLam(:, sortIdx(li)) = cPath;
        end

        % Accumulate held-out squared error for each lambda.
        for li = 1:nLam
            res = yva - Xva * cByLam(:, li);
            cvLossPerLambda(li) = cvLossPerLambda(li) + sum(res.^2);
        end
    end

    [~, bestIdx] = min(cvLossPerLambda);
    lambdaBest = lambdaGrid(bestIdx);

    % Final fit on the full seen pool with the chosen lambda.
    c = fit_ridge_cheb(Xall, yall, lambdaBest, cWarm, chebOrder);
end

function [r2, lossSSE] = eval_on_pairs(pairsEval, yField, c)
%EVAL_ON_PAIRS Pool all points across the given pairs, return R^2 and SSE.
    chebOrder = numel(c) - 1;
    [X, y] = stack_pairs(pairsEval, yField, chebOrder);
    if isempty(y)
        r2 = NaN; lossSSE = NaN; return
    end
    yhat = X * c;
    res = y - yhat;
    lossSSE = sum(res.^2);
    sse_tot = sum((y - mean(y)).^2);
    if sse_tot <= eps
        r2 = NaN;
    else
        r2 = 1 - lossSSE / sse_tot;
    end
end

function H = init_history(nSteps, chebOrder, nLam)
%INIT_HISTORY Pre-allocate the diagnostics container.
    H = struct();
    H.step           = zeros(nSteps, 1);
    H.nSeen          = zeros(nSteps, 1);
    H.lambda_dU      = nan(nSteps, 1);
    H.lambda_E       = nan(nSteps, 1);
    H.coef_dU        = nan(chebOrder + 1, nSteps);
    H.coef_E         = nan(chebOrder + 1, nSteps);
    H.seenR2_dU      = nan(nSteps, 1);
    H.seenR2_E       = nan(nSteps, 1);
    H.seenLoss_dU    = nan(nSteps, 1);
    H.seenLoss_E     = nan(nSteps, 1);
    H.unseenR2_dU    = nan(nSteps, 1);
    H.unseenR2_E     = nan(nSteps, 1);
    H.unseenLoss_dU  = nan(nSteps, 1);
    H.unseenLoss_E   = nan(nSteps, 1);
    H.cvLoss_dU      = nan(nLam, nSteps);
    H.cvLoss_E       = nan(nLam, nSteps);
end

function plot_progression(H, yField, yLabelTex, baseColor, fontSize, figTitle)
%PLOT_PROGRESSION Three-panel summary: unseen R^2, coefficients, lambda.
    isSSdU = strcmp(yField, "ratioSSdU");
    if isSSdU
        lambdaTrace = H.lambda_dU;
        coefTrace   = H.coef_dU;
        unseenR2    = H.unseenR2_dU;
        seenR2      = H.seenR2_dU;
    else
        lambdaTrace = H.lambda_E;
        coefTrace   = H.coef_E;
        unseenR2    = H.unseenR2_E;
        seenR2      = H.seenR2_E;
    end

    nSeen = H.nSeen;

    fig = figure('Name', figTitle);

    subplot(3, 1, 1); hold on
    plot(nSeen, seenR2,   '-', 'LineWidth', 1.6, 'Color', [baseColor, 0.55]);
    plot(nSeen, unseenR2, '-', 'LineWidth', 2.2, 'Color', baseColor);
    ylabel(sprintf("$R^2$ for %s", yLabelTex));
    legend({"seen", "unseen"}, "Location", "best");
    grid on; box on
    if exist('set_font_size', 'file') == 2; set_font_size(fontSize); end

    subplot(3, 1, 2); hold on
    nC = size(coefTrace, 1);
    cmap = lines(nC);
    for j = 1:nC
        plot(nSeen, coefTrace(j, :), '-', 'LineWidth', 1.6, 'Color', cmap(j, :));
    end
    ylabel("Chebyshev coefficients");
    legend(arrayfun(@(k) sprintf("c_%d", k - 1), 1:nC, ...
        "UniformOutput", false), "Location", "best");
    grid on; box on
    if exist('set_font_size', 'file') == 2; set_font_size(fontSize); end

    subplot(3, 1, 3);
    semilogy(nSeen, lambdaTrace, '-', 'LineWidth', 2.0, 'Color', baseColor);
    ylabel("$\lambda^\star$");
    xlabel("Number of seen simulation pairs");
    grid on; box on
    if exist('set_font_size', 'file') == 2; set_font_size(fontSize); end

    sgtitle(figTitle, "Interpreter", "none");
    if exist('set_fig_size', 'file') == 2
        set_fig_size(900, 850);
    end
end

function plot_r2_comparison(H, colSSdU, colSSE, fontSize)
%PLOT_R2_COMPARISON Seen and unseen R^2 for both surrogates on one figure.
    nSeen = H.nSeen;

    figure('Name', 'R^2 comparison -- SSdU vs SSE'); hold on
    plot(nSeen, H.seenR2_dU,   '-', 'LineWidth', 1.6, 'Color', [colSSdU, 0.45]);
    plot(nSeen, H.unseenR2_dU, '-', 'LineWidth', 2.2, 'Color', colSSdU);
    plot(nSeen, H.seenR2_E,    '-', 'LineWidth', 1.6, 'Color', [colSSE,  0.45]);
    plot(nSeen, H.unseenR2_E,  '-', 'LineWidth', 2.2, 'Color', colSSE);
    yline(0, ':', 'Color', [0.4 0.4 0.4]);
    legend({"$J_{TV}$ seen", "$J_{TV}$ unseen", ...
            "$J_{track}$ seen", "$J_{track}$ unseen"}, ...
        "Location", "best");
    xlabel("Number of seen simulation pairs");
    ylabel("$R^2$");
    title("Surrogate generalisation: $J_{TV}$ vs $J_{track}$", ...
        "Interpreter", "latex");
    grid on; box on
    if exist('set_font_size', 'file') == 2; set_font_size(fontSize); end
    if exist('set_fig_size', 'file') == 2;  set_fig_size(900, 400);  end
end

function plot_final_fit(zGrid, Xgrid, pairs, H, yField, yLabelTex, ...
    baseColor, accentColor, fontSize, figTitle)
%PLOT_FINAL_FIT Overlay all pair trajectories and the final fitted curve.
    isSSdU = strcmp(yField, "ratioSSdU");
    if isSSdU
        cFinal = H.coef_dU(:, end);
    else
        cFinal = H.coef_E(:,  end);
    end

    yhat = Xgrid * cFinal;

    figure('Name', figTitle); hold on
    for i = 1:numel(pairs)
        plot(pairs(i).z, pairs(i).(yField), '-', ...
            'Color', [baseColor, 0.18], 'LineWidth', 0.8);
    end
    plot(zGrid, yhat, '--', 'LineWidth', 2.4, 'Color', accentColor);
    yline(1, ':', 'Color', [0 0 0]);
    xlabel("Fidelity $z$");
    ylabel(sprintf("Normalised %s ratio", yLabelTex));
    xlim([0, 1]); ylim([-0.1, 1.3]);
    title(figTitle, "Interpreter", "none");
    grid on; box on
    if exist('set_font_size', 'file') == 2; set_font_size(fontSize); end
    if exist('set_fig_size', 'file') == 2; set_fig_size(900, 600); end
end