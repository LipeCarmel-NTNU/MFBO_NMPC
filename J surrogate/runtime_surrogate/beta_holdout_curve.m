%% Beta runtime surrogate tested on complete runs, against training size
%
% Model
%   J(z) = J(1) * f(z), with f the regularised incomplete beta function
%       f(z) = I_z(a, b),   a > 0,   b > 0,
%   fitted by profile least squares on log J. For each training run i the total
%   J_i(1) is unknown, so log J_i(z) = c_i + log f(z) with c_i free. The
%   least-squares c_i is the mean of u_i(z) = log J_i(z) - log f(z) over that
%   run's samples, and substituting it leaves
%       cost(a, b) = mean_iz ( u_i(z) - mean_z u_i(z) )^2.
%   The mean rather than the sum keeps the data term independent of the number of
%   samples, so one L2 strength is comparable across training sizes and folds.
%
% Training and test sets
%   The training pool is the truncated runs of one BO campaign, in timestamp
%   order. It starts as the DOE, the first nDoe files, and grows by nAdd runs at
%   a time drawn without replacement from the rest. The ordering is redrawn nRep
%   times and the plotted curves are the mean over those repetitions with a band
%   spanning them, so a curve reflects training size and not one arbitrary
%   ordering. At the last size the training set is the whole pool and every
%   ordering coincides, which is why the band closes there.
%
%   The test set is a separate folder of runs carried to the full horizon, so
%   z = 1 is present and J_i(1) is measured. Nothing about it changes as the
%   training set grows. Two consequences follow. The shape is observed directly,
%   f_i(z) = J_i(z) / J_i(1), with no constant left to profile out, so the fitted
%   f is compared against a measurement rather than against a ratio of two
%   fidelities. And the extrapolation the controller performs is reproduced
%   exactly: truncate at zTrunc, carry the cost observed there to the finish,
%       log Jhat_i(1) = log J_i(zTrunc) - log f(zTrunc),
%   and compare against the measured log J_i(1).
%
% Reported quantities
%   R2_fit  on the training runs used in the fit. Numerator and denominator are
%           both centred per run, so it is the share of the within-run variance
%           of log J that the shared shape explains. It is the only one of the
%           three that a truncated run can support.
%   R2_f    on the test samples, f_i(z) against I_z(a, b), a direct check of the
%           shape over the whole horizon.
%   R2_J    on the test runs, log Jhat_i(1) against log J_i(1). This is the
%           deployment quantity. Its spread comes from the run-to-run spread in
%           the total cost, which is wide, so read it together with R2_f rather
%           than on its own.
%
%   The log ratio log( J_i(1) / J_i(zTrunc) ) is not scored. With zTrunc fixed
%   the surrogate predicts -log f(zTrunc) for every run, a constant, and R2
%   against a constant prediction cannot exceed zero whatever the fit does.
%
% Regularisation and its selection
%   The fitted parameters are pulled towards a = b = 1, the identity f(z) = z:
%       cost_lambda(a, b) = cost(a, b) + lambda * ( (a - 1)^2 + (b - 1)^2 ).
%   lambda is chosen inside the training set alone, by kFold-fold cross-validation
%   over runs and not over samples, because the profiled per-run constant makes a
%   run the indivisible unit. The partition is redrawn nCV times and the held-out
%   score is the unregularised cost on the held-out runs, which is the only score
%   a truncated run admits. The lambda minimising the average is then used to fit
%   the whole training set.

clear; clc; close all;

%% Configuration
horizonHours  = 10.0;                 % T in z = t / T
runNames      = ["run1", "run2"];     % one BO campaign per block
testNames     = ["final_fidelity_same_noise/run1_full_f1_same_noise", ...
                 "final_fidelity_same_noise/run2_full_f1_same_noise"];
zMin          = 0.1;                  % samples below this fidelity are dropped
zTrunc        = 0.75;                 % truncation the extrapolation starts from
nDoe          = 20;                   % first nDoe files of a run are the DOE
nAdd          = 20;                   % runs added to the fit per step
nRep          = 10;                   % random orderings averaged over
seed          = 0;                    % fixes the orderings and the CV partitions

kFold         = 10;                   % folds of the inner cross-validation
nCV           = 5;                    % fold partitions redrawn per training set
lamGrid       = [0, logspace(-3, 2, 4)];  % L2 strengths tried

targetNames   = ["J_track", "J_TV"];  % column order of the figures
targetFields  = ["partial_SSE", "partial_SSdU"];
targetOffsets = [0, 1];               % index in out.T of the first partial entry

ab0  = [1; 1];                        % warm start, f(z) = z
abLb = [0.05; 0.05];
abUb = [25; 25];

%% Paths
scriptDir   = fileparts(mfilename("fullpath"));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(fullfile(projectRoot, "dependencies")));

set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");

%% Work plan, so the progress report knows the total before the first fit
% The file count alone fixes the training sizes, and reading the directory is
% cheap, so the schedule of every block is known before any file is loaded. One
% unit of work is one (repetition, training size) pair. Its cost grows with the
% training size, so the units are weighted by that size.
nRuns = numel(runNames);
nTgt  = numel(targetNames);

plan = repmat(struct('nPool', 0, 'sizes', [], 'nStep', 0), nRuns, 1);
for r = 1:nRuns
    files         = list_files(fullfile(projectRoot, "results", runNames(r)), "out_*.mat");
    plan(r).nPool = numel(files);
    plan(r).sizes = unique([nDoe:nAdd:plan(r).nPool, plan(r).nPool]);
    plan(r).nStep = numel(plan(r).sizes);
end

totalUnits = nTgt * nRep * sum([plan.nStep]);
totalWork  = nTgt * nRep * sum(arrayfun(@(p) sum(p.sizes), plan));
doneUnits  = 0;
doneWork   = 0;
tStart     = tic;

fprintf("Fitting %d units over %d blocks, %d-fold CV repeated %d times on %d " + ...
    "values of lambda\n", totalUnits, nRuns * nTgt, kFold, nCV, numel(lamGrid));

%% Learning curve for each BO run and target
L = struct([]);

for r = 1:nRuns
    poolDir = fullfile(projectRoot, "results", runNames(r));
    testDir = fullfile(projectRoot, "results", testNames(r));

    for t = 1:nTgt
        trainPairs = load_pairs(poolDir, "out_*.mat", targetFields(t), ...
            targetOffsets(t), horizonHours);
        testPairs  = load_pairs(testDir, "out_*.mat", targetFields(t), ...
            targetOffsets(t), horizonHours);

        nPool = numel(trainPairs);
        fd    = pool_samples(trainPairs, zMin);
        te    = test_set(testPairs, zMin, zTrunc);

        sizes = plan(r).sizes;
        nStep = plan(r).nStep;
        label = sprintf("%-5s %-8s", runNames(r), targetNames(t));

        r2Fit  = nan(nStep, nRep);
        r2F    = nan(nStep, nRep);
        r2J    = nan(nStep, nRep);
        lamRep = nan(nStep, nRep);

        for rep = 1:nRep
            rng(seed + rep);
            order = [1:nDoe, nDoe + randperm(nPool - nDoe)];
            for s = 1:nStep
                inFit = false(nPool, 1);
                inFit(order(1:sizes(s))) = true;

                sel = select_runs(fd, inFit);
                lam = select_lambda(sel, lamGrid, kFold, nCV, ab0, abLb, abUb, ...
                    seed + rep + 1000 * s);
                ab  = fit_beta(sel, ab0, abLb, abUb, lam);

                r2Fit(s, rep)  = r2_profiled(sel, ab);
                r2F(s, rep)    = r2_shape(te, ab);
                r2J(s, rep)    = r2_total_cost(te, ab, zTrunc);
                lamRep(s, rep) = lam;

                doneUnits = doneUnits + 1;
                doneWork  = doneWork + sizes(s);
                report_progress(doneUnits, totalUnits, doneWork, totalWork, ...
                    tStart, label, ...
                    [r2Fit(s, rep), r2F(s, rep), r2J(s, rep)], ...
                    (rep == nRep) && (s == nStep));
            end
        end

        % The last step trains on the whole pool, so its fit does not depend on
        % the ordering. It is repeated once here to hold the parity plots.
        selAll = select_runs(fd, true(nPool, 1));
        lamAll = select_lambda(selAll, lamGrid, kFold, nCV, ab0, abLb, abUb, seed);
        abAll  = fit_beta(selAll, ab0, abLb, abUb, lamAll);

        L(r, t).pct    = 100 * sizes(:) / nPool;
        L(r, t).nFit   = sizes(:);
        L(r, t).r2Fit  = band(r2Fit);
        L(r, t).r2F    = band(r2F);
        L(r, t).r2J    = band(r2J);
        L(r, t).lam    = grid_median(lamRep);
        L(r, t).ab     = abAll;
        L(r, t).lamAll = lamAll;
        L(r, t).te     = te;

        fprintf("\n%s %s: pool %d runs, %d usable in the fit, test %d complete runs\n", ...
            runNames(r), targetNames(t), nPool, sum(fd.usable), te.nRun);
        fprintf("  full pool fit: a %.4f, b %.4f, lambda %g\n", ...
            abAll(1), abAll(2), lamAll);
        fprintf("  full pool R2: fit %.4f, f(z) %.4f, log J(1) %.4f\n", ...
            r2_profiled(selAll, abAll), r2_shape(te, abAll), ...
            r2_total_cost(te, abAll, zTrunc));
        print_curve(L(r, t));
    end
end

%% Figure: learning curves, one column per target
cFit  = [0.45 0.45 0.45];
cF    = [0.00 0.45 0.70];
cJ    = [0.84 0.37 0.00];
cLam  = [0.35 0.35 0.35];
cZero = [0.85 0.10 0.10];

% A selected lambda of zero has no place on a logarithmic axis, so the steps that
% chose the unregularised fit are marked by a thin red vertical line instead and
% the curve skips them.
lamPos = lamGrid(lamGrid > 0);

for r = 1:nRuns
    figure('Name', sprintf('Held-out learning curve, %s', runNames(r)));

    for t = 1:nTgt
        subplot(3, nTgt, t); hold on
        hI = plot_band(L(r, t).pct, L(r, t).r2Fit, cFit);
        hF = plot_band(L(r, t).pct, L(r, t).r2F, cF);
        yline(0, ':', 'Color', [0.5 0.5 0.5]);
        xlim([0 100]);
        xlabel("Share of the pool used in the fit [\%]");
        ylabel("$R^2$");
        title(tex_name(targetNames(t)));
        legend([hI, hF], ["in sample, $\log J$ shape", "test, $f(z)$"], ...
            "Location", "southeast");
        grid on; box on

        subplot(3, nTgt, nTgt + t); hold on
        plot_band(L(r, t).pct, L(r, t).r2J, cJ);
        yline(0, ':', 'Color', [0.5 0.5 0.5]);
        xlim([0 100]);
        xlabel("Share of the pool used in the fit [\%]");
        ylabel(sprintf("$R^2$, $\\log J(1)$ from $z = %.2f$", zTrunc));
        grid on; box on

        subplot(3, nTgt, 2 * nTgt + t); hold on
        isZero = L(r, t).lam <= 0;
        for k = find(isZero(:)).'
            xline(L(r, t).pct(k), '-', 'Color', cZero, 'LineWidth', 0.5);
        end
        plot(L(r, t).pct(~isZero), L(r, t).lam(~isZero), '-d', 'Color', cLam, ...
            'MarkerFaceColor', cLam, 'MarkerSize', 4, 'LineWidth', 1.6);
        set(gca, 'YScale', 'log');
        xlim([0 100]); ylim([min(lamPos) / 2, max(lamPos) * 2]);
        xlabel("Share of the pool used in the fit [\%]");
        ylabel("$\lambda$ chosen by CV");
        grid on; box on
    end

    sgtitle(sprintf("%s against %d complete runs, %d-fold CV repeated %d times", ...
        tex_name(runNames(r)), L(r, 1).te.nRun, kFold, nCV));
end

%% Figure: parity of the full pool fit
for r = 1:nRuns
    figure('Name', sprintf('Parity, %s', runNames(r)));

    subplot(1, 3, 1); hold on
    hs  = gobjects(nTgt, 1);
    lbl = strings(nTgt, 1);
    col = [cF; cJ];
    for t = 1:nTgt
        te   = L(r, t).te;
        fHat = beta_ratio(te.z, L(r, t).ab);
        hs(t) = scatter(fHat, te.f, 6, col(t, :), 'filled', ...
            'MarkerFaceAlpha', 0.15);
        lbl(t) = sprintf("%s, $R^2 = %.3f$", tex_name(targetNames(t)), ...
            r2_shape(te, L(r, t).ab));
    end
    plot([0 1], [0 1], 'k--');
    axis square; xlim([0 1]); ylim([0 1]);
    xlabel("$I_z(a, b)$ predicted");
    ylabel("$J(z) / J(1)$ measured");
    title("Shape $f(z)$");
    legend(hs, lbl, "Location", "southeast");
    grid on; box on

    for t = 1:nTgt
        subplot(1, 3, 1 + t); hold on
        te    = L(r, t).te;
        logJ  = log(te.jFull);
        logHat = log(te.jTrunc) - log(beta_ratio(zTrunc, L(r, t).ab));
        lo    = min([logJ; logHat]);
        hi    = max([logJ; logHat]);
        pad   = 0.05 * (hi - lo) + eps;
        plot([lo - pad, hi + pad], [lo - pad, hi + pad], 'k--');
        scatter(logHat, logJ, 36, col(t, :), 'filled');
        axis square; xlim([lo - pad, hi + pad]); ylim([lo - pad, hi + pad]);
        xlabel(sprintf("$\\log \\hat{J}(1)$ from $z = %.2f$", zTrunc));
        ylabel("$\log J(1)$ measured");
        title(sprintf("%s, $R^2 = %.3f$", tex_name(targetNames(t)), ...
            r2_total_cost(te, L(r, t).ab, zTrunc)));
        grid on; box on
    end

    sgtitle(sprintf("%s: full pool fit against the complete runs", ...
        tex_name(runNames(r))));
end

%% Figure: measured shapes of the complete runs against the fitted f
% Every complete run carries its own f_i(z) = J_i(z) / J_i(1), and the fit offers
% one curve for all of them. The spread between the grey lines is the part of the
% test error that no choice of a and b can remove, so this panel says what R2_f
% is measuring. All the curves meet at z = 1 by construction.
zPlot = linspace(0, 1, 400).';

for r = 1:nRuns
    figure('Name', sprintf('Shape overlay, %s', runNames(r)));

    for t = 1:nTgt
        subplot(1, nTgt, t); hold on
        te  = L(r, t).te;
        ids = unique(te.idx);
        hM  = gobjects(1);
        for k = 1:numel(ids)
            m  = te.idx == ids(k);
            hM = plot(te.z(m), te.f(m), '-', 'Color', [0.65 0.65 0.65], ...
                'LineWidth', 0.8);
        end
        hP = plot(zPlot, beta_ratio(zPlot, L(r, t).ab), '-', 'Color', cF, ...
            'LineWidth', 2.2);
        hT = xline(zTrunc, '--', 'Color', cJ, 'LineWidth', 1.2);
        xlim([0 1]); ylim([0 1]);
        xlabel("$z$");
        ylabel("$J(z) / J(1)$");
        title(sprintf("%s, $a = %.3f$, $b = %.3f$, $R^2 = %.3f$", ...
            tex_name(targetNames(t)), L(r, t).ab(1), L(r, t).ab(2), ...
            r2_shape(te, L(r, t).ab)));
        legend([hM, hP, hT], [sprintf("measured, %d runs", te.nRun), ...
            "fitted $I_z(a, b)$", sprintf("$z = %.2f$", zTrunc)], ...
            "Location", "southeast");
        grid on; box on
    end

    sgtitle(sprintf("%s: shape of the complete runs, full pool fit", ...
        tex_name(runNames(r))));
end

%% Local functions: model and fit

function f = beta_ratio(z, ab)
%BETA_RATIO f(z) = I_z(a, b), the regularised incomplete beta function.
% Increasing in z with f(0) = 0 and f(1) = 1 for every positive a and b.
% Every argument is clipped to the domain betainc accepts before the call. The
% bounds passed to fmincon hold at the solution but not at every trial point a
% line search proposes, and betainc raises an error outside its domain instead of
% returning a value the search could reject. The clip limits on (a, b) are wide
% enough to contain any parameter box the script would use, so they act only on
% infeasible trial points. The two-argument min and max omit NaN, so a non-finite
% entry collapses onto the nearest limit and the call still returns.
    abMin = 1e-6;
    abMax = 1e4;
    a = min(max(ab(1), abMin), abMax);
    b = min(max(ab(2), abMin), abMax);
    f = betainc(min(max(z(:), 0), 1), a, b);
end

function fd = pool_samples(pairs, zMin)
%POOL_SAMPLES Flat arrays of every usable sample of the pool, tagged by run.
% One flat vector lets the cost call betainc once per evaluation instead of once
% per run. A run needs two samples to say anything after its own constant is
% profiled out, so shorter ones are marked unusable.
    n      = numel(pairs);
    z      = cell(n, 1);
    logJ   = cell(n, 1);
    idx    = cell(n, 1);
    usable = false(n, 1);
    for i = 1:n
        keep = pairs(i).z >= zMin & pairs(i).J > 0;
        if sum(keep) < 2
            continue
        end
        usable(i) = true;
        z{i}      = pairs(i).z(keep);
        logJ{i}   = log(pairs(i).J(keep));
        idx{i}    = repmat(i, sum(keep), 1);
    end
    fd.z      = vertcat(z{:});
    fd.logJ   = vertcat(logJ{:});
    fd.idx    = vertcat(idx{:});
    fd.usable = usable;
end

function sel = select_runs(fd, inFit)
%SELECT_RUNS Samples of the selected runs, with a compact per-run group index.
% The group index and the group sizes are built once per training set so that the
% cost does not repeat them at every trial (a, b).
    inFit = inFit(:) & fd.usable;
    take  = inFit(fd.idx);
    sel.z    = fd.z(take);
    sel.logJ = fd.logJ(take);
    [~, ~, sel.g] = unique(fd.idx(take));
    sel.cnt  = accumarray(sel.g, 1);
    sel.nG   = numel(sel.cnt);
end

function sub = subset_groups(sel, keep)
%SUBSET_GROUPS Samples of the selected groups, with the group index rebuilt.
% keep is a logical mask over the compact group index of sel. Splitting by group
% and never by sample is what the per-run constant requires: a run either informs
% the fit or is scored against it.
    keep = keep(:);
    take = keep(sel.g);
    sub.z    = sel.z(take);
    sub.logJ = sel.logJ(take);
    [~, ~, sub.g] = unique(sel.g(take));
    sub.cnt  = accumarray(sub.g, 1);
    sub.nG   = numel(sub.cnt);
end

function loss = beta_loss(ab, sel)
%BETA_LOSS Mean squared deviation of u = log J - log f from its per-run mean.
% Centring on the per-run mean is what profiling out the unknown log J_i(1)
% amounts to, so the loss holds no run total. Averaging over samples instead of
% summing makes the value comparable between training sets of different size,
% which one shared grid of L2 strengths needs. The floor of 1e-12 keeps the
% logarithm defined at a trial point where f underflows.
    loss = mean(profile_residual(ab, sel).^2);
end

function res = profile_residual(ab, sel)
%PROFILE_RESIDUAL u = log J - log f, centred on its per-run mean.
    u    = sel.logJ - log(max(beta_ratio(sel.z, ab), 1e-12));
    uBar = accumarray(sel.g, u) ./ sel.cnt;
    res  = u - uBar(sel.g);
end

function cost = beta_cost(ab, sel, lam)
%BETA_COST Data loss plus an L2 pull of (a, b) towards (1, 1).
% The penalty centre is the identity f(z) = z, so shrinkage moves the surrogate
% towards a cost that accumulates at a constant rate rather than towards zero.
    cost = beta_loss(ab, sel) + lam * sum((ab(:) - 1).^2);
end

function ab = fit_beta(sel, ab0, lb, ub, lam)
%FIT_BETA Minimise beta_cost over the two shape parameters.
% Only bounds constrain the problem. The start point a = b = 1 is f(z) = z.
    opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
        'MaxFunctionEvaluations', 2000, 'MaxIterations', 400);
    ab = fmincon(@(p) beta_cost(p, sel, lam), ab0, [], [], [], [], lb, ub, [], opts);
end

function [lamBest, cvLoss] = select_lambda(sel, lamGrid, kFold, nCV, ab0, lb, ub, seed)
%SELECT_LAMBDA Pick the L2 strength by repeated k-fold cross-validation over runs.
% The runs of the training set are split into k folds, each fold in turn is scored
% by the unregularised beta_loss under parameters fitted on the other folds, and
% the partition is redrawn nCV times. The grid is walked from weak to strong
% shrinkage with the previous solution as the start point, which keeps the number
% of fmincon iterations small along the path. The returned lambda is the minimiser
% of the averaged held-out loss.
    nLam    = numel(lamGrid);
    k       = min(kFold, sel.nG);
    cvLoss  = zeros(nLam, 1);
    nScored = 0;

    rng(seed);
    for c = 1:nCV
        fold = mod(randperm(sel.nG), k).' + 1;
        for f = 1:k
            held  = (fold == f);
            train = ~held;
            if ~any(held) || ~any(train)
                continue
            end
            sTr = subset_groups(sel, train);
            sHe = subset_groups(sel, held);
            ab  = ab0;
            for iL = 1:nLam
                ab = fit_beta(sTr, ab, lb, ub, lamGrid(iL));
                cvLoss(iL) = cvLoss(iL) + beta_loss(ab, sHe);
            end
            nScored = nScored + 1;
        end
    end

    cvLoss = cvLoss / max(nScored, 1);
    [~, iBest] = min(cvLoss);
    lamBest = lamGrid(iBest);
end

%% Local functions: test set and scores

function te = test_set(pairs, zMin, zTrunc)
%TEST_SET Complete runs, with the shape measured and the truncation point read.
% Every run of this folder reaches z = 1, so J_i(1) is measured and the shape
% f_i(z) = J_i(z) / J_i(1) needs no constant profiled out. jTrunc is the cost
% observed at the truncation the controller would apply.
    n     = numel(pairs);
    z     = cell(n, 1);
    f     = cell(n, 1);
    idx   = cell(n, 1);
    jFull = nan(n, 1);
    jTrunc = nan(n, 1);
    keepRun = false(n, 1);

    for i = 1:n
        jEnd = pairs(i).J(end);
        jLo  = interp1(pairs(i).z, pairs(i).J, zTrunc, "linear");
        keep = pairs(i).z >= zMin & pairs(i).J > 0;
        if ~(isfinite(jEnd) && jEnd > 0 && isfinite(jLo) && jLo > 0) || sum(keep) < 2
            continue
        end
        if abs(pairs(i).zEnd - 1) > 1e-6
            warning("test run %d ends at z = %.4f, not 1", i, pairs(i).zEnd);
        end
        keepRun(i) = true;
        z{i}       = pairs(i).z(keep);
        f{i}       = pairs(i).J(keep) / jEnd;
        idx{i}     = repmat(i, sum(keep), 1);
        jFull(i)   = jEnd;
        jTrunc(i)  = jLo;
    end

    te.z      = vertcat(z{:});
    te.f      = vertcat(f{:});
    te.idx    = vertcat(idx{:});
    te.jFull  = jFull(keepRun);
    te.jTrunc = jTrunc(keepRun);
    te.nRun   = sum(keepRun);
end

function r2 = r2_profiled(sel, ab)
%R2_PROFILED Share of the within-run variance of log J explained by the shape.
% Both terms are centred on the per-run mean, so the unknown run total cancels
% from each and the comparison is against a flat log J, which is the model that
% carries no shape at all.
    res   = profile_residual(ab, sel);
    logJ  = sel.logJ;
    jBar  = accumarray(sel.g, logJ) ./ sel.cnt;
    ssTot = sum((logJ - jBar(sel.g)).^2);
    r2    = ratio_r2(sum(res.^2), ssTot);
end

function r2 = r2_shape(te, ab)
%R2_SHAPE Fitted f against the measured J(z) / J(1) over the test samples.
    fHat = beta_ratio(te.z, ab);
    r2   = r_squared(te.f, fHat);
end

function r2 = r2_total_cost(te, ab, zTrunc)
%R2_TOTAL_COST Extrapolation of the truncated cost to the finish, on a log scale.
% log Jhat(1) = log J(zTrunc) - log f(zTrunc), the operation the controller
% performs when it stops a run early and needs its total.
    logHat = log(te.jTrunc) - log(beta_ratio(zTrunc, ab));
    r2     = r_squared(log(te.jFull), logHat);
end

function r2 = r_squared(y, yHat)
%R_SQUARED Coefficient of determination of yHat against y.
    ok = isfinite(y) & isfinite(yHat);
    y  = y(ok); yHat = yHat(ok);
    if numel(y) < 3
        r2 = NaN;
        return
    end
    r2 = ratio_r2(sum((y - yHat).^2), sum((y - mean(y)).^2));
end

function r2 = ratio_r2(ssRes, ssTot)
%RATIO_R2 1 - ssRes / ssTot, with a guard on a degenerate denominator.
    if ~(ssTot > 0)
        r2 = NaN;
    else
        r2 = 1 - ssRes / ssTot;
    end
end

function s = band(x)
%BAND Mean and extremes over the repetitions, as one struct per curve.
    s.mid = mean(x, 2, "omitnan");
    s.lo  = min(x, [], 2, "omitnan");
    s.hi  = max(x, [], 2, "omitnan");
end

function m = grid_median(x)
%GRID_MEDIAN Lower median along the second dimension, kept on the sampled grid.
% lambda comes from a finite grid, so the ordinary median of an even number of
% repetitions would report a strength that was never tried, and a geometric mean
% would collapse to zero as soon as one repetition selected the unregularised
% fit. The lower of the two central values is returned, which is a value the
% cross-validation actually chose.
    x = sort(x, 2);
    k = max(ceil(size(x, 2) / 2), 1);
    m = x(:, k);
end

%% Local functions: data

function files = list_files(folder, pattern)
%LIST_FILES Return the files of one folder that match the pattern, sorted.
% The names carry a timestamp, so the sorted order is the order of the runs and
% the first nDoe of them are the DOE.
    files = dir(fullfile(folder, pattern));
    if isempty(files)
        error("No %s files in %s", pattern, folder);
    end
    [~, order] = sort({files.name});
    files = files(order);
end

function pair = load_one_pair(filePath, partialField, offset, horizonHours)
%LOAD_ONE_PAIR Read one file into a cumulative cost curve.
% The two cases of the file are added together. offset is the index in out.T
% that holds the first entry of the partial array.
    data = load(filePath, "out");
    out  = data.out;

    time = double(out.T(:));
    p1   = double(out.case(1).(partialField));
    p2   = double(out.case(2).(partialField));
    p1   = p1(:);
    p2   = p2(:);
    n    = min([numel(p1), numel(p2), numel(time) - offset]);

    pair.z    = time(offset + (1:n)) / horizonHours;
    pair.J    = cumsum(p1(1:n)) + cumsum(p2(1:n));
    pair.zEnd = pair.z(end);
end

function pairs = load_pairs(folder, pattern, partialField, offset, horizonHours)
%LOAD_PAIRS Read every matching file of one folder into an array of pairs.
    files = list_files(folder, pattern);
    pairs = struct('z', {}, 'J', {}, 'zEnd', {});
    for i = 1:numel(files)
        pairs(i) = load_one_pair(fullfile(files(i).folder, files(i).name), ...
            partialField, offset, horizonHours);
    end
end

%% Local functions: report

function h = plot_band(x, s, colour)
%PLOT_BAND Mean curve with the spread over the repetitions shaded behind it.
% The band closes at the last training size because every ordering then holds the
% whole pool and the fits coincide.
    fill([x; flipud(x)], [s.lo; flipud(s.hi)], colour, ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none');
    h = plot(x, s.mid, '-o', 'Color', colour, 'MarkerFaceColor', colour, ...
        'MarkerSize', 4, 'LineWidth', 1.6);
end

function s = tex_name(name)
%TEX_NAME Escape the underscores of a name for the latex interpreter.
    s = strrep(name, "_", "\_");
end

function print_curve(c)
%PRINT_CURVE Print the learning curve as a table, means over the repetitions.
    fprintf("    %5s %6s %9s %9s %9s %10s\n", ...
        "pct", "nFit", "R2_fit", "R2_f", "R2_J", "lambda");
    for k = 1:numel(c.pct)
        fprintf("    %5.1f %6d %9.4f %9.4f %9.4f %10.2e\n", ...
            c.pct(k), c.nFit(k), c.r2Fit.mid(k), c.r2F.mid(k), c.r2J.mid(k), ...
            c.lam(k));
    end
end

function report_progress(done, total, doneWork, totalWork, tStart, label, r2, isLast)
%REPORT_PROGRESS Overwrite one console line with the scores and the time left.
% One unit is a lambda selection by cross-validation followed by the fit on the
% whole training set, and s/unit is the mean over the units finished so far. The
% remaining time comes from the weighted work rather than from the unit count,
% because a unit at the last training size fits several times more samples than
% one at the first. The line is rewritten in place with backspaces and padded to
% the previous width so that a shorter message leaves nothing behind; isLast
% closes it with a newline before the block prints its table.
    persistent nBack
    if isempty(nBack) || done == 1
        nBack = 0;
    end

    elapsed = toc(tStart);
    perUnit = elapsed / max(done, 1);

    msg = sprintf("  %s %s %5.1f%%  R2 fit %6.3f  f %6.3f  J %6.3f  %5.2f s/unit  ETA %s", ...
        label, progress_bar(doneWork / totalWork, 16), ...
        100 * doneWork / totalWork, r2(1), r2(2), r2(3), perUnit, ...
        clock_time(elapsed * (totalWork - doneWork) / max(doneWork, eps)));
    if numel(msg) < nBack
        msg = [msg, blanks(nBack - numel(msg))];
    end

    fprintf("%s%s", repmat(char(8), 1, nBack), msg);
    nBack = numel(msg);
    if isLast
        fprintf("\n");
        nBack = 0;
    end
end

function bar = progress_bar(frac, width)
%PROGRESS_BAR Fixed-width bar of the share done.
    n   = round(min(max(frac, 0), 1) * width);
    bar = "[" + string(repmat('=', 1, n)) + string(repmat('.', 1, width - n)) + "]";
end

function s = clock_time(sec)
%CLOCK_TIME Seconds as mm:ss, or h:mm:ss past an hour.
    if ~isfinite(sec)
        s = "--:--";
        return
    end
    sec = max(sec, 0);
    h   = floor(sec / 3600);
    m   = floor(mod(sec, 3600) / 60);
    sc  = floor(mod(sec, 60));
    if h > 0
        s = sprintf("%d:%02d:%02d", h, m, sc);
    else
        s = sprintf("%02d:%02d", m, sc);
    end
end
