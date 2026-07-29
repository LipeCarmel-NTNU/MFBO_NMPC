%% Learning curve of the beta runtime surrogate: out-of-sample R2 against training size
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
% Regularisation and its selection
%   The fitted parameters are pulled towards a = b = 1, the identity f(z) = z:
%       cost_lambda(a, b) = cost(a, b) + lambda * ( (a - 1)^2 + (b - 1)^2 ).
%   lambda is chosen inside the training set alone, by kFold-fold cross-validation
%   over runs and not over samples, because the profiled per-run constant makes a
%   run the indivisible unit. The partition is redrawn nCV times and the held-out
%   score is the unregularised cost evaluated on the held-out runs, averaged over
%   folds and redraws. The lambda minimising that average is then used to fit the
%   whole training set. Every repetition of the outer ordering carries its own
%   training set and so its own selection, and the bottom row of the figure plots
%   the median of those selections against training size.
%
% Out-of-sample quantity
%   Each run in the pool gets one evaluation point, drawn once from a fixed seed
%   and held constant across the whole experiment. The high fidelity is
%       z ~ U(2 zMin, zEnd_i)
%   and the low one follows zTestRule. Both lie inside the run, so the ratio of
%   its own cumulative costs is measured and needs no J_i(1). The comparison is
%   made on the logarithm of that ratio, the same space the loss uses:
%       y    = log( J_i(z) / J_i(z_test) )   measured, treated as ground truth
%       yhat = log( f(z) / f(z_test) )       surrogate, from the fitted a and b
%   This is what a truncated run asks the surrogate for: carry a cost observed at
%   z_test forward to a higher fidelity. RMSE on y then reads directly as a
%   fractional error on the ratio. A run that stops below 2 zMin has no admissible
%   z and is left out of the evaluation, though its samples still enter the fit.
%
%   The logarithm matters for partial_SSdU. Its cumulative sum starts near zero
%   whenever the controller moves little early in the horizon, so J_i(z_test) can
%   be many orders of magnitude below J_i(z) and the raw ratio then reaches the
%   hundreds. A few such points dominate both sums in R2 and leave it near zero
%   whatever the fit does. minShare guards the same failure directly by requiring
%       J_i(z_test) >= minShare * J_i(z),
%   which caps the ratio at 1 / minShare, and the count of points dropped this way
%   is reported per target.
%
%   R2 = 1 - sum (y - yhat)^2 / sum (y - mean y)^2 over a set of runs. The
%   out-of-sample curve uses the held-out runs, the total curve uses every
%   evaluable run in the pool. Read it as the answer to one question: across a
%   pool of runs with different theta, does knowing the fidelity pair predict the
%   ratio? A converged fit together with R2 near zero says the run-to-run spread
%   in the shape of the cost accumulation is larger than the part that any single
%   shared f can explain, which is a statement about the shared-f premise and not
%   about the choice of f. std_y is printed alongside so the denominator is
%   visible rather than implied.
%
% Choice of zTestRule
%   "half" sets z_test = z / 2. The fidelity ratio z / z_test is then the same at
%   every evaluation point, and because f is close to a power law over the
%   observed range, y is close to the constant 2^a. The denominator of R2 is then
%   set by measurement noise instead of by any spread in the quantity being
%   predicted, and R2 stays near 0.5 whatever the fit does. On a synthetic pool of
%   121 runs generated from f = I_z(0.7, 1.4), inserting those exact parameters
%   gave R2 = 0.496 with y spanning 1.40 to 1.65.
%   "uniform" draws z_test ~ U(zMin, z), so the fidelity ratio varies from point
%   to point and y carries real spread. The same synthetic check gave R2 = 0.998
%   at the true parameters, with a lower RMSE than the "half" rule, 0.022 against
%   0.034. It is the default for that reason. The console table prints the RMSE
%   of y - yhat under either rule, which no normalisation touches.
%
% Design
%   The training set starts as the DOE, the first nDoe files of the run in
%   timestamp order, and grows by nAdd runs at a time drawn without replacement
%   from the rest of the pool. The ordering is redrawn nRep times and the plotted
%   curves are means over those repetitions, so the curve reflects training size
%   and not one arbitrary ordering. The evaluation points are fixed, so every
%   change along the x axis comes from which runs entered the fit. The last step
%   leaves minEval runs held out, which is few enough that the out-of-sample R2
%   there is noisy; that is a property of the estimator on a small sample.

clear; clc; close all;

%% Configuration
horizonHours  = 10.0;                 % T in z = t / T
runNames      = ["run1", "run2"];     % one experiment per BO run
zMin          = 0.1;                  % samples below this fidelity are dropped
nDoe          = 20;                   % first nDoe files of a run are the DOE
nAdd          = 50;                    % runs added to the fit per step
nRep          = 1;                   % random orderings averaged over
seed          = 0;                    % fixes the z draws and the orderings
minEval       = 6;                    % held-out runs left at the last step
zTestRule     = "uniform";            % "uniform" for z_test ~ U(zMin, z), or "half" for z / 2
minShare      = 0.02;                 % require J(z_test) >= minShare * J(z), capping the ratio at 50
kFold         = 10;                   % folds of the inner cross-validation, split over runs
nCV           = 1;                    % fold partitions redrawn per training set
lamGrid       = [0, logspace(-2, 2, 3)];  % L2 strengths tried, shrinking (a, b) towards (1, 1)

targetNames   = ["J_track", "J_TV"];  % panel order of the figure
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
% training size, so the units are weighted by that size and the time left is
% extrapolated from the weighted share done.
nRuns = numel(runNames);
nTgt  = numel(targetNames);

plan = repmat(struct('nPool', 0, 'sizes', [], 'nStep', 0), nRuns, 1);
for r = 1:nRuns
    files         = list_files(fullfile(projectRoot, "results", runNames(r)), "out_*.mat");
    plan(r).nPool = numel(files);
    maxFit        = plan(r).nPool - minEval;
    plan(r).sizes = unique([nDoe:nAdd:maxFit, maxFit]);
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

    for t = 1:nTgt
        pairs = load_pairs(poolDir, "out_*.mat", targetFields(t), ...
            targetOffsets(t), horizonHours);
        nPool = numel(pairs);

        fd = pool_samples(pairs, zMin);
        rng(seed);
        ev = eval_points(pairs, zMin, zTestRule, minShare);

        sizes = plan(r).sizes;
        nStep = plan(r).nStep;
        label = sprintf("%-6s %-8s", runNames(r), targetNames(t));

        r2Oos  = nan(nStep, nRep);
        r2Tot  = nan(nStep, nRep);
        rmse   = nan(nStep, nRep);
        aRep   = nan(nStep, nRep);
        bRep   = nan(nStep, nRep);
        lamRep = nan(nStep, nRep);

        for rep = 1:nRep
            rng(seed + rep);
            order = [1:nDoe, nDoe + randperm(nPool - nDoe)];
            for s = 1:nStep
                inFit = false(nPool, 1);
                inFit(order(1:sizes(s))) = true;

                sel  = select_runs(fd, inFit);
                lam  = select_lambda(sel, lamGrid, kFold, nCV, ab0, abLb, abUb, ...
                    seed + rep + 1000 * s);
                ab   = fit_beta(sel, ab0, abLb, abUb, lam);
                yHat = predict_log_ratio(ev, ab);

                held           = ev.usable & ~inFit;
                r2Oos(s, rep)  = r_squared(ev.y(held), yHat(held), minEval);
                r2Tot(s, rep)  = r_squared(ev.y(ev.usable), yHat(ev.usable), minEval);
                rmse(s, rep)   = rms_error(ev.y(held), yHat(held), minEval);
                aRep(s, rep)   = ab(1);
                bRep(s, rep)   = ab(2);
                lamRep(s, rep) = lam;

                doneUnits = doneUnits + 1;
                doneWork  = doneWork + sizes(s);
                isLast    = (rep == nRep) && (s == nStep);
                report_progress(doneUnits, totalUnits, doneWork, totalWork, ...
                    tStart, label, isLast);
            end
        end

        L(r, t).pct   = 100 * sizes(:) / nPool;
        L(r, t).nFit  = sizes(:);
        L(r, t).nHeld = nPool - sizes(:);
        L(r, t).r2Oos = mean(r2Oos, 2, "omitnan");
        L(r, t).r2Tot = mean(r2Tot, 2, "omitnan");
        L(r, t).rmse  = mean(rmse,  2, "omitnan");
        L(r, t).a     = mean(aRep,  2, "omitnan");
        L(r, t).b     = mean(bRep,  2, "omitnan");
        L(r, t).lam   = grid_median(lamRep);
        L(r, t).stdY  = std(ev.y(ev.usable));

        fprintf("\n%s %s: pool %d runs, %d usable in the fit, %d evaluable\n", ...
            runNames(r), targetNames(t), nPool, sum(fd.usable), sum(ev.usable));
        fprintf("  dropped %d for zEnd below %.2f, %d for J(z_test) below %g of J(z)\n", ...
            ev.nShort, 2 * zMin, ev.nThin, minShare);
        fprintf("  y = log ratio: spread std %.4f over %d points\n", ...
            L(r, t).stdY, sum(ev.usable));
        print_curve(L(r, t));
    end
end

%% Figure: one column per target, rows for R2, shape parameters and lambda
cOos = [0.84 0.37 0.00];
cTot = [0.00 0.45 0.70];
cA   = [0.20 0.55 0.85];
cB   = [0.80 0.40 0.60];
cLam = [0.35 0.35 0.35];
cZero = [0.85 0.10 0.10];

% A selected lambda of zero has no place on a logarithmic axis, so the steps that
% chose the unregularised fit are marked by a thin red vertical line instead and
% the curve skips them.
lamPos = lamGrid(lamGrid > 0);

for r = 1:nRuns
    figure('Name', sprintf('Learning curve, %s', runNames(r)));

    for t = 1:nTgt
        subplot(3, nTgt, t); hold on
        hT = plot(L(r, t).pct, L(r, t).r2Tot, '-o', 'Color', cTot, ...
            'MarkerFaceColor', cTot, 'MarkerSize', 4, 'LineWidth', 1.6);
        hO = plot(L(r, t).pct, L(r, t).r2Oos, '-s', 'Color', cOos, ...
            'MarkerFaceColor', cOos, 'MarkerSize', 4, 'LineWidth', 1.8);
        xlim([0 100]);
        xlabel("Share of the pool used in the fit [\%]");
        ylabel("$R^2$");
        title(tex_name(targetNames(t)));
        legend([hT, hO], ["total", "out of sample"], "Location", "southeast");
        grid on; box on

        subplot(3, nTgt, nTgt + t); hold on
        hA = plot(L(r, t).pct, L(r, t).a, '-', 'Color', cA, 'LineWidth', 1.8);
        hB = plot(L(r, t).pct, L(r, t).b, '-', 'Color', cB, 'LineWidth', 1.8);
        yline(1, ':', 'Color', [0.5 0.5 0.5]);
        xlim([0 100]);
        xlabel("Share of the pool used in the fit [\%]");
        ylabel("Shape parameter");
        legend([hA, hB], ["$a$", "$b$"], "Location", "best");
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

    sgtitle(sprintf("%s: out-of-sample fit of $f(z)/f(z_{\\mathrm{test}})$ " + ...
        "against training size, $z_{\\mathrm{test}}$ rule %s, " + ...
        "%d-fold CV repeated %d times", ...
        tex_name(runNames(r)), zTestRule, kFold, nCV));
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
    u    = sel.logJ - log(max(beta_ratio(sel.z, ab), 1e-12));
    uBar = accumarray(sel.g, u) ./ sel.cnt;
    loss = mean((u - uBar(sel.g)).^2);
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
        'MaxFunctionEvaluations', 2000, 'MaxIterations', 400, 'FiniteDifferenceType','central');
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

%% Local functions: evaluation

function ev = eval_points(pairs, zMin, rule, minShare)
%EVAL_POINTS One evaluation point per run, drawn once and then held fixed.
% z is uniform on [2 zMin, zEnd]. Under "uniform" the low fidelity is uniform on
% [zMin, z], under "half" it is z / 2. Both fidelities sit inside the run, so y is
% the logarithm of a measured ratio of its own cumulative costs.
% Two exclusions, each counted. A run stopping below 2 zMin has no admissible z.
% A point whose J(z_test) falls below minShare of J(z) is dropped, which bounds
% the ratio at 1 / minShare: the cumulative partial_SSdU starts near zero when the
% controller moves little early on, and without this guard a handful of such
% points set both sums in R2 on their own.
    n         = numel(pairs);
    ev.z      = nan(n, 1);
    ev.zTest  = nan(n, 1);
    ev.y      = nan(n, 1);
    ev.usable = false(n, 1);
    ev.nShort = 0;
    ev.nThin  = 0;
    for i = 1:n
        if pairs(i).zEnd < 2 * zMin
            ev.nShort = ev.nShort + 1;
            continue
        end
        ev.z(i) = 2 * zMin + (pairs(i).zEnd - 2 * zMin) * rand();
        switch rule
            case "uniform"
                ev.zTest(i) = zMin + (ev.z(i) - zMin) * rand();
            case "half"
                ev.zTest(i) = ev.z(i) / 2;
            otherwise
                error("zTestRule must be ""uniform"" or ""half"", got ""%s""", rule);
        end
        jHi = interp1(pairs(i).z, pairs(i).J, ev.z(i),     "linear");
        jLo = interp1(pairs(i).z, pairs(i).J, ev.zTest(i), "linear");
        if ~(isfinite(jHi) && isfinite(jLo) && jHi > 0 && jLo > 0)
            ev.nShort = ev.nShort + 1;
            continue
        end
        if jLo < minShare * jHi
            ev.nThin = ev.nThin + 1;
            continue
        end
        ev.y(i)      = log(jHi) - log(jLo);
        ev.usable(i) = true;
    end
end

function yHat = predict_log_ratio(ev, ab)
%PREDICT_LOG_RATIO Surrogate log f(z) - log f(z_test) at every usable point.
% Unusable runs keep NaN and are never passed to betainc, which rejects NaN.
    yHat = nan(numel(ev.z), 1);
    u    = ev.usable;
    yHat(u) = log(beta_ratio(ev.z(u), ab)) - log(beta_ratio(ev.zTest(u), ab));
end

function r2 = r_squared(y, yHat, minEval)
%R_SQUARED Coefficient of determination of yHat against y.
% NaN when fewer than minEval points remain, because the denominator is then set
% by a handful of draws.
    ok = isfinite(y) & isfinite(yHat);
    y  = y(ok); yHat = yHat(ok);
    if numel(y) < minEval
        r2 = NaN;
        return
    end
    ssTot = sum((y - mean(y)).^2);
    if ssTot <= 0
        r2 = NaN;
        return
    end
    r2 = 1 - sum((y - yHat).^2) / ssTot;
end

function e = rms_error(y, yHat, minEval)
%RMS_ERROR Root mean squared error of yHat against y, with no normalisation.
% y is a log ratio, so e reads directly as a fractional error on the ratio.
    ok = isfinite(y) & isfinite(yHat);
    y  = y(ok); yHat = yHat(ok);
    if numel(y) < minEval
        e = NaN;
        return
    end
    e = sqrt(mean((y - yHat).^2));
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

function m = grid_median(x)
%GRID_MEDIAN Lower median along the second dimension, kept on the sampled grid.
% lambda comes from a finite grid, so the ordinary median of an even number of
% repetitions would report a strength that was never tried, and a geometric mean
% would collapse to zero as soon as one repetition selected the unregularised
% fit. The lower of the two central values is returned, which is a value the
% cross-validation actually chose. NaN sorts last and so never becomes the
% median unless most repetitions failed.
    x = sort(x, 2);
    k = max(ceil(size(x, 2) / 2), 1);
    m = x(:, k);
end

function report_progress(done, total, doneWork, totalWork, tStart, label, isLast)
%REPORT_PROGRESS Overwrite one console line with the rate and the time left.
% One unit is a lambda selection by cross-validation followed by the fit on the
% whole training set, and s/unit is the mean over the units finished so far. The
% remaining time comes from the weighted work rather than from the unit count,
% because a unit at the last training size fits several times more samples than
% one at the first. The line is rewritten in place with backspaces and padded to
% the previous width so that a shorter message leaves nothing behind; isLast
% closes it with a newline before the block prints its table.
    fprintf('\n\n')
    persistent nBack
    if isempty(nBack) || done == 1
        nBack = 0;
    end

    elapsed = toc(tStart);
    perUnit = elapsed / max(done, 1);
    eta     = elapsed * (totalWork - doneWork) / max(doneWork, eps);

    msg = sprintf("  %s %s %5.1f%%  %5d/%d  %6.2f s/unit  elapsed %s  ETA %s", ...
        label, progress_bar(doneWork / totalWork, 24), ...
        100 * doneWork / totalWork, done, total, perUnit, ...
        clock_time(elapsed), clock_time(eta));
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

function s = tex_name(name)
%TEX_NAME Escape the underscores of a name for the latex interpreter.
    s = strrep(name, "_", "\_");
end

function print_curve(c)
%PRINT_CURVE Print the learning curve as a table.
% held is the number of runs left out of the fit, which the out-of-sample R2 is
% computed over and which explains the noise in the last rows. lambda is the
% median over the repetitions of the value selected by cross-validation, and a
% zero there means most repetitions kept the unregularised fit.
    fprintf("    %5s %6s %6s %9s %9s %9s %8s %8s %8s %10s\n", ...
        "pct", "nFit", "held", "R2_total", "R2_oos", "RMSE_oos", "std_y", ...
        "a", "b", "lambda");
    for k = 1:numel(c.pct)
        fprintf("    %5.1f %6d %6d %9.4f %9.4f %9.4f %8.4f %8.4f %8.4f %10.2e\n", ...
            c.pct(k), c.nFit(k), c.nHeld(k), c.r2Tot(k), c.r2Oos(k), ...
            c.rmse(k), c.stdY, c.a(k), c.b(k), c.lam(k));
    end
end
