%% Fit alternative runtime surrogates f(z) on the DOE, test them on the full-horizon reruns
%
% Model
%   J(z) = J(1) * f(z). The total J(1) is unknown for a truncated run, so only
%   ratios of the same run are observable:
%       J(zHi) / J(zLo) = f(zHi) / f(zLo).
%   Every candidate f below therefore has to satisfy f(1) = 1 to fix the scale,
%   and f(0) = 0 so that the ratio is well posed. The fit reads out.T and the two
%   partial_* arrays only. It never reads out.SSE or out.SSdU, because those two
%   fields hold partial / frac(f) from the curve hard-coded in main_BO.m.
%
% Candidates
%   cheb_deltalog  Chebyshev polynomial of order chebOrder in x = 2z - 1, fitted
%                  on the residual of the log ratio. This reproduces
%                  deltalog_surrogate.m and is the reference row of every table.
%   cheb_ratio     The same polynomial fitted on the relative residual of the raw
%                  ratio, with no logarithm anywhere in the objective. The
%                  residual stays scale free, so pairs at high and low z remain
%                  comparable and the only change from the reference is the
%                  transform itself.
%   expsat         f(z) = expm1(-lambda z) / expm1(-lambda), one parameter,
%                  warm started at lambda = 0.3. The normalisation gives
%                  f(0) = 0 and f(1) = 1 for every lambda, and f is increasing
%                  for every real lambda: lambda > 0 is concave, lambda < 0 is
%                  convex, and lambda -> 0 recovers f(z) = z.
%   monotone       f(z) = int_0^z exp(p(s)) ds / int_0^1 exp(p(s)) ds with p a
%                  Chebyshev polynomial of order monoOrder. The integrand is
%                  positive, so f is strictly increasing, and the normalisation
%                  makes the endpoints exact. Both properties hold for every
%                  parameter value, so the running maximum and the cap of
%                  deploy_model never bind and the shipped surrogate equals the
%                  raw fit. A constant added to p cancels in the ratio of the two
%                  integrals, so the constant term is held at zero and only
%                  monoOrder coefficients are free.
%   beta_cdf       f(z) = I_z(a, b), the regularised incomplete beta function.
%                  Two parameters, endpoints exact and increasing by
%                  construction, and the shape covers concave, convex and
%                  sigmoid without the free ends of a polynomial. a = b = 1 is
%                  f(z) = z and is the warm start.
%
% Observations and weights
%   Every DOE pair gives the grid zMin:zStep:zEnd, and every ordered combination
%   of two grid points of the same pair gives one observation, so a grid of n
%   points gives n (n - 1) / 2 of them. The DOE draws the fidelity from a Sobol
%   sample of [0, 1], so a short run reaches only the low grid points and the low
%   combinations appear in many more pairs than the high ones. One bin is one
%   combination (zLo, zHi), and each observation carries the weight 1 / count of
%   its own bin, so every combination has the same total weight.
%
% Test
%   The files in results/final_fidelity_same_noise are the frontier controllers
%   rerun with theta(1) = 1. Those runs cover the whole horizon, so J(end) is the
%   measured total and f_measured(z) = J(z) / J(1) is exact at every z. Nothing
%   is extrapolated on the test side. The score is the error of the total that
%   the surrogate would report from a run truncated at z:
%       Jhat(1) = J(z) / deploy_model(z),   error = Jhat(1) / J(1) - 1.
%
% Output
%   One test table per run, target and candidate, then one summary table per run
%   and target holding the median absolute error of every candidate side by side,
%   then one figure per target with the measured curves in gray and one coloured
%   line per candidate.

clear; clc; close all;

%% Configuration
horizonHours  = 10.0;                 % T in z = t / T
runNames      = ["run1", "run2"];     % one fit per run
nDoe          = 20;                   % first nDoe files of a run are the DOE
zMin          = 0.1;                  % first point of the grid
zStep         = 0.05;                 % step of the grid
nDeploy       = 100;                  % samples inside deploy_model
zReport       = 0.1:0.1:0.9;          % fidelities scored in the test
targetNames   = ["J_TV", "J_track"];
targetFields  = ["partial_SSdU", "partial_SSE"];
targetOffsets = [1, 0];               % index in out.T of the first partial entry

chebOrder = 4;                        % order of the Chebyshev candidates
monoOrder = 3;                        % free coefficients of the monotone candidate
nQuad     = 401;                      % quadrature nodes of the monotone candidate
lambda0   = 0.3;                      % warm start of expsat

% Candidates to fit. Comment a name out to skip it.
modelNames = ["cheb_deltalog", "cheb_ratio", "expsat", "monotone", "beta_cdf"];

%% Paths
scriptDir   = fileparts(mfilename("fullpath"));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(fullfile(projectRoot, "dependencies")));

set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");

%% Model registry
models = build_models(chebOrder, monoOrder, nQuad, lambda0);
[keep, order] = ismember(modelNames, [models.name]);
if ~all(keep)
    error("Unknown model name: %s", strjoin(modelNames(~keep), ", "));
end
models = models(order);

%% Fit every candidate on the DOE pairs of each run, test on that run's reruns
nRuns  = numel(runNames);
nTgt   = numel(targetNames);
nMdl   = numel(models);
S      = struct([]);

for r = 1:nRuns
    doeDir  = fullfile(projectRoot, "results", runNames(r));
    testDir = fullfile(projectRoot, "results", "final_fidelity_same_noise", ...
        runNames(r) + "_full_f1_same_noise");

    for t = 1:nTgt
        allPairs  = load_pairs(doeDir, "out_*.mat", targetFields(t), ...
            targetOffsets(t), horizonHours);
        doePairs  = allPairs(1:nDoe);
        testPairs = load_pairs(testDir, "out_full_*.mat", targetFields(t), ...
            targetOffsets(t), horizonHours);

        obs = pair_observations(doePairs, zMin, zStep);
        if isempty(obs.dLog)
            error("No observations: no DOE run of %s reaches z = %.2f", ...
                runNames(r), zMin + zStep);
        end
        w = bin_weights(obs.iLo, obs.iHi);

        fprintf("\n%s %s: fit on %d DOE pairs, test on %d full-horizon reruns\n", ...
            runNames(r), targetNames(t), numel(doePairs), numel(testPairs));
        counts = bin_counts(obs.iLo, obs.iHi);
        fprintf("  observations %d in %d bins, counts per bin from %d to %d\n", ...
            numel(obs.dLog), numel(counts), min(counts), max(counts));
        warn_if_truncated(testPairs);

        for m = 1:nMdl
            tFit    = tic;                % local handle, fmincon may call tic
            theta   = fit_model(models(m), obs, w);
            fitTime = toc(tFit);

            S(r, t, m).theta     = theta;
            S(r, t, m).testPairs = testPairs;

            fprintf("\n  [%s] %s, fitted in %.2f s\n", ...
                models(m).name, models(m).label, fitTime);
            fprintf("    parameters %s\n", ...
                strjoin(compose("%+.6f", theta(:)).', "  "));
            print_test_table(testPairs, models(m), theta, zReport, nDeploy);
        end

        print_summary_table(S, r, t, models, zReport, nDeploy);
    end
end

%% Figure: measured f(z) against every candidate
colors = lines(nMdl);
for t = 1:nTgt
    figure('Name', sprintf('Runtime surrogate candidates, %s', targetNames(t)));
    zPlot = linspace(0, 1, 201).';
    for r = 1:nRuns
        subplot(1, nRuns, r); hold on

        testPairs = S(r, t, 1).testPairs;
        hData     = gobjects(0, 1);
        for i = 1:numel(testPairs)
            [zData, fData] = measured_ratio(testPairs(i));
            hData(1) = plot(zData, fData, '-', 'Color', [0.6 0.6 0.6 0.6], ...
                'LineWidth', 1.0);
        end

        hMdl = gobjects(nMdl, 1);
        for m = 1:nMdl
            g = deploy_model(zPlot, models(m), S(r, t, m).theta, nDeploy);
            hMdl(m) = plot(zPlot, g, '-', 'LineWidth', 2.0, ...
                'Color', colors(m, :));
        end

        yline(1, ':', 'Color', [0.4 0.4 0.4]);
        xlim([0 1]); ylim([0 1.15]);
        xlabel("Fidelity $z$"); ylabel("$f(z)$");
        title(sprintf("%s, %s", runNames(r), targetNames(t)), "Interpreter", "none");
        if r == 1
            legend([hMdl; hData], ...
                [[models.name].'; repmat("measured", numel(hData), 1)], ...
                "Location", "southeast", "Interpreter", "none");
        end
        grid on; box on
    end
    sgtitle(sprintf("Deployed surrogate of each candidate, %s", targetNames(t)), ...
        "Interpreter", "none");
end

%% Local functions: model registry

function models = build_models(chebOrder, monoOrder, nQuad, lambda0)
%BUILD_MODELS Registry of candidate forms for f(z).
% Each entry holds the warm start theta0, the bounds, the linear constraints
% used by fmincon, an eval handle f = eval(z, theta), the name of the loss, and
% the flag monotone. monotone is true when the form is nondecreasing and has
% exact endpoints for every parameter value, in which case deploy_model needs no
% running maximum.
    models = struct('name', {}, 'label', {}, 'loss', {}, 'theta0', {}, ...
        'lb', {}, 'ub', {}, 'A', {}, 'b', {}, 'Aeq', {}, 'beq', {}, ...
        'eval', {}, 'monotone', {});

    % Chebyshev polynomial. The endpoints are linear equalities on the
    % coefficients, and positivity is imposed on a grid so that the log of the
    % model stays defined.
    zPos     = linspace(0.05, 1, 40);
    chebEval = @(z, c) cheb_features(z, chebOrder) * c(:);
    chebC0   = zeros(chebOrder + 1, 1);
    chebC0(1) = 0.5;
    chebC0(2) = 0.5;                  % start at f(z) = z, which is feasible

    k = 1;
    models(k).name     = "cheb_deltalog";
    models(k).label    = sprintf("Chebyshev order %d, log-ratio residual", chebOrder);
    models(k).loss     = "deltalog";
    models(k).theta0   = chebC0;
    models(k).lb       = [];
    models(k).ub       = [];
    models(k).A        = -cheb_features(zPos, chebOrder);
    models(k).b        = -1e-4 * ones(numel(zPos), 1);
    models(k).Aeq      = [cheb_features(0, chebOrder); cheb_features(1, chebOrder)];
    models(k).beq      = [0; 1];
    models(k).eval     = chebEval;
    models(k).monotone = false;

    k = k + 1;
    models(k)          = models(k - 1);
    models(k).name     = "cheb_ratio";
    models(k).label    = sprintf("Chebyshev order %d, relative ratio residual", chebOrder);
    models(k).loss     = "ratio";

    % Saturating exponential, normalised to f(1) = 1.
    k = k + 1;
    models(k).name     = "expsat";
    models(k).label    = "expm1(-lambda z) / expm1(-lambda)";
    models(k).loss     = "deltalog";
    models(k).theta0   = lambda0;
    models(k).lb       = -30;
    models(k).ub       = 30;
    models(k).A        = [];
    models(k).b        = [];
    models(k).Aeq      = [];
    models(k).beq      = [];
    models(k).eval     = @(z, th) exp_sat(z, th);
    models(k).monotone = true;

    % Normalised integral of a positive integrand. The constant term of the
    % exponent cancels in the normalisation and is held at zero.
    k = k + 1;
    models(k).name     = "monotone";
    models(k).label    = sprintf("normalised integral of exp(Chebyshev order %d)", monoOrder);
    models(k).loss     = "deltalog";
    models(k).theta0   = zeros(monoOrder, 1);
    models(k).lb       = -8 * ones(monoOrder, 1);
    models(k).ub       =  8 * ones(monoOrder, 1);
    models(k).A        = [];
    models(k).b        = [];
    models(k).Aeq      = [];
    models(k).beq      = [];
    models(k).eval     = @(z, a) monotone_integral(z, a, nQuad);
    models(k).monotone = true;

    % Regularised incomplete beta function.
    k = k + 1;
    models(k).name     = "beta_cdf";
    models(k).label    = "regularised incomplete beta I_z(a, b)";
    models(k).loss     = "deltalog";
    models(k).theta0   = [1; 1];
    models(k).lb       = [0.05; 0.05];
    models(k).ub       = [25; 25];
    models(k).A        = [];
    models(k).b        = [];
    models(k).Aeq      = [];
    models(k).beq      = [];
    models(k).eval     = @(z, ab) beta_cdf(z, ab);
    models(k).monotone = true;
end

function X = cheb_features(z, order)
%CHEB_FEATURES Chebyshev design matrix on the map z in [0,1] to x in [-1,1].
% Column k holds T_k(2z - 1). The recurrence is T_0 = 1, T_1 = x and
% T_{k+1} = 2 x T_k - T_{k-1}.
    z = z(:);
    x = 2 * z - 1;
    X = zeros(numel(z), order + 1);
    X(:, 1) = 1;
    if order >= 1
        X(:, 2) = x;
    end
    for k = 2:order
        X(:, k+1) = 2 * x .* X(:, k) - X(:, k-1);
    end
end

function f = exp_sat(z, lambda)
%EXP_SAT Saturating exponential normalised to f(0) = 0 and f(1) = 1.
% expm1 keeps the small-lambda quotient accurate. At lambda exactly zero the
% quotient is 0 / 0, and the limit is f(z) = z.
    z = z(:);
    if abs(lambda) < 1e-8
        f = z;
        return
    end
    f = expm1(-lambda * z) / expm1(-lambda);
end

function f = monotone_integral(z, a, nQuad)
%MONOTONE_INTEGRAL Normalised integral of a positive integrand.
% The integrand is exp(p(s)) with p a Chebyshev polynomial whose constant term
% is zero. The cumulative trapezoid on a fixed grid of nQuad nodes gives an
% increasing curve, and dividing by its last value sets f(1) = 1. Subtracting
% max(p) before the exponential keeps the integrand inside the double range at a
% trial point with a large exponent, and the constant factor it removes cancels
% in the normalisation. The returned f is the piecewise-linear interpolant of the
% trapezoid rule, not the exact integral, and the fit and the deployed surrogate
% use that same approximation.
    z     = z(:);
    a     = a(:);
    sGrid = linspace(0, 1, nQuad).';
    p     = cheb_features(sGrid, numel(a)) * [0; a];
    g     = exp(p - max(p));          % shift cancels in the normalisation
    F     = cumtrapz(sGrid, g);
    F     = F / F(end);
    f     = interp1(sGrid, F, min(max(z, 0), 1), "linear");
end

function f = beta_cdf(z, ab)
%BETA_CDF Regularised incomplete beta function I_z(a, b).
% Increasing in z with I_0 = 0 and I_1 = 1 for every positive a and b.
    z = min(max(z(:), 0), 1);
    f = betainc(z, ab(1), ab(2));
end

%% Local functions: data

function files = list_files(folder, pattern)
%LIST_FILES Return the files of one folder that match the pattern, sorted.
% The names carry a timestamp, so the sorted order is the order of the runs.
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

function zGrid = obs_grid(zEnd, zMin, zStep)
%OBS_GRID Return the grid zMin:zStep:zEnd as a column.
% The grid is empty when the run stops before zMin + zStep.
    zGrid = (zMin:zStep:zEnd).';
    if numel(zGrid) < 2
        zGrid = [];
    end
end

function obs = observations_one_pair(pair, zMin, zStep)
%OBSERVATIONS_ONE_PAIR Ratios of every ordered pair of grid points of one run.
% The loops take every zLo on the grid with every zHi above it. The fields are
% columns of equal length. ratio holds J(zHi) / J(zLo) and dLog holds its
% logarithm, so a loss can use either transform. iLo and iHi hold the grid
% indices of the two points, and the weights use them later.
    obs   = struct('zLo', [], 'zHi', [], 'ratio', [], 'dLog', [], ...
        'iLo', [], 'iHi', []);
    zGrid = obs_grid(pair.zEnd, zMin, zStep);
    if isempty(zGrid)
        return
    end
    J = interp1(pair.z, pair.J, zGrid, "linear");

    for a = 1:numel(zGrid) - 1
        for b = a + 1:numel(zGrid)
            obs.zLo   = [obs.zLo;   zGrid(a)];
            obs.zHi   = [obs.zHi;   zGrid(b)];
            obs.ratio = [obs.ratio; J(b) / J(a)];
            obs.dLog  = [obs.dLog;  log(J(b)) - log(J(a))];
            obs.iLo   = [obs.iLo;   a];
            obs.iHi   = [obs.iHi;   b];
        end
    end
end

function obs = pair_observations(pairs, zMin, zStep)
%PAIR_OBSERVATIONS Stack the ratios of every pair into one set.
    obs = struct('zLo', [], 'zHi', [], 'ratio', [], 'dLog', [], ...
        'iLo', [], 'iHi', []);
    for i = 1:numel(pairs)
        one = observations_one_pair(pairs(i), zMin, zStep);
        obs.zLo   = [obs.zLo;   one.zLo];
        obs.zHi   = [obs.zHi;   one.zHi];
        obs.ratio = [obs.ratio; one.ratio];
        obs.dLog  = [obs.dLog;  one.dLog];
        obs.iLo   = [obs.iLo;   one.iLo];
        obs.iHi   = [obs.iHi;   one.iHi];
    end
end

function id = bin_id(iLo, iHi)
%BIN_ID Number every distinct combination (iLo, iHi) from 1 upward.
% Two observations share an id when they use the same two grid points. They then
% come from different simulation pairs.
    [~, ~, id] = unique([iLo, iHi], "rows");
end

function counts = bin_counts(iLo, iHi)
%BIN_COUNTS Number of observations in each bin, in the order of the bin ids.
    counts = accumarray(bin_id(iLo, iHi), 1);
end

function w = bin_weights(iLo, iHi)
%BIN_WEIGHTS Weight of each observation, one over the count of its own bin.
% Every bin then carries the same total weight.
    id     = bin_id(iLo, iHi);
    counts = accumarray(id, 1);
    w      = 1 ./ counts(id);
end

%% Local functions: fit

function f = model_value(z, model, theta)
%MODEL_VALUE Raw value of one candidate at z, with no envelope and no cap.
    f = model.eval(z, theta);
    f = f(:);
end

function cost = model_cost(theta, model, obs, w)
%MODEL_COST Weighted mean squared residual of one candidate.
% deltalog compares the logarithm of the measured ratio with the logarithm of
% the model ratio. ratio compares the two ratios themselves, divided by the
% model ratio so that the residual stays scale free. The floor of 1e-6 keeps
% both quotients defined at a trial point that the solver may place outside the
% feasible set. It does not clip the model.
    fLo = max(model_value(obs.zLo, model, theta), 1e-6);
    fHi = max(model_value(obs.zHi, model, theta), 1e-6);

    switch model.loss
        case "deltalog"
            r = obs.dLog - (log(fHi) - log(fLo));
        case "ratio"
            r = obs.ratio ./ (fHi ./ fLo) - 1;
        otherwise
            error("Unknown loss %s", model.loss);
    end
    cost = sum(w .* r.^2) / sum(w);
end

function theta = fit_model(model, obs, w)
%FIT_MODEL Minimise model_cost under the constraints of the candidate.
% fmincon uses finite differences. Every candidate except expsat starts from the
% parameter value that gives f(z) = z, which is feasible for all of them. expsat
% starts at lambda = lambda0.
    opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
        'MaxFunctionEvaluations', 20000, 'MaxIterations', 2000);
    theta = fmincon(@(th) model_cost(th, model, obs, w), model.theta0, ...
        model.A, model.b, model.Aeq, model.beq, model.lb, model.ub, [], opts);
end

function g = deploy_model(z, model, theta, nSamples)
%DEPLOY_MODEL Surrogate value as it ships.
% A candidate flagged monotone is already nondecreasing with exact endpoints, so
% the value is used directly. For the others the function samples nSamples
% points on [0, z] and takes the largest value, which makes the curve
% nondecreasing in z, then caps the result at 1 so that f(1) = 1 holds.
    z = z(:);
    if model.monotone
        g = min(1, max(model_value(z, model, theta), 0));
        return
    end
    g = zeros(numel(z), 1);
    for i = 1:numel(z)
        s    = linspace(0, z(i), nSamples).';
        g(i) = min(1, max(model_value(s, model, theta)));
    end
end

%% Local functions: test and report

function [z, f] = measured_ratio(pair)
%MEASURED_RATIO Measured f(z) of one full-horizon run.
% J(end) is the measured total, so J(z) / J(end) is f(z) with no model and no
% extrapolation. The result is only meaningful when the run reached z = 1.
    z = pair.z;
    f = pair.J / pair.J(end);
end

function f = measured_ratio_at(pair, z)
%MEASURED_RATIO_AT Measured f at one fidelity, by interpolation of the curve.
    f = interp1(pair.z, pair.J, z, "linear") / pair.J(end);
end

function e = total_error(pair, z, model, theta, nSamples)
%TOTAL_ERROR Error of the total that the surrogate reports from a run cut at z.
% The estimate is Jhat(1) = J(z) / deploy_model(z). The error is
% Jhat(1) / J(1) - 1, which equals f_measured(z) / deploy_model(z) - 1.
    e = measured_ratio_at(pair, z) / deploy_model(z, model, theta, nSamples) - 1;
end

function [medAbs, bias] = test_scores(pairs, model, theta, zReport, nSamples)
%TEST_SCORES Median absolute error and mean error at each reported fidelity.
    medAbs = zeros(numel(zReport), 1);
    bias   = zeros(numel(zReport), 1);
    for k = 1:numel(zReport)
        e = zeros(numel(pairs), 1);
        for i = 1:numel(pairs)
            e(i) = total_error(pairs(i), zReport(k), model, theta, nSamples);
        end
        medAbs(k) = median(abs(e));
        bias(k)   = mean(e);
    end
end

function warn_if_truncated(pairs)
%WARN_IF_TRUNCATED Report any test run that did not reach the full horizon.
% The measured f(z) of such a run is wrong, because J(end) is not the total.
    for i = 1:numel(pairs)
        if pairs(i).zEnd < 0.999
            warning("Test run %d stops at z = %.3f and is not a full horizon.", ...
                i, pairs(i).zEnd);
        end
    end
end

function print_test_table(pairs, model, theta, zReport, nSamples)
%PRINT_TEST_TABLE Print the measured f, the model, and the error at each z.
% f_meas is the median over the test runs. The error column is the median
% absolute error of the predicted total, and bias is its mean.
    [medAbs, bias] = test_scores(pairs, model, theta, zReport, nSamples);
    fprintf("      %5s %6s %8s %8s %8s %9s %9s\n", ...
        "z", "n", "f_meas", "f_fit", "deploy", "median|e|", "bias");
    for k = 1:numel(zReport)
        z     = zReport(k);
        fMeas = zeros(numel(pairs), 1);
        for i = 1:numel(pairs)
            fMeas(i) = measured_ratio_at(pairs(i), z);
        end
        fprintf("      %5.2f %6d %8.4f %8.4f %8.4f %9.4f %+9.4f\n", ...
            z, numel(pairs), median(fMeas), model_value(z, model, theta), ...
            deploy_model(z, model, theta, nSamples), medAbs(k), bias(k));
    end
end

function print_summary_table(S, r, t, models, zReport, nSamples)
%PRINT_SUMMARY_TABLE Median absolute error of every candidate side by side.
% One row per candidate, one column per reported fidelity, and a last column
% holding the mean over the reported fidelities.
    nMdl = numel(models);
    fprintf("\n  median|e| by candidate\n");
    fprintf("    %-16s", "z");
    fprintf("%8.2f", zReport);
    fprintf("%9s\n", "mean");
    for m = 1:nMdl
        medAbs = test_scores(S(r, t, m).testPairs, models(m), ...
            S(r, t, m).theta, zReport, nSamples);
        fprintf("    %-16s", models(m).name);
        fprintf("%8.4f", medAbs);
        fprintf("%9.4f\n", mean(medAbs));
    end
end
