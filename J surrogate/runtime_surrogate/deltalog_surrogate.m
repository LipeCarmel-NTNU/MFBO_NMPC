%% Fit the runtime surrogate f(z) on the DOE, test it on the full-horizon reruns
%
% Model
%   J(z) = J(1) * f(z). Take two fidelities of the same simulation. The unknown
%   total J(1) cancels in the difference of the logarithms:
%       log J(z2) - log J(z1) = log f(z2) - log f(z1).
%   A truncated run is therefore usable on its own. The fit reads out.T and the
%   two partial_* arrays only. It never reads out.SSE or out.SSdU, because those
%   two fields hold partial / frac(f) from the curve hard-coded in main_BO.m.
%
% Fit
%   f(z) is a Chebyshev polynomial of order chebOrder in x = 2z - 1. Every DOE
%   pair gives the grid zMin:zStep:zEnd. The fit uses every ordered combination
%   of two grid points of the same pair, so a grid of n points gives
%   n (n - 1) / 2 observations. The constraints are f(0) = 0 and f(1) = 1. No DOE
%   run reaches z = 1. The fit uses the raw polynomial and does not clip to
%   [0, 1].
%
% Weights
%   The DOE draws the fidelity from a Sobol sample of [0, 1]. A short run reaches
%   only the low grid points, so the low combinations appear in many pairs and
%   the high combinations appear in few. One bin is one combination (zLo, zHi).
%   Each observation gets the weight 1 / count of its own bin, so every
%   combination carries the same total weight.
%
% Test
%   The files in results/final_fidelity_same_noise are the frontier controllers
%   rerun with theta(1) = 1. Those runs cover the whole horizon, so J(end) is the
%   measured total and
%       f_measured(z) = J(z) / J(1)
%   is exact at every z. Nothing is extrapolated on the test side.
%   The score is the error of the total that the surrogate would report from a
%   run truncated at z:
%       Jhat(1) = J(z) / deploy_model(z),   error = Jhat(1) / J(1) - 1.
%
% Figure
%   One panel per run and target. Gray lines are the measured f(z) of each
%   full-horizon rerun. The orange line is the raw fit. The dashed black line is
%   deploy_model, which is the raw fit after the running maximum and the cap.

clear; clc; close all;

%% Configuration
horizonHours  = 10.0;                 % T in z = t / T
runNames      = ["run1", "run2"];     % one fit per run
nDoe          = 20;                   % first nDoe files of a run are the DOE
zMin          = 0.1;                 % first point of the log grid
zStep         = 0.05;                 % step of the log grid
chebOrder     = 4;                    % order of the polynomial
nDeploy       = 100;                  % samples inside deploy_model
zReport       = 0.1:0.1:0.9;          % fidelities scored in the test
targetNames   = ["J_TV", "J_track"];
targetFields  = ["partial_SSdU", "partial_SSE"];
targetOffsets = [1, 0];               % index in out.T of the first partial entry

%% Paths
scriptDir   = fileparts(mfilename("fullpath"));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(fullfile(projectRoot, "dependencies")));

set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");

%% Fit on the DOE pairs of each run, test on that run's full-horizon reruns
nRuns = numel(runNames);
nTgt  = numel(targetNames);
S     = struct([]);

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

        obs = log_ratios(doePairs, zMin, zStep);
        w   = bin_weights(obs.iLo, obs.iHi);
        c   = fit_coefficients(obs, w, chebOrder);

        S(r, t).c         = c;
        S(r, t).testPairs = testPairs;

        fprintf("\n%s %s: fit on %d DOE pairs, test on %d full-horizon reruns\n", ...
            runNames(r), targetNames(t), numel(doePairs), numel(testPairs));
        counts = bin_counts(obs.iLo, obs.iHi);
        fprintf("  observations %d in %d bins, counts per bin from %d to %d\n", ...
            numel(obs.dLog), numel(counts), min(counts), max(counts));
        fprintf("  coefficients %s\n", sprintf("% .6f", c));
        warn_if_truncated(testPairs);
        print_test_table(testPairs, c, zReport, nDeploy);
    end
end

%% Figure: measured f(z) against the fit
figure('Name', 'Measured f(z) and the fit');
zPlot = linspace(0, 1, 201).';
panel = 0;
for t = 1:nTgt
    for r = 1:nRuns
        panel = panel + 1;
        subplot(nTgt, nRuns, panel); hold on

        testPairs = S(r, t).testPairs;
        for i = 1:numel(testPairs)
            [zData, fData] = measured_ratio(testPairs(i));
            hData = plot(zData, fData, '-', 'Color', [0.6 0.6 0.6 0.6], ...
                'LineWidth', 1.0);
        end
        hFit = plot(zPlot, predict_ratio(zPlot, S(r, t).c), '-', ...
            'LineWidth', 2.4, 'Color', [0.84 0.37 0.00]);
        hDep = plot(zPlot, deploy_model(zPlot, S(r, t).c, nDeploy), '--', ...
            'LineWidth', 1.6, 'Color', [0 0 0]);

        yline(1, ':', 'Color', [0.4 0.4 0.4]);
        xlim([0 1]); ylim([0 1.15]);
        xlabel("Fidelity $z$"); ylabel("$f(z)$");
        title(sprintf("%s, %s", runNames(r), targetNames(t)), "Interpreter", "none");
        if panel == 1
            legend([hFit, hDep, hData], ["raw fit", "deploy model", "measured"], ...
                "Location", "southeast");
        end
        grid on; box on
    end
end
sgtitle("Fit on the DOE, measured f(z) from the full-horizon reruns", ...
    "Interpreter", "none");

%% Local functions

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

function grid = log_grid(zEnd, zMin, zStep)
%LOG_GRID Return the grid zMin:zStep:zEnd as a column.
% The grid is empty when the run stops before zMin + zStep.
    grid = (zMin:zStep:zEnd).';
    if numel(grid) < 2
        grid = [];
    end
end

function obs = log_ratios_one_pair(pair, zMin, zStep)
%LOG_RATIOS_ONE_PAIR Log differences of every ordered pair of grid points.
% The loops take every zLo on the grid with every zHi above it. The fields are
% columns of equal length. dLog holds log J(zHi) - log J(zLo). iLo and iHi hold
% the grid indices of the two points, and the weights use them later.
    obs  = struct('zLo', [], 'zHi', [], 'dLog', [], 'iLo', [], 'iHi', []);
    grid = log_grid(pair.zEnd, zMin, zStep);
    if isempty(grid)
        return
    end
    logJ = log(interp1(pair.z, pair.J, grid, "linear"));

    for a = 1:numel(grid) - 1
        for b = a + 1:numel(grid)
            obs.zLo  = [obs.zLo;  grid(a)];
            obs.zHi  = [obs.zHi;  grid(b)];
            obs.dLog = [obs.dLog; logJ(b) - logJ(a)];
            obs.iLo  = [obs.iLo;  a];
            obs.iHi  = [obs.iHi;  b];
        end
    end
end

function obs = log_ratios(pairs, zMin, zStep)
%LOG_RATIOS Stack the log differences of every pair into one set.
    obs = struct('zLo', [], 'zHi', [], 'dLog', [], 'iLo', [], 'iHi', []);
    for i = 1:numel(pairs)
        one = log_ratios_one_pair(pairs(i), zMin, zStep);
        obs.zLo  = [obs.zLo;  one.zLo];
        obs.zHi  = [obs.zHi;  one.zHi];
        obs.dLog = [obs.dLog; one.dLog];
        obs.iLo  = [obs.iLo;  one.iLo];
        obs.iHi  = [obs.iHi;  one.iHi];
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

function f = predict_ratio(z, c)
%PREDICT_RATIO Value of the raw polynomial at z.
% No running maximum and no cap. The fit and the figure use this form.
    f = cheb_features(z, numel(c) - 1) * c;
end

function cost = delta_log_cost(c, obs, w)
%DELTA_LOG_COST Weighted mean squared residual of the log differences.
% The floor of 1e-6 keeps the logarithm defined at a trial point that the solver
% may place outside the feasible set. It does not clip the model.
    fLo  = max(predict_ratio(obs.zLo, c), 1e-6);
    fHi  = max(predict_ratio(obs.zHi, c), 1e-6);
    r    = obs.dLog - (log(fHi) - log(fLo));
    cost = sum(w .* r.^2) / sum(w);
end

function c = fit_coefficients(obs, w, order)
%FIT_COEFFICIENTS Minimize delta_log_cost under the endpoint conditions.
% The equalities are f(0) = 0 and f(1) = 1. The inequalities keep f(z) positive
% on a grid of [0.05, 1], which the logarithm needs. The start point is
% f(z) = z, which satisfies every constraint. fmincon uses finite differences.
    Aeq = [cheb_features(0, order); cheb_features(1, order)];
    beq = [0; 1];
    A   = -cheb_features(linspace(0.05, 1, 40), order);
    b   = -1e-4 * ones(40, 1);

    c0    = zeros(order + 1, 1);
    c0(1) = 0.5;
    c0(2) = 0.5;

    opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
        'MaxFunctionEvaluations', 20000, 'MaxIterations', 2000);
    c = fmincon(@(cc) delta_log_cost(cc, obs, w), c0, A, b, Aeq, beq, ...
        [], [], [], opts);
end

function g = deploy_model(z, c, nSamples)
%DEPLOY_MODEL Surrogate value as it ships.
% The function samples nSamples points on [0, z], takes the largest value of the
% polynomial, and caps the result at 1. The largest value makes the curve
% nondecreasing in z. The cap holds f(1) = 1.
    z = z(:);
    g = zeros(numel(z), 1);
    for i = 1:numel(z)
        s    = linspace(0, z(i), nSamples).';
        g(i) = min(1, max(predict_ratio(s, c)));
    end
end

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

function e = total_error(pair, z, c, nSamples)
%TOTAL_ERROR Error of the total that the surrogate reports from a run cut at z.
% The estimate is Jhat(1) = J(z) / deploy_model(z). The error is
% Jhat(1) / J(1) - 1, which equals f_measured(z) / deploy_model(z) - 1.
    e = measured_ratio_at(pair, z) / deploy_model(z, c, nSamples) - 1;
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

function print_test_table(pairs, c, zReport, nSamples)
%PRINT_TEST_TABLE Print the measured f, the model, and the error at each z.
% f_meas is the median over the test runs. The error column is the median
% absolute error of the predicted total, and bias is its mean.
    fprintf("      %5s %6s %8s %8s %8s %9s %9s\n", ...
        "z", "n", "f_meas", "f_fit", "deploy", "median|e|", "bias");
    for z = zReport
        fMeas = zeros(numel(pairs), 1);
        e     = zeros(numel(pairs), 1);
        for i = 1:numel(pairs)
            fMeas(i) = measured_ratio_at(pairs(i), z);
            e(i)     = total_error(pairs(i), z, c, nSamples);
        end
        fprintf("      %5.2f %6d %8.4f %8.4f %8.4f %9.4f %+9.4f\n", ...
            z, numel(pairs), median(fMeas), predict_ratio(z, c), ...
            deploy_model(z, c, nSamples), median(abs(e)), mean(e));
    end
end
