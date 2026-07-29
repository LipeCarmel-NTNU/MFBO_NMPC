%% Fit the beta runtime surrogate f(z) on the DOE, test it on the full-horizon reruns
%
% Model
%   J(z) = J(1) * f(z), with f the regularised incomplete beta function
%       f(z) = I_z(a, b),   a > 0,   b > 0.
%   For every a and b this gives f(0) = 0, f(1) = 1 and f'(z) > 0 on (0, 1), so
%   the fitted curve is already nondecreasing and already normalised. It ships as
%   it is fitted: there is no running maximum, no cap and no separate deploy
%   model. a = b = 1 is f(z) = z, a constant cost rate.
%
% Fit
%   Take logarithms of J(z) = J(1) f(z):
%       log J_i(z) = c_i + log f(z),    c_i = log J_i(1),
%   where i indexes the DOE runs. A truncated run never reaches z = 1, so its
%   total is unknown and c_i is a free constant, one per run. Write
%       u_i(z) = log J_i(z) - log f(z).
%   For a trial (a, b) the least-squares value of c_i is the mean of u_i over
%   that run's samples. Substituting it removes every c_i and leaves
%       cost(a, b) = sum_i sum_z ( u_i(z) - mean_z u_i(z) )^2,
%   the squared deviations of u_i from its own mean. Minimising that over (a, b)
%   is the whole fit. It uses every stored sample of every DOE run above zMin,
%   and forms no ratios, no pairs, no bins and no interpolation.
%
% Test
%   The files in results/final_fidelity_same_noise are the frontier controllers
%   rerun with theta(1) = 1. Those runs cover the whole horizon, so J(end) is the
%   measured total and f_measured(z) = J(z) / J(1) is exact at every z. Nothing
%   is extrapolated on the test side. The score is the error of the total that
%   the surrogate would report from a run truncated at z:
%       Jhat(1) = J(z) / f(z),   error = Jhat(1) / J(1) - 1.
%   The error goes to zero as z goes to 1 for both sides by construction, so the
%   low-z rows of the table carry the information.

clear; clc; close all;

%% Configuration
horizonHours  = 10.0;                 % T in z = t / T
runNames      = ["run1", "run2"];     % one fit per run
nDoe          = 20;                   % first nDoe files of a run are the DOE
zMin          = 0.1;                  % samples below this fidelity are dropped
zReport       = 0.1:0.1:0.9;          % fidelities scored in the test
targetNames   = ["J_TV", "J_track"];
targetFields  = ["partial_SSdU", "partial_SSE"];
targetOffsets = [1, 0];               % index in out.T of the first partial entry

ab0 = [1; 1];                         % warm start, f(z) = z
abLb = [0.05; 0.05];
abUb = [25; 25];

%% Paths
scriptDir   = fileparts(mfilename("fullpath"));
projectRoot = fileparts(fileparts(scriptDir));
addpath(genpath(fullfile(projectRoot, "dependencies")));

set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");

%% Fit on the DOE of each run, test on that run's full-horizon reruns
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

        samples = fit_samples(doePairs, zMin);
        if numel(samples) < 2
            error("Only %d usable DOE runs of %s reach z = %.2f", ...
                numel(samples), runNames(r), zMin);
        end

        ab        = fit_beta(samples, ab0, abLb, abUb);
        [~, nObs] = beta_cost(ab, samples);

        S(r, t).ab        = ab;
        S(r, t).testPairs = testPairs;

        fprintf("\n%s %s\n", runNames(r), targetNames(t));
        fprintf("  fit on %d of %d DOE runs, %d samples above z = %.2f\n", ...
            numel(samples), numel(doePairs), nObs, zMin);
        fprintf("  a = %.4f, b = %.4f\n", ab(1), ab(2));
        warn_if_truncated(testPairs);
        print_test_table(testPairs, ab, zReport);
    end
end

%% Figure: measured f(z) against the fit
figure('Name', 'Beta runtime surrogate');
zPlot = linspace(0, 1, 201).';
panel = 0;
for t = 1:nTgt
    for r = 1:nRuns
        panel = panel + 1;
        subplot(nTgt, nRuns, panel); hold on

        testPairs = S(r, t).testPairs;
        hData     = gobjects(0, 1);
        for i = 1:numel(testPairs)
            [zData, fData] = measured_ratio(testPairs(i));
            hData(1) = plot(zData, fData, '-', 'Color', [0.6 0.6 0.6 0.6], ...
                'LineWidth', 1.0);
        end
        hFit = plot(zPlot, beta_ratio(zPlot, S(r, t).ab), '-', ...
            'LineWidth', 2.4, 'Color', [0.84 0.37 0.00]);

        yline(1, ':', 'Color', [0.4 0.4 0.4]);
        xlim([0 1]); ylim([0 1.15]);
        xlabel("Fidelity $z$"); ylabel("$f(z)$");
        title(sprintf("%s, %s, $a = %.2f$, $b = %.2f$", runNames(r), ...
            targetNames(t), S(r, t).ab(1), S(r, t).ab(2)));
        if panel == 1 && ~isempty(hData)
            legend([hFit; hData], ["beta fit", "measured"], ...
                "Location", "southeast");
        end
        grid on; box on
    end
end
sgtitle("Fit on the DOE, measured f(z) from the full-horizon reruns", ...
    "Interpreter", "none");

%% Local functions: model and fit

function f = beta_ratio(z, ab)
%BETA_RATIO f(z) = I_z(a, b), the regularised incomplete beta function.
% Increasing in z with f(0) = 0 and f(1) = 1 for every positive a and b.
    f = betainc(min(max(z(:), 0), 1), ab(1), ab(2));
end

function samples = fit_samples(pairs, zMin)
%FIT_SAMPLES Stored samples of each DOE run above zMin, with log of the cost.
% A run needs two usable samples to say anything once its own constant is
% profiled out, so shorter ones are dropped.
    samples = struct('z', {}, 'logJ', {});
    for i = 1:numel(pairs)
        keep = pairs(i).z >= zMin & pairs(i).J > 0;
        if sum(keep) < 2
            continue
        end
        samples(end + 1).z  = pairs(i).z(keep);   %#ok<SAGROW>
        samples(end).logJ   = log(pairs(i).J(keep));
    end
end

function [cost, nObs] = beta_cost(ab, samples)
%BETA_COST Squared deviations of u_i = log J_i - log f from its per-run mean.
% Centring on the mean is what profiling out the unknown log J_i(1) amounts to,
% so the cost is free of every run total. The floor of 1e-12 keeps the logarithm
% defined at a trial point where f underflows.
    cost = 0;
    nObs = 0;
    for i = 1:numel(samples)
        u    = samples(i).logJ - log(max(beta_ratio(samples(i).z, ab), 1e-12));
        cost = cost + sum((u - mean(u)).^2);
        nObs = nObs + numel(u);
    end
end

function ab = fit_beta(samples, ab0, lb, ub)
%FIT_BETA Minimise beta_cost over the two shape parameters.
% Only bounds constrain the problem. The start point a = b = 1 is f(z) = z.
    opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
        'MaxFunctionEvaluations', 5000, 'MaxIterations', 1000);
    ab = fmincon(@(p) beta_cost(p, samples), ab0, [], [], [], [], ...
        lb, ub, [], opts);
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

function print_test_table(pairs, ab, zReport)
%PRINT_TEST_TABLE Print the measured f, the fit, and the error at each z.
% f_meas is the median over the test runs. The error columns hold the median
% absolute error of the predicted total and its mean.
    fprintf("      %5s %6s %8s %8s %9s %9s\n", ...
        "z", "n", "f_meas", "f_fit", "median|e|", "bias");
    for z = zReport
        fFit  = beta_ratio(z, ab);
        fMeas = zeros(numel(pairs), 1);
        e     = zeros(numel(pairs), 1);
        for i = 1:numel(pairs)
            fMeas(i) = measured_ratio_at(pairs(i), z);
            e(i)     = fMeas(i) / fFit - 1;
        end
        fprintf("      %5.2f %6d %8.4f %8.4f %9.4f %+9.4f\n", ...
            z, numel(pairs), median(fMeas), fFit, median(abs(e)), mean(e));
    end
end
