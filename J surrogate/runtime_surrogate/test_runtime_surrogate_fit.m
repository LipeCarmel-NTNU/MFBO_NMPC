%% Pairwise log-ratio fitting of normalized fidelity curve
%
% Goal:
%   We want to fit the normalized curve
%
%       J(z) / J(1) = f(z, theta)
%
%   but for truncated simulations z_i < 1, J(1) is unavailable.
%
% True normalized model:
%
%       f_true(z) = (1 - exp(-z)) / (1 - exp(-1))
%
%   which corresponds to theta = 1.
%
% Artificial data model:
%
%       J_i(z) = k_i * f_true(z) * noise
%
%   where k_i > 0 is random and unknown for each simulation/run.
%
%   Because k_i cancels in ratios,
%
%       J_i(z_a) / J_i(z_b) = f(z_a) / f(z_b),
%
%   we can fit theta without knowing k_i or J_i(1).
%
% Fitted model:
%
%       f(z, theta) = (1 - exp(-theta*z)) / (1 - exp(-theta)),
%       theta > 0.
%
%   The constraint theta > 0 is enforced by optimizing
%
%       theta = exp(eta),
%
%   over unconstrained eta.

clear; clc; close all;

rng(1);

%% Synthetic data settings

nRuns = 40;              % number of independent truncated simulations
nZmax = 20;              % maximum number of z samples per run
zMin = 1/nZmax;          % avoid z = 0 because log(J(0)) = -Inf
noiseSigma = 0.03;       % multiplicative log-normal noise level

% Candidate z grid.
% This mimics oversampling near the beginning of the simulation.
zGrid = linspace(zMin, 1.0, nZmax);

% True normalized function. The true theta is 1.
fTrue = @(z) (1 - exp(-z)) ./ (1 - exp(-1));

% Model to be fitted.
fModel = @(z, theta) (1 - exp(-theta*z)) ./ (1 - exp(-theta));

%% Generate truncated artificial datasets

data = struct([]);

for i = 1:nRuns

    % Each run ends at a random truncated fidelity z_i <= 1.
    % Many runs do not reach full fidelity.
    ZEnd_i = randi([1, nZmax]);
    zEnd = zGrid(ZEnd_i);

    % Available z values for this run.
    zi = zGrid(zGrid <= zEnd);

    % Unknown positive scale for this run.
    % This plays the role of the unavailable J_i(1).
    ki = exp(1.0*randn());

    % Multiplicative noise.
    noise = exp(noiseSigma*randn(size(zi)));

    % Observed unnormalized cost.
    Ji = ki * fTrue(zi) .* noise;

    data(i).z = zi(:);
    data(i).J = Ji(:);
    data(i).k = ki;
    data(i).zEnd = max(zi);
end

%% Construct pairwise log-ratio residual data
%
% For each run i and every valid pair (z_a, z_b), define
%
%   r_iab(theta)
%       = log(J_i(z_a)) - log(J_i(z_b))
%         - [log(f(z_a,theta)) - log(f(z_b,theta))].
%
% Since
%
%   J_i(z) = k_i f_true(z),
%
% the unknown scale k_i cancels:
%
%   log(J_i(z_a)) - log(J_i(z_b))
%       = log(f_true(z_a)) - log(f_true(z_b)).
%
% Therefore the residual depends only on the shape of f, not on k_i.

za_all = [];
zb_all = [];
yPair_all = [];
runIdx_all = [];

for i = 1:nRuns
    z = data(i).z;
    J = data(i).J;
    n = numel(z);

    for a = 1:n-1
        for b = a+1:n

            za = z(a);
            zb = z(b);

            % Observed log-ratio.
            yPair = log(J(a)) - log(J(b));

            za_all(end+1,1) = za; %#ok<SAGROW>
            zb_all(end+1,1) = zb; %#ok<SAGROW>
            yPair_all(end+1,1) = yPair; %#ok<SAGROW>
            runIdx_all(end+1,1) = i; %#ok<SAGROW>
        end
    end
end

nPairs = numel(yPair_all);

%% Coverage-based weights
%
% If early z intervals are oversampled, pairwise residuals involving those
% intervals can dominate the objective. To reduce this bias, define a
% coverage count C(z), equal to the number of pair intervals covering each
% bin in z.
%
% For pair interval I_ab = [z_a, z_b], use
%
%   w_ab = 1 / average_C_over_Iab.
%
% Hence, pairs covering frequently represented regions receive lower weight.

nBins = 20;
binEdges = linspace(0, 1, nBins+1);
coverage = zeros(nBins,1);

% First count how many pair intervals cover each bin.
for p = 1:nPairs
    lo = min(za_all(p), zb_all(p));
    hi = max(za_all(p), zb_all(p));

    for k = 1:nBins
        binLo = binEdges(k);
        binHi = binEdges(k+1);

        overlaps = (hi > binLo) && (lo < binHi);

        if overlaps
            coverage(k) = coverage(k) + 1;
        end
    end
end

% Now assign one inverse-coverage weight per pair.
w = zeros(nPairs,1);

for p = 1:nPairs
    lo = min(za_all(p), zb_all(p));
    hi = max(za_all(p), zb_all(p));

    coveredBins = false(nBins,1);

    for k = 1:nBins
        binLo = binEdges(k);
        binHi = binEdges(k+1);

        overlaps = (hi > binLo) && (lo < binHi);

        if overlaps
            coveredBins(k) = true;
        end
    end

    meanCoverage = mean(coverage(coveredBins));

    % Avoid division by zero, though meanCoverage should be positive.
    w(p) = 1 / max(meanCoverage, eps);
end

% Optional normalization: keep the average weight near 1.
% This does not change the minimizer, but improves numerical scaling.
w = w / mean(w);

%% Weighted least-squares objective
%
% Optimize eta, where theta = exp(eta), so theta is always positive.

objective = @(eta) pairwiseObjective( ...
    eta, za_all, zb_all, yPair_all, w, fModel);

eta0 = log(0.5);   % initial guess theta = 0.5

opts = optimset( ...
    'Display', 'iter', ...
    'TolX', 1e-10, ...
    'TolFun', 1e-12, ...
    'MaxIter', 1000, ...
    'MaxFunEvals', 5000);

etaHat = fminsearch(objective, eta0, opts);
thetaHat = exp(etaHat);

fprintf('\nEstimated theta = %.6f\n', thetaHat);
fprintf('True theta      = %.6f\n', 1.0);

%% Compare fitted and true normalized curves

zFine = linspace(0, 1, 300).';
fTrueFine = fTrue(zFine);
fFitFine = fModel(zFine, thetaHat);

figure;
plot(zFine, fTrueFine, 'k-', 'LineWidth', 2); hold on;
plot(zFine, fFitFine, 'r--', 'LineWidth', 2);

xlabel('z');
ylabel('f(z)');
legend('True f(z), \theta = 1', ...
       sprintf('Fitted f(z), \\theta = %.4f', thetaHat), ...
       'Location', 'southeast');
grid on;
title('Normalized fidelity curve fit');

%% Plot artificial unnormalized data
%
% The raw J_i(z) curves have different vertical scales k_i.
% The fitting does not require knowing those scales.

figure; hold on;

for i = 1:nRuns
    plot(data(i).z, data(i).J, '.-', 'LineWidth', 0.75);
end

xlabel('z');
ylabel('J_i(z)');
grid on;
title('Artificial truncated data with unknown positive scales k_i');

%% Plot coverage function

binCenters = 0.5*(binEdges(1:end-1) + binEdges(2:end));

figure;
bar(binCenters, coverage, 1.0);
xlabel('z');
ylabel('coverage C(z)');
grid on;
title('Pair-interval coverage used for weighting');

%% Local function

function sse = pairwiseObjective(eta, za, zb, yPair, w, fModel)
    % Convert unconstrained eta to positive theta.
    theta = exp(eta);

    % Model-predicted log-ratio:
    %
    %   log(f(za,theta)) - log(f(zb,theta)).
    %
    % f is strictly positive for theta > 0 and z > 0.
    logRatioModel = log(fModel(za, theta)) - log(fModel(zb, theta));

    % Pairwise residual.
    r = yPair - logRatioModel;

    % Weighted least-squares objective.
    sse = sum(w .* r.^2);
end