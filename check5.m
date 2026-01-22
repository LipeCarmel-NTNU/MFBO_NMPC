%% Fit model: y = Tmin + (Tplus - Tmin) * (1 + ap*p.^ep + am*m.^em) .* (1 - exp(-Iter_max/N))
% Requires vectors Iter_max, y, m, p (same length). Iter_max integer.
clear; close all; clc
rng(1)

%% LOADING
results_csv = fullfile("results","results.csv");
T = readtable(results_csv, 'TextType','string');

% Response
y = T.runtime_s;

% Theta matrix
theta_cols = startsWith(T.Properties.VariableNames, "theta_");
Theta = T{:, theta_cols};

% Decode (nx=3, nu=3)
nx = 3; nu = 3;
expected_theta_len = 1 + 2 + nx + nu + nu;
if size(Theta,2) ~= expected_theta_len
    error("Expected %d theta columns, found %d.", expected_theta_len, size(Theta,2));
end

Iter_max = Theta(:,1);
theta_p  = Theta(:,2);
theta_m  = Theta(:,3);

m = theta_m + 1;
p = theta_p + m;

q1   = Theta(:,4);  q2   = Theta(:,5);  q3   = Theta(:,6);
ru1  = Theta(:,7);  ru2  = Theta(:,8);  ru3  = Theta(:,9);
rdu1 = Theta(:,10); rdu2 = Theta(:,11); rdu3 = Theta(:,12);

%% --- Model definition ---
% theta = [Tplus, Tmin, N, ap, ep, am, em]
model = @(theta, X) ...
    theta(2) + ...
    (theta(1)) .* ...
    (1 + theta(4).*X(:,2).^theta(5) + theta(6).*X(:,3).^theta(7)) .* ...
    (1 - exp(-X(:,1)./theta(3)));

X = [Iter_max, p, m];

%% --- Initial guess ---
Tplus0 = 5.8652;
Tmin0  = 0.4066;
N0     = 6.5001;

ap0 = 0.0;  ep0 = 1.0;
am0 = 0.0;  em0 = 1.0;

theta0 = [Tplus0, Tmin0, N0, ap0, ep0, am0, em0];

%% --- Bounds ---
lb = [0, 0, 1e-1, 0, 0, 0, 0];
ub = [Inf, 1e3, Inf, Inf, 10, Inf, 10];

%% --- Regularisation setup ---
lambda = 1e-2;          % regularisation weight
idx_reg = 3:7;          

residual_fun = @(theta) [
    model(theta, X) - y;
    sqrt(lambda) * theta(idx_reg).'
];

%% --- Fit (regularised nonlinear least squares) ---
opts = optimoptions('lsqnonlin', ...
    'Display','iter', ...
    'MaxFunctionEvaluations',5e4, ...
    'MaxIterations',2e3);

theta = lsqnonlin(residual_fun, theta0, lb, ub, opts);

Tplus = theta(1);
Tmin  = theta(2);
N     = theta(3);
ap    = theta(4);
ep    = theta(5);
am    = theta(6);
em    = theta(7);

fprintf('\nFit:\n  Tplus = %.6g\n  Tmin = %.6g\n  N    = %.6g\n  ap   = %.6g\n  ep   = %.6g\n  am   = %.6g\n  em   = %.6g\n', ...
    Tplus, Tmin, N, ap, ep, am, em);

%% --- Goodness of fit ---
yhat = model(theta, X);

SS_res = sum((y - yhat).^2);
SS_tot = sum((y - mean(y)).^2);
R2 = 1 - SS_res/SS_tot;

fprintf('\nGoodness of fit:\n  R^2 = %.6f\n', R2);

%% --- y vs yhat ---
figure; hold on; grid on
scatter(y, yhat, 16, 'filled')

ymin = min([y; yhat]);
ymax = max([y; yhat]);
plot([ymin ymax], [ymin ymax], 'k--', 'LineWidth', 1.5)

xlabel('Observed y')
ylabel('Predicted \hat{y}')
title(sprintf('y vs \\hat{y} (R^2 = %.4f)', R2))
axis equal
xlim([ymin ymax])
ylim([ymin ymax])
