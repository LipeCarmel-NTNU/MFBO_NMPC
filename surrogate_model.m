% fit_runtime_from_results_csv.m
% Fit runtime_s as a function of f (=theta_1), p and m from results.csv
%
% Model:
%   t = (m/10)^alfa * (p/30)^beta * g( a0 + a1*x )
%   x = 2*f - 1
%   g(y) = y + softplus_k(-y)   (smooth positivity mapping)
%
% Notes:
% - Reads results/results.csv
% - Skips the first data row after loading
% - Uses fminunc (quasi-newton)
% - Converts runtime_s (seconds) to hours once

clear; close all; clc
rng(1)

%% --- Load and parse ---
results_csv = fullfile("results_bckp","results.csv");
T = readtable(results_csv, 'TextType','string');

% Response: runtime in hours (convert once)
t = double(T.runtime_s) / 3600;

% Predictors
f       = double(T.theta_1);   % fraction (0,1)
theta_p = double(T.theta_2);
theta_m = double(T.theta_3);

% Decode m,p (same convention as surrogate scripts)
m = theta_m + 1;
p = theta_p + m;

% Affine coordinate in [-1, 1]
x = 2*f - 1;

%% --- Optimisation options ---
opts = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'MaxIterations', 2000, ...
    'OptimalityTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'Display','iter');

%% --- Fit time model ---
% th = [alfa; beta; a0; a1]
theta0 = [ ...
    0.1; ...   % alfa
    0.1; ...   % beta
    1.0; ...   % a0
    1.0  ...   % a1
];

[theta_hat, fval] = fminunc(@(th) obj_time(th, x, m, p, t), theta0, opts);

alfa = theta_hat(1);
beta = theta_hat(2);
a    = theta_hat(3:4);

%% --- Predictions and diagnostics ---
t_hat = time_model(x, m, p, alfa, beta, a);

SS_res = sum((t - t_hat).^2);
SS_tot = sum((t - mean(t)).^2);
R2 = 1 - SS_res/SS_tot;

fprintf('\n===== Runtime model fit (skip first row) =====\n');
fprintf('alfa = %.8g\n', alfa);
fprintf('beta = %.8g\n', beta);
fprintf('a_time (a0,a1):\n');
fprintf('% .6e  % .6e\n', a(1), a(2));
fprintf('R^2 = %.6f\n', R2);

figure;
plot(f, t, 'o'); hold on
plot(f, t_hat, 'x');
grid on
xlabel('f = theta_1');
ylabel('runtime (h)');
legend('data','model','Location','best');

figure;
plot(t, t_hat, 'o', 'MarkerSize',10); grid on
hold on
mx = max([t; t_hat]);
plot([0 mx], [0 mx], 'k-')
xlabel('runtime (h)');
ylabel('runtime\_hat (h)');
title(sprintf('Parity plot, R^2 = %.4f', R2));

%% -----------------------------
% Shape term as a function of f
% -----------------------------
f_plot = linspace(0, 1, 200).';
x_plot = 2*f_plot - 1;

% Raw affine term in x
y_aff = affine_in_x(x_plot, a);

% Positivity mapping used in the time model
k = 20;
softplus_k = @(z) log1p(exp(k*z)) / k;
shape_pos = y_aff + softplus_k(-y_aff);

figure;
plot(f_plot, shape_pos, 'LineWidth', 1.5); hold on
plot(f_plot, y_aff, '--', 'LineWidth', 1.2);
grid on
xlabel('f');
ylabel('Value');
legend('Positive-mapped shape term', 'Raw affine term', 'Location','best');
title('Shape contribution to runtime as a function of f');

%% =============================
% Objective
% =============================
function J = obj_time(th, x, m, p, t)
    alfa = th(1);
    beta = th(2);
    a    = th(3:end);

    t_hat = time_model(x, m, p, alfa, beta, a);
    r = t - t_hat;

    lambda = 1e-1;  % regularisation weight

    fit_cost = r.'*r;
    reg_cost = lambda*(alfa^2 + beta^2 + sumsqr(a(2:end)));  % do not penalise offset
    J = fit_cost + reg_cost;

    if ~isfinite(J)
        J = realmax;
    end
end

%% =============================
% Models
% =============================
function y = affine_in_x(x, a)
    a = a(:);
    if numel(a) ~= 2
        error('affine_in_x expects 2 coefficients [a0; a1].');
    end
    y = a(1) + a(2)*x;
end

function t_hat = time_model(x, m, p, alfa, beta, a)
    k = 20;
    softplus_k = @(z) log1p(exp(k*z)) / k;

    y = affine_in_x(x, a);
    shape_pos = y + softplus_k(-y);  % smooth positivity mapping

    scale = (m/10).^(alfa^2) .* (p/30).^(beta^2);
    t_hat = scale .* shape_pos;
end
