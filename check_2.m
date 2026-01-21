clear; close all; clc

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

q1 = Theta(:,4);  q2 = Theta(:,5);  q3 = Theta(:,6);
ru1 = Theta(:,7); ru2 = Theta(:,8); ru3 = Theta(:,9);
rdu1 = Theta(:,10); rdu2 = Theta(:,11); rdu3 = Theta(:,12);

%% Clean rows used for fitting (only those needed by the model)
ok = isfinite(y) & isfinite(Iter_max) & isfinite(m) & isfinite(p);
y = y(ok);
Iter_max = Iter_max(ok);
m = m(ok);
p = p(ok);

%% Nonlinear model fit with fmincon
% Model:
%   yhat = (T_base + c10*10^(cm*m + cp*p))*(1 - exp(-Iter_max/(2*N)))

x0 = [0;  ...   % T_base
      1; ... % c10
      0.11; ...                 % cm
      0.168; ...                 % cp
      200; ... % N
      0]; % intercept

lb = [0,    0,   -5,  -5,  1e-2, 0];
ub = [Inf, Inf,   5,   5,  1e6, inf];

obj = @(x) sse_runtime_model(x, Iter_max, m, p, y);

opts = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 5e4, ...
    'MaxIterations', 2e5);

[xhat, fval, exitflag, output] = fmincon(obj, x0, [], [], [], [], lb, ub, [], opts);

% Unpack
T_base = xhat(1);
c10    = xhat(2);
cm     = xhat(3);
cp     = xhat(4);
N      = xhat(5);

% Predictions + diagnostics
yhat = runtime_model(xhat, Iter_max, m, p);

rmse = sqrt(mean((y - yhat).^2));
r2   = 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);

fprintf('\nNonlinear runtime model fit (fmincon)\n');
fprintf('  exitflag: %d\n', exitflag);
fprintf('  SSE:      %.6g\n', fval);
fprintf('  RMSE:     %.6g s\n', rmse);
fprintf('  R^2:      %.6g\n', r2);
fprintf('  T_base:   %.6g\n', T_base);
fprintf('  c10:      %.6g\n', c10);
fprintf('  cm:       %.6g\n', cm);
fprintf('  cp:       %.6g\n', cp);
fprintf('  N:        %.6g\n', N);

figure('Color','w');
plot(y, yhat, '.', 'MarkerSize', 12); grid on; box on
xlabel('Measured runtime\_s');
ylabel('Predicted runtime\_s');
title('Nonlinear model: measured vs predicted');

%% Local functions
function yhat = runtime_model(x, Iter_max, m, p)
    T_base = x(1);
    c10    = x(2);
    cm     = x(3);
    cp     = x(4);
    N      = x(5);
    inter      = x(6);

    % Stable exponent handling
    expo = cm.*m + cp.*p;
    expo = max(min(expo, 100), -100);           % avoids 10^expo overflow
    amp  = T_base + c10 .* 10.^expo;

    z = Iter_max ./ (2*N);
    z = max(min(z, 700), -700);                 % avoids exp overflow
    rise = 1 - exp(-z);

    yhat = amp .* rise + inter;
end

function J = sse_runtime_model(x, Iter_max, m, p, y)
    yhat = runtime_model(x, Iter_max, m, p);
    r = y - yhat;
    J = r.'*r;
    if ~isfinite(J); J = realmax; end
end
