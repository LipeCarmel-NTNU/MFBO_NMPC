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

sq  = q1  + q2  + q3;
sru = ru1 + ru2 + ru3;
sdu = rdu1 + rdu2 + rdu3;

%% Clean rows used for fitting (only those needed by the model)
ok = isfinite(y) & isfinite(Iter_max) & isfinite(m) & isfinite(p);
y = y(ok);
y = log10(y);
Iter_max = Iter_max(ok);
m = m(ok);
p = p(ok);

%% Nonlinear model fit with fmincon
% Model:
%   yhat = (T_base + c10*10^(cm*m + cp*p))*(1 - exp(-Iter_max/(2*N)))

x0 = [0;  ...                       % T_base
      1; ...                        % cInter
      0.11; ...                     % cm
      0.168; ...                    % cp
      0; ...                        % N
      0];                           % 

lb = [0,    0,   -5,  -5,  0,     0];
ub = [Inf, Inf,   5,   5,  200, inf];

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

%%
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
plot(10.^y, 10.^yhat, '.', 'MarkerSize', 12); grid on; box on
xlabel('Measured runtime\_s');
ylabel('Predicted runtime\_s');
title('Nonlinear model: measured vs predicted');

%% Residuals vs each predictor (one plot per variable)
r = 10.^y - 10.^yhat;   % residuals in log10(runtime)

vars = {Iter_max, m, p, sq, sru};
names = {'Iter_{max}', 'm', 'p'};

for i = 1:numel(vars)
    figure('Color','w'); hold on; grid on; box on
    plot(vars{i}, r, 'o', 'MarkerSize', 12)
    yline(0, 'k-')
    xlabel(names{i});
    ylabel('Residual  (log_{10}(runtime_s) - log_{10}(\hat{runtime}_s))');
    title(sprintf('Residuals vs %s', names{i}));
end


%% Local functions
function yhat = runtime_model(x, Iter_max, m, p)
    T_base = x(1);
    cInter    = x(2);
    cm     = x(3);
    cp     = x(4);
    % Stable exponent handling
    yhat = T_base + cm.*m + cp.*p;
end

function J = sse_runtime_model(x, Iter_max, m, p, y)

    yhat = runtime_model(x, Iter_max, m, p);
    r = y - yhat;

    % Data fidelity term
    J_data = r.' * r;

    % Small Tikhonov regularisation
    % Do not regularise T_base aggressively
    lambda = 1e-6;
    w = [0, ones(1, numel(x) - 1)];     % weights per parameter
    J_reg = lambda * sum(w(:) .* (x(:).^2));

    J = J_data + J_reg;

    if ~isfinite(J)
        J = realmax;
    end
end
