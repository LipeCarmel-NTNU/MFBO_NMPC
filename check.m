
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

%%
% Predictors (theta-derived, using log10-parameterisation plus horizons)
Xnames = ["Iter_{max}","m","p","q1","q2","q3","ru1","ru2","ru3","rdu1","rdu2","rdu3"];
X = [Iter_max, m, p, q1, q2, q3, ru1, ru2, ru3, rdu1, rdu2, rdu3];

% Remove any non-finite rows (robust to partial logs)
ok = all(isfinite(X),2) & isfinite(y);
X = X(ok,:);
y = y(ok);
% y = log10(y);
%% Ridge regression with K-fold CV (via lasso with Alpha=0)
K = 10;
[B, FitInfo] = lasso(X, y, ...
    'Alpha', 0.5, ...
    'Standardize', true, ...
    'Intercept', true, ...
    'CV', K);

idx = FitInfo.IndexMinMSE;

beta      = B(:,idx);
intercept = FitInfo.Intercept(idx);
lambda    = FitInfo.Lambda(idx);

yhat = intercept + X*beta;

% Diagnostics
rmse = sqrt(mean((y - yhat).^2));
r2   = 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);

fprintf('Ridge regression (Alpha=0)\n');
fprintf('  K-fold CV: %d\n', K);
fprintf('  Lambda*:  %.6g\n', lambda);
fprintf('  RMSE:     %.6g s\n', rmse);
fprintf('  R^2:      %.6g\n', r2);

% Coefficients table (note: predictors were standardised inside lasso)
coef_tbl = table(Xnames(:), beta, 'VariableNames', {'Predictor','Beta'});
disp(coef_tbl);

%% Plots: runtime_s as a function of each predictor (one figure per predictor)
for j = 1:size(X,2)
    figure('Color','w');
    plot(X(:,j), y, '.', 'MarkerSize', 12); grid on; box on
    xlabel(Xnames(j), 'Interpreter','tex');
    ylabel('runtime\_s');

    % Optional: add linear trend line for visual reference (not a model claim)
    hold on
    pfit = polyfit(X(:,j), y, 1);
    xg = linspace(min(X(:,j)), max(X(:,j)), 200);
    yg = polyval(pfit, xg);
    plot(xg, yg, '-', 'LineWidth', 2);
    hold off

    title(sprintf('runtime\\_s vs %s', Xnames(j)), 'Interpreter','tex');
end
