%% Fit model: y = Tmin + (Tmean - Tmin) * (1 + ap*p.^ep + am*m.^em) .* (1 - exp(-Iter_max/N))
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

q1 = Theta(:,4);  q2 = Theta(:,5);  q3 = Theta(:,6);
ru1 = Theta(:,7); ru2 = Theta(:,8); ru3 = Theta(:,9);
rdu1 = Theta(:,10); rdu2 = Theta(:,11); rdu3 = Theta(:,12);


%% --- Model definition ---
% theta = [Tmean, Tmin, N, ap, ep, am, em]
model = @(theta, X) ...
    theta(2) + ...
    (theta(1)-theta(2)) .* (1 + theta(4).*X(:,2).^theta(5) + theta(6).*X(:,3).^theta(7)) .* ...
    (1 - exp(-X(:,1)./theta(3)));

X = [Iter_max, p, m];

%% --- Initial guess (use your previous fit as a stable start) ---
Tmean0 = 5.8652;
Tmin0 = 0.4066;
N0    = 6.5001;

% Start with mild modulation; exponents ~1
ap0 = 0.0;  ep0 = 1.0;
am0 = 0.0;  em0 = 1.0;

theta0 = [Tmean0, Tmin0, N0, ap0, ep0, am0, em0];

%% --- Bounds (keep the problem well-posed) ---
% Ensure N>0. Exponents >=0 to avoid blow-ups at p=0 or m=0.
lb = [0, 0, 1e-1, 0, 0, 0, 0];
ub = [ Inf,  1e3,  Inf,  Inf, 10,  Inf, 10];

%% --- Fit (nonlinear least squares) ---
opts = optimoptions('lsqcurvefit', ...
    'Display','iter', ...
    'MaxFunctionEvaluations',5e4, ...
    'MaxIterations',2e3);

theta = lsqcurvefit(model, theta0, X, y, lb, ub, opts);

Tmean = theta(1);
Tmin = theta(2);
N    = theta(3);
ap   = theta(4);
ep   = theta(5);
am   = theta(6);
em   = theta(7);

fprintf('\nFit:\n  Tmean = %.6g\n  Tmin = %.6g\n  N    = %.6g\n  ap   = %.6g\n  ep   = %.6g\n  am   = %.6g\n  em   = %.6g\n', ...
    Tmean, Tmin, N, ap, ep, am, em);

% %% --- Plot: raw data, mean per Iter_max, and fitted mean curve at representative (p,m) ---
% % Empirical mean per Iter_max (ignores p,m; shown for reference only)
% [uIter, ~, idx] = unique(Iter_max);
% y_mean = accumarray(idx, y, [], @mean);
% 
% % Choose representative p,m for a 1D fitted curve (medians)
% p_med = median(p);
% m_med = median(m);
% 
% xfit = linspace(min(uIter), max(uIter), 300).';
% Xfit = [xfit, p_med*ones(size(xfit)), m_med*ones(size(xfit))];
% yfit = model(theta, Xfit);
% 
% figure; hold on; grid on
% scatter(Iter_max, y, 12, 'filled')
% plot(uIter, y_mean, 'o-', 'LineWidth', 1.5)
% plot(xfit, yfit, '--', 'LineWidth', 2)
% 
% xlabel('Iter_{max}')
% ylabel('y')
% legend('Raw data','Mean per Iter_{max}','Fit at median (p,m)','Location','best')
% title(sprintf('Fit: Tmean=%.3g Tmin=%.3g N=%.3g ap=%.3g ep=%.3g am=%.3g em=%.3g', ...
%     Tmean, Tmin, N, ap, ep, am, em));

%% --- 4) Goodness of fit: y vs yhat and R^2 ---
yhat = model(theta, X);

% R-squared
SS_res = sum((y - yhat).^2);
SS_tot = sum((y - mean(y)).^2);
R2 = 1 - SS_res/SS_tot;

fprintf('\nGoodness of fit:\n  R^2 = %.6f\n', R2);

% y vs yhat plot
figure; hold on; grid on
scatter(y, yhat, 16, 'filled')

% 45-degree reference line
ymin = min([y; yhat]);
ymax = max([y; yhat]);
plot([ymin ymax], [ymin ymax], 'k--', 'LineWidth', 1.5)

xlabel('Observed y')
ylabel('Predicted \hat{y}')
title(sprintf('y vs \\hat{y} (R^2 = %.4f)', R2))
axis equal
xlim([ymin ymax])
ylim([ymin ymax])
