%% Plot y vs Iter_max, mean per Iter_max, and fit exponential saturation model
% y_mean(Iter) = (Tmax - Tmin)*(1 - exp(-Iter/N)) + Tmin

%clear; close all; clc

%% --- INPUT DATA (assumed already in workspace) ---
% Iter_max : integer vector
% y        : response vector

Iter_max = Iter_max(:);
y        = y(:);

%% --- 1) Mean y for each unique Iter_max ---
[uIter, ~, idx] = unique(Iter_max);
y_mean = accumarray(idx, y, [], @mean);

%% --- 2) Fit exponential saturation model ---
model = @(p,x) (p(1) - p(2)) .* (1 - exp(-x ./ p(3))) + p(2);
% p = [Tmax, Tmin, N]

Tmin0 = min(y_mean);
Tmax0 = max(y_mean);
N0    = mean(uIter);

p0 = [Tmax0, Tmin0, N0];

opts = optimoptions('lsqcurvefit','Display','off');
p = lsqcurvefit(model, p0, uIter, y_mean, [], [], opts);

Tmax = p(1);
Tmin = p(2);
N    = p(3);

%% --- 3) Plot raw data, means, and fitted curve ---
xfit = linspace(min(uIter), max(uIter), 300);
yfit = model(p, xfit);

figure; hold on; grid on
scatter(Iter_max, y, 12, 'filled')
plot(uIter, y_mean, 'o-', 'LineWidth', 1.5)
plot(xfit, yfit, '--', 'LineWidth', 2)

xlabel('Iter_{max}')
ylabel('y')
legend('Raw data','Mean per Iter_{max}','Exponential fit','Location','best')
title(sprintf('Fit: Tmax=%.3g, Tmin=%.3g, N=%.3g', Tmax, Tmin, N))
