% fit_surrogates_all.m
% Fit on all surrogate_data_<i>.mat files:
%   SSdU: Cheb5 in x on grid f = (1:N)/N  (length N)
%   SSE : Cheb5 in x on grid f = (1:N)/N  (length N)
%   time: t = (m/10)^alfa*(p/30)^beta*softmax(cheb3(x,c)), on grid f (length N)

clear; close all; clc
rng(1)
addpath(genpath("dependencies"));

set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");
fontSize = 18;
plotColors = good_colors(2);

files = dir("surrogate_data_*.mat");
if isempty(files)
    error('No files matching surrogate_data_*.mat in the current folder.');
end

% -----------------------------
% Aggregate all data
% -----------------------------
x_all     = [];
SSdU_all  = [];
SSE_all   = [];
t_all     = [];
m_all     = [];
p_all     = [];
file_id   = [];

for kf = 1:numel(files)
    S = load(fullfile(files(kf).folder, files(kf).name));
    if ~isfield(S, "out")
        error("File %s does not contain variable 'out'.", files(kf).name);
    end
    out = S.out;

    % Grid
    N  = length(out.case(1).SSE);
    f  = (1:N).'/N;
    x  = 2*f - 1;

    % Decode m,p
    m = out.theta(3) + 1;
    p = out.theta(2) + m;

    % Signals
    time = cumsum(out.case(1).RUNTIME) + cumsum(out.case(2).RUNTIME);

    SSdU = [0; cumsum(out.case(1).SSdU)] + [0; cumsum(out.case(2).SSdU)];
    SSdU = SSdU / SSdU(end);

    SSE  = cumsum(out.case(1).SSE) + cumsum(out.case(2).SSE);
    SSE  = SSE / SSE(end);

    % Sanity checks (your intended construction)
    if numel(SSdU) ~= N
        error("File %s: expected SSdU length N, got %d (N=%d).", files(kf).name, numel(SSdU), N);
    end
    if numel(SSE) ~= N || numel(time) ~= N
        error("File %s: expected SSE and time length N, got SSE=%d, time=%d (N=%d).", ...
            files(kf).name, numel(SSE), numel(time), N);
    end

    % Append
    x_all    = [x_all;    x];
    SSdU_all = [SSdU_all; SSdU];
    SSE_all  = [SSE_all;  SSE];
    t_all    = [t_all;    time];
    m_all    = [m_all;    repmat(m, N, 1)];
    p_all    = [p_all;    repmat(p, N, 1)];
    file_id  = [file_id;  repmat(kf, N, 1)];
end

% -----------------------------
% Optimisation options
% -----------------------------
opts = optimoptions('fminunc', ...
    'Algorithm','quasi-newton', ...
    'Display','off', ...
    'MaxIterations', 2000, ...
    'OptimalityTolerance', 1e-10, ...
    'StepTolerance', 1e-12);

% -----------------------------
% Fit SSdU and SSE with Cheb5
% -----------------------------
c5_0   = zeros(6,1);
c_SSdU = fminunc(@(c) obj_Cheb5(c, x_all, SSdU_all), c5_0, opts);
c_SSE  = fminunc(@(c) obj_Cheb5(c, x_all, SSE_all ), c5_0, opts);

% Goodness-of-fit on the aggregated training set
SSdU_hat_all = clamp01(Cheb5(x_all, c_SSdU));
SSE_hat_all  = clamp01(Cheb5(x_all, c_SSE));
R2_SSdU = compute_r2(SSdU_all, SSdU_hat_all);
R2_SSE  = compute_r2(SSE_all,  SSE_hat_all);

% Check rng(1) results. If non destructive changes are made to the code, then
% the results should be the same within rouding error
c_SSdU_check = [6.442657e-01   4.682368e-01  -1.455242e-01   2.457724e-02   2.198219e-02  -1.598959e-02]';
c_SSE_check = [7.827330e-01   3.771709e-01  -2.433549e-01   1.091471e-01  -2.846259e-02   1.358572e-03]';
exp_SSdU = floor(log10(abs(c_SSdU_check)));
exp_SSE  = floor(log10(abs(c_SSE_check)));

% Absolute rounding tolerance per entry (within significant digits)
significant_digits = 6;
tol_SSdU = 0.5 * 10.^(exp_SSdU - significant_digits);
tol_SSE  = 0.5 * 10.^(exp_SSE  - significant_digits);

if all(abs(c_SSdU - c_SSdU_check) <= tol_SSdU) && ...
   all(abs(c_SSE  - c_SSE_check ) <= tol_SSE)

    fprintf("Chebyshev coefficient check passed (within printed rounding tolerance).\n");

else
    error("Chebyshev coefficient check failed: deviation exceeds printed rounding tolerance.");
end

% -----------------------------
% Fit time model
%   t = (m/10)^alfa * (p/30)^beta * softmax(cheb3(x,c))
% theta = [alfa; beta; c0; c1; c2; c3]
% -----------------------------
theta0_time = [1.0; 1.7; 1e3; 1e3; 1e2; 1e1];

theta_time = fminunc(@(th) obj_time(th, x_all, m_all, p_all, t_all), theta0_time, opts);

alfa = theta_time(1);
beta = theta_time(2);
cT   = theta_time(3:6);

% -----------------------------
% Print coefficients
% -----------------------------
fprintf('Cheb5 coeffs SSdU (c0..c5):\n');
fprintf('% .6e  ', c_SSdU); fprintf('\n\n');

fprintf('Cheb5 coeffs SSE  (c0..c5):\n');
fprintf('% .6e  ', c_SSE); fprintf('\n\n');

fprintf('Time params:\n');
fprintf('alfa = %.6g, beta = %.6g\n', alfa, beta);
fprintf('c_time (c0..c3):\n');
fprintf('% .6e  ', cT); fprintf('\n');

fprintf('\nR^2 on aggregated fit data:\n');
fprintf('R^2(SSdU) = %.6f\n', R2_SSdU);
fprintf('R^2(SSE)  = %.6f\n', R2_SSE);

% Export numerical summary for SSdU/SSE surrogate coefficients.
numDir = fullfile("results", "numerical results");
if ~isfolder(numDir)
    mkdir(numDir);
end
numTxtPath = fullfile(numDir, "surrogate_cheb_coeffs.txt");
fid = fopen(numTxtPath, "w");
if fid == -1
    error("Unable to open file for writing: %s", numTxtPath);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, "=== Surrogate Coefficients (Cheb5) ===\n");
fprintf(fid, "SSdU c0..c5:\n");
fprintf(fid, "% .10e\n", c_SSdU);
fprintf(fid, "\nSSE c0..c5:\n");
fprintf(fid, "% .10e\n", c_SSE);
fprintf(fid, "\nR^2:\n");
fprintf(fid, "R^2(SSdU) = %.10f\n", R2_SSdU);
fprintf(fid, "R^2(SSE)  = %.10f\n", R2_SSE);

%% -----------------------------
% Quick diagnostics plots (optional)
% -----------------------------
S = load(fullfile(files(1).folder, files(1).name));
out = S.out;
N   = length(out.case(1).SSE);

f   = (1:N).'/N;  x = 2*f - 1;

m = out.theta(3) + 1;
p = out.theta(2) + m;

time = cumsum(out.case(1).RUNTIME) + cumsum(out.case(2).RUNTIME);
SSdU = [0; cumsum(out.case(1).SSdU)] + [0; cumsum(out.case(2).SSdU)];
SSdU = SSdU / SSdU(end);
SSE  = cumsum(out.case(1).SSE) + cumsum(out.case(2).SSE);
SSE  = SSE / SSE(end);

SSdU_hat = clamp01(Cheb5(x, c_SSdU));
SSE_hat  = clamp01(Cheb5(x, c_SSE));
t_hat    = time_model(x, m, p, alfa, beta, cT); %#ok<NASGU>

figure; plot(f, SSdU, '-', 'LineWidth', 2.0, 'Color', plotColors(1,:)); hold on
plot(f, SSdU_hat, '--', 'LineWidth', 2.0, 'Color', plotColors(2,:))
xlim([0, 1]);
ylim([0, 1.005]);
xlabel('Fidelity $z$ (dimensionless)'); ylabel('$J_{\mathrm{TV}}$');

grid off; box off
set_font_size(fontSize);
set_fig_size(920, 520);
format_tick(1, 1);

figure; plot(f, SSE, '-', 'LineWidth', 2.0, 'Color', plotColors(1,:)); hold on
plot(f, SSE_hat, '--', 'LineWidth', 2.0, 'Color', plotColors(2,:))
xlim([0, 1]);
ylim([0, 1.005]);
xlabel('Fidelity $z$ (dimensionless)'); ylabel('$J_{\mathrm{track}}$');
grid off; box off
set_font_size(fontSize);
set_fig_size(920, 520);
format_tick(1, 1);

% figure; plot(f, time, 'LineWidth', 1.5); hold on
% plot(f, t_hat, '--', 'LineWidth', 1.5)
% xlabel('f'); ylabel('Time'); grid on
% legend('Data','Time model','Location','best')

% figure; plot(f, time./t_hat); grid on
% xlabel('f'); ylabel('time / time_hat')

% =============================
% Objectives
% =============================
function J = obj_Cheb5(c, x, y)
    r = y - Cheb5(x, c);
    J = r.'*r;
    if ~isfinite(J); J = realmax; end
end

function J = obj_time(th, x, m, p, t)
    alfa = th(1);
    beta = th(2);
    cT   = th(3:6);

    t_hat = time_model(x, m, p, alfa, beta, cT);

    r = log(t) - log(t_hat);
    J = r.'*r;
    if ~isfinite(J); J = realmax; end
end

% =============================
% Models
% =============================
function y = Cheb5(x, c)
    c = c(:);
    if numel(c) ~= 6
        error('Cheb5 expects 6 coefficients (c0..c5).');
    end
    T0 = ones(size(x));
    T1 = x;
    T2 = 2*x.^2 - 1;
    T3 = 4*x.^3 - 3*x;
    T4 = 8*x.^4 - 8*x.^2 + 1;
    T5 = 16*x.^5 - 20*x.^3 + 5*x;

    y = c(1)*T0 + c(2)*T1 + c(3)*T2 + c(4)*T3 + c(5)*T4 + c(6)*T5;
end

function y = Cheb3(x, c)
    c = c(:);
    if numel(c) ~= 4
        error('Cheb3 expects 4 coefficients (c0..c3).');
    end
    T0 = ones(size(x));
    T1 = x;
    T2 = 2*x.^2 - 1;
    T3 = 4*x.^3 - 3*x;

    y = c(1)*T0 + c(2)*T1 + c(3)*T2 + c(4)*T3;
end

function t_hat = time_model(x, m, p, alfa, beta, cT)
    k = 20;
    softplus = @(z) log1p(exp(k*z)) / k;

    g = Cheb3(x, cT);
    s = 1 + softplus(g - 1);

    scale = (m/10).^alfa .* (p/30).^beta;
    t_hat = scale .* s;

    t_hat = max(t_hat, 1);
end

function y = clamp01(y)
    y = min(max(y, 0), 1);
end

function r2 = compute_r2(y, y_hat)
    y = y(:);
    y_hat = y_hat(:);
    sse_res = sum((y - y_hat).^2);
    sse_tot = sum((y - mean(y)).^2);
    if sse_tot <= eps
        r2 = NaN;
    else
        r2 = 1 - sse_res / sse_tot;
    end
end
