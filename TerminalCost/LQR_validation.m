%% Monte Carlo validation of fixed LQR tuning (opt_var)
clear; close all; clc
rng(1)

get_par
model   = @(x, u) dilution_reduced(0, x, u, par);
ode_opt = odeset('NonNegative', 2:3, 'RelTol', 1e-3, 'AbsTol', 1e-3);

% Fixed tuning (already in log10-space as expected by build_LQR)
opt_var = [-2.8127  0.5583  2.5095 -0.4645 -0.9051];

% Simulation settings
Ts = 1/60;
Tf = 20;

% Monte Carlo settings
N_lin = 100;      % number of random linearisation points (random steady states)
N_ic  = 10;       % number of random initial conditions per linearisation
dev_frac = 0.50;  % up to 50% deviation from steady-state (non-negative)

% Nominal ranges for random operating points (50% around nominal V=1, X=10)
V_nom = 1.0;
X_nom = 10.0;
V_rng = [0.5, 1.5];
X_rng = [2, 20];

% Storage
SSE      = nan(N_lin, N_ic);
E1_gt_End = false(N_lin, N_ic);

% Store some details for the worst case
worst.SSE = -inf;
worst.i_lin = NaN;
worst.i_ic  = NaN;
worst.xss = [];
worst.uss = [];
worst.A = [];
worst.B = [];
worst.K = [];
worst.x0 = [];
worst.T = [];
worst.Y = [];
worst.U = [];

fprintf('Running Monte Carlo: %d linearisations x %d ICs each = %d sims\n', ...
    N_lin, N_ic, N_lin*N_ic);

for i = 1:N_lin
    disp(i)
    % --- Random operating point -> steady state -> linearisation ---
    V = V_rng(1) + (V_rng(2)-V_rng(1))*rand;
    X = X_rng(1) + (X_rng(2)-X_rng(1))*rand;

    try
        [xss, uss] = find_ss(V, X, par, model, ode_opt);
        [A, B]     = linearize(xss, uss, model);
    catch ME
        warning('Skipping linearisation %d due to failure: %s', i, ME.message);
        continue
    end

    % Optional controllability check (skip uncontrollable)
    try
        Ctrb = ctrb(A,B);
        unco = size(A,1) - rank(Ctrb);
        if unco > 0
            warning('Skipping linearisation %d (uncontrollable: %d uncontrollable states).', i, unco);
            continue
        end
    catch
        % if ctrb/rank fails, continue anyway
    end

    % LQR gain using fixed tuning at this (A,B)
    try
        [K, ~, ~] = build_LQR(opt_var, A, B, Ts);
    catch ME
        warning('Skipping linearisation %d (dlqr/build failed): %s', i, ME.message);
        continue
    end

    % --- Multiple random initial conditions around xss ---
    for j = 1:N_ic
        x0 = sample_ic_nonneg(xss(:), dev_frac);

        try
            [Yode, Tode, Uode] = LQR_simulation(@(t,x,u) model(x,u), Ts, Tf, x0, K, xss, uss, ode_opt);
        catch ME
            warning('Simulation failed at lin %d ic %d: %s', i, j, ME.message);
            continue
        end

        e = Yode - xss(:).';                 % deviation from operating point
        e2 = sum(e.^2, 2);                   % squared error per time step
        SSE(i,j) = sum(e2);                  % sum over time
        E1_gt_End(i,j) = (e2(1) > e2(end)); 

        if SSE(i,j) > worst.SSE
            worst.SSE  = SSE(i,j);
            worst.i_lin = i;
            worst.i_ic  = j;
            worst.xss = xss(:);
            worst.uss = uss(:);
            worst.A = A; worst.B = B;
            worst.K = K;
            worst.x0 = x0(:);
            worst.T = Tode;
            worst.Y = Yode;
            worst.U = Uode;
        end
    end
end

% Summary statistics (ignoring NaNs)
valid = ~isnan(SSE);
n_valid = nnz(valid);
fprintf('\nValid simulations: %d / %d\n', n_valid, N_lin*N_ic);

SSE_vec = SSE(valid);
fprintf('SSE: mean = %.4g, median = %.4g, min = %.4g, max = %.4g\n', ...
    mean(SSE_vec), median(SSE_vec), min(SSE_vec), max(SSE_vec));

% Count of runs where e2(1) > e2(end)
dec_mask = E1_gt_End & valid;
fprintf('Runs with e^2(1) > e^2(end): %d / %d\n', nnz(dec_mask), n_valid);

% List a few "decreasing-error" cases (indices)
[idx_lin, idx_ic] = find(dec_mask);
n_show = min(20, numel(idx_lin));
if n_show > 0
    fprintf('\nFirst %d cases where e^2(1) > e^2(end):\n', n_show);
    for k = 1:n_show
        fprintf('  lin %3d, ic %2d, SSE = %.4g\n', idx_lin(k), idx_ic(k), SSE(idx_lin(k), idx_ic(k)));
    end
end

% Plot worst case trajectory and inputs
if isfinite(worst.SSE)
    figure('Color','w'); tiledlayout(4,1,'TileSpacing','tight','Padding','tight')

    % states
    for s = 1:3
        nexttile
        plot(worst.T, worst.Y(:,s), 'LineWidth', 1.8); hold on
        yline(worst.xss(s), '--', 'LineWidth', 1.2);
        grid on; box on
        ylabel(sprintf('x_%d', s))
        if s == 1
            title(sprintf('Worst case: lin %d, ic %d, SSE = %.4g', worst.i_lin, worst.i_ic, worst.SSE))
            legend({'state','x_{ss}'}, 'Location','best')
        end
    end

    % inputs
    nexttile
    stairs(worst.T, worst.U, 'LineWidth', 1.4);
    grid on; box on
    xlabel('Time (h)')
    ylabel('u')
else
    warning('No successful simulations to plot.')
end

%% ---- helper: sample initial condition with bounded deviation and non-negativity
function x0 = sample_ic_nonneg(xss, frac)
    % Sample each state independently in [xss*(1-frac), xss*(1+frac)],
    % then clip to >= 0. Also handles xss=0 robustly.
    xss = xss(:);
    lb = xss .* (1 - frac);
    ub = xss .* (1 + frac);

    % If any xss is near zero, use an absolute small range instead of zero-width.
    tiny = 1e-3;
    zeroish = abs(xss) < tiny;
    ub(zeroish) = 2;

    x0 = lb + (ub - lb).*rand(size(xss));
    x0 = max(x0, 0);
end

