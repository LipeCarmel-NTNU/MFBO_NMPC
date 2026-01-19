%% Monte Carlo validation of fixed incremental (xu,du) LQR tuning (log10 weights)
clear; close all; clc
rng(1)

get_par
model   = @(x, u) dilution_reduced(0, x, u, par);
ode_opt = odeset('NonNegative', 2:3, 'RelTol', 1e-3, 'AbsTol', 1e-3);

% optimal tuning (log10w = [q(1:nx-1), r1(1:nu), r2(1:nu)])
lqr_tuning = [-1.9980    0.0003    1.4849    0.5267   -0.9742    0.0425    0.1074   -0.1175];

% Simulation settings
Ts = 1/60;
Tf = 40;

% Monte Carlo settings
N_lin = 100;      % number of random linearisation points (random steady states)
N_ic  = 20;       % number of random initial conditions per linearisation
dev_frac = 0.50;  % up to 50% deviation from steady-state (non-negative)

% Nominal ranges for random operating points
V_rng = [0.5, 1.5];
X_rng = [2, 20];

% Storage
SSE       = nan(N_lin, N_ic);
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
        nx = size(A,1);
        nu = size(B,2);
    catch ME
        warning('Skipping linearisation %d due to failure: %s', i, ME.message);
        continue
    end

    % Validate tuning length against (nx-1)+2*nu
    expected_len = (nx-1) + 2*nu;
    if numel(lqr_tuning) ~= expected_len
        error('lqr_tuning length mismatch: expected %d (= (nx-1)+2*nu), got %d.', expected_len, numel(lqr_tuning));
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
    end

    % Build incremental model and LQR gain using fixed tuning at this (A,B)
    try
        [Ai, Bi] = incremental(A, B, Ts);
        K = build_LQR_full(lqr_tuning, Ai, Bi, nx, nu);
    catch ME
        warning('Skipping linearisation %d (build/dlqr failed): %s', i, ME.message);
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

        e  = Yode - xss(:).';
        e2 = sum(e.^2, 2);
        SSE(i,j) = sum(e2);
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

            if isfinite(worst.SSE)
                figure(1); pause(2)
                clf;
                tiledlayout(4,1,'TileSpacing','tight','Padding','tight')
            
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
            
                nexttile
                stairs(worst.T, worst.U, 'LineWidth', 1.4);
                grid on; box on
                xlabel('Time (h)')
                ylabel('u')
            end
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
dec_mask = not(E1_gt_End & valid);
nonzero = nnz(dec_mask);
fprintf('Runs with e^2(1) <= e^2(end): %d / %d\n', nonzero, n_valid);
if nonzero > 0
    warning('Found e^2(1) > e^2(end)')
end

%% ---- helper: sample initial condition with bounded deviation and non-negativity
function x0 = sample_ic_nonneg(xss, frac)
    xss = xss(:);
    lb = xss .* (1 - frac);
    ub = xss .* (1 + frac);

    tiny = 1e-3;
    zeroish = abs(xss) < tiny;
    ub(zeroish) = 2;

    x0 = lb + (ub - lb).*rand(size(xss));
    x0 = max(x0, 0);
end