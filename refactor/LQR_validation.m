%% Monte Carlo LQR screening over multiple steady states

clear; close all; clc

%%
get_par
model = @(t, x, u) dilution_reduced(t, x, u, par);
ode_opt = odeset('NonNegative', 1:2, 'RelTol', 1e-9, 'AbsTol', 1e-9);

% opt_tuning = [-0.3628    4.5735   -3.9818    1.9878   -2.7478];
opt_tuning = [-0.3628    4.5735   -1.9818    3.9878   -0.7478];

N_ss   = 10;    % number of steady states
N_runs = 10;     % Monte Carlo runs per steady state

Tf  = 20;        % simulation horizon
dt  = 0.05;      % evaluation/sample step
T   = 0:dt:Tf;   % evaluation grid

% How to sample steady states (edit to your feasible region)
V_range = [0.5, 2.0];
X_range = [5.0, 30.0];

% Saturation / feasibility (optional). If you do not want saturation, set to [].
u_min =  0.0*ones(3, 1);
u_max =  0.4*ones(3, 1);

% Random IC perturbation (relative, applied elementwise to xss)
rel_dev_max = 0.50; % up to 50% deviation

rng(1); % reproducibility

%% Storage
all = struct();
all.ss(N_ss) = struct();

%% Main loop over steady states
for iSS = 1:N_ss
    % --- sample a steady state "target" (xsp, usp) via random V, X (edit if needed)
    Vsp = V_range(1) + (V_range(2)-V_range(1))*rand();
    Xsp = X_range(1) + (X_range(2)-X_range(1))*rand();

    % Find steady state and linearise about it
    [xsp, usp] = find_ss(Vsp, Xsp, par, @(x,u) model(0,x,u), ode_opt);
    nx = numel(xsp);
    nu = numel(usp);

    [A, B] = linearize(xsp, usp, @(x,u) model(0,x,u));

    % Controllability sanity check (warn only)
    Cc = ctrb(A,B);
    if rank(Cc) < nx
        error('SS %d: Linearisation not fully controllable (rank %d < %d).', iSS, rank(Cc), nx);
    end

    % --- Design LQR using opt_tuning (LQR weights, NOT evaluation Q,R)
    [Q_lqr, R_lqr] = build_lqr_weights_from_opt_tuning(opt_tuning, nx);

    % LQR gain
    try
        K = lqr(A, B, Q_lqr, R_lqr);
    catch ME
        error('SS %d: lqr failed (%s).', iSS, ME.message);
    end

    % --- Evaluation weights (fixed, distinct from LQR weights)
    Qeval = eye(nx);
    Reval = 0.1*eye(nu);

    % --- Monte Carlo runs for this SS
    J_runs = nan(N_runs,1);
    stage_first = nan(N_runs,1);
    stage_last  = nan(N_runs,1);
    runs = struct();
    runs(N_runs) = struct();

    for r = 1:N_runs
        % Random initial state within +/- rel_dev_max, clipped to non-negative
        rel = (2*rand(nx,1)-1) * rel_dev_max;
        x0  = max(0, xsp(:) .* (1 + rel));

        % Simulate with sampled-data state feedback:
        % u(t_k) = usp - K*(x(t_k)-xsp), held constant on [t_k, t_{k+1})
        [X, U] = simulate_sampled_feedback(model, ode_opt, T, x0, xsp, usp, K, u_min, u_max);

        % Compute e and delta_U on the sampling grid
        e = (X - xsp(:).');          % size: [Nt x nx]
        dU = diff(U - usp(:).');         % size: [Nt x nu]

        % Stage cost per k (on grid)
        Nt = numel(T);
        stage = zeros(Nt,1);
        for k = 1:Nt-1
            ek  = e(k,:);    % 1 x nx
            duk = dU(k,:);   % 1 x nu
            stage(k) = ek*Qeval*ek' + duk*Reval*duk';
        end
        J = sum(stage);
        % Store
        J_runs(r) = J;
        stage_first(r) = stage(1);
        stage_last(r)  = stage(end);

        runs(r).x0 = x0;
        runs(r).X = X;
        runs(r).U = U;
        runs(r).e = e;
        runs(r).dU = dU;
        runs(r).stage = stage;
        runs(r).J = J;
    end

    % Summary for this steady state
    all.ss(iSS).valid = true;
    all.ss(iSS).Vsp = Vsp;
    all.ss(iSS).Xsp = Xsp;
    all.ss(iSS).xsp = xsp;
    all.ss(iSS).usp = usp;
    all.ss(iSS).A = A;
    all.ss(iSS).B = B;
    all.ss(iSS).K = K;
    all.ss(iSS).Q_lqr = Q_lqr;
    all.ss(iSS).R_lqr = R_lqr;

    all.ss(iSS).Qeval = Qeval;
    all.ss(iSS).Reval = Reval;

    all.ss(iSS).J_runs = J_runs;
    all.ss(iSS).J_mean = mean(J_runs, 'omitnan');
    all.ss(iSS).J_max  = max(J_runs, [], 'omitnan');

    all.ss(iSS).stage_first = stage_first;
    all.ss(iSS).stage_last  = stage_last;

    % Flag runs where stage(1) < stage(end)
    all.ss(iSS).increasing_stage_runs = find(stage_first < stage_last);

    all.ss(iSS).runs = runs;
end

%% Identify the worst-performing steady state (by max cost across its 10 runs)
valid_idx = find(arrayfun(@(s) isfield(s,'valid') && s.valid, all.ss));
if isempty(valid_idx)
    error('No valid steady states were simulated (all LQR designs failed).')
end

Jmax_perSS = arrayfun(@(s) s.J_max, all.ss(valid_idx));
[~, worst_loc] = max(Jmax_perSS);
iWorst = valid_idx(worst_loc);

% Identify the worst run within that steady state (by J)
[~, rWorst] = max(all.ss(iWorst).J_runs);

%% Report: worst and flagged cases
fprintf('\n=== Summary over %d steady states (valid: %d) ===\n', N_ss, numel(valid_idx));
fprintf('Worst SS index: %d\n', iWorst);
fprintf('  Vsp=%.4g, Xsp=%.4g\n', all.ss(iWorst).Vsp, all.ss(iWorst).Xsp);
fprintf('  J_mean=%.6g, J_max=%.6g\n', all.ss(iWorst).J_mean, all.ss(iWorst).J_max);
fprintf('  Worst run within worst SS: run %d (J=%.6g)\n', rWorst, all.ss(iWorst).runs(rWorst).J);

% Any cases where stage(1) < stage(end)
flagged_ss = [];
for k = valid_idx(:).'
    if ~isempty(all.ss(k).increasing_stage_runs)
        flagged_ss(end+1) = k;
    end
end

fprintf('\nSteady states with at least one run where stage(1) < stage(end): %d\n', numel(flagged_ss));
if ~isempty(flagged_ss)
    fprintf('Indices: '); fprintf('%d ', flagged_ss); fprintf('\n');
end

%% Visualise the worst-performing one
S = all.ss(iWorst);
R = S.runs(rWorst);

figure('Color','w'); 
plot(T, R.X, 'LineWidth', 1.2);
xlabel('Time'); ylabel('States'); grid on;
title(sprintf('Worst SS %d, run %d: state trajectories', iWorst, rWorst));

figure('Color','w');
plot(T, R.U, 'LineWidth', 1.2);
xlabel('Time'); ylabel('Inputs'); grid on;
title(sprintf('Worst SS %d, run %d: inputs', iWorst, rWorst));

figure('Color','w');
plot(T, R.stage, 'LineWidth', 1.2);
xlabel('Time'); ylabel('Stage cost: e^T Q_e e + dU^T R_e dU'); grid on;
title(sprintf('Worst SS %d, run %d: stage cost (k=1 vs end)', iWorst, rWorst));
xline(T(1)); xline(T(end));

%% If there are flagged cases, plot them compactly (up to first 5 SS)
max_plot = 5;
if ~isempty(flagged_ss)
    nplot = min(max_plot, numel(flagged_ss));
    for ii = 1:nplot
        kSS = flagged_ss(ii);
        rr  = all.ss(kSS).increasing_stage_runs(1); % show first flagged run
        Rf  = all.ss(kSS).runs(rr);

        figure('Color','w');
        plot(T, Rf.stage, 'LineWidth', 1.2);
        xlabel('Time'); ylabel('Stage cost'); grid on;
        title(sprintf('Flagged SS %d, run %d: stage(1)=%.3g < stage(end)=%.3g', ...
            kSS, rr, all.ss(kSS).stage_first(rr), all.ss(kSS).stage_last(rr)));
    end
end

%% ---------------- Local functions ----------------

function [Q_lqr, R_lqr] = build_lqr_weights_from_opt_tuning(var, nx)
    var = 10.^var;
    q = [1 var(1:nx-1)];
    r = var(nx:end);

    Q_lqr = diag(q);
    R_lqr = diag(r);

end

function [X, U] = simulate_sampled_feedback(model, ode_opt, T, x0, xsp, usp, K, u_min, u_max)
% Sampled-data simulation on grid T with ZOH inputs:
% u_k = usp - K*(x_k - xsp), applied over [T(k), T(k+1)).
%
% Outputs:
% X: [Nt x nx] states on grid
% U: [Nt x nu] inputs on grid (u at each sample)

    Nt = numel(T);
    nx = numel(x0);
    nu = numel(usp);

    X = zeros(Nt, nx);
    U = zeros(Nt, nu);

    xk = x0(:);
    X(1,:) = xk.';

    for k = 1:Nt
        % Compute control at sample k
        uk = usp(:) - K*(xk - xsp(:));

        % Optional saturation
        if ~isempty(u_min); uk = max(uk, u_min(:)); end
        if ~isempty(u_max); uk = min(uk, u_max(:)); end

        U(k,:) = uk.';

        if k == Nt
            break
        end

        % Integrate to next sample with ZOH input
        tspan = [T(k), T(k+1)];
        rhs = @(t,x) model(t, x, uk);
        [~, xtraj] = ode45(rhs, tspan, xk, ode_opt);
        xk = xtraj(end,:).';

        % Enforce non-negativity at the grid explicitly (numerical guard)
        xk = max(xk, 0);

        X(k+1,:) = xk.';
    end
end


function [xss, uss] = find_ss(V, X, par, model, ode_opt)
    % Starvation:
    % Solve for low flow rate steady-state
    Fin = X * par.Y_XSinv * par.kd / par.Sin; % approximate solution

    var = [Fin 0]; % solve ss for Fin and S
    var_opt = fsolve(@(var) model([V, X, var(2)], [var(1), 0, var(1)]), var);
    Fin = var_opt(1);
    S = var_opt(2);

    uss = [Fin 0 Fin];
    u = uss;
    TSPAN = [0 10];
    [t, x] = ode45(@(t,x) model(x, u), TSPAN, [V, X, S], ode_opt); % check
    xss = x(end, :);
end
