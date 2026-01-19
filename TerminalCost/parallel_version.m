%% Monte Carlo validation with parallel computing
clear; close all; clc
rng(1)

get_par
model   = @(x, u) dilution_reduced(0, x, u, par);
ode_opt = odeset('NonNegative', 2:3, 'RelTol', 1e-3, 'AbsTol', 1e-3);

lqr_tuning = [-1.9980 0.0003 1.4849 0.5267 -0.9742 0.0425 0.1074 -0.1175];

Ts = 1/60;
Tf = 40;

N_lin = 100;
N_ic  = 20;
dev_frac = 0.50;

V_rng = [0.5, 1.5];
X_rng = [2, 20];

%%
% Initialize parallel pool with explicit worker count to avoid RAM bottleneck
N_workers = 2;  % Set to 4 to balance speedup vs RAM usage

poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool('Threads', N_workers);
    fprintf('Created new parallel pool with %d workers\n', N_workers);
elseif poolobj.NumWorkers ~= N_workers
    delete(poolobj);
    poolobj = parpool('Threads', N_workers);
    fprintf('Restarted parallel pool with %d workers\n', N_workers);
else
    fprintf('Using existing parallel pool with %d workers\n', poolobj.NumWorkers);
end
%%
nx = 3;
nu = 3;
[Asym, Bsym, X, U] = symbolic_ss(model, nx, nu);
Afun = matlabFunction(Asym, 'Vars', {X, U});
Bfun = matlabFunction(Bsym, 'Vars', {X, U});

% Pre-allocate storage
SSE       = nan(N_lin, N_ic);
E1_gt_End = false(N_lin, N_ic);

% Worst case tracking (collect all, find worst after)
worst_candidates = cell(N_lin, 1);

fprintf('Running Monte Carlo: %d linearisations x %d ICs each = %d sims\n', ...
    N_lin, N_ic, N_lin*N_ic);

tic;
parfor i = 1:N_lin
    % Local storage for this worker
    local_SSE = nan(1, N_ic);
    local_E1_gt_End = false(1, N_ic);
    local_worst = struct('SSE', -inf, 'i_lin', i, 'i_ic', NaN, ...
                         'xss', [], 'uss', [], 'A', [], 'B', [], ...
                         'K', [], 'x0', [], 'T', [], 'Y', [], 'U', []);
    
    % Random operating point
    V = V_rng(1) + (V_rng(2)-V_rng(1))*rand;
    X = X_rng(1) + (X_rng(2)-X_rng(1))*rand;
    
    try
        [xss, uss] = find_ss(V, X, par, model, ode_opt);
        A = Afun(xss', uss');
        B = Bfun(xss', uss');
        nx = size(A,1);
        nu = size(B,2);
    catch ME
        warning('Skipping linearisation %d: %s', i, ME.message);
        continue
    end
    
    expected_len = (nx-1) + 2*nu;
    if numel(lqr_tuning) ~= expected_len
        error('lqr_tuning length mismatch: expected %d, got %d.', expected_len, numel(lqr_tuning));
    end
    
    % Controllability check
    try
        Ctrb = ctrb(A,B);
        unco = size(A,1) - rank(Ctrb);
        if unco > 0
            warning('Skipping linearisation %d (uncontrollable: %d states).', i, unco);
            continue
        end
    catch
        continue
    end
    
    % Build LQR gain
    try
        [Ai, Bi] = incremental(A, B, Ts);
        K = build_LQR_full(lqr_tuning, Ai, Bi, nx, nu);
    catch ME
        warning('Skipping linearisation %d (LQR failed): %s', i, ME.message);
        continue
    end
    
    % Inner loop: multiple ICs (kept serial within each worker)
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
        local_SSE(j) = sum(e2);
        local_E1_gt_End(j) = (e2(1) > e2(end));
        
        % Track worst case for this linearization
        if local_SSE(j) > local_worst.SSE
            local_worst.SSE  = local_SSE(j);
            local_worst.i_ic  = j;
            local_worst.xss = xss(:);
            local_worst.uss = uss(:);
            local_worst.A = A;
            local_worst.B = B;
            local_worst.K = K;
            local_worst.x0 = x0(:);
            local_worst.T = Tode;
            local_worst.Y = Yode;
            local_worst.U = Uode;
        end
    end
    
    % Store results from this linearization
    SSE(i, :) = local_SSE;
    E1_gt_End(i, :) = local_E1_gt_End;
    worst_candidates{i} = local_worst;
end
elapsed = toc;

fprintf('\nParallel execution time: %.2f seconds\n', elapsed);

% Find global worst case
worst = struct('SSE', -inf);
for i = 1:N_lin
    if ~isempty(worst_candidates{i}) && worst_candidates{i}.SSE > worst.SSE
        worst = worst_candidates{i};
    end
end

% Summary statistics
valid = ~isnan(SSE);
n_valid = nnz(valid);
fprintf('\nValid simulations: %d / %d\n', n_valid, N_lin*N_ic);

SSE_vec = SSE(valid);
fprintf('SSE: mean = %.4g, median = %.4g, min = %.4g, max = %.4g\n', ...
    mean(SSE_vec), median(SSE_vec), min(SSE_vec), max(SSE_vec));

dec_mask = not(E1_gt_End & valid);
nonzero = nnz(dec_mask);
fprintf('Runs with e^2(1) <= e^2(end): %d / %d\n', nonzero, n_valid);
if nonzero > 0
    warning('Found e^2(1) > e^2(end)')
end

% Plot worst case
if isfinite(worst.SSE)
    figure(1);
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

%% Helper function
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