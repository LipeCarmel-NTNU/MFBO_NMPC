%% NMPC test runs based on timestamp lists (full horizon, no noise)

clear all; close all; clc;

% Add folders and subfolders to path (ensures project functions are visible)
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))

% Cleanup lock files from a previous interrupted run (best-effort)
delete_if_exists('.lock')
delete_if_exists('matlab.lock')

rng(1)

% Do once at the start
USE_PARALLEL = true;

% Configure a parallel pool. If USE_PARALLEL==false, the script ensures no pool.
p = gcp('nocreate');
if USE_PARALLEL
    NumWorkers = 2;   % Choose based on machine capacity (cores/RAM).
    if isempty(p) || p.NumWorkers ~= NumWorkers
        if ~isempty(p)
            delete(p);
        end
        parpool('Processes', NumWorkers);
    end
else
    NumWorkers = 1;
    if ~isempty(p)
        delete(p)
    end
end

% =========================================================================
% List mode: run stored controllers by timestamp
% =========================================================================
cfg_run = struct();
cfg_run.mode              = "list";                 % list only
cfg_run.results_root      = "results";
cfg_run.sigma_y           = [0 0 0];        % no measurement noise
cfg_run.NumWorkers        = NumWorkers;

% One-time initialisation (shared across all theta evaluations):
base = nmpc_init_base(cfg_run.sigma_y);

% Timestamp lists (Pareto-optimal controllers)
t1 = [
    "20260201_232106";
    "20260201_111557";
    "20260201_223337";
    "20260131_151035";
    "20260201_111928";
    "20260201_192807";
    "20260201_120223";
    "20260201_153030";
    "20260201_154920";
    "20260201_193032"
];

t2 = [
    "20260210_134744";
    "20260211_133012";
    "20260210_174745";
    "20260210_153818";
    "20260210_154010";
    "20260211_134235";
    "20260210_171107";
    "20260210_180826";
    "20260210_151703";
    "20260211_122653"
];

run_timestamp_list(cfg_run, base, "run1", t1);
run_timestamp_list(cfg_run, base, "run2", t2);


%% FUNCTIONS

% Run a fixed list of controllers from a previous run folder with f forced to 1.0 and no noise
function run_timestamp_list(cfg_run, base, run_label, timestamps)
    src_dir = fullfile(cfg_run.results_root, run_label);
    out_dir = fullfile(cfg_run.results_root, run_label + "_full_f1_no_noise");
    ensure_dir(out_dir);
    results_csv = fullfile(out_dir, "results_full.csv");
    theta_len = 1 + 2 + base.nx + 2*base.nu;
    init_results_csv(results_csv, theta_len);

    for i = 1:numel(timestamps)
        ts = timestamps(i);
        mat_path = fullfile(src_dir, "out_" + ts + ".mat");
        if ~isfile(mat_path)
            warning("Missing controller file: %s", mat_path);
            continue
        end

        S = load(mat_path, "out");
        if ~isfield(S, "out")
            warning("No 'out' struct in %s", mat_path);
            continue
        end

        theta = S.out.theta;
        theta(1) = 1; % full simulation horizon

        out = simulate_nmpc(base, theta);
        SSE       = out.SSE;
        SSdU      = out.SSdU;
        runtime_s = out.runtime_s;
        J = SSE + 1e4 * SSdU;

        append_results_row(results_csv, char(ts), SSE, SSdU, J, runtime_s, theta);
        save(fullfile(out_dir, "out_full_" + ts + ".mat"), "ts", "out", "theta", "cfg_run", "base");
    end
end

% =========================================================================
% Logging / IO helpers
% =========================================================================
function delete_if_exists(p)
    if exist(p,'file') == 2
        delete(p);
    end
end

function init_results_csv(results_csv, theta_len)

    % If CSV already exists, keep it (do not rewrite header)
    if isfile(results_csv)
        return
    end

    fid = fopen(results_csv, "w");
    if fid < 0
        error("Could not open results file for writing: %s", results_csv);
    end

    % Header: fixed fields + theta_1..theta_K
    fprintf(fid, "timestamp,SSE,SSdU,J,runtime_s");
    for k = 1:theta_len
        fprintf(fid, ",theta_%d", k);
    end
    fprintf(fid, "\n");
    fclose(fid);
end

function append_results_row(results_csv, ts, SSE, SSdU, J, runtime_s, theta)

    fid = fopen(results_csv, "a");
    if fid < 0
        error("Could not open results file for appending: %s", results_csv);
    end

    % Use %.17g for reproducible-ish float text without excessive length
    fprintf(fid, "%s,%.17g,%.17g,%.17g,%.17g", ts, SSE, SSdU, J, runtime_s);

    theta = theta(:).';
    for k = 1:numel(theta)
        fprintf(fid, ",%.17g", theta(k));
    end
    fprintf(fid, "\n");
    fclose(fid);
end

function ensure_dir(p)
    if ~isfolder(p)
        mkdir(p);
    end
end

% =========================================================================
% NMPC code
% =========================================================================
function base = nmpc_init_base(sigma_y)
    %NMPC_INIT_BASE One-time initialisation shared across all theta evaluations.
    %
    % Creates a "base" struct reused across runs:
    %   - dt, max tf, max horizon arrays
    %   - model/plant handles
    %   - setpoints and steady-state
    %   - LQR data for terminal cost construction
    %   - a fixed noise trajectory base.noise
    %   - Chebyshev coefficients used to estimate full-horizon costs from truncated runs

    if nargin < 1 || isempty(sigma_y)
        sigma_y = [0.001 0.1 0.1];
    end

    current_dir = fileparts(mfilename('fullpath'));
    addpath(genpath(current_dir))

    % Enforce non-negativity for states 2 and 3 during ODE integration.
    base.ode_opt = odeset('NonNegative', [2 3]);

    % Loads/sets parameters into variable "par" (then stored in base.par).
    get_par
    base.par = par;

    base.dt = 1/60;   % hours (1 minute sampling)

    base.tf_max = 10;                          % hours (corresponds to f=1)
    base.N_max  = ceil(base.tf_max/base.dt) + 1;
    base.T_max  = (0:base.N_max-1) * base.dt;
    base.tspan  = [0 base.dt];                 % one-step integration window

    % Model/plant (here plant = model)
    % dilution_reduced(t, x, u, par) is expected to return dx/dt (or next-state derivative)
    base.model = @(x, u) dilution_reduced(0, x, u(:)', base.par);
    base.plant = base.model;

    base.nx = 3;
    base.nu = 3;

    % Setpoints (interpreted by find_ss and used as NMPC targets)
    base.V_sp = 1;
    base.X_sp = 20;

    % find_ss is expected to compute steady-state (xss, uss) for given setpoints.
    [xss, uss] = find_ss(base.V_sp, base.X_sp, base.par, base.model, base.ode_opt);
    base.xsp = xss;
    base.usp = uss;

    % LQR_data.mat is expected to provide data needed to build terminal cost P.
    S = load('LQR_data.mat', 'LQR_data');
    base.LQR_data = S.LQR_data;

    base.sigma_y = sigma_y;

    % Fixed noise realization used for all theta evaluations (reproducible / comparable)
    base.noise = randn(base.N_max, base.nx) .* base.sigma_y;

    base.optimizer_max_iter = 100;

    % Chebyshev fraction models (c0..c5)
    % These approximate the fraction of total SSE / SSdU accumulated by time f.
    % Used to extrapolate partial costs from truncated runs to a full-horizon estimate.
    base.cheb_c_SSdU = [ ...
         6.442657e-01 ...
         4.682368e-01 ...
        -1.455242e-01 ...
         2.457724e-02 ...
         2.198219e-02 ...
        -1.598959e-02 ...
    ].';

    base.cheb_c_SSE = [ ...
         7.827330e-01 ...
         3.771709e-01 ...
        -2.433549e-01 ...
         1.091471e-01 ...
        -2.846259e-02 ...
         1.358572e-03 ...
    ].';
end

function out = simulate_nmpc(base, theta)
    %SIMULATE_NMPC Evaluate one theta using preinitialised base.
    %
    % Uses per-run tf = 10*f (hours), where f is theta(1) clamped to [0,1].
    % If f < 1, the run is shorter than the base.tf_max (=10h).
    %
    % Surrogate extrapolation of totals:
    %   SSE_total_est  = SSE_partial / frac_SSE(f)
    %   SSdU_total_est = SSdU_partial / frac_SSdU(f)
    % where frac_* is Cheb5(2*f-1, c_*) clamped and lower-bounded (>=0.01).

    cfg = decode_theta(theta, base.nx, base.nu);
    disp('Run cfg:')
    disp(cfg)
    
    % Determine whether a parallel pool exists. If no pool, disable UseParallel in optimizer.
    pool_empty = isempty(gcp('nocreate'));

    tf = 10 * cfg.f;                      % hours
    N = ceil(tf/base.dt) + 1;

    % Use precomputed max horizon arrays; slice to run horizon
    T     = base.T_max(1:N);
    noise = base.noise(1:N, :);

    % Compute surrogate fractions based on f (map f in [0,1] to x in [-1,1] via 2*f-1)
    frac_SSE  = min(Cheb5(2*cfg.f - 1, base.cheb_c_SSE),  1);
    frac_SSdU = min(Cheb5(2*cfg.f - 1, base.cheb_c_SSdU), 1);

    % Prevent extreme blow-ups from tiny fractions (floor at 0.01 => at most 100x scaling)
    frac_SSE  = max(frac_SSE,  0.01);
    frac_SSdU = max(frac_SSdU, 0.01);

    % NMPC_terminal is expected to construct an NMPC controller object configured for:
    %   - model, state dimension, input dimension
    %   - .solve(y_meas, u_prev): returns control action
    %   - .constraints(): applies bounds/constraints to the optimization
    NMPC = NMPC_terminal(base.model, base.nx, base.nu);
    NMPC.optimizer_options.MaxIterations = base.optimizer_max_iter;
    NMPC.optimizer_options.UseParallel = ~pool_empty;
    NMPC.optimizer_options.Display = 'final-detailed';

    NMPC.Ts  = base.dt;
    NMPC.p   = cfg.p;
    NMPC.m   = cfg.m;

    NMPC.Q   = cfg.Q;
    NMPC.Ru  = cfg.Ru;
    NMPC.Rdu = cfg.Rdu;

    NMPC.constraints();

    NMPC.xsp = base.xsp;
    NMPC.usp = base.usp;

    % construct_P is expected to produce a terminal cost matrix P from LQR_data and the weights.
    NMPC.P = construct_P(base.LQR_data, NMPC.Q, NMPC.Ru, NMPC.Rdu);

    out = struct();
    out.theta = theta(:).';
    out.cfg   = cfg;
    out.tf    = tf;
    out.N     = N;
    out.T     = T;

    % Aggregated (across two cases) *estimated full-horizon* totals
    out.SSE = 0;
    out.SSdU = 0;
    out.runtime_s = 0;

    % Evaluate theta on two different initial conditions to improve robustness
    for case_id = 1:2
        if case_id == 1
            x0 = [1.0, 10, 0];
        else
            x0 = [1.1, 25, 5];
        end

        % Initial control guess to warm-start NMPC solve
        uk = zeros(1, base.nu);

        % Preallocate logging arrays
        Y       = zeros(N, base.nx);
        Y_meas  = zeros(N, base.nx);
        Ysp     = repmat(NMPC.xsp(1:base.nx), N, 1);
        U       = zeros(N, base.nu);
        RUNTIME = zeros(N, 1);

        xk = x0;

        timer = tic;

        for i = 1:N
            iter_timer = tic;

            % True plant state
            Y(i,:) = xk;

            % Measurement = state + Gaussian noise; then clipped to nonnegative
            yk_meas = xk + noise(i,:);
            yk_meas(yk_meas < 0) = 0;
            Y_meas(i,:) = yk_meas;

            % NMPC solve: uses measured output and previous control for warm start
            uk = NMPC.solve(yk_meas(:)', uk(:)');
            U(i,:) = uk;

            % Verbose printing (can slow down long runs; kept as-is by request)
            fprintf('\nMeasurement: \n')
            disp(yk_meas)

            fprintf('\nControl action: \n')
            disp(uk)

            % Plant simulation: one-step integrate using ode45
            if N > 1 && i < N
                [~, y] = ode45(@(t,x) base.plant(x, uk), base.tspan, xk, base.ode_opt);
                xk = y(end,:);
            end

            RUNTIME(i) = toc(iter_timer);

            % Progress estimate (elapsed vs estimated total)
            elapsed_s = toc(timer);
            elapsed_min = elapsed_s/60;
            progress = i/max(N,1);
            total_min = elapsed_min/max(progress, eps);
            remaining_min = total_min - elapsed_min;

            fprintf('Case %d: %.1f %% | elapsed %.2f min | total %.2f min | left %.2f min\n', ...
                case_id, 100*progress, elapsed_min, total_min, remaining_min);
        end

        runtime_s = toc(timer);

        % Tracking error relative to setpoint; state 1 is weighted 10x (E.*[10 1 1])
        E = Y - Ysp;
        E = E.*[10 1 1];

        % Partial SSE accumulated over this truncated horizon
        SSE_partial_vec = sum(E.^2, 2);

        % Partial SSdU accumulated over this truncated horizon (sum of squared input increments)
        dU = diff(U, 1, 1);
        SSdU_partial_vec = sum(dU.^2, 2);

        SSE_partial  = sum(SSE_partial_vec);
        SSdU_partial = sum(SSdU_partial_vec);

        % Surrogate extrapolation to estimated full-horizon totals
        SSE_est  = SSE_partial  / frac_SSE;
        SSdU_est = SSdU_partial / frac_SSdU;

        J_partial = SSE_partial + 1e4 * SSdU_partial;
        J_est     = SSE_est     + 1e4 * SSdU_est;

        % Store full trajectories and both partial + estimated totals for this case
        out.case(case_id).case_id    = case_id;
        out.case(case_id).x0         = x0;

        out.case(case_id).Y          = Y;
        out.case(case_id).Y_meas     = Y_meas;
        out.case(case_id).Ysp        = Ysp;
        out.case(case_id).U          = U;

        out.case(case_id).noise      = noise;
        out.case(case_id).dt         = base.dt;
        out.case(case_id).tf         = tf;

        % Existing fields hold the approximated (full-horizon) totals
        out.case(case_id).SSE        = SSE_est;
        out.case(case_id).SSdU       = SSdU_est;
        out.case(case_id).cost_total = J_est;

        % Partial values from the truncated simulation
        out.case(case_id).partial_SSE        = SSE_partial_vec;
        out.case(case_id).partial_SSdU       = SSdU_partial_vec;
        out.case(case_id).partial_cost_total = J_partial;

        out.case(case_id).RUNTIME    = RUNTIME;
        out.case(case_id).runtime_s  = runtime_s;

        out.case(case_id).frac_SSE   = frac_SSE;
        out.case(case_id).frac_SSdU  = frac_SSdU;

        % Aggregate across cases
        out.SSE       = out.SSE  + SSE_est;
        out.SSdU      = out.SSdU + SSdU_est;
        out.runtime_s = out.runtime_s + runtime_s;
    end
end

function cfg = decode_theta(theta, nx, nu)
    %DECODE_THETA Decode theta vector into NMPC configuration.
    %
    % theta = [f, theta_p, theta_m, q_exp(1:nx), ru_exp(1:nu), rdu_exp(1:nu)]
    % where q_exp, ru_exp, rdu_exp are log10(diagonal entries).
    %
    % Horizons are encoded as:
    %   m = theta_m + 1
    %   p = theta_p + m
    %
    % f is clamped into [0,1] (safety against invalid external inputs).

    theta = theta(:).';
    expected_len = 1 + 2 + nx + nu + nu;
    if numel(theta) ~= expected_len
        error('theta must have length %d.', expected_len)
    end

    k = 1;
    cfg = struct();

    cfg.f = theta(k); k = k + 1;

    theta_p = theta(k); k = k + 1;
    theta_m = theta(k); k = k + 1;

    cfg.m = theta_m + 1;
    cfg.p = theta_p + cfg.m;

    q_exp   = theta(k:k+nx-1); k = k + nx;
    ru_exp  = theta(k:k+nu-1); k = k + nu;
    rdu_exp = theta(k:k+nu-1);

    cfg.Q   = diag(10.^q_exp);
    cfg.Ru  = diag(10.^ru_exp);
    cfg.Rdu = diag(10.^rdu_exp);

    cfg.f = max(0, min(1, cfg.f));
end

function y = Cheb5(x, c)
    %CHEB5 Evaluate Chebyshev polynomial of degree 5:
    % y = sum_{k=0..5} c_k * T_k(x), where x is expected in [-1,1].
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



