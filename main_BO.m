%% NMPC simulation runner + logger (external polling or DOE sweep)
%
% This script evaluates an NMPC controller across many parameter vectors
% ("theta") and logs performance metrics for later analysis.
%
% What it produces:
%   1) A CSV file at cfg_run.results_csv with one row per theta run:
%        timestamp, SSE, SSdU, J, runtime_s, theta_1..theta_K
%      where J = SSE + 1e4*SSdU (consistent with simulate_nmpc()).
%   2) A MAT file per run in cfg_run.out_dir named out_<timestamp>.mat,
%      containing the full simulation trajectories and metadata.
%
% Two runtime modes:
%   - mode="external": continuously polls inbox/theta.txt for a new theta
%     vector (last non-empty line), then runs + logs each new signature.
%   - mode="doe": uses a predefined DOE matrix (ThetaDOE) and runs all rows.
%
% Important: "surrogate scaling" for truncated simulations
%   The simulation length is tf = 10*f hours where f = theta(1) clamped to [0,1].
%   If f < 1, the run is intentionally shorter than the "full horizon".
%   The code then estimates full-horizon SSE and SSdU by dividing the partial
%   sums by fitted Chebyshev-based fractions (frac_SSE, frac_SSdU). This is a
%   surrogate / extrapolation, not a physical guarantee.
%
% Operational guidance:
%   - Parallel workers: NumWorkers is set to 8 when USE_PARALLEL=true.
%     Pick a value compatible with your machine (CPU cores, RAM). Too many
%     workers can slow down due to overhead or memory pressure.
%   - Common warnings/errors to watch:
%       * "Failed to read theta" (external mode): theta.txt missing/locked/empty.
%       * "Stale theta signature": theta file has not changed since last run.
%       * Errors from decode_theta: wrong theta length/format.
%       * ODE/NMPC errors: integration failures, optimizer failures, etc.
%   - Lock files:
%       matlab.lock is created/deleted around each run in external mode.
%       NOTE: this script does not check the lock before reading theta; it is
%       mainly a signal to other processes (if they choose to honor it).
%
% Timestamp note:
%   timestamp_compact() actually uses Europe/Oslo timezone.
%
% Reproducibility note:
%   rng(1) fixes randomness. Noise is generated once in base.noise and then
%   reused across all theta evaluations (fair comparisons, but only one noise draw).
%

clear all; close all; clc;

% Add folders and subfolders to path (ensures project functions are visible)
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))

% Cleanup lock files from a previous interrupted run (best-effort)
delete_if_exists('.lock')
delete_if_exists('matlab.lock')
%%
rng(1)

% Do once at the start
USE_PARALLEL = true;

% Configure a parallel pool. If USE_PARALLEL==false, the script ensures no pool.
p = gcp('nocreate');
if USE_PARALLEL
    NumWorkers = 31;   % Choose based on machine capacity (cores/RAM).
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
% Runtime mode
%   mode = "external"  -> theta is read from a txt file (polled)
%   mode = "doe"       -> theta is taken from a predefined matrix
% =========================================================================
cfg_run = struct();
cfg_run.mode              = "external";                 % "external" | "doe"
cfg_run.theta_txt         = fullfile("inbox","theta.txt");
cfg_run.poll_s            = 2.0;                   % pause between polls if theta is stale/unreadable
cfg_run.results_csv       = fullfile("results","results.csv");   % <- CSV summary per theta
cfg_run.out_dir           = fullfile("results");   % saves results/out_<timestamp>.mat
cfg_run.sigma_y           = [0.001 0.1 0.1];        % measurement noise std dev (per state)
cfg_run.NumWorkers        = NumWorkers;

% DOE matrix: each row is one theta (only used when cfg_run.mode == "doe")
ThetaDOE = theta_doe_generator();
cfg_run.ThetaDOE = ThetaDOE;

% One-time initialisation (shared across all theta evaluations):
% Builds model/plant, setpoints, steady-state, LQR data, max horizon arrays,
% a fixed noise realization, and Chebyshev coefficients used for surrogate scaling.
base = nmpc_init_base(cfg_run.sigma_y);

% Ensure output locations exist
ensure_parent_dir(cfg_run.results_csv);
ensure_dir(cfg_run.out_dir);

% Initialise results header if missing (true CSV header)
% theta_len corresponds to theta layout = [f, theta_p, theta_m, log10(Qdiag), log10(Ru), log10(Rdu)]
theta_len = 1 + 2 + base.nx + 2*base.nu;
init_results_csv(cfg_run.results_csv, theta_len);

switch cfg_run.mode
    case "external"
        run_external_theta_loop(cfg_run, base);

    case "doe"
        run_doe(cfg_run, base);

    otherwise
        error('Unknown mode "%s". Use "external" or "doe".', cfg_run.mode);
end

%% FUNCTIONS

% Main run and save method: evaluates one theta, logs CSV metrics, saves .mat payload
function [] = run_and_log(cfg_run, base, theta)
    ts = timestamp_compact(); % Oslo time

    % Core evaluation: runs 2 initial-condition cases and returns aggregated SSE/SSdU estimates.
    out = simulate_nmpc(base, theta);

    SSE       = out.SSE;       % sum of per-case estimated SSE totals (full-horizon estimate)
    SSdU      = out.SSdU;      % sum of per-case estimated SSdU totals (full-horizon estimate)
    runtime_s = out.runtime_s; % sum of per-case wall-clock runtimes

    % Consistent with simulate_nmpc() where J = SSE + 1e4*SSdU per case
    J = SSE + 1e4 * SSdU;

    % Append a single row to results CSV (timestamp + metrics + raw theta values)
    append_results_row(cfg_run.results_csv, ts, SSE, SSdU, J, runtime_s, theta);

    % Save full output structure for offline inspection / plotting
    mat_path = fullfile(cfg_run.out_dir, "out_" + ts + ".mat");
    save(mat_path, "ts", "out", "cfg_run", "base");
end

% =========================================================================
% External theta mode
% =========================================================================
function run_external_theta_loop(cfg_run, base)

    % Simple lock file convention: created before running, deleted after.
    % NOTE: This script does NOT check for lock existence before reading theta.
    lock = 'matlab.lock';
    unlock = @() delete_if_exists(lock);

    % signature prevents re-running the same theta.txt contents repeatedly
    last_signature = "";
    while true
        % Reads last non-empty line of theta.txt and parses floats
        [theta, signature, ok] = read_theta_from_txt(cfg_run.theta_txt);

        if ~ok
            % Display current time in Oslo timezone (useful when running unattended)
            osloTimeZone = 'Europe/Oslo';
            currentTimeOslo = datetime('now', 'TimeZone', osloTimeZone);
            w = ['Failed to read theta. Time: ', datestr(currentTimeOslo)];
            warning(w)
            pause(cfg_run.poll_s);
            continue
        end

        % Stale signature means the file (mtime + last line) is unchanged
        if signature == last_signature
            disp('Stale theta signature')
            pause(cfg_run.poll_s);
            continue
        end

        last_signature = signature;

        % Create a lock file as an external signal "MATLAB is busy running"
        fid = fopen(lock,'w');
        if fid < 0
            warning('Failed to create lock (pwd = %s)', pwd);
            keyboard
        end
        fclose(fid);

        % Run and log. Ensure lock is removed on failure.
        try
            run_and_log(cfg_run, base, theta);
        catch ME
            unlock();
            rethrow(ME);
        end

        unlock();
    end
end

% =========================================================================
% DOE mode
% =========================================================================
function run_doe(cfg_run, base)

    Theta = cfg_run.ThetaDOE;
    if isempty(Theta)
        error("ThetaDOE is empty.");
    end

    for i = 1:size(Theta,1)
        clc
        theta = Theta(i,:);
        run_and_log(cfg_run, base, theta)
    end
end

function ThetaDOE = theta_doe_generator()
%THETA_DOE_GENERATOR DOE for theta with targeted weight variants and horizon sweep.
%
% theta layout (exponents are base-10 logs):
%   theta = [f, theta_p, theta_m, log10(q_diag), log10(ru_diag), log10(rdu_diag)]
%
%   - f in [0,1] controls simulation length tf = 10*f hours.
%   - theta_m and theta_p encode horizons as:
%       m = theta_m + 1
%       p = theta_p + m
%   - q_diag, ru_diag, rdu_diag are converted back via diag(10.^exponent).
%
% NOTE: f_set is currently (1:1)/100 -> only a single value f=0.01 (6 minutes).
%       If you intended a sweep, adjust f_set outside this script (execution unchanged here).

    nx = 3;
    nu = 3;

    f_set = (1 : 1)./100;

    % Baseline diagonal weights (not log-space yet)
    q0   = [10 1 1];
    ru0  = [2 2 1];
    rdu0 = [100 100 10];

    % Variants of state penalty weights (Q diagonal candidates)
    Qdiag_A = [
        q0
        q0 .* [2 1 1]
        q0 .* [1 2 1]
        q0 .* [1 1 2]
    ];

    % Variants of delta-u penalties (Rdu diagonal candidates)
    Rdu_A = [
        rdu0
        10*rdu0
    ];

    % Candidate horizons (m and p). Only pairs with p>m are accepted.
    m_set = [3 6 12];
    p_set = [20 40 60];

    % Count valid (im, ip) pairs
    n_valid = 0;
    for im = m_set
        for ip = p_set
            if ip > im
                n_valid = n_valid + 1;
            end
        end
    end

    % Total DOE size = (#f) * (#valid horizons) * (#Q variants) * (#Rdu variants)
    nTheta = numel(f_set) * n_valid * size(Qdiag_A,1) * size(Rdu_A,1);
    ThetaDOE = zeros(nTheta, 1 + 2 + nx + 2*nu);

    % Ru is kept fixed across DOE (converted to log10 exponents once)
    log_ru0 = log10(ru0);

    row = 0;

    for im = m_set
        for ip = p_set
            if ip > im
                % Encode horizons into theta fields:
                %  theta_m = m-1
                %  theta_p = p-m
                theta_m = im - 1;
                theta_p = ip - im;

                for iQ = 1:size(Qdiag_A,1)
                    q = Qdiag_A(iQ,:);
                    log_q = log10(q);

                    for iRdu = 1:size(Rdu_A,1)
                        rdu = Rdu_A(iRdu,:);
                        log_rdu = log10(rdu);

                        for f = f_set
                            row = row + 1;
                            ThetaDOE(row,:) = [
                                f, ...
                                theta_p, ...
                                theta_m, ...
                                log_q, ...
                                log_ru0, ...
                                log_rdu
                            ];
                        end
                    end
                end
            end
        end
    end

    % Safety: trim if count mismatched for any reason
    if row ~= nTheta
        ThetaDOE = ThetaDOE(1:row,:);
    end

    % Randomize DOE order (so incomplete sweeps are still informative)
    [~, I] = sort(rand(size(ThetaDOE, 1), 1));
    ThetaDOE = ThetaDOE(I, :);
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

function [theta, signature, ok] = read_theta_from_txt(theta_txt)
    % Reads theta from the *last non-empty line* of theta_txt.
    % Returns ok=false if file missing/unreadable/empty/unparseable.

    theta = [];
    signature = "";
    ok = false;

    if ~isfile(theta_txt)
        return
    end

    try
        raw = fileread(theta_txt);
    catch
        return
    end

    lines = splitlines(string(raw));
    lines = strip(lines);
    lines = lines(lines ~= "");

    if isempty(lines)
        return
    end

    last_line = lines(end);
    vals = sscanf(last_line, "%f").';
    if isempty(vals)
        return
    end

    theta = vals;
    % Signature uses file modification time + content of last line.
    signature = compute_signature(theta_txt, last_line);
    ok = true;
end

function sig = compute_signature(theta_txt, last_line)
    % Signature changes if either:
    %   - file modification time changes, or
    %   - last line content changes
    % Used to avoid re-running the same theta repeatedly.

    d = dir(theta_txt);
    if isempty(d)
        sig = "missing";
        return
    end

    sig = string(d.datenum) + "|" + string(last_line);
end

function ts = timestamp_compact()
    % NOTE: Timezone is Europe/Oslo, not UTC.
    t = datetime("now","TimeZone",'Europe/Oslo');
    ts = char(datestr(t, "yyyymmdd_HHMMSS"));
end

function ensure_dir(p)
    if ~isfolder(p)
        mkdir(p);
    end
end

function ensure_parent_dir(filepath)
    [parent,~,~] = fileparts(filepath);
    if parent ~= "" && ~isfolder(parent)
        mkdir(parent);
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
    % Used to extrapolate partial costs from truncated simulations to a full-horizon estimate.
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

function plot_simulation(out, case_id)
    %PLOT_SIMULATION Plot states and inputs for a given simulation case.
    % Utility function (not called by the main script). Use after loading a saved out_*.mat.

    Y   = out.case(case_id).Y;
    Ysp = out.case(case_id).Ysp;
    U   = out.case(case_id).U;
    dt  = out.case(case_id).dt;

    N = size(Y,1);
    T = 0 : dt : (N - 1)*dt;

    figure(1);
    clf

    subplot(3,1,1);
    plot(T(1:N), Y(1:N,1), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
    plot(T(1:N), Ysp(1:N,1), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
    grid on; box on;
    xlabel('Time (h)');
    ylabel('State 1');
    legend('Location','best');
    hold off;

    subplot(3,1,2);
    plot(T(1:N), Y(1:N,2), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
    plot(T(1:N), Ysp(1:N,2), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
    grid on; box on;
    xlabel('Time (h)');
    ylabel('State 2');
    legend('Location','best');
    hold off;

    subplot(3,1,3);
    plot(T(1:N), Y(1:N,3), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
    plot(T(1:N), Ysp(1:N,3), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
    grid on; box on;
    xlabel('Time (h)');
    ylabel('State 3');
    legend('Location','best');
    hold off;

    ax = findall(gcf, 'type', 'axes');
    for j = 1:numel(ax)
        ax(j).FontSize = 15;
        ax(j).XLabel.FontSize = 15;
        ax(j).YLabel.FontSize = 15;
    end

    figure(2);
    clf
    plot(T(1:N), U(1:N,:), 'LineWidth', 2);
    grid on; box on;
    xlabel('Time (h)');
    ylabel('Inputs');
end
