clear all; close all; clc;

% Add folders and subfolders to path
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))

delete_if_exists('.lock')
delete_if_exists('matlab.lock')
%%
rng(1)

% Do once at the start
USE_PARALLEL = true;

p = gcp('nocreate');
if USE_PARALLEL
    NumWorkers = 8;
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
cfg_run.mode              = "doe";                 % "external" | "doe"
cfg_run.theta_txt         = fullfile("inbox","theta.txt");
cfg_run.poll_s            = 2.0;                   % pause between polls if theta is stale
cfg_run.results_csv       = fullfile("results","results.csv");   % <- CSV
cfg_run.out_dir           = fullfile("results");   % saves results/out_<timestamp>.mat
cfg_run.sigma_y           = [0.001 0.1 0.1];
cfg_run.NumWorkers        = NumWorkers;

% DOE matrix: each row is one theta (only used when cfg_run.mode == "doe")
ThetaDOE = theta_doe_generator_v2();
cfg_run.ThetaDOE = ThetaDOE;

% One-time initialisation (shared across all theta evaluations)
base = nmpc_init_base(cfg_run.sigma_y);

% Ensure output locations exist
ensure_parent_dir(cfg_run.results_csv);
ensure_dir(cfg_run.out_dir);

% Initialise results header if missing (true CSV header)
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

% Main run and save method
function [] = run_and_log(cfg_run, base, theta)
    ts = timestamp_utc_compact();

    out = simulate_nmpc(base, theta);

    SSE       = out.SSE;
    SSdU      = out.SSdU;
    runtime_s = out.runtime_s;

    % Consistent with simulate_nmpc() where J = SSE + 1e4*SSdU per case
    J = SSE + 1e4 * SSdU;

    append_results_row(cfg_run.results_csv, ts, SSE, SSdU, J, runtime_s, theta);

    mat_path = fullfile(cfg_run.out_dir, "out_" + ts + ".mat");
    save(mat_path, "ts", "out", "cfg_run", "base");
end

% =========================================================================
% External theta mode
% =========================================================================
function run_external_theta_loop(cfg_run, base)

    lock = 'matlab.lock';
    unlock = @() delete_if_exists(lock);

    last_signature = "";
    while true
        [theta, signature, ok] = read_theta_from_txt(cfg_run.theta_txt);

        if ~ok
            warning('Failed to read theta')
            pause(cfg_run.poll_s);
            continue
        end

        if signature == last_signature
            disp('Stale theta signature')
            pause(cfg_run.poll_s);
            continue
        end

        last_signature = signature;

        fid = fopen(lock,'w');
        if fid < 0
            warning('Failed to create lock (pwd = %s)', pwd);
            keyboard
        end
        fclose(fid);

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

function ThetaDOE = theta_doe_generator_v2()
%THETA_DOE_GENERATOR_V2 DOE for theta with targeted weight variants and horizon sweep.
%
% theta layout:
%   theta = [f, theta_p, theta_m, log10(q_diag), log10(ru_diag), log10(rdu_diag)]
%
% f in [0,1], and simulation runs with tf = 10*f (hours).
%
% Horizons:
%   m = theta_m + 1
%   p = theta_p + m
%
% Nominal:
%   Q  = diag([10 1 1])
%   Ru = diag([2 2 1])
%   Rdu = diag([100 100 10])
%
% Horizon set:
%   m_set = [3 6 12]
%   p_set = [20 40 60]
%   only ip > im
%
% f set:
%   f_set = [0.1 0.2 0.4 0.6 0.8 1.0]

    nx = 3;
    nu = 3;

    f_set = (1 : 10)./100;

    q0   = [10 1 1];
    ru0  = [2 2 1];
    rdu0 = [100 100 10];

    Qdiag_A = [
        q0
        q0 .* [2 1 1]
        q0 .* [1 2 1]
        q0 .* [1 1 2]
    ];

    Rdu_A = [
        rdu0
        10*rdu0
    ];

    m_set = [3 6 12];
    p_set = [20 40 60];

    n_valid = 0;
    for im = m_set
        for ip = p_set
            if ip > im
                n_valid = n_valid + 1;
            end
        end
    end

    nTheta = numel(f_set) * n_valid * size(Qdiag_A,1) * size(Rdu_A,1);
    ThetaDOE = zeros(nTheta, 1 + 2 + nx + 2*nu);

    log_ru0 = log10(ru0);

    row = 0;

    for im = m_set
        for ip = p_set
            if ip > im
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

    if row ~= nTheta
        ThetaDOE = ThetaDOE(1:row,:);
    end

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

    if isfile(results_csv)
        return
    end

    fid = fopen(results_csv, "w");
    if fid < 0
        error("Could not open results file for writing: %s", results_csv);
    end

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

    fprintf(fid, "%s,%.17g,%.17g,%.17g,%.17g", ts, SSE, SSdU, J, runtime_s);

    theta = theta(:).';
    for k = 1:numel(theta)
        fprintf(fid, ",%.17g", theta(k));
    end
    fprintf(fid, "\n");
    fclose(fid);
end

function [theta, signature, ok] = read_theta_from_txt(theta_txt)

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
    signature = compute_signature(theta_txt, last_line);
    ok = true;
end

function sig = compute_signature(theta_txt, last_line)

    d = dir(theta_txt);
    if isempty(d)
        sig = "missing";
        return
    end

    sig = string(d.datenum) + "|" + string(last_line);
end

function ts = timestamp_utc_compact()

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

    if nargin < 1 || isempty(sigma_y)
        sigma_y = [0.001 0.1 0.1];
    end

    current_dir = fileparts(mfilename('fullpath'));
    addpath(genpath(current_dir))

    base.ode_opt = odeset('NonNegative', [2 3]);

    get_par
    base.par = par;

    base.dt = 1/60;

    base.tf_max = 10;                          % hours (f=1)
    base.N_max  = ceil(base.tf_max/base.dt) + 1;
    base.T_max  = (0:base.N_max-1) * base.dt;
    base.tspan  = [0 base.dt];

    base.model = @(x, u) dilution_reduced(0, x, u(:)', base.par);
    base.plant = base.model;

    base.nx = 3;
    base.nu = 3;

    base.V_sp = 1;
    base.X_sp = 20;

    [xss, uss] = find_ss(base.V_sp, base.X_sp, base.par, base.model, base.ode_opt);
    base.xsp = xss;
    base.usp = uss;

    S = load('LQR_data.mat', 'LQR_data');
    base.LQR_data = S.LQR_data;

    base.sigma_y = sigma_y;

    base.noise = randn(base.N_max, base.nx) .* base.sigma_y;

    base.optimizer_max_iter = 100;
end

function out = simulate_nmpc(base, theta)
    %SIMULATE_NMPC Evaluate one theta using preinitialised base.
    % Uses per-run tf = 10*f, where f is theta(1).

    cfg = decode_theta(theta, base.nx, base.nu);
    pool_empty = isempty(gcp('nocreate'));

    tf = 10 * cfg.f;                      % hours
    N = ceil(tf/base.dt) + 1;

    T     = base.T_max(1:N);
    noise = base.noise(1:N, :);

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

    NMPC.P = construct_P(base.LQR_data, NMPC.Q, NMPC.Ru, NMPC.Rdu);

    out = struct();
    out.theta = theta(:).';
    out.cfg   = cfg;
    out.tf    = tf;
    out.N     = N;
    out.T     = T;

    out.SSE = 0;
    out.SSdU = 0;
    out.runtime_s = 0;

    for case_id = 1:2
        if case_id == 1
            x0 = [1.0, 10, 0];
        else
            x0 = [1.1, 25, 5];
        end

        uk = zeros(1, base.nu);

        Y       = zeros(N, base.nx);
        Y_meas  = zeros(N, base.nx);
        Ysp     = repmat(NMPC.xsp(1:base.nx), N, 1);
        U       = zeros(N, base.nu);
        RUNTIME = zeros(N, 1);

        xk = x0;

        timer = tic;

        for i = 1:N
            iter_timer = tic;

            Y(i,:) = xk;

            yk_meas = xk + noise(i,:);
            yk_meas(yk_meas < 0) = 0;
            Y_meas(i,:) = yk_meas;

            uk = NMPC.solve(yk_meas(:)', uk(:)');
            U(i,:) = uk;

            fprintf('\nMeasurement: \n')
            disp(yk_meas)

            fprintf('\nControl action: \n')
            disp(uk)

            if N > 1 && i < N
                [~, y] = ode45(@(t,x) base.plant(x, uk), base.tspan, xk, base.ode_opt);
                xk = y(end,:);
            end

            RUNTIME(i) = toc(iter_timer);

            elapsed_s = toc(timer);
            elapsed_min = elapsed_s/60;
            progress = i/max(N,1);
            total_min = elapsed_min/max(progress, eps);
            remaining_min = total_min - elapsed_min;

            fprintf('Case %d: %.1f %% | elapsed %.2f min | total %.2f min | left %.2f min\n', ...
                case_id, 100*progress, elapsed_min, total_min, remaining_min);
        end

        runtime_s = toc(timer);

        E = Y - Ysp;
        E = E.*[10 1 1];
        SSE = sum(E.^2, 2);

        dU   = diff(U, 1, 1);
        SSdU = sum(dU.^2, 2);

        J = sum(SSE) + 1e4 * sum(SSdU);

        out.case(case_id).case_id    = case_id;
        out.case(case_id).x0         = x0;

        out.case(case_id).Y          = Y;
        out.case(case_id).Y_meas     = Y_meas;
        out.case(case_id).Ysp        = Ysp;
        out.case(case_id).U          = U;

        out.case(case_id).noise      = noise;
        out.case(case_id).dt         = base.dt;
        out.case(case_id).tf         = tf;

        out.case(case_id).cost_total = J;

        out.case(case_id).SSE        = SSE;
        out.case(case_id).SSdU       = SSdU;
        out.case(case_id).RUNTIME    = RUNTIME;

        out.case(case_id).runtime_s  = runtime_s;

        out.SSE       = out.SSE  + sum(SSE);
        out.SSdU      = out.SSdU + sum(SSdU);
        out.runtime_s = out.runtime_s + runtime_s;
    end
end

function cfg = decode_theta(theta, nx, nu)
    %DECODE_THETA theta = [f, theta_p, theta_m, q(1:nx), ru(1:nu), rdu(1:nu)]
    % f in [0,1], and tf = 10*f (hours)

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

function plot_simulation(out, case_id)
    %PLOT_SIMULATION Plot states and inputs for a given simulation case.

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
