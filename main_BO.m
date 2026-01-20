clear all; close all; clc;

rng(1)
theta = [3 2 1 0 0.3 0.3 0.3 0 2 2 1];
base = nmpc_init_base([0.001 0.1 0.1]);   % one-time

out1 = nmpc_eval_theta(base, theta);


function base = nmpc_init_base(sigma_y)
%NMPC_INIT_BASE One-time initialisation shared across all theta evaluations.

    if nargin < 1 || isempty(sigma_y)
        sigma_y = [0.001 0.1 0.1];
    end

    tf = 1;

    current_dir = fileparts(mfilename('fullpath'));
    addpath(genpath(current_dir))

    base.ode_opt = odeset('NonNegative', [2 3]);

    get_par
    base.par = par;

    base.dt = 1/60;
    base.tf = tf;

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

    base.N = ceil(base.tf/base.dt) + 1;
    base.T = (0:base.N-1) * base.dt;
    base.tspan = [0 base.dt];

    base.sigma_y = sigma_y;

    base.noise = randn(base.N, base.nx) .* base.sigma_y;
end


function out = nmpc_eval_theta(base, theta)
%NMPC_EVAL_THETA Evaluate one theta using preinitialised base.
%
% Returns one struct with fields:
%   out.case(1) and out.case(2), each containing trajectories + costs.

    cfg = decode_theta(theta, base.nx, base.nu);

    NMPC = NMPC_terminal(base.model, base.nx, base.nu);
    NMPC.optimizer_options.UseParallel = true;

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
    out.T     = base.T;

    for case_id = 1:2
        if case_id == 1
            x0 = [1.0, 10, 0];
        else
            x0 = [1.1, 25, 5];
        end

        uk = zeros(1, base.nu);

        Y      = zeros(base.N, base.nx);
        Y_meas = zeros(base.N, base.nx);
        Y_sp   = repmat(NMPC.xsp(1:base.nx), base.N, 1);
        U      = zeros(base.N, base.nu);

        timer = tic;
        xk = x0;

        for i = 1:base.N
            Y(i,:) = xk;

            yk_meas = xk + base.noise(i,:);
            yk_meas(yk_meas < 0) = 0;
            Y_meas(i,:) = yk_meas;

            uk = NMPC.solve(yk_meas(:)', uk(:)');
            U(i,:) = uk;

            [~, y] = ode45(@(t,x) base.plant(x, uk), base.tspan, xk, base.ode_opt);
            xk = y(end,:);
        end
        runtime_s = toc(timer);

        E = Y - Y_sp;
        SSE = sum(E(:).^2);

        dU  = diff(U, 1, 1);
        SSR = sum(dU(:).^2);

        J = SSE + 1e4 * SSR;

        out.case(case_id).case_id    = case_id;
        out.case(case_id).x0         = x0;

        out.case(case_id).Y          = Y;
        out.case(case_id).Y_meas     = Y_meas;
        out.case(case_id).Y_sp       = Y_sp;
        out.case(case_id).U          = U;

        out.case(case_id).noise      = base.noise;

        out.case(case_id).cost_total = J;
        out.case(case_id).cost_SSE   = SSE;
        out.case(case_id).cost_SSR   = SSR;

        out.case(case_id).runtime_s  = runtime_s;
    end
end


function cfg = decode_theta(theta, nx, nu)
%DECODE_THETA theta = [theta_p, theta_m, q(1:nx), ru(1:nu), rdu(1:nu)]
%
% Horizons:
%   m = theta_m + 1
%   p = theta_p + m

    theta = theta(:).';
    expected_len = 2 + nx + nu + nu;
    if numel(theta) ~= expected_len
        error('theta must have length %d.', expected_len)
    end

    k = 1;
    theta_p = theta(k); k = k + 1;
    theta_m = theta(k); k = k + 1;

    cfg = struct();
    cfg.m = theta_m + 1;
    cfg.p = theta_p + cfg.m;

    q_exp   = theta(k:k+nx-1); k = k + nx;
    ru_exp  = theta(k:k+nu-1); k = k + nu;
    rdu_exp = theta(k:k+nu-1);

    cfg.Q   = diag(10.^q_exp);
    cfg.Ru  = diag(10.^ru_exp);
    cfg.Rdu = diag(10.^rdu_exp);
end


function P = construct_P(LQR_data, Q, R1, R2)
    Sx  = LQR_data.Sx;
    Su  = LQR_data.Su;
    K   = LQR_data.K;
    Acl = LQR_data.Acl;

    Qbar = (Sx.'*Q*Sx) ...
         + (Su - K).'*R1*(Su - K) ...
         + K.'*R2*K;

    P = dlyap(Acl', Qbar);
end
