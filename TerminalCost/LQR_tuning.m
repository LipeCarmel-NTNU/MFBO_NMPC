clear all; close all; clc
rng(1)

%%
get_par

model = @(x, u)  dilution_reduced(0, x, u, par);
ode_opt = odeset('NonNegative', 2:3, 'RelTol', 1e-3, 'AbsTol', 1e-3);

%% Linearization conditions
V = 1;
X = 10;

% Check for controllability
[xss, uss] = find_ss(V, X, par, model, ode_opt);
[A, B] = linearize(xss, uss, model);
controllability = ctrb(A,B);
unco = length(A) - rank(controllability);

%% Construct the LQR

% Targets
Xsp = 15;
Vsp = 1;
[xsp, usp] = find_ss(Vsp, Xsp, par, model, ode_opt);

% Sampling
Ts = 1/60;

% Initial conditions
x0 = xss;
x0(1) = 1.1;
x0(2) = 5;
x0(3) = 1;

Tf = 20;

%% Tuning (log10 weights), diagonal by construction:
% var = [ q(1:n-1), r1(1:m), r2(1:m) ]
% where Q = diag([1, 10.^q]), R1 = diag(10.^r1), R2 = diag(10.^r2)
nx = size(A,1);
nu = size(B,2);
% var0 = zeros(1, (nx-1) + 2*nu);

var0 =[0, 0, -1, 0, -1, 0, 0, 0];
% var0 = [-2.8127  0.5583,  -1 0 -1,   2.5095  -0.4645  -0.9051]; % Initially optimized with one R
%var0 = [-1.9567    0.0003    1.3962    0.4969   -0.8722    0.0568    0.1182   -0.0894]; % var0 was [0, 0, -1, 0, -1, 0, 0, 0]

% Check initial
% J0 = cost_LQR(var0, A, B, model, Ts, Tf, x0, xsp, usp, ode_opt);
J0 = obj_func(var0, A, B, model, Ts, Tf, xsp, usp, ode_opt)

% Optimize
opt_options = optimoptions('fminunc', 'Display','iter-detailed','UseParallel',true, 'MaxFunctionEvaluations', 2000);
[optvar, J] = fminunc(@(var) obj_func(var, A, B, model, Ts, Tf, xsp, usp, ode_opt), var0, opt_options);
disp(optvar)
%% Simulation
tuning_par = optvar;
% tuning_par = var0;

[Ai, Bi] = incremental(A,B,Ts);
[K, Qz, R, N] = build_LQR_full(tuning_par, Ai, Bi, nx, nu);
[Y, T, U] = LQR_simulation(@(t,x,u) model(x, u), Ts, Tf, x0, K, xsp, usp, ode_opt);

t_interp = 0 : Ts : T(end);
X_interp = interp1(T, Y(:, 2), t_interp);
Jx = trapz(t_interp, Xsp - X_interp);

colors = good_colors(3);
figure(2); clf
tiledlayout(3,1,'TileSpacing','tight','Padding','tight')

for i = 1:3
    nexttile
    plot(T, Y(:,i), '-', 'LineWidth', 3, 'Color', colors(i,:)); hold on
    yline(xsp(i), '--', 'LineWidth', 3, 'Color', 'k')
    grid on; box on
    ylabel(sprintf('x_%d', i))
    if i == 1
        title('States')
        legend({'state','setpoint'}, 'Location','best', 'Interpreter','latex')
    end
    if i == 3
        xlabel('Time (h)')
    end
end

set_font_size()

figure(1);
stairs(T, U)
figure(2)
function J = obj_func(var, A, B, model, Ts, Tf, xsp, usp, ode_opt)
    % disp(var)
    x0 = xsp;
    x0(1) = 1.1;
    x0(2) = 5;
    x0(3) = 1;
    J1 = cost_LQR(var, A, B, model, Ts, Tf, x0, xsp, usp, ode_opt);

    x0 = xsp;
    x0(2) = x0(2) + 2;
    J2 = cost_LQR(var, A, B, model, Ts, Tf, x0, xsp, usp, ode_opt);

    J = (J1 + J2)/2;
end
function J = cost_LQR(var, A, B, model, Ts, Tf, x0, xsp, usp, ode_opt)
    nx = size(A,1);
    nu = size(B,2);

    [Ai, Bi] = incremental(A,B,Ts);
    K = build_LQR_full(var, Ai, Bi, nx, nu);

    [Yode, ~, Uode] = LQR_simulation(@(t,x,u) model(x, u), Ts, Tf, x0, K, xsp, usp, ode_opt);

    r = Yode - xsp;
    r = r.*[10 1 1];
    % sr = sum(r.^2);
    sr = sum(abs(r));
    ssr = sum(sr);
    reg = 1e4*sum(sumsqr(diff(Uode)));
    reg2 = 1e2*sumsqr(var);
    J = ssr + reg + reg2;
    % disp(J)
    % disp(reg2)
end

function [Y, T, U] = LQR_simulation(system, Ts, Tf, y0, K, yss, uss, ode_opt)

    % Steps
    num_sim = ceil(Tf/Ts);

    % Dimensions
    n = length(yss);
    m = length(uss);

    % Preallocate
    Y = zeros(num_sim, n);
    U = zeros(num_sim, m);
    T = zeros(num_sim, 1);

    % Initialise
    Y(1,:) = y0(:).';
    T(1)   = 0;

    y_current = y0(:);
    t_current = 0;

    % Initialise previous input at setpoint (piecewise-constant hold)
    u_prev = uss(:);

    U(1,:) = u_prev(:).';

    for k = 2:num_sim
        % Deviations
        x_tilde      = y_current - yss(:);
        u_tilde_prev = u_prev    - uss(:);

        % Augmented deviation state
        z = [x_tilde; u_tilde_prev];

        % Incremental control in deviation coordinates
        du_tilde = -K * z;

        % Absolute input update
        uk = u_prev + du_tilde;

        % Constraints
        uk = max(uk, zeros(m,1));
        uk = min(uk, 0.4*ones(m,1));

        % ODE over the interval with piecewise-constant uk
        ode_current = @(t, y) system(t, y, uk);
        [t_span, y_span] = ode45(ode_current, [t_current, t_current + Ts], y_current, ode_opt);

        % Update state/time
        y_current = y_span(end,:).';
        y_current(y_current < 0) = 0;
        t_current = t_span(end);

        % Store
        Y(k,:) = y_current.';
        U(k,:) = uk.';
        T(k)   = t_current;

        % Update previous input for next increment
        u_prev = uk;
    end
end

function [Ai, Bi] = incremental(A,B,Ts)
    % Assumes C = I, D = 0;

    nx = size(A,1);
    nu = size(B,2);

    sysc = ss(A,B,eye(nx),zeros(nx,nu));
    sysd = c2d(sysc, Ts, 'zoh');
    Ad = sysd.A;  Bd = sysd.B;

    % Augmentation with u_{k-1}
    Ai = [Ad, Bd;
        zeros(nu,nx), eye(nu)];
    Bi = [Bd;
        eye(nu)];
end

function [K, Qz, R, N] = build_LQR_full(log10w, Ai, Bi, nx, nu)
    % Stage cost in (z, du):
    %   x'Qx + u'R1u + du'R2du, with u = u_{k-1} + du
    %
    % Expanded form:
    %   z'Qz z + 2 z' N du + du' R du
    %
    % with Qz = blkdiag(Q, R1), N = [0; R1], R = R1 + R2.

    w = 10.^log10w(:);

    nq = nx - 1; % Q(1,1) = 1

    q  = [1; w(1:nq)];
    r1 = w(nq + (1:nu));
    r2 = w(nq + nu + (1:nu));

    Q  = diag(q);
    R1 = diag(r1);
    R2 = diag(r2);

    Qz = blkdiag(Q, R1);
    R  = R1 + R2;
    N  = [zeros(nx,nu);
        R1];

    try
        K = dlqr(Ai, Bi, Qz, R, N);
    catch ME
        warning('Crashed')
        keyboard
        throw(ME)
    end
end
