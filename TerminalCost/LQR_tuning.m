clear all; close all; clc
rng(1)

%% Model
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

var0 =[0, 0, 0, 0, 0, 0, 0, 0];
%var0 = [-1.9980    0.0003    1.4849    0.5267   -0.9742    0.0425    0.1074   -0.1175];

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

%% FUNCTIONS

function J = obj_func(var, A, B, model, Ts, Tf, xsp, usp, ode_opt)
    % disp(var)

    % Case 1
    x0 = xsp;
    x0(1) = 1.1;
    x0(2) = 5;
    x0(3) = 1;
    J1 = cost_LQR(var, A, B, model, Ts, Tf, x0, xsp, usp, ode_opt);

    % Case 2
    x0 = xsp;
    x0(2) = x0(2) + 2;
    J2 = cost_LQR(var, A, B, model, Ts, Tf, x0, xsp, usp, ode_opt);

    % Mean
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
    sr = sum(r.^2);
    ssr = sum(sr);
    reg = 1e4*sum(sumsqr(diff(Uode)));
    reg2 = 1e2*sumsqr(var);
    J = ssr + reg + reg2;
    % disp(J)
    % disp(reg2)
end
