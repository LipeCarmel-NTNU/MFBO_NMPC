%% Terminal cost for incremental (delta-u) LQR with x'Qx + u'R1u + du'R2du
clear; close all; clc
rng(1)

get_par
model   = @(x,u) dilution_reduced(0, x, u, par);
ode_opt = odeset('NonNegative', 2:3, 'RelTol', 1e-3, 'AbsTol', 1e-3);

% Setpoint used to compute (xss, uss)
Vsp = 1;
Xsp = 20;

% Sample time
Ts = 1/60;

% Evaluation weights (used ONLY for terminal matrix + validation cost)
Q_eval  = diag([10 1 1]);
R1_eval = diag([1 2 1]);       % penalty on absolute input deviation u_k - u_ss
R2_eval = diag([1 1 1e-1]);    % penalty on delta-u

lqr_tuning = [-1.9980    0.0003    1.4849    0.5267   -0.9742    0.0425    0.1074   -0.1175];

% Terminal matrix for the fixed incremental feedback du = -K z
[P, K, Ai, Bi, xss, uss, LQR_data] = terminal_P_xu_du( ...
    Vsp, Xsp, par, model, ode_opt, Ts, lqr_tuning, Q_eval, R1_eval, R2_eval);

%% Verification by direct simulation (linear incremental closed-loop)

% Initial deviation
x0 = xss; x0(2) = 2; x0(3) = 1;
u_prev = uss;

z0 = [x0(:) - xss(:);
      u_prev(:) - uss(:)];

n  = numel(xss);
m  = numel(uss);
Sx = [eye(n), zeros(n,m)];
Su = [zeros(m,n), eye(m)];     % maps z -> (u_{k-1} - u_ss)

Acl = Ai - Bi*K;

N = 1000;
J_num = 0;
z = z0;

for k = 1:N
    du = -K*z;

    x_tilde = Sx*z;
    u_tilde = Su*z + du;       % u_k - u_ss = (u_{k-1}-u_ss) + du_k

    J_num = J_num ...
        + x_tilde.'*Q_eval*x_tilde ...
        + u_tilde.'*R1_eval*u_tilde ...
        + du.'*R2_eval*du;

    z = Acl*z;
end

J_P = z0.' * P * z0;

fprintf('J_num                 = %.12g\n', J_num);
fprintf('J_P (infinite-horizon) = %.12g\n', J_P);
fprintf('abs error              = %.3g\n', abs(J_num - J_P));
fprintf('rel error              = %.3g\n', abs(J_num - J_P)/max(1,abs(J_P)));

save('LQR_data.mat', "LQR_data")

%% FUNCTIONS

function [P, K, Ai, Bi, xss, uss, LQR_data] = terminal_P_xu_du( ...
        Vsp, Xsp, par, model, ode_opt, Ts, log10w, Q_eval, R1_eval, R2_eval)

    [xss, uss] = find_ss(Vsp, Xsp, par, model, ode_opt);

    [A, B] = linearize(xss, uss, model);
    [Ai, Bi] = incremental(A,B,Ts);

    nx = size(A,1);
    nu = size(B,2);

    % Design incremental LQR gain for stage cost:
    %   x'Qx + u'R1u + du'R2du with u = u_{k-1} + du
    [K, Qz, R, N] = build_LQR_full(log10w, Ai, Bi, nx, nu);

    % Closed-loop for the fixed policy du = -K z
    Acl = Ai - Bi*K;

    % Evaluation quadratic under du = -K z:
    %   x = Sx z
    %   u = (Su - K) z
    Sx = [eye(nx), zeros(nx,nu)];
    Su = [zeros(nu,nx), eye(nu)];

    LQR_data.Sx = Sx;
    LQR_data.K   = K;
    LQR_data.Acl = Acl;
    LQR_data.Su  = Su;
    P = construct_P(LQR_data, Q_eval, R1_eval, R2_eval);
end
