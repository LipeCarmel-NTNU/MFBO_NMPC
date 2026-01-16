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
R1_eval = diag([1 1 1]);       % penalty on absolute input deviation u_k - u_ss
R2_eval = diag([1 1 1e-1]);    % penalty on delta-u

% LQR design tuning (log10 weights), diagonal by construction:
%   log10w = [ q(1:n-1), r1(1:m), r2(1:m) ]
% where Q = diag([1, 10.^q]), R1 = diag(10.^r1), R2 = diag(10.^r2)
%
% For n=3, m=3 this must have length 8.
% If you want to mimic your previous "only du-penalty" design, choose r1 very small (e.g. -6).
lqr_tuning = [-2.8127  0.5583,  -6 -6 -6,   2.5095  -0.4645  -0.9051];

% Terminal matrix for the fixed incremental feedback du = -K z
[P, K, Ai, Bi, xss, uss] = terminal_P_xu_du( ...
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

%% ========================= FUNCTIONS =========================

function [P, K, Ai, Bi, xss, uss] = terminal_P_xu_du( ...
        Vsp, Xsp, par, model, ode_opt, Ts, log10w, Q_eval, R1_eval, R2_eval)

    [xss, uss] = find_ss(Vsp, Xsp, par, model, ode_opt);

    [A, B] = linearize(xss, uss, model);
    n = size(A,1);
    m = size(B,2);

    sysc = ss(A,B,eye(n),zeros(n,m));
    sysd = c2d(sysc, Ts, 'zoh');
    Ad = sysd.A;  Bd = sysd.B;

    % Augmentation with u_{k-1}
    Ai = [Ad, Bd;
          zeros(m,n), eye(m)];
    Bi = [Bd;
          eye(m)];

    % Design incremental LQR gain for stage cost:
    %   x'Qx + u'R1u + du'R2du with u = u_{k-1} + du
    [K, Qz, R, N] = build_LQR_xu_du(log10w, Ai, Bi, n, m);

    % Closed-loop for the fixed policy du = -K z
    Acl = Ai - Bi*K;

    % Evaluation quadratic under du = -K z:
    %   x = Sx z
    %   u = (Su - K) z
    Sx = [eye(n), zeros(n,m)];
    Su = [zeros(m,n), eye(m)];

    Qbar = (Sx.'*Q_eval*Sx) ...
         + (Su - K).'*R1_eval*(Su - K) ...
         + K.'*R2_eval*K;

    % Infinite-horizon evaluation matrix
    P = dlyap(Acl', Qbar);
end

function [K, Qz, R, N] = build_LQR_xu_du(log10w, Ai, Bi, n, m)
    % Stage cost in (z, du):
    %   x'Qx + u'R1u + du'R2du, with u = u_{k-1} + du
    %
    % Expanded form:
    %   z'Qz z + 2 z' N du + du' R du
    %
    % with Qz = blkdiag(Q, R1), N = [0; R1], R = R1 + R2.

    w = 10.^log10w(:);

    nq = n - 1;
    if numel(w) ~= (nq + 2*m)
        error('log10w must have length (n-1)+2m = %d (got %d).', nq + 2*m, numel(w));
    end

    q  = w(1:nq);
    r1 = w(nq + (1:m));
    r2 = w(nq + m + (1:m));

    Q  = diag([1; q(:)]);
    R1 = diag(r1(:));
    R2 = diag(r2(:));

    Qz = blkdiag(Q, R1);
    R  = R1 + R2;
    N  = [zeros(n,m);
          R1];

    K = dlqr(Ai, Bi, Qz, R, N);
end
