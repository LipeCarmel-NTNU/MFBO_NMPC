%% Terminal cost for incremental (delta-u) LQR evaluated with alternative weights
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

% LQR tuning variables (log10-parameterisation used by build_LQR)
lqr_tuning = [-2.8127  0.5583  2.5095  -0.4645  -0.9051];

% Evaluation weights for terminal matrix (stage cost in discrete time)
%   ell_k = (x_k-xss)'Q_eval(x_k-xss) + (dU_k)'R_eval(dU_k)
Q_eval = diag([10 1 1]);
R_eval = diag([1  1 1e-1]);

% Terminal matrix for the fixed incremental feedback dU = -K z
[P, K, Ai, Bi, xss, uss] = terminal_P(Vsp, Xsp, par, model, ode_opt, Ts, lqr_tuning, Q_eval, R_eval);

%% Verification of terminal matrix P by direct simulation (incremental closed loop)

% Initial deviation
x0 = xss; x0(2)=2; x0(3)=1;
u_prev = uss;

z0 = [x0(:) - xss(:);
      u_prev(:) - uss(:)];

Acl = Ai - Bi*K;

N = 1000;
J_num = 0;
z = z0;

for k = 1:N
    dU = -K*z;
    x_tilde = z(1:numel(xss));

    J_num = J_num + x_tilde.'*Q_eval*x_tilde + dU.'*R_eval*dU;

    z = Acl*z;
end

J_P      = z0.' * P * z0;

fprintf('J_num                 = %.12g\n', J_num);
fprintf('J_P (infinite-horizon) = %.12g\n', J_P);

fprintf('abs(J_num - (J_P-tail)) = %.3g\n', abs(J_num - J_P));
fprintf('rel error               = %.3g\n', abs(J_num - J_P)/max(1,abs(J_P)));

function [P, K, Ai, Bi, xss, uss] = terminal_P(Vsp, Xsp, par, model, ode_opt, Ts, lqr_tuning, Q_eval, R_eval)
    % terminal_P  Infinite-horizon terminal matrix for a fixed incremental LQR policy.

    % Steady-state at setpoint
    [xss, uss] = find_ss(Vsp, Xsp, par, model, ode_opt);

    % Linearise plant at (xss, uss)
    [A, B] = linearize(xss, uss, model);
    nu = size(B,2);

    % Incremental LQR gain and corresponding discrete-time lifted model
    [K, Qi, Ri, Ai, Bi] = build_LQR(lqr_tuning, A, B, Ts);

    Acl = Ai - Bi*K;
    Qz   = blkdiag(Q_eval, zeros(nu));
    Qbar = Qz + K' * R_eval * K;

    % Discrete-time infinite-horizon evaluation matrix
    P = dlyap(Acl', Qbar);
end

function [K, Qi, Ri, Ai, Bi] = build_LQR(log10w, A, B, Ts)
    w = 10.^log10w;
    q = w(1:2);
    r = w(3:end);

    Q = diag([1 q]);
    R = diag(r);

    % LQR
    n = size(A,1);
    m = size(B,2);

    % Discretise
    sysc = ss(A,B,eye(n),zeros(n,m));
    sysd = c2d(sysc, Ts, 'zoh');
    Ad = sysd.A;  Bd = sysd.B;

    % incremental form
    Ai = [Ad, Bd;
        zeros(m,n), eye(m)];
    Bi = [Bd;
        eye(m)];


    Qi = blkdiag(Q, zeros(m));  % Q only on the natural states
    Ri = R;                     % R acts on delta u explicitly

    try
        K = dlqr(Ai, Bi, Qi, Ri);
    catch
        warning('something went wrong. Probably infinite cost somewhere')
        keyboard
    end

    % The following yields an interesting structure, but cant be trusted
    % Ai = [A, B;
    %       zeros(m,n), zeros(m)];
    % Bi = [B;
    %       eye(m)];
    % Qi = blkdiag(Q, zeros(m));   % penalise x only
    % Ri = R*Ts^2;                 % The input is dudt and scalling matters
    % [K2,S,e] = lqrd(Ai,Bi,Qi,Ri,Ts);

end