%% Exact ZOH discretisation in incremental coordinates: enforce (approximately) constant T*P across sampling times
% z_k = [x_k; u_{k-1}],  v_k = Delta u_k,  z_{k+1} = Ai z_k + Bi v_k
%
% Methodology identical to your original:
%   - Compute Pref at Tref by dlqr(Ai_ref,Bi_ref,Qref,Rref)
%   - Pd_ref = Tref*Pref
%   - For T=alpha*Tref:
%       Ptarget = (Tref/T)*Pref
%       R(T)=alpha*Rref
%       Q(T) from DARE identity so that Ptarget satisfies DARE for (Ai,Bi,Q,R)
%       Solve dlqr(Ai,Bi,Q,R) and verify T*Pdlqr ≈ Pd_ref
%   - Simulate for constant physical horizon Tsim (N = ceil(Tsim/T))
%
clear; close all; clc
rng(3)

%% 1) Random continuous-time controllable system (A,B) (continuous time)
nx = 5;
nu = 2;

max_tries = 200;
for ktry = 1:max_tries
    A = randn(nx) - 0.5*eye(nx);
    B = randn(nx,nu);
    if rank(ctrb(A,B)) == nx
        break
    end
end
assert(rank(ctrb(A,B)) == nx, "Failed to generate controllable (A,B).");

%% 2) Reference sample time and baseline discrete LQR at Tref (incremental form)
Tref = 0.20;

sysc = ss(A,B,eye(nx),zeros(nx,nu));
    sysd = c2d(sysc, Ts, 'zoh');
    Ad = sysd.A;  Bd = sysd.B;
    
[A_ref, B_ref] = incremental(A, B, Tref);
nzi = nx + nu;

% Reference weights in incremental coordinates:
%   stage cost = z'Qref z + du'Rref du
% Choose simple, coordinate-agnostic baseline weights:
Qx_ref = eye(nx);
Ru_ref = eye(nu);
Qref   = blkdiag(Qx_ref, Ru_ref);

Rref   = eye(nu);         % penalty on Delta u

[Kref, Pref, ~] = dlqr(A_ref, B_ref, Qref, Rref);

Pd_ref = Tref * Pref;     % "continuous-time value proxy" in incremental coordinates

assert(is_symmetric(Pref), "Pref should be symmetric.");

%% 3) Sweep alphas, enforce Ptarget(T)=(Tref/T)*Pref, compute Q(T), solve transformed LQR
alphas = logspace(log10(1.0), log10(0.01), 30);

Nsteps_ref = 300;
Tsim = Nsteps_ref * Tref;              % constant physical simulation time

x0 = randn(nx,1);
u_prev0 = zeros(nu,1);
z0 = [x0; u_prev0];

minEigQ        = nan(size(alphas));
minEigM        = nan(size(alphas));
condM          = nan(size(alphas));
rhoAcl         = nan(size(alphas));
Perr_to_target = nan(size(alphas));
Perr_TP        = nan(size(alphas));
dareResFro     = nan(size(alphas));
Jsim           = nan(size(alphas));
Jriccati       = nan(size(alphas));

% "structure" diagnostics: coupling block sizes in Q and P
Q_xu_norm      = nan(size(alphas));    % ||Q(1:nx,nx+1:end)||_F
P_xu_norm      = nan(size(alphas));    % ||P(1:nx,nx+1:end)||_F

dlqrOK         = false(size(alphas));
QpsdGoal       = false(size(alphas));
clStableGoal   = false(size(alphas));

for i = 1:numel(alphas)
    alpha = alphas(i);
    T = alpha * Tref;

    [Ai, Bi] = incremental(A, B, T);

    % Your imposed scaling on the *input* weight (here: Delta u penalty)
    R = alpha * Rref;

    % Target Riccati scaling: enforce constant T*P
    Ptarget = (Tref / T) * Pref;
    Ptarget = (Ptarget + Ptarget')/2;

    % Must-hold: M(T) = R + B'PB is PD for alpha>0 if Rref is PD
    S = Bi' * Ptarget * Bi;      % PSD
    M = R + S;                   % PD for alpha>0
    M = (M + M')/2;

    evM = eig(M);
    minEigM(i) = min(evM);
    condM(i) = cond(M);
    assert(minEigM(i) > 0, "Derivation requires M(T) invertible (PD).");

    % Construct Q(T) so that Ptarget satisfies the DARE identity
    Q = Ptarget - Ai' * Ptarget * Ai + Ai' * Ptarget * Bi * (M \ (Bi' * Ptarget * Ai));
    Q = (Q + Q')/2;

    eQ = eig(Q);
    minEigQ(i) = min(eQ);
    QpsdGoal(i) = (minEigQ(i) >= -1e-10);

    % "Structure" diagnostics: coupling between x and u_{k-1}
    Q_xu_norm(i) = norm(Q(1:nx, nx+1:end), 'fro');
    P_xu_norm(i) = norm(Ptarget(1:nx, nx+1:end), 'fro');

    % Must-hold: DARE identity residual should be ~0
    [~, dare_res] = dare_residual(Ai, Bi, Q, R, Ptarget);
    dareResFro(i) = norm(dare_res, 'fro');
    assert(dareResFro(i) < 1e-9, "DARE identity residual too large; check numerics.");

    % Goal-level: attempt dlqr
    try
        [Ktr, Pdlqr, ~] = dlqr(Ai, Bi, Q, R);
        Pdlqr = (Pdlqr + Pdlqr')/2;
        dlqrOK(i) = true;

        Perr_to_target(i) = norm(Pdlqr - Ptarget, 'fro');
        Perr_TP(i) = norm(T * Pdlqr - Pd_ref, 'fro');

        Acl = Ai - Bi * Ktr;
        rhoAcl(i) = max(abs(eig(Acl)));
        clStableGoal(i) = (rhoAcl(i) < 1 - 1e-10);

        % Constant physical time horizon
        N = ceil(Tsim / T);
        [Jsim(i), zN] = simulate_cost_incremental(Ai, Bi, Ktr, Q, R, z0, N);
        Jriccati(i) = z0' * Pdlqr * z0;

        % Optional exact identity check for truncation:
        % Jriccati should equal Jsim + zN'*Pdlqr*zN (up to numerical error)
        % gap = Jriccati(i) - (Jsim(i) + zN'*Pdlqr*zN);
        % if abs(gap) > 1e-7, fprintf("alpha=%.4g, truncation identity gap=%g\n",alpha,gap); end

    catch ME
        dlqrOK(i) = false;
        fprintf("alpha=%.4g: dlqr failed (%s)\n", alpha, ME.message);
    end
end

%% 4) Plots
figure;
semilogx(alphas, minEigQ, 'o-'); grid on
xlabel('\alpha = T/T_{ref}'); ylabel('min eig(Q(\alpha))')
title('Goal-level: Q(\alpha) \succeq 0 (not guaranteed)')

figure;
semilogx(alphas, minEigM, 'o-', alphas, condM, 'x-'); grid on
xlabel('\alpha = T/T_{ref}'); ylabel('min eig(M) and cond(M)')
legend('min eig(M)','cond(M)','Location','best')
title('Must-hold: M \succ 0; conditioning matters numerically')

figure;
ok = dlqrOK;
semilogx(alphas(ok), Perr_TP(ok), 'o-'); grid on
xlabel('\alpha = T/T_{ref}'); ylabel('||T P_{dlqr}(T) - T_{ref} P_{ref}||_F')
title('Goal-level: how well T*P is preserved (incremental coordinates)')

figure;
semilogx(alphas(ok), Jsim(ok), 'o-', alphas(ok), Jriccati(ok), 'x-'); grid on
xlabel('\alpha = T/T_{ref}'); ylabel('Cost')
legend('Simulated finite-horizon', 'z_0^T P_{dlqr} z_0', 'Location','best')
title('Consistency: finite-horizon cost vs Riccati value (incremental coordinates)')

figure;
semilogx(alphas, Q_xu_norm, 'o-', alphas, P_xu_norm, 'x-'); grid on
xlabel('\alpha = T/T_{ref}'); ylabel('Frobenius norm')
legend('||Q_{x,u_{prev}}||_F', '||Ptarget_{x,u_{prev}}||_F', 'Location','best')
title('Incremental structure: x–u_{prev} coupling lives in Q and/or P')

%% 5) Table
fprintf("\n%-10s %-12s %-12s %-10s %-12s %-10s %-8s %-8s %-10s %-12s\n", ...
    "alpha","minEig(Q)","minEig(M)","cond(M)","||TP-Pd||","rho(Acl)","Q_PSD","dlqrOK","||Q_xu||","||P_xu||");
for i = 1:numel(alphas)
    if dlqrOK(i)
        fprintf("%-10.4g %-12.4g %-12.4g %-10.3g %-12.3g %-10.6g %-8d %-8d %-10.3g %-12.3g\n", ...
            alphas(i), minEigQ(i), minEigM(i), condM(i), Perr_TP(i), rhoAcl(i), QpsdGoal(i), dlqrOK(i), Q_xu_norm(i), P_xu_norm(i));
    else
        fprintf("%-10.4g %-12.4g %-12.4g %-10.3g %-12s %-10s %-8d %-8d %-10.3g %-12.3g\n", ...
            alphas(i), minEigQ(i), minEigM(i), condM(i), "NaN", "NaN", QpsdGoal(i), dlqrOK(i), Q_xu_norm(i), P_xu_norm(i));
    end
end

%% -------- Local functions --------

function [Ai, Bi] = incremental(A,B,Ts)
    % Assumes C = I, D = 0;
    nx = size(A,1);
    nu = size(B,2);

    sysc = ss(A,B,eye(nx),zeros(nx,nu));
    sysd = c2d(sysc, Ts, 'zoh');
    Ad = sysd.A;  Bd = sysd.B;

    Ai = [Ad, Bd;
          zeros(nu,nx), eye(nu)];
    Bi = [Bd;
          eye(nu)];
end

function tf = is_symmetric(X)
    tf = norm(X - X', 'fro') <= 1e-12 * max(1, norm(X,'fro'));
end

function [K, res] = dare_residual(A, B, Q, R, P)
    S = B' * P * B;
    M = (R + S);
    M = (M + M')/2;
    K = (M \ (B' * P * A));
    res = P - (A' * P * A - A' * P * B * K + Q);
    res = (res + res')/2;
end

function [J, zN] = simulate_cost_incremental(A, B, K, Q, R, z0, N)
    z = z0;
    J = 0;
    for k = 1:N
        du = -K * z;
        J = J + z' * Q * z + du' * R * du;
        z = A*z + B*du;
    end
    zN = z;
end
