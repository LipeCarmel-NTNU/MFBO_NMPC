%% Demonstration: LQR designed with (Q_lqr, R_lqr), evaluated with (Q_eval, R_eval)
% Linear, discrete-time
clear; close all; clc

%% Plant (choose any stabilisable pair)
A = [ 0.95  0.10  0.00;
     -0.05  0.90  0.10;
      0.00 -0.08  0.92];

B = [0.10 0.00;
     0.00 0.08;
     0.05 0.02];

n = size(A,1); m = size(B,2);

% Sanity: controllability rank (for reference)
Co = ctrb(A,B);
fprintf('ctrb rank = %d (n=%d)\n', rank(Co), n);

%% 1) LQR design weights (tuning)
Q_lqr = diag([1 1e0 1e2]);
R_lqr = diag([1 1e-2]);

% Discrete-time LQR gain
[K, P_lqr, eigAcl] = dlqr(A,B,Q_lqr,R_lqr);

fprintf('max |eig(A-BK)| = %.6f\n', max(abs(real(eigAcl))));

%% 2) Evaluation weights (different from tuning)
Q_eval = diag(10*rand(3, 1));     % change state weighting
R_eval = diag(rand(2, 1));     % change input weighting

% Closed-loop matrix
Acl = A - B*K;

% Stage cost under fixed K is x' (Q_eval + K' R_eval K) x
Qbar = Q_eval + K'*R_eval*K;

%% 3) "Terminal"/infinite-horizon evaluation cost via Lyapunov
% For stable Acl, the infinite-horizon cost for u=-Kx is:
% J(x0) = sum_{k=0}^\infty (xk'Q_eval xk + uk'R_eval uk) = x0' P_eval x0
% where P_eval solves: P = Acl' P Acl + Qbar
P_eval = dlyap(Acl', Qbar);

%% 4) Numerical approximation of the same infinite-horizon cost
% Use a long horizon sum to approximate the infinite series.
x0 = [10; -8; 5];

N  = 5000;                 % increase if needed
J_num = 0;
xk = x0;

for k = 1:N
    uk = -K*xk;
    J_num = J_num + (xk'*Q_eval*xk + uk'*R_eval*uk);
    xk = Acl*xk;
end

J_lyap = x0'*P_eval*x0;

fprintf('\nEVALUATION (fixed K):\n');
fprintf('J_lyap = %.12g\n', J_lyap);
fprintf('J_num  = %.12g   (N=%d)\n', J_num, N);
fprintf('abs err = %.3g\n', abs(J_num - J_lyap));
fprintf('rel err = %.3g\n', abs(J_num - J_lyap)/max(1,abs(J_lyap)));

%% 5) Show that P_lqr (from dlqr) is NOT the evaluation cost matrix in general
J_wrong = x0'*P_lqr*x0;
fprintf('\nIf you (incorrectly) reuse P_lqr as evaluation terminal cost:\n');
fprintf('J_using_P_lqr = %.12g\n', J_wrong);
fprintf('bias vs correct J_lyap = %.3g\n', J_wrong - J_lyap);

%% 6) Optional: convergence diagnostic of the numerical sum
% The residual of the discrete Lyapunov equation should be ~0:
res = norm(P_eval - (Acl'*P_eval*Acl + Qbar), 'fro');
fprintf('\nLyapunov residual (Frobenius): %.3g\n', res);
