%% EKF-EM estimation of Q (process noise) given fixed R
% Two cases:
%   Case A: NO model mismatch (true Y_XS = model Y_XS)
%   Case B: WITH model mismatch (true Y_XS changes, filter model keeps nominal)
%
% Measurements:
%   Every minute: [V; X; CO2]
%   Every hour:   S
%
% Notes:
% - Nonlinear state transition -> EKF-based EM (approximate E-step via EKF + RTS smoother).
% - M-step updates Q using smoothed moments and linearization.
% - Here we estimate Q as DIAGONAL only (common in practice).
%

clear; clc; close all;

%% -------------------------
% User settings
% --------------------------
dt   = 1;          % minutes
Tend = 2500;       % minutes
t    = 0:dt:Tend;
N    = numel(t);

rng(4);

p_nom = get_p0();     % model used by EKF and EM (nominal)
R_fast = diag([1e-5, 2e-2, 2e-4]);  % measurement noise for [V,X,CO2]
R_S    = (0.30)^2;                 % substrate measurement variance (hourly)

% Initial conditions
x0_true = [1.00; 0.9; 10.0; 0.10];
x0_est  = [1.02; 1.0; 10.5; 0.15];

P0 = diag([1e-4, 5e-2, 2.0, 5e-2]);

% True Q used to generate data (unknown to EM)
Q_true = diag([1e-7, 2e-5, 1e-4, 2e-5]);

% EM settings
nEM = 100*3; % number of EM iterations
Q_init = diag([1e-4, 1e-2, 1e-2, 1e-2]); % deliberately imperfect starting Q
q_floor = 1e-12;  % keep Q positive
r_floor = 1e-12;  % keep R positive
R_S_init = R_S;     % same initial guess as the true one, it should decrease because of the EM optimization step

%% -------------------------
% Generate synthetic data: Case A (no mismatch)
% --------------------------
fprintf('Generating Case A (no mismatch)...\n');
p_true_A = p_nom;
[x_true_A, y_fast_A, yS_A] = simulate_truth_and_measurements(x0_true, p_true_A, Q_true, R_fast, R_S, t, dt, false);

%% -------------------------
% Generate synthetic data: Case B (with mismatch in Y_XS)
% True system changes Y_XS after t_mis minutes, but EKF model remains nominal.
% --------------------------
fprintf('Generating Case B (with mismatch in Y_XS)...\n');
p_true_B = p_nom;
t_mis = 1;      % minutes
YXS_new = 0.45;% 0.55;%;    % true yield after mismatch time
[x_true_B, y_fast_B, yS_B] = simulate_truth_and_measurements(x0_true, p_true_B, Q_true, R_fast, R_S, t, dt, true, t_mis, YXS_new);

%% -------------------------
% Run EM to estimate Q (Case A)
% --------------------------
case_type = 'onlyQ';
fprintf('Running EM for Case A only Q...\n');
[Qhat_hist_A, x_s_A] = em_estimate_Q( ...
    y_fast_A, yS_A, x0_est, P0, Q_init, R_fast, R_S, p_nom, t, dt, nEM, q_floor);

% case_type = 'QandR';
% fprintf('Running EM for Case A, Q and R...\n');
% [Qhat_hist_A, Rhat_fast_hist_A, Rhat_S_hist_A, x_s_A] = em_estimate_QR( ...
%     y_fast_A, yS_A, x0_est, P0, Q_init, R_fast, R_S, p_nom, t, dt, nEM, q_floor, r_floor);

case_type = 'QandR_S';
% fprintf('Running EM for Case A, Q and R_S...\n');
% [Qhat_hist_A, Rhat_S_hist_A, x_s_A] = em_estimate_Q_and_RS( ...
%     y_fast_A, yS_A, x0_est, P0, Q_init, R_fast, R_S_init, p_nom, t, dt, nEM, q_floor, r_floor);

%% -------------------------
% Run EM to estimate Q (Case B)
% --------------------------

% fprintf('Running EM for Case B only Q...\n');
% [Qhat_hist_B, x_s_B] = em_estimate_Q( ...
%     y_fast_B, yS_B, x0_est, P0, Q_init, R_fast, R_S, p_nom, t, dt, nEM, q_floor);

% fprintf('Running EM for Case B only Q...\n');
% [Qhat_hist_B, Rhat_fast_hist_B, Rhat_S_hist_B, x_s_B] = em_estimate_QR( ...
%     y_fast_B, yS_B, x0_est, P0, Q_init, R_fast, R_S, p_nom, t, dt, nEM, q_floor, r_floor);


fprintf('Running EM for Case A, Q and R_S...\n');
[Qhat_hist_B, Rhat_S_hist_B, x_s_B] = em_estimate_Q_and_RS( ...
    y_fast_B, yS_B, x0_est, P0, Q_init, R_fast, R_S_init, p_nom, t, dt, nEM, q_floor, r_floor);

%% 
% run the EKF
Q_final_A  = Qhat_hist_A(:,:,end);
if exist('Rhat_S_hist_A')
    RS_final_A = Rhat_S_hist_A(end);
else
    RS_final_A = R_S;
end

Q_final_B  = Qhat_hist_B(:,:,end);
RS_final_B = Rhat_S_hist_B(end);

% --- run EKF forward pass with final tuning (Case A) ---
[xf_A, Pf_A, xp_A, Pp_A, Fk_A] = ekf_forward_multirate( ...
    y_fast_A, yS_A, x0_est, P0, Q_final_A, R_fast, RS_final_A, p_nom, t, dt);

% --- run EKF forward pass with final tuning (Case B) ---
[xf_B, Pf_B, xp_B, Pp_B, Fk_B] = ekf_forward_multirate( ...
    y_fast_B, yS_B, x0_est, P0, Q_final_B, R_fast, RS_final_B, p_nom, t, dt);


%% -------------------------
% Plot EM convergence of Q diagonals
% --------------------------
qA = squeeze(Qhat_hist_A(1:4,1:4,:));
qB = squeeze(Qhat_hist_B(1:4,1:4,:));

diagA = zeros(4,nEM);
diagB = zeros(4,nEM);
for it = 1:nEM
    diagA(:,it) = diag(qA(:,:,it));
    diagB(:,it) = diag(qB(:,:,it));
end

names = {'V','X','S','CO2'};

figure;
for i = 1:4
    plot(1:nEM, diagA(i,:), 'LineWidth', 1.5); hold on;
end
yline(diag(Q_true(1:4,1:4)),'k:'); % true diagonals (dotted)
xlabel('EM iteration'); ylabel('Q diagonal value');
title('Case A (no mismatch): EM estimates of diag(Q)');
legend([strcat('Q_{',names,'}') 'Q true'], 'Location','best');
grid on;

figure;
for i = 1:4
    plot(1:nEM, diagB(i,:), 'LineWidth', 1.5); hold on;
end
yline(diag(Q_true(1:4,1:4)),'k:'); % true diagonals (for comparison)
xlabel('EM iteration'); ylabel('Q diagonal value');
title('Case B (with mismatch): EM estimates of diag(Q)');
legend([strcat('Q_{',names,'}') 'Q true'], 'Location','best');
grid on;

%% -------------------------
% Plot EM convergence of R_fast diagonals
% --------------------------
if strcmp(case_type,'QandR')
    RfastA = Rhat_fast_hist_A;   % 3x3xnEM
    RfastB = Rhat_fast_hist_B;
    
    diagRfastA = zeros(3,nEM);
    diagRfastB = zeros(3,nEM);
    
    for it = 1:nEM
        diagRfastA(:,it) = diag(RfastA(:,:,it));
        diagRfastB(:,it) = diag(RfastB(:,:,it));
    end
    
    meas_names = {'V','X','CO2'};
    
    figure;
    for i = 1:3
        plot(1:nEM, diagRfastA(i,:), 'LineWidth', 1.5); hold on;
    end
    yline(diag(R_fast), 'k:');   % true used for simulation (since you set it)
    xlabel('EM iteration'); ylabel('R_{fast} diagonal value');
    title('Case A: EM estimates of diag(R_{fast})');
    legend([strcat('R_{',meas_names,'}') 'R true'], 'Location','best');
    grid on;
    
    figure;
    for i = 1:3
        plot(1:nEM, diagRfastB(i,:), 'LineWidth', 1.5); hold on;
    end
    yline(diag(R_fast), 'k:');
    xlabel('EM iteration'); ylabel('R_{fast} diagonal value');
    title('Case B: EM estimates of diag(R_{fast})');
    legend([strcat('R_{',meas_names,'}') 'R true'], 'Location','best');
    grid on;
end

%% -------------------------
% Plot EM convergence of R_S (hourly substrate variance)
% --------------------------
figure;
if exist('Rhat_S_hist_A')
    plot(1:nEM, Rhat_S_hist_A, 'LineWidth', 1.5); hold on;
    yline(R_S, 'k:');
    xlabel('EM iteration'); ylabel('R_S value');
    title('Case A: EM estimates of R_S (substrate measurement variance)');
    legend('R_S estimate','R_S true','Location','best');
    grid on;
end

figure;
plot(1:nEM, Rhat_S_hist_B, 'LineWidth', 1.5); hold on;
yline(R_S, 'k:');
xlabel('EM iteration'); ylabel('R_S value');
title('Case B: EM estimates of R_S (substrate measurement variance)');
legend('R_S estimate','R_S true','Location','best');
grid on;

%% -------------------------
% Compare smoothed state estimates vs true (last EM iter)
% --------------------------
figName = sprintf("Case A - No Mismatch");   % format as needed
figure('Name',figName);

subplot(4,1,1);
h1 = plot(t, x_true_A(1,:), ...
          t, xf_A(1,:),  ...
          t, y_fast_A(1,:),'.', 'LineWidth', 1.2);
grid on; ylabel('V');
title(figName);
legend(h1, {'True state','Estimate','Measurement'}, 'Location','best');

subplot(4,1,2);
h2 = plot(t, x_true_A(2,:), ...
          t, xf_A(2,:),    ...
          t, y_fast_A(2,:),'.', 'LineWidth', 1.2);
grid on; ylabel('X');
legend(h2, {'True state','Estimate','Measurement'}, 'Location','best');

subplot(4,1,3);
h3 = plot(t, x_true_A(3,:), ...
          t, xf_A(3,:),   ...
          t, yS_A,         'o', 'LineWidth', 1.2);
grid on; ylabel('S');
legend(h3, {'True state','Estimate','Measurement'}, 'Location','best');

subplot(4,1,4);
h4 = plot(t, x_true_A(4,:),  ...
          t, xf_A(4,:),   ...
          t, y_fast_A(3,:),'.', 'LineWidth', 1.2);
grid on; ylabel('CO2'); xlabel('Time (min)');
legend(h4, {'True state','Estimate','Measurement'}, 'Location','best');

% mismatch
figName = sprintf("Case B - Mismatch Y_XS=%.2f", YXS_new);   % format as needed
figure('Name',figName);

subplot(4,1,1);
h1 = plot(t, x_true_B(1,:),  ...
          t, xf_B(1,:),    ...
          t, y_fast_B(1,:),'.', 'LineWidth', 1.2);
grid on; ylabel('V');
title(figName);
legend(h1, {'True state','Estimate','Measurement'}, 'Location','best');

subplot(4,1,2);
h2 = plot(t, x_true_B(2,:),  ...
          t, xf_B(2,:),    ...
          t, y_fast_B(2,:),'.', 'LineWidth', 1.2);
grid on; ylabel('X');
legend(h2, {'True state','Estimate','Measurement'}, 'Location','best');

subplot(4,1,3);
h3 = plot(t, x_true_B(3,:),  ...
          t, xf_B(3,:),    ...
          t, yS_B,         'o', 'LineWidth', 1.2);
grid on; ylabel('S');
legend(h3, {'True state','Estimate','Measurement'}, 'Location','best');

subplot(4,1,4);
h4 = plot(t, x_true_B(4,:),  ...
          t, xf_B(4,:),    ...
          t, y_fast_B(3,:),'.', 'LineWidth', 1.2);
grid on; ylabel('CO2'); xlabel('Time (min)');
legend(h4, {'True state','Estimate','Measurement'}, 'Location','best');


%% ========================================================================
% FUNCTIONS
% ========================================================================
function [Qhist, RShist, x_s_last] = em_estimate_Q_and_RS( ...
    y_fast, yS, x0_est, P0, Q0, R_fast_fixed, R_S0, p, t, dt, nEM, q_floor, r_floor)
% EKF-EM estimation of diagonal Q and scalar R_S, with fixed R_fast.
%
% Inputs:
%   y_fast         3xN measurements [V; X; CO2] each minute
%   yS             1xN substrate measurement (NaN when missing, e.g. hourly)
%   x0_est, P0     initial state mean/cov for EKF
%   Q0             initial process noise covariance (4x4)
%   R_fast_fixed   fixed measurement noise covariance for fast sensors (3x3)
%   R_S0           initial variance for substrate measurement (scalar)
%   p, t, dt       model params/time
%   nEM            number of EM iterations
%   q_floor        minimum allowed diagonal Q element
%   r_floor        minimum allowed R_S
%
% Outputs:
%   Qhist          4x4xnEM history of Q estimates (diagonal)
%   RShist         1xnEM history of R_S estimates
%   x_s_last       4xN smoothed state (last EM iteration)

    N = size(y_fast,2);

    % Measurement matrices
    H_fast = [1 0 0 0;
              0 1 0 0;
              0 0 0 1];
    H_S    = [0 0 1 0];

    Q   = Q0;
    R_S = R_S0;

    Qhist  = zeros(4,4,nEM);
    RShist = zeros(1,nEM);

    for it = 1:nEM

        % =========================
        % E-step: EKF + RTS smoother
        % =========================
        [xf, Pf, xp, Pp, Fk] = ekf_forward_multirate( ...
            y_fast, yS, x0_est, P0, Q, R_fast_fixed, R_S, p, t, dt);

        [xs, Ps, Pcs] = rts_smoother(xf, Pf, xp, Pp, Fk);

        % =========================
        % M-step: update Q (diagonal)
        % =========================
        Qnew = zeros(4);

        for k = 1:N-1
            xk_s    = xs(:,k);
            xkp1_s  = xs(:,k+1);
            
            if k >1500 && k < 1800
                p.Fin = 0.1/60;
                p.Fout = 0.1/60;
            else
                p.Fin = 0.;
                p.Fout = 0.;
            end

            % mean prediction from smoothed state point
            xkp1_pred = rk4_step(xk_s, dt, p);

            % linearization
            A = jacobian_A(xk_s, p);
            F = eye(4) + dt*A;

            % smoothed covariances
            Pk     = Ps(:,:,k);
            Pk1    = Ps(:,:,k+1);
            Pk_k1  = Pcs(:,:,k);

            % expected moments 
            Ex1x1 = Pk1 + xkp1_s*xkp1_s';
            Exfx  = (Pk_k1')*F' + xkp1_s*xkp1_pred';
            Efx1  = F*Pk_k1 + xkp1_pred*xkp1_s';
            Eff   = (F*Pk*F') + xkp1_pred*xkp1_pred';

            Qk = Ex1x1 - Exfx - Efx1 + Eff;
            Qnew = Qnew + Qk;
        end

        Qnew = Qnew/(N-1);

        % enforce diagonal-only + floor
        Q = diag(max(diag(Qnew), q_floor));

        % =========================
        % M-step: update R_S (scalar)
        % =========================
        idxS = find(~isnan(yS));  % indices where slow measurement exists

        if ~isempty(idxS)
            RS_acc = 0;
            for ii = 1:numel(idxS)
                k  = idxS(ii);

                % residual using smoothed state
                rs = yS(k) - H_S*xs(:,k);

                % E[(y - Hx)^2] = residual^2 + H P H'
                RS_acc = RS_acc + (rs^2 + H_S*Ps(:,:,k)*H_S');
            end

            R_S = max(RS_acc/numel(idxS), r_floor);
        end

        % store history
        Qhist(:,:,it)  = Q;
        RShist(it)     = R_S;

        fprintf('EM %3d/%3d: diag(Q)=[%g %g %g %g], R_S=%g\n', ...
            it, nEM, diag(Q), R_S);
    end

    x_s_last = xs;
end


function [Qhist, Rhist_fast, Rhist_S, x_s_last] = em_estimate_QR( ...
    y_fast, yS, x0_est, P0, Q0, R0_fast, R0_S, p, t, dt, nEM, q_floor, r_floor)

    N = size(y_fast,2);

    Q = Q0;
    R_fast = R0_fast;
    R_S = R0_S;

    Qhist      = zeros(4,4,nEM);
    Rhist_fast = zeros(3,3,nEM);
    Rhist_S    = zeros(1,nEM);

    % Measurement matrices
    H_fast = [1 0 0 0;
              0 1 0 0;
              0 0 0 1];
    H_S = [0 0 1 0];

    for it = 1:nEM
        % ----- E-step -----
        [xf, Pf, xp, Pp, Fk] = ekf_forward_multirate(y_fast, yS, x0_est, P0, Q, R_fast, R_S, p, t, dt);
        [xs, Ps, Pcs]        = rts_smoother(xf, Pf, xp, Pp, Fk);

        % ----- M-step: Q (your current logic) -----
        Qnew = zeros(4);
        for k = 1:N-1
            xk_s   = xs(:,k);
            xkp1_s = xs(:,k+1);

            xkp1_pred = rk4_step(xk_s, dt, p);

            A = jacobian_A(xk_s, p);
            F = eye(4) + dt*A;

            Pk     = Ps(:,:,k);
            Pk1    = Ps(:,:,k+1);
            Pk_k1  = Pcs(:,:,k);

            Ex1x1 = Pk1 + xkp1_s*xkp1_s';
            Exfx  = (Pk_k1')*F' + xkp1_s*xkp1_pred';
            Efx1  = F*Pk_k1 + xkp1_pred*xkp1_s';
            Eff   = (F*Pk*F') + xkp1_pred*xkp1_pred';

            Qk = Ex1x1 - Exfx - Efx1 + Eff;
            Qnew = Qnew + Qk;
        end
        Qnew = Qnew/(N-1);
        Q = diag(max(diag(Qnew), q_floor));

        % ----- M-step: R_fast -----
        Rfast_acc = zeros(3);
        for k = 1:N
            rk = y_fast(:,k) - H_fast*xs(:,k);                 % smoothed residual
            Rfast_acc = Rfast_acc + (rk*rk' + H_fast*Ps(:,:,k)*H_fast');
        end
        Rfast_new = Rfast_acc / N;

        % diagonal-only (recommended unless you really need correlations)
        R_fast = diag(max(diag(Rfast_new), r_floor));

        % ----- M-step: R_S (hourly only) -----
        idxS = find(~isnan(yS));  % time indices where substrate measured
        if ~isempty(idxS)
            RS_acc = 0;
            for ii = 1:numel(idxS)
                k = idxS(ii);
                rs = yS(k) - H_S*xs(:,k);
                RS_acc = RS_acc + (rs^2 + H_S*Ps(:,:,k)*H_S');
            end
            R_S = max(RS_acc/numel(idxS), r_floor);
        end

        % store
        Qhist(:,:,it)       = Q;
        Rhist_fast(:,:,it)  = R_fast;
        Rhist_S(it)         = R_S;

        fprintf('EM %3d/%3d: diag(Q)=[%g %g %g %g], diag(Rfast)=[%g %g %g], R_S=%g\n', ...
            it, nEM, diag(Q), diag(R_fast), R_S);
    end

    x_s_last = xs;
end

function [Qhist, x_s_last] = em_estimate_Q(y_fast, yS, x0_est, P0, Q0, R_fast, R_S, p, t, dt, nEM, q_floor)
% EKF-EM to estimate diagonal Q given fixed R.
%
% E-step: EKF (multi-rate) + RTS smoother to get x_s, P_s, and P_{k,k+1}^s
% M-step: update Q using smoothed moments and linearization around x_s(k)

    N = size(y_fast,2);
    Q = Q0;
    Qhist = zeros(4,4,nEM);

    for it = 1:nEM
        % ----- E-step: EKF forward pass -----
        [xf, Pf, xp, Pp, Fk] = ekf_forward_multirate(y_fast, yS, x0_est, P0, Q, R_fast, R_S, p, t, dt);

        % ----- RTS smoother -----
        [xs, Ps, Pcs] = rts_smoother(xf, Pf, xp, Pp, Fk);

        % ----- M-step: update Q -----
        Qnew = zeros(4);

        for k = 1:N-1
            % Linearize f around smoothed state xs(:,k)
            xk_s = xs(:,k);
            xkp1_s = xs(:,k+1);

            if k >1500 && k < 1800
                p.Fin = 0.1/60;
                p.Fout = 0.1/60;
            else
                p.Fin = 0.;
                p.Fout = 0.;
            end
            % Predicted mean from smoothed point
            xkp1_pred = rk4_step(xk_s, dt, p);

            % Discrete Jacobian F_k approximated as I + dt*A(xk_s)
            A = jacobian_A(xk_s, p);
            F = eye(4) + dt*A;

            % Smoothed covariances
            Pk  = Ps(:,:,k);
            Pk1 = Ps(:,:,k+1);
            Pk_k1 = Pcs(:,:,k);  % approx smoothed cross-cov Cov(x_k, x_{k+1})

            % Compute expected residual covariance:
            % E[(x_{k+1}-f(x_k))(...)^T] with f linearized at xk_s:
            % f(x_k) ≈ xkp1_pred + F (x_k - xk_s)
            %
            % Qk = E[x_{k+1}x_{k+1}^T] - E[x_{k+1}f^T] - E[f x_{k+1}^T] + E[ff^T]

            Ex1x1 = Pk1 + xkp1_s*xkp1_s';
            Exfx  = (Pk_k1')*F' + xkp1_s*xkp1_pred';   % E[x_{k+1} f(x_k)^T]
            Efx1  = F*Pk_k1 + xkp1_pred*xkp1_s';       % E[f(x_k) x_{k+1}^T]
            Eff   = (F*Pk*F') + xkp1_pred*xkp1_pred';  % E[f f^T]

            Qk = Ex1x1 - Exfx - Efx1 + Eff;

            Qnew = Qnew + Qk;
        end

        Qnew = Qnew / (N-1);

        % Enforce diagonal-only estimate (common practical choice)
        Q = diag(max(diag(Qnew), q_floor));

        Qhist(:,:,it) = Q;

        fprintf('  EM %2d/%2d: diag(Q) = [%g, %g, %g, %g]\n', it, nEM, diag(Q));
    end

    x_s_last = xs;
end

function [xf, Pf, xp, Pp, Fk] = ekf_forward_multirate(y_fast, yS, x0, P0, Qd, R_fast, R_S, p, t, dt)
% EKF forward pass with:
% - fast measurements y_fast(:,k) = [V; X; CO2] at every k
% - optional hourly substrate yS(k) (NaN otherwise)
%
% Returns filtered means/covs (xf,Pf), predicted (xp,Pp), and Fk for RTS.

    N = size(y_fast,2);

    H_fast = [1 0 0 0;
              0 1 0 0;
              0 0 0 1];
    H_S = [0 0 1 0];

    xf = zeros(4,N);
    Pf = zeros(4,4,N);
    xp = zeros(4,N);
    Pp = zeros(4,4,N);
    Fk = zeros(4,4,N-1);

    xf(:,1) = x0;
    Pf(:,:,1) = P0;

    I = eye(4);

    for k = 2:N
        % Predict
        x_prev = xf(:,k-1);

        if k >1500 && k < 1800
            p.Fin = 0.1/60;
            p.Fout = 0.1/60;
        else
            p.Fin = 0.;
            p.Fout = 0.;
        end

        x_pred = rk4_step(x_prev, dt, p);

        A = jacobian_A(x_prev, p);
        F = I + dt*A;

        P_pred = F*Pf(:,:,k-1)*F' + Qd;

        xp(:,k) = x_pred;
        Pp(:,:,k) = P_pred;
        Fk(:,:,k-1) = F;

        % Update with fast measurements
        yk = y_fast(:,k);
        ypred = H_fast*x_pred;

        S = H_fast*P_pred*H_fast' + R_fast;
        K = P_pred*H_fast'/S;

        x_upd = x_pred + K*(yk - ypred);
        P_upd = (I - K*H_fast)*P_pred;

        % Optional update with hourly substrate
        if ~isnan(yS(k))
            hs = H_S;
            ys = yS(k);

            Ss = hs*P_upd*hs' + R_S;
            Ks = (P_upd*hs')/Ss;

            x_upd = x_upd + Ks*(ys - hs*x_upd);
            P_upd = (I - Ks*hs)*P_upd;
        end

        xf(:,k) = max(x_upd, [0.2;0;0;0]);
        Pf(:,:,k) = P_upd;
    end

    % Define xp(:,1), Pp(:,:,1) for completeness
    xp(:,1) = xf(:,1);
    Pp(:,:,1) = Pf(:,:,1);
end

function [xs, Ps, Pcs] = rts_smoother(xf, Pf, xp, Pp, Fk)
% RTS smoother for time-varying linearizations.
% Outputs:
%   xs(:,k), Ps(:,:,k) smoothed means/covs
%   Pcs(:,:,k) approx smoothed cross-cov Cov(x_k, x_{k+1}) for k=1..N-1
%
% NOTE: cross-cov formula here is a common practical approximation used in EKF-EM.
% It is adequate for estimating diagonal Q in many applications.

    [n, N] = size(xf);
    xs = zeros(n,N);
    Ps = zeros(n,n,N);
    Pcs = zeros(n,n,N-1);

    xs(:,N) = xf(:,N);
    Ps(:,:,N) = Pf(:,:,N);

    for k = N-1:-1:1
        F = Fk(:,:,k);
        Ppred_next = Pp(:,:,k+1);

        % Smoother gain
        J = Pf(:,:,k)*F'/(Ppred_next);

        xs(:,k) = xf(:,k) + J*(xs(:,k+1) - xp(:,k+1));
        Ps(:,:,k) = Pf(:,:,k) + J*(Ps(:,:,k+1) - Ppred_next)*J';

        % Approx cross-covariance Cov(x_k, x_{k+1})
        % Practical approximation:
        Pcs(:,:,k) = J*Ps(:,:,k+1);
    end
end

function [x_true, y_fast, yS] = simulate_truth_and_measurements(x0_true, p_true, Q_true, R_fast, R_S, t, dt, useMismatch, t_mis, YXS_new)
% Simulate truth + measurements.
% If useMismatch=true: p_true.Y_XS changes after t_mis minutes (true model only).

    N = numel(t);

    H_fast = [1 0 0 0;
              0 1 0 0;
              0 0 0 1];

    x_true = zeros(4,N);
    y_fast = zeros(3,N);
    yS = nan(1,N);

    x_true(:,1) = x0_true;

    for k = 2:N
        if useMismatch && t(k) >= t_mis
            p_true.Y_XS = YXS_new;
        end
        
        if k >1500 && k < 1800
            p_true.Fin = 0.1/60;
            p_true.Fout = 0.1/60;
        else
            p_true.Fin = 0.;
            p_true.Fout = 0.;
        end
        x_pred = rk4_step(x_true(:,k-1), dt, p_true);
        w = mvnrnd([0 0 0 0], Q_true)'; % process noise
        x_true(:,k) = max(x_pred + w, [0.2;0;0;0]);
    end

    for k = 1:N
        vfast = mvnrnd([0 0 0], R_fast)';
        y_fast(:,k) = max(H_fast*x_true(:,k) + vfast, 0.0);

        if mod(t(k), 60) == 0
            yS(k) = max(x_true(3,k) + sqrt(R_S)*randn, 0.0);
        end
    end
end

function p = get_p0()
    p.mu_max   = 0.1944/60;   % 1/min
    p.Ks       = 0.2;         % g/L
    p.k_d      = 0.006/60;    % 1/min
    p.Y_XS     = 0.5;         % gX/gS
    p.Y_XCO2   = 0.8;         % gX/gCO2
    p.Qgas     = 2/60;        % 1/min
    p.Fin      = 0.00;        % L/min
    p.Fout     = 0.0;         % L/min
    p.Sin      = 70.0;        % g/L
end

function xnext = rk4_step(x, dt, p)
    k1 = f_cstr(x, p);
    k2 = f_cstr(x + 0.5*dt*k1, p);
    k3 = f_cstr(x + 0.5*dt*k2, p);
    k4 = f_cstr(x + dt*k3, p);
    xnext = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function dx = f_cstr(x, p)
    V   = max(x(1), 1e-6);
    X   = x(2);
    S   = x(3);
    CO2 = x(4);

    monod = S/(p.Ks + S + 1e-12);

    dV   = p.Fin - p.Fout;
    dX   = -X*(p.Fin/V) + p.mu_max*monod*X - p.k_d*X;
    dS   = (p.Sin - S)*(p.Fin/V) - p.mu_max*monod*X/p.Y_XS;
    dCO2 = -CO2*p.Qgas + X*monod*p.mu_max/p.Y_XCO2;

    dx = [dV; dX; dS; dCO2];
end

function A = jacobian_A(x, p)
    V   = max(x(1), 1e-6);
    X   = x(2);
    S   = x(3);

    denom = (p.Ks + S + 1e-12);
    monod = S/denom;
    dmonod_dS = p.Ks/(denom^2);

    A = zeros(4,4);

    % dV/dx = 0

    % dX
    A(2,1) = X*p.Fin/(V^2);
    A(2,2) = -p.Fin/V + p.mu_max*monod - p.k_d;
    A(2,3) = p.mu_max*X*dmonod_dS;

    % dS
    A(3,1) = -(p.Sin - S)*p.Fin/(V^2);
    A(3,2) = -(p.mu_max*monod)/p.Y_XS;
    A(3,3) = -(p.Fin/V) - (p.mu_max*X*dmonod_dS)/p.Y_XS;

    % dCO2
    A(4,2) = (p.mu_max*monod)/p.Y_XCO2;
    A(4,3) = (p.mu_max*X*dmonod_dS)/p.Y_XCO2;
    A(4,4) = -p.Qgas;
end