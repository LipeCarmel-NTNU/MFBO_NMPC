%% EKF-EM estimation of Q (process noise) given fixed R
% Two cases:
%   Case A: NO model mismatch (true Y_XS = model Y_XS)
%   Case B: WITH model mismatch (true Y_XS changes, filter model keeps nominal)
%
% Measurements:
%   Every minute: [V; X; CO2]
%   Every hour:   S

% clear; clc; close all;


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

% True Q used to generate data (unknown to filter)
Q_true = diag([1e-7, 2e-5, 1e-4, 2e-5]);

Q_final = diag([7.82705e-08 5.46446e-05 0.000423615 2.41746e-05]); % Q after EM to use in the filter

R_S_final = 0.048575;     % same initial guess as the true one, it should decrease because of the EM optimization step


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


xf_cont = zeros(4,N);
Pf_cont = zeros(4,4,N);

for k=1:N
    if k >1500 && k < 1800
        p_nom.Fin = 0.1/60;
        p_nom.Fout = 0.1/60;
    else
        p_nom.Fin = 0.;
        p_nom.Fout = 0.;
    end
    [x0_est, P0] = ekf_forward_multirate(y_fast_B(:,k), yS_B(:,k), x0_est, P0, Q_final, R_fast, R_S_final, p_nom, t, dt, k);
    xf_cont(:,k) = x0_est;
    Pf_cont(:,:,k) = P0;
end

figure()
subplot(4,1,1)
plot(t,xf_cont(1,:))
hold on
plot(t,x_true_B(1,:))
plot(t,y_fast_B(1,:),'.')
grid on; ylabel('V'); xlabel('Time (min)');
legend( {'Estimate','True state','Measurement'}, 'Location','best');
subplot(4,1,2)
plot(t,xf_cont(2,:))
hold on
plot(t,x_true_B(2,:))
plot(t,y_fast_B(2,:),'.')
grid on; ylabel('X'); xlabel('Time (min)');
legend( {'Estimate','True state','Measurement'}, 'Location','best');
subplot(4,1,3)
plot(t,xf_cont(3,:))
hold on
plot(t,x_true_B(3,:))
plot(t,yS_B(:),'o')
grid on; ylabel('S'); xlabel('Time (min)');
legend( {'Estimate','True state','Measurement'}, 'Location','best');
subplot(4,1,4)
plot(t,xf_cont(4,:))
hold on
plot(t,x_true_B(4,:))
plot(t,y_fast_B(3,:),'.')
grid on; ylabel('CO2'); xlabel('Time (min)');
legend( {'Estimate','True state','Measurement'}, 'Location','best');

function [xf, Pf] = ekf_forward_multirate(y_fast, yS, x0, P0, Qd, R_fast, R_S, p, t, dt, k)
% EKF forward pass with:
% - fast measurements y_fast = [V; X; CO2] 3X1
% - optional hourly substrate yS scalar (NaN otherwise)
% - p contains parameters and inputs athe the current time
%
% Returns filtered means/covs (xf,Pf), predicted (xp,Pp).

    

    H_fast = [1 0 0 0;
              0 1 0 0;
              0 0 0 1];
    H_S = [0 0 1 0];

    I = eye(4);

    if k == 1
        xf = x0;
        Pf = P0;
    else 
        % Predict
        x_prev = x0;

        
        x_pred = rk4_step(x_prev, dt, p);
        x_pred(x_pred<0) = 0;
        
        A = jacobian_A(x_prev, p);
        F = I + dt*A;

        P_pred = F*P0*F' + Qd;


        % Update with fast measurements
        yk = y_fast(:);
        ypred = H_fast*x_pred;

        S = H_fast*P_pred*H_fast' + R_fast;
        K = P_pred*H_fast'/S;

        x_upd = x_pred + K*(yk - ypred);
        P_upd = (I - K*H_fast)*P_pred;

        % Optional update with hourly substrate
        if ~isnan(yS)
            hs = H_S;
            ys = yS;

            Ss = hs*P_upd*hs' + R_S;
            Ks = (P_upd*hs')/Ss;

            x_upd = x_upd + Ks*(ys - hs*x_upd);
            P_upd = (I - Ks*hs)*P_upd;
        end

        xf = max(x_upd, [0.2;0;0;0]);
        Pf = P_upd;
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
    V   = max(x(1), 1e-3);
    X   = x(2);
    S   = x(3);

    denom = p.Ks + S;
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