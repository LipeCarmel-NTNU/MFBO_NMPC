%% Plant
clear; clc; close all;

Ts    = 0.1;   % sampling time [s]
tau_p = 2.0;   % plant time constant [s]

% Continuous-time model: xdot = -x + u
[A, B] = first_order_discrete(Ts, tau_p);

x0    = 10.0;
u0    = 0.0;   % this is u_{k-1} at k = 0
r     = 0.0;   %#ok<NASGU>  % reference, included for completeness

%% LQR synthesis
Q   = 1.0;
Rdu = 0.01;

[K, S, Abar, Bbar, Qbar] = synthesize_delta_u_lqr(A, B, Q, Rdu); %#ok<NASGU>

fprintf('Discrete augmented matrices:\n');
disp('Abar ='); disp(Abar);
disp('Bbar ='); disp(Bbar);
fprintf('LQR gain K = [%.12f  %.12f]\n', K(1), K(2));

%% Filtered output definition
alpha_f = 0.5;        % filter

%% Define N = 1000
N       = 1000;
T       = 500;      % finite horizon for Monte Carlo cost
sigma_x = 5;      % measurement noise std

rng(1, 'twister');
noise_sequences = sigma_x * randn(T, N);   % identical noise realisations for both cases

%% Nominal simulation (plot)
% LQR cost-to-go at z0                 : 201.430859834329112
% Closed-loop simulated cost           : 201.430859834329340
% Absolute difference                  : 2.274e-13
nom = simulate_nominal(A, B, K, Q, Rdu, x0, u0, T);

z0       = [x0; u0];
J_lqr    = z0.' * S * z0;
J_closed = nom.J;

figure('Name', 'Nominal closed-loop simulation', 'Color', 'w');
tiledlayout(3,1);

nexttile;
plot(0:T, nom.x, 'LineWidth', 1.5);
grid on;
xlabel('k');
ylabel('x_k');
title('Nominal state');

nexttile;
stairs(0:T-1, nom.u, 'LineWidth', 1.5);
grid on;
xlabel('k');
ylabel('u_k');
title('Nominal input');

nexttile;
stairs(0:T-1, nom.du, 'LineWidth', 1.5);
grid on;
xlabel('k');
ylabel('\Delta u_k');
title('Nominal input increments');

fprintf('\nNominal case\n');
fprintf('LQR cost-to-go at z0                 : %.15f\n', J_lqr);
fprintf('Closed-loop simulated cost           : %.15f\n', J_closed);
fprintf('Absolute difference                  : %.3e\n', abs(J_lqr - J_closed));

%% disp(####)
disp(repmat('#',1,70))

%% N simulations without filtering
% Noisy case without filtering
% Mean closed-loop cost                : 10319.217595
%  (expected closed-loop cost)
% Std. dev. closed-loop cost           : 635.438108
J_mc_unfiltered = simulate_monte_carlo_unfiltered(A, B, K, Q, Rdu, x0, u0, noise_sequences);

figure('Name', 'Cost distributions', 'Color', 'w');
tiledlayout(2,1);

nexttile;
histogram(J_mc_unfiltered, 40);
hold on;
xline(J_closed, 'r', 'LineWidth', 1.5, 'Label', 'nominal', ...
    'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Closed-loop cost');
ylabel('Count');
title(sprintf('Noisy case without filtering, N = %d', N));

fprintf('Noisy case without filtering\n');
fprintf('Mean closed-loop cost                : %.6f\n', mean(J_mc_unfiltered));
fprintf('Std. dev. closed-loop cost           : %.6f\n', std(J_mc_unfiltered));

%% disp(####)
disp(repmat('#',1,70))

%% Repeat with filtering
% Noisy case with filtering
% Mean closed-loop cost                : 8207.916163
%  (expected closed-loop cost with a filter)
% Std. dev. closed-loop cost           : 538.097251
J_mc_filtered = simulate_monte_carlo_filtered(A, B, K, Q, Rdu, x0, u0, alpha_f, noise_sequences);

nexttile;
histogram(J_mc_filtered, 40);
hold on;
xline(J_closed, 'r', 'LineWidth', 1.5, 'Label', 'nominal', ...
    'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Closed-loop cost');
ylabel('Count');

fprintf('Noisy case with filtering\n');
fprintf('Mean closed-loop cost                : %.6f\n', mean(J_mc_filtered));
fprintf('Std. dev. closed-loop cost           : %.6f\n', std(J_mc_filtered));

%% FUNCTIONS
function [A, B] = first_order_discrete(Ts, tau)
    % Continuous-time plant: xdot = -(1/tau)x + (1/tau)u
    A = exp(-Ts/tau);
    B = 1 - A;
end

function [K, S, Abar, Bbar, Qbar] = synthesize_delta_u_lqr(A, B, Q, Rdu)
    % Augmented model:
    % z_k   = [x_k; u_{k-1}]
    % z_{k+1} = Abar*z_k + Bbar*delta_u_k
    Abar = [A, B;
            0, 1];
    Bbar = [B;
            1];
    Qbar = [Q, 0;
            0, 0];

    [K, S, ~] = dlqr(Abar, Bbar, Qbar, Rdu);
end

function sim = simulate_nominal(A, B, K, Q, Rdu, x0, u0, T)
    x  = zeros(T+1,1);
    u  = zeros(T,1);
    du = zeros(T,1);

    x(1)   = x0;
    u_prev = u0;
    J      = 0.0;

    for k = 1:T
        zk    = [x(k); u_prev];
        du(k) = -K * zk;
        u(k)  = u_prev + du(k);

        J = J + Q*x(k)^2 + Rdu*du(k)^2;

        x(k+1) = A*x(k) + B*u(k);
        u_prev = u(k);
    end

    sim.x  = x;
    sim.u  = u;
    sim.du = du;
    sim.J  = J;
end

function J_all = simulate_monte_carlo_unfiltered(A, B, K, Q, Rdu, x0, u0, noise_sequences)
    [T, N] = size(noise_sequences);
    J_all  = zeros(N,1);

    for i = 1:N
        x      = x0;
        u_prev = u0;
        J      = 0.0;

        for k = 1:T
            x_meas = x + noise_sequences(k,i);

            z_meas = [x_meas; u_prev];
            du     = -K * z_meas;
            u      = u_prev + du;

            J = J + Q*x^2 + Rdu*du^2;

            x      = A*x + B*u;
            u_prev = u;
        end

        J_all(i) = J;
    end
end

function J_all = simulate_monte_carlo_filtered(A, B, K, Q, Rdu, x0, u0, alpha_f, noise_sequences)
    [T, N] = size(noise_sequences);
    J_all  = zeros(N,1);

    for i = 1:N
        x          = x0;
        u_prev   = u0;
        J          = 0.0;

        for k = 1:T
            x_meas = x + noise_sequences(k,i);

            z_meas = [x_meas; u_prev];
            du     = -K * z_meas;
            u_cmd  = u_prev + du;

            % First-order low-pass on the control signal
            u = alpha_f*u_prev + (1 - alpha_f)*u_cmd;

            J = J + Q*x^2 + Rdu*du^2;

            x          = A*x + B*u;
            u_prev   = u;
        end

        J_all(i) = J;
    end
end