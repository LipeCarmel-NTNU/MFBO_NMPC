%% Demo of the tangent/ suite on a self-contained nonlinear system
% Damped pendulum with torque input:
%   x1dot = x2
%   x2dot = -sin(x1) - 0.3 x2 + u
% Output map: y = x1 (angle). The suite is exercised end to end:
% trim -> linearize -> plant -> LQR synthesis -> weight tuning -> terminal P.

clear; close all; clc
rng(1)

f = @(x, u) [x(2); -sin(x(1)) - 0.3*x(2) + u(1)];
h = @(x) x(1);

%% 1) Given yss, find (xss, uss)
% yss = pi/4 constrains the angle; uss = NaN marks the input as free
yss = pi/4;
[xss, uss, trim_info] = trim_point(f, yss, [0; 0], NaN, h = h);
fprintf('xss = [%g %g], uss = %g (residual %.2e)\n', xss, uss, trim_info.resnorm)

%% 2) Given (xss, uss), find the state space
[A, B, C] = linearize_fd(f, xss, uss, h = h);
fprintf('Controllability rank deficiency: %d\n', size(A,1) - rank(ctrb(A, B)))

%% 3) Build design plant and synthesize LQR (incremental here; 'positional' also works)
Ts = 0.05;
plant = build_plant(A, B, Ts, 'incremental');

w0.Q  = eye(2);
w0.R1 = 1;
w0.R2 = 0.1;
K0 = lqr_synth(plant, w0);

%% 4) Optimize the LQR weights against the nonlinear closed loop
scenarios = { [0; 0], [pi/2; -0.5] }; % initial conditions to tune over

obj = @(var) tuning_cost(var, f, plant, xss, uss, scenarios, ...
    Tf = 10, y_weight = [10 1], reg_du = 1e2, reg_var = 1e-1, ...
    u_min = -2, u_max = 2);

var0 = zeros(1, plant.nx - 1 + 2*plant.nu);
solver_options = optimoptions('fminunc', 'Display', 'final', 'MaxFunctionEvaluations', 500);
[var_opt, J_opt] = lqr_tune(obj, var0, solver_options = solver_options);
fprintf('Tuned var = [%s], J = %.4g\n', num2str(var_opt, '%.3f '), J_opt)

K_opt = lqr_synth(plant, log10_weights(var_opt, plant));

%% 5) Terminal cost under NMPC evaluation weights
eval_w.Q  = diag([10 1]);
eval_w.R1 = 1;
eval_w.R2 = 0.1;
[P, P_info] = terminal_cost(plant, K_opt, eval_w);
fprintf('Terminal P validation rel. error: %.2e\n', P_info.rel_err)

%% 6) Compare initial vs tuned gain on the nonlinear plant
x0 = [pi/2; -0.5];
Tf = 10;
[Y0, T0, U0] = lqr_sim(f, plant, K0,    x0, xss, uss, Tf, u_min = -2, u_max = 2);
[Y1, T1, U1] = lqr_sim(f, plant, K_opt, x0, xss, uss, Tf, u_min = -2, u_max = 2);

figure(1); clf
tiledlayout(3, 1, 'TileSpacing', 'tight', 'Padding', 'tight')

nexttile
plot(T0, Y0(:,1), 'LineWidth', 2); hold on
plot(T1, Y1(:,1), 'LineWidth', 2)
yline(xss(1), '--k', 'LineWidth', 1.5)
ylabel('x_1'); legend({'initial', 'tuned'}, 'Location', 'best'); grid on

nexttile
plot(T0, Y0(:,2), 'LineWidth', 2); hold on
plot(T1, Y1(:,2), 'LineWidth', 2)
yline(xss(2), '--k', 'LineWidth', 1.5)
ylabel('x_2'); grid on

nexttile
stairs(T0, U0, 'LineWidth', 2); hold on
stairs(T1, U1, 'LineWidth', 2)
yline(uss, '--k', 'LineWidth', 1.5)
ylabel('u'); xlabel('Time'); grid on
