% simulate_linear_2x2_all_options
% Small closed-loop NMPC example with all optional features enabled.

clear
clc

here = fileparts(mfilename('fullpath'));
addpath(here)

if isempty(which('fmincon'))
    error('simulate_linear_2x2_all_options:fmincon', ...
        'fmincon is required for this NMPC example.')
end

%% Stable 2x2 continuous-time plant: xdot = A*x + B*u
A = [-0.8  0.0;
      0.0 -1.1];
B = [1.0 0.2;
     0.1 0.8];
xdot = @(x, u) (A * x(:) + B * u(:)).';

Ts = 0.1;
p  = 6;
m  = 3;

x_target = [1.0 0.5];
u_target = steady_state_input(A, B, x_target);

opts = optimoptions('fmincon', ...
    'Display', 'off', ...
    'Algorithm', 'sqp', ...
    'MaxIterations', 300, ...
    'MaxFunEvals', Inf, ...
    'StepTolerance', 1e-9, ...
    'OptimalityTolerance', 1e-7, ...
    'ScaleProblem', true);

nmpc = NMPC( ...
    xdot = xdot, ...
    nx = 2, ...
    nu = 2, ...
    Ts = Ts, ...
    p = p, ...
    m = m, ...
    x_sp = [0 0], ...
    u_sp = [0 0], ...
    Q = diag([25 20]), ...
    R_u = diag([0.05 0.05])/10, ...
    Xmin = [ 0.0 0.0], ...
    Xmax = [ 0.8 0.4], ...
    umin = [-2 -2], ...
    umax = [ 2  2], ...
    P = diag([80 60]), ...
    x_scale = [1.5 1.0], ...
    u_scale = [1.2 1.0], ...
    soft_mask = [true true], ...
    rho_L2 = [1 1]*1e5, ...
    R_du = diag([0.3 0.3]), ...
    dumax = [0.25 0.25], ...
    optimizer_options = opts, ...
    log_enabled = true);

%% Closed-loop simulation with a setpoint step
n_steps = 60;
step_at = 6;

x = zeros(n_steps + 1, nmpc.nx);
u = zeros(n_steps,     nmpc.nu);
du = zeros(n_steps,    nmpc.nu);
slack_max = zeros(n_steps, nmpc.n_soft);
flag = zeros(n_steps, 1);

x(1, :) = [1 0];
u_prev = [0 0];

for k = 1:n_steps
    if k == step_at
        nmpc.x_sp = x_target;
        nmpc.u_sp = u_target;
    end

    [uk, ~, ~, info] = nmpc.solve(x(k, :), u_prev);
    u(k, :) = uk;
    du(k, :) = uk - u_prev;
    flag(k) = info.flag;

    if ~isempty(nmpc.latest_wopt) && nmpc.n_soft > 0
        [~, ~, s] = nmpc.unpack_phys(nmpc.latest_wopt);
        slack_max(k, :) = max(s, [], 1);
    end

    x(k + 1, :) = nmpc.step(x(k, :), uk);
    u_prev = uk;
end

t = (0:n_steps).' * Ts;
tu = (0:n_steps - 1).' * Ts;
x_ref = zeros(n_steps + 1, 2);
x_ref(step_at + 1:end, :) = repmat(x_target, n_steps - step_at + 1, 1);

fprintf('Final state:      [% .4f, % .4f]\n', x(end, 1), x(end, 2));
fprintf('Target state:     [% .4f, % .4f]\n', x_target(1), x_target(2));
fprintf('Max abs(du):      [% .4f, % .4f]\n', max(abs(du(:, 1))), max(abs(du(:, 2))));
fprintf('Max pred. slack:  [% .4f, % .4f]\n', max(slack_max(:, 1)), max(slack_max(:, 2)));
fprintf('Solver flags: min=%g, max=%g\n', min(flag), max(flag));
fprintf('Logged samples:   %d\n', numel(nmpc.log));

figure('Name', 'NMPC 2x2 all options')
tiledlayout(3, 1)

nexttile
plot(t, x, 'LineWidth', 1.2)
hold on
plot(t, x_ref, '--', 'LineWidth', 1.0)
grid on
ylabel('x')
legend('x_1', 'x_2', 'x_{sp,1}', 'x_{sp,2}', 'Location', 'best')

nexttile
stairs(tu, u, 'LineWidth', 1.2)
grid on
ylabel('u')
legend('u_1', 'u_2', 'Location', 'best')

nexttile
stairs(tu, abs(du), 'LineWidth', 1.2)
hold on
yline(nmpc.dumax(1), '--')
grid on
ylabel('|du|')
xlabel('time')
legend('|du_1|', '|du_2|', 'dumax', 'Location', 'best')

function u_ss = steady_state_input(A, B, x_ss)
    u_ss = (B \ (-A * x_ss(:))).';
end
