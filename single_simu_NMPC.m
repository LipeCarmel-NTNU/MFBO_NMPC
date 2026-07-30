clear all; close all; clc

% Add folders and subfolders to path
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))

USE_PARALLEL = true;
pool = gcp('nocreate');
if USE_PARALLEL
    NumWorkers = 8;
    if isempty(pool) || pool.NumWorkers ~= NumWorkers
        if ~isempty(pool)
            delete(pool);
        end
        parpool('Processes', NumWorkers);
    end
else
    NumWorkers = 1;
    delete(pool)
end


%% Settings
ode_opt = odeset('NonNegative', [2 3]);

%% Get model parameters
get_par

%% System sampling rate
dt = 1/60;

%% Model definition
model = @(x, u) dilution_reduced(0, x, u(:)', par);
plant = model;

% The controller integrates xdot(x, u) with a row state and expects a 1-by-nx
% return, while model follows the orientation of its input.
xdot = @(x, u) reshape(model(x(:), u), 1, []);

%% Initial conditions
tf = 5/60;                % h

V_0   = 100;
X_0   = 15;
S_0   = 0;
x0_plant = [V_0, X_0, S_0];

nx = 3;
nu = 3;

%% Setpoints
V_sp   = 1;
X_sp   = 20;
[xss, uss] = find_ss(V_sp, X_sp, par, model, ode_opt);

%% Weights
Q_V = 10;
Q_X = 1;
Q_S = 2;
Q   = diag([Q_V, Q_X, Q_S]);
R_u  = diag([2 2 1]);
R_du = diag([100 100 10]);

%% Terminal cost
% P is augmented, (nx+nu)-by-(nx+nu), on z = [x_end - x_term; u_last - u_term],
% and construct_P linearises about the setpoint, so the terminal reference is
% the setpoint itself.
load('LQR_data.mat')
P = construct_P(LQR_data, Q, R_u, R_du);

%% MPC definition
% The biomass and sugar bounds are relaxed by an L1 slack, so a prediction that
% starts outside one stays feasible and pays rho_L1 per unit of violation. The
% volume bounds stay hard because the model divides by V.
% Control actions are subjected to non-negativity constraints.
nmpc = NMPC( ...
    xdot = xdot, nx = nx, nu = nu, Ts = dt, ...
    p = 5, m = 3, ...
    x_sp = xss, u_sp = uss, ...
    Q = Q, R_u = R_u, R_du = R_du, ...
    P = P, x_term = xss, u_term = uss, ...
    Xmin = [0.5 0 -0.1], Xmax = [2 50 20], ...
    umin = zeros(1, nu), umax = 0.4 * ones(1, nu), ...
    x_scale = [1 20 1], u_scale = 0.4 * ones(1, nu), ...
    soft_mask = [false true true], rho_L1 = 1e3);

nmpc.optimizer_options.UseParallel = USE_PARALLEL;
nmpc.optimizer_options.FiniteDifferenceType = 'central';

% To calculate control actions use nmpc.solve(x, u) where x and u are row
% vectors representing the current state and previous action.
% See:
help NMPC.solve

%% Initial input
uk = zeros(1, nu);

%% Simulation setup
tspan = [0 dt];             % integration interval per control step
N = ceil(tf/dt) + 1;        % number of samples including initial condition

U = zeros(N, nu);
Y = zeros(N, nx);
Y_sp = zeros(N, nx);

U(1, :) = uk;

RUNTIME = zeros(N, 1);

%% Simulation
for i = 1 : N
    timer = tic;

    fprintf('Simulated: %.1f %% \n', i/N*100)
    fprintf('Time elapsed: %.1f minutes \n', toc(timer)/60)

    Y(i, :) = x0_plant;
    if i > 1
        x0_plant(1) = 1;
    end

    %% Calculate control action
    Y_sp(i,:) = nmpc.x_sp(1:nx);
    uk = nmpc.solve(x0_plant(:)', uk(:)');

    U(i, :) = uk;

    fprintf('\nControl action: \n')
    disp(uk)

    %% Plant propagation
    [~, y] = ode45(@(t,x) plant(x, uk), tspan, x0_plant, ode_opt);
    x0_plant = y(end, :);

    RUNTIME(i) = toc(timer);

end
E = Y - Y_sp;
E = E.*[10 1 1];
SSE = sum(E.^2, 2);

dU  = diff(U, 1, 1);
SSdU = sum(dU.^2, 2);

%% Results
i = length(Y);
T = 0 : dt : (i-1)*dt;

figure(1);
clf

% State 1
subplot(3,1,1);
plot(T, Y(1:i,1), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
plot(T, Y_sp(1:i,1), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
grid on; box on;
xlabel('Time (h)');
ylabel('State 1');
legend('Location','best');
hold off;

% State 2
subplot(3,1,2);
plot(T, Y(1:i,2), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
plot(T, Y_sp(1:i,2), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
grid on; box on;
xlabel('Time (h)');
ylabel('State 2');
legend('Location','best');
hold off;

% State 3
subplot(3,1,3);
plot(T, Y(1:i,3), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
plot(T, Y_sp(1:i,3), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
grid on; box on;
xlabel('Time (h)');
ylabel('State 3');
legend('Location','best');
hold off;

ax = findall(gcf, 'type', 'axes');
for j = 1:length(ax)
    ax(j).FontSize = 15;
    ax(j).XLabel.FontSize = 15;
    ax(j).YLabel.FontSize = 15;
end

figure(2)
plot(T, U, 'LineWidth',2)
grid on; box on;
xlabel('Time (h)');
ylabel('Inputs');
