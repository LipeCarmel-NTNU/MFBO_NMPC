clear all; close all; clc

% Add folders and subfolders to path
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))

%% Settings
ode_opt = odeset('NonNegative', [2 3]);

%% Get model parameters
get_par

%% System sampling rate
dt = 1/60;

%% Model definition
model = @(x, u) dilution_reduced(0, x, u(:)', par);
plant = model;

%% MPC definition

% Control actions are subjected to non-negativity constraints.
nx = 3;
nu = 3;
NMPC = NMPC_terminal(model, nx, nu);
NMPC.optimizer_options.UseParallel = true;

% Sampling time in hours
NMPC.Ts = dt;

% Update optimisation horizons
NMPC.p = 20;
NMPC.m = 3;

% Weights
Q_V = 10;
Q_X = 1;
Q_S = 2;
NMPC.Q = diag([Q_V, Q_X, Q_S]);

% Update constraint vectors used by the optimiser
NMPC.constraints();

% To calculate control actions use MPC.solve(x, u) where x and u are row
% vectors representing the current state and previous action.
% See:
help NMPC.solve

%% Initial conditions
tf = 10;                % h

V_0   = 1;
X_0   = 10;
S_0   = 0;
x0_plant = [V_0, X_0, S_0];

%% Setpoints
V_sp   = 1;
X_sp   = 20;
[xss, uss] = find_ss(V_sp, X_sp, par, model, ode_opt);
NMPC.xsp = xss;
NMPC.usp = uss;

lqr_tuning = [-1.9980    0.0003    1.4849    0.5267   -0.9742    0.0425    0.1074   -0.1175];

%% Terminal cost
load('LQR_data.mat')
P = construct_P(LQR_data, NMPC.Q, NMPC.Ru, NMPC.Rdu);
NMPC.P = P;

%% Initial input
uk = [0, 0, 0];

nu = length(uk);
nx = length(x0_plant);
NMPC.nx = nx;
NMPC.nu = nu;

%% Simulation setup
tspan = [0 dt];             % integration interval per control step
N = ceil(tf/dt) + 1;        % number of samples including initial condition

U = zeros(N, nu);
Y = zeros(N, nx);
Y_sp = zeros(N, nx);

U(1, :) = uk;


%% Simulation
timer = tic;
for i = 1 : N
    fprintf('Simulated: %.1f %% \n', i/N*100)
    fprintf('Time elapsed: %.1f minutes \n', toc(timer)/60)

    Y(i, :) = x0_plant;

    %% Calculate control action

    Y_sp(i,:) = NMPC.xsp(1:3);
    uk = NMPC.solve(x0_plant(:)', uk(:)');

    U(i, :) = uk;

    fprintf('\nControl action: \n')
    disp(uk)

    %% Plant propagation
    [~, y] = ode45(@(t,x) plant(x, uk), tspan, x0_plant, ode_opt);
    x0_plant = y(end, :);
end

%% Results
i = length(Y);
T = 0 : dt : (i-1)*dt;

figure(1);
clf

% --- State 1 ---
subplot(3,1,1);
plot(T, Y(1:i,1), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
plot(T, Y_sp(1:i,1), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
grid on; box on;
xlabel('Time (h)');
ylabel('State 1');
legend('Location','best');
hold off;

% --- State 2 ---
subplot(3,1,2);
plot(T, Y(1:i,2), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
plot(T, Y_sp(1:i,2), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
grid on; box on;
xlabel('Time (h)');
ylabel('State 2');
legend('Location','best');
hold off;

% --- State 3 ---
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


function P = construct_P(LQR_data, Q, R1, R2)

    Sx  = LQR_data.Sx;
    Su  = LQR_data.Su;
    K   = LQR_data.K;
    Acl = LQR_data.Acl;

    Qbar = (Sx.'*Q*Sx) ...
         + (Su - K).'*R1*(Su - K) ...
         + K.'*R2*K;

    % Infinite-horizon evaluation matrix
    P = dlyap(Acl', Qbar);
end