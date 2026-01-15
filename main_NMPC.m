clear all; close all; clc

% Add folders and subfolders to path
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))

%% Get model parameters
get_par

%% System sampling rate
dt = 6/60;

%% Model definition
nmpc_model = @(x, u) dilution_Monod_model(0, x, [2, u(:)'], par);
plant = nmpc_model;

%% MPC definition

% MPC_mono_dilution object requires a 4 x 3 process model to be constructed
% Default setpoints, tuning parameters, and constraints are predefined.
% Control actions are subjected to non-negativity constraints.
MPC = MPC_mono_dilution(nmpc_model);

% Sampling time in hours
MPC.Ts = dt;

% Update optimisation horizons
MPC.p = 5;
MPC.m = 2;

% Weights
Q_V = 10;
Q_X = 1;
Q_S = 2;
Q_CO2 = 0;
MPC.Q = diag([Q_V, Q_X, Q_S, Q_CO2]);

% Update constraint vectors used by the optimiser
MPC.constraints();

% Penalise slack on the dilution-related term
MPC.dilution = 1e7;

% To calculate control actions use MPC.solve(x, u) where x and u are row
% vectors representing the current state and previous action.
% See:
help MPC.solve

%% Initial conditions
tf = 10;                % h

V_0   = 1;
X_0   = 10;
S_0   = 0;
CO2_0 = 0.0049;
x0_plant = [V_0, X_0, S_0, CO2_0];

%% Setpoints
V_sp   = 1;
X_sp   = 20;
S_sp   = 0;
CO2_sp = 0;
MPC.xsp = [V_sp, X_sp, S_sp, CO2_sp];

%% Initial input
uk = [0, 0, 0];

nu = length(uk);
nx = length(x0_plant);

%% Simulation setup
tspan = [0 dt];             % integration interval per control step
N = ceil(tf/dt) + 1;        % number of samples including initial condition

U = zeros(N, nu);
Y = zeros(N, nx);
Y_sp = zeros(N, nx - 1);    % last state is not tracked here

U(1, :) = uk;

opts = odeset('NonNegative', [2 3]);

%% Simulation
timer = tic;
for i = 1 : N
    fprintf('Simulated: %.1f %% \n', i/N*100)
    fprintf('Time elapsed: %.1f minutes \n', toc(timer)/60)

    Y(i, :) = x0_plant;

    %% Calculate control action
    Kp = 0.8/par.Y_XS;
    Ssp = Kp*(MPC.xsp(2) - x0_plant(2));
    Ssp = min([3 Ssp]);
    MPC.xsp(3) = Ssp;

    if x0_plant(2) > MPC.xsp(2)
        MPC.dilution = 0;
    else
        MPC.dilution = 1e7;
    end

    Y_sp(i,:) = MPC.xsp(1:3);
    uk = MPC.solve(x0_plant(:)', uk(:)');

    U(i, :) = uk;

    fprintf('\nControl action: \n')
    disp(uk)

    %% Plant propagation
    [~, y] = ode45(@(t,x) plant(x, uk), tspan, x0_plant, opts);
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
