clear all; close all; clc

% Add folders and subfolders to path
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))

USE_PARALLEL = true;
p = gcp('nocreate');
if USE_PARALLEL
    NumWorkers = 8;
    if isempty(p) || p.NumWorkers ~= NumWorkers
        if ~isempty(p)
            delete(p);
        end
        parpool('Processes', NumWorkers);
    end
else
    NumWorkers = 1;
    delete(p)
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

%% MPC definition

% Control actions are subjected to non-negativity constraints.
nx = 3;
nu = 3;
NMPC = NMPC_terminal(model, nx, nu);
NMPC.optimizer_options.UseParallel = USE_PARALLEL;

% Sampling time in hours
NMPC.Ts = dt;

% Update optimisation horizons
NMPC.p = 5;
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
tf = 5/60;                % h

V_0   = 100;
X_0   = 15;
S_0   = 0;
x0_plant = [V_0, X_0, S_0];

%% Setpoints
V_sp   = 1;
X_sp   = 20;
[xss, uss] = find_ss(V_sp, X_sp, par, model, ode_opt);
NMPC.xsp = xss;
NMPC.usp = uss;

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

RUNTIME = zeros(N, 1);

%% Simulationx
for i = 1 : N
    timer = tic;

    fprintf('Simulated: %.1f %% \n', i/N*100)
    fprintf('Time elapsed: %.1f minutes \n', toc(timer)/60)

    Y(i, :) = x0_plant;
    % if i > 1
    %     x0_plant(1) = 100;
    % end
    if i > 1
        x0_plant(1) = 1;
    end
    %% Calculate control action

    Y_sp(i,:) = NMPC.xsp(1:3);
    uk = NMPC.solve(x0_plant(:)', uk(:)');

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


