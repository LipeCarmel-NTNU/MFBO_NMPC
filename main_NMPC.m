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

% Properties can be overwritten
MPC.Ts = dt; % Sampling time in hours

% Properties that modify the constraint matrices of the optimization
% require the constraints() method to be enforced
MPC.p = 5;
MPC.m = 2;
MPC.Q_S = 2;
MPC.constraints();

MPC.dilution = 1e7; % dilution is not allowedl

% Parallel computation can be faster on most machines
% if not(MPC.optimizer_options.UseParallel)
%     % setup parallel
%     c = parcluster('Processes');
%     delete(c.Jobs)
%     delete(gcp('nocreate'))
%     parpool(5, 'IdleTimeout', inf);
% 
%     MPC.optimizer_options.UseParallel = true;
%     MPC.optimizer_options.MaxIterations = 1000;
% end

% To calculate control actions use MPC.solve(x, u) where x and u are row
% vectors representing the current state and previous action.
% See:
help MPC.solve

%% Initial conditions

% Duration of the simulation
tf = 20;                % h

V_0 = 1;
X_0 = 2;
S_0 = 0;
CO2_0 = 0.0049;

x0_plant = [V_0, X_0, S_0, CO2_0];

MPC.Xsp = 20;
MPC.Vsp = 1;
MPC.Ssp = 0;


uk = [0, 0, 0];

nu = length(uk);
nx = length(x0_plant);

%% Simulation setup
tspan= [0 dt]; % from t to t+dt, or just 0 to dt because the ODE is time invariant

N = ceil(tf/dt) + 1;   % number of time samples + initial condition

% Initialization
U = zeros(N, nu);
Y = zeros(N, nx);
Y_sp = zeros(N, nx - 1); % CO2 is open-loop


U(1, :) = uk;


opts = odeset('NonNegative', [2 3]);

%% Simulation
timer = tic;
for i = 1 : N
    fprintf('Simulated: %.1f %% \n', i/N*100)
    fprintf('Time elapsed: %.1f minutes \n', toc(timer)/60)

    Y(i, :) = x0_plant;

    %% Calculate control action
    % Cascade
    Kp = 0.8/par.Y_XS;
    Ssp = Kp*(MPC.Xsp - x0_plant(2));
    Ssp = min([3 Ssp]);
    MPC.Ssp = Ssp;

    if x0_plant(2) > MPC.Xsp
        MPC.dilution = 0; % dilution is allowed
    else
        % Toggle slacked contraint
        MPC.dilution = 1e7;
    end

    Y_sp(i,:) = [MPC.Vsp, MPC.Xsp, MPC.Ssp];
    uk = MPC.solve(x0_plant(:)', uk(:)');

    U(i, :) = uk;

    fprintf('\nControl action: \n')
    disp(uk)

    %% Apply control action to the process and obtain y k+1

    % Plant
    [T,y] = ode45(@(t,x) plant(x, uk), tspan, x0_plant, opts);

    %% ---------------------------- k+1 -----------------------------------
    x0_plant = y(end, :);

end

%% Results
i = length(Y);
T = 0 : dt : (i-1)*dt;

figure(1);
clf

% --- Volume ---
subplot(3,1,1);
plot(T, Y(1:i,1), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
plot(T, Y_sp(1:i,1), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
grid on; box on;
xlabel('Time (h)');
ylabel('V (L)');
legend('Location','best');
hold off;

% --- Biomass concentration ---
subplot(3,1,2);
plot(T, Y(1:i,2), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
plot(T, Y_sp(1:i,2), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
grid on; box on;
xlabel('Time (h)');
ylabel('X (g/L)');
legend('Location','best');
hold off;

% --- Substrate concentration ---
subplot(3,1,3);
plot(T, Y(1:i,3), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
plot(T, Y_sp(1:i,3), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
grid on; box on;
xlabel('Time (h)');
ylabel('S (g/L)');
legend('Location','best');
hold off;


ax = findall(gcf, 'type', 'axes');

for i = 1:length(ax)
    ax(i).FontSize = 15;
    ax(i).XLabel.FontSize = 15;
    ax(i).YLabel.FontSize = 15;
end
%%
figure(2)
plot(T, U, 'LineWidth',2)
