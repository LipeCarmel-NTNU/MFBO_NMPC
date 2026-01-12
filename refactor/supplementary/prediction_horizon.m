% This script is used to understand the prediction horizon dependency of
% Fin when substrate error is low
clear all; close all; clc

%% Get model parameters
get_par

%% System sampling rate
dt = 1/60;

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
MPC.p = 60;
MPC.m = 3;
MPC.constraints();

MPC.dilution = 1e7; % dilution is not allowed


% Parallel computation can be faster on most machines
if not(MPC.optimizer_options.UseParallel)
    % setup parallel
    delete(gcp('nocreate'))
    c = parcluster('Processes');
    parpool(2, 'IdleTimeout', inf);

    MPC.optimizer_options.UseParallel = true;
    MPC.optimizer_options.MaxIterations = 1000;
end
%% Test
MPC.Xsp = 20;
MPC.Ssp = 3;
MPC.m = 3;

p = [3, 6, 10, 15, 20, 30, 40, 60, 100];
x_init = [1 2 2.5 0.3];
u_init = [1 0 1]*0.4;
Fin = p;
for i = 1 : length(p)
MPC.p = p(i);
MPC.constraints();
[uk, x, u] = solve_first_iter(MPC, x_init, u_init);
Fin(i) = uk(1);
end

%%
colors = good_colors(1);
figure(1)
plot(p, Fin, 'o-','Color', colors, 'LineWidth', 2);
xlabel('Prediction Horizon p', 'Interpreter','latex');
ylabel('$F_{in}$ (L/h)', 'Interpreter','latex');
grid on;
box on;
set_fig_size
set_font_size
save_figure('supplementary/Fin_p.pdf')