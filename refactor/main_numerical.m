clear all; close all; clc

% Add folders and subfolders to path
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))

%% Get model parameters
get_par

%% Load results struct to prevent erasing results if the code is partially
% executed
load results.mat


%% Steady-state solutions
plot_steady_state(par)
save_figure('results\steady state.pdf')

%% Numerical integration
disp('Tesing RK45...')
rng("default")
median_ratio = comparing_RK45(par);
results.median_ratio = median_ratio;
save results results



