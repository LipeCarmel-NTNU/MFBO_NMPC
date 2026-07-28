%% NMPC test runs from timestamp lists (full horizon, no measurement noise)
%
% Reruns stored controllers with f forced to 1, so the cost is measured over
% the full horizon and no surrogate extrapolation is involved.

clear all; close all; clc;

current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))

delete_if_exists('.lock')
delete_if_exists('matlab.lock')

rng(1)

USE_PARALLEL = true;
NumWorkers = 2;
configure_pool(USE_PARALLEL, NumWorkers);

%% Run configuration
cfg_run = struct();
cfg_run.source_root = "results";
cfg_run.output_root = fullfile("results", "test_run");
cfg_run.log_path = fullfile("SIMULATIONS_LOG.txt");
cfg_run.sigma_y = [0 0 0];                  % no measurement noise
cfg_run.NumWorkers = NumWorkers;

base = nmpc_base(sigma_y = cfg_run.sigma_y);

%% Timestamp lists (Pareto-optimal controllers)
t1 = [
    "20260201_232106";
    "20260202_000115";
    "20260201_111557";
    "20260201_223337";
    "20260131_151035";
    "20260201_111928";
    "20260201_192807";
    "20260201_120223";
    "20260201_153030";
    "20260201_154920";
    "20260201_193032"
];

t2 = [
    "20260210_134744";
    "20260211_133012";
    "20260210_174745";
    "20260210_153818";
    "20260210_154010";
    "20260211_134235";
    "20260210_171107";
    "20260210_180826";
    "20260210_151703";
    "20260211_122653"
];

rerun_stored_controllers(cfg_run, base, "run1", t1, out_suffix = "_full_f1_no_noise");
rerun_stored_controllers(cfg_run, base, "run2", t2, out_suffix = "_full_f1_no_noise");
