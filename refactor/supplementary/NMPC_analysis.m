addpath(genpath(pwd))
get_par
Kp = 0.8 / par.Y_XS;

% Load simulation data
load results\simulation_data.mat
if isfield(simulation_data, 'run')
    n = numel(simulation_data.run);
else
    warning('simulation_data.run is missing or empty. Nothing to process.');
    return
end

% Preallocate results table
results = table('Size', [0, 8], ...
    'VariableTypes', {'int64', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'idx', 'p', 'm', 'Q_S', 'Ts', 'Jx', 'Js', 'TV'});

% Common scalars
Jx_ideal = simulation_data.Jx_ideal;
Jx_adj   = simulation_data.Jx_adj;
Xsp      = simulation_data.Xsp;


idx = [];
for i = 1:n
    if contains(simulation_data.run(i).controller.type, 'MPC')
    idx = [idx i];
        % Extract trajectories and sampling
        Y  = simulation_data.run(i).results.Y;
        U  = simulation_data.run(i).results.U;
        T  = simulation_data.run(i).results.T;
        Ts = simulation_data.run(i).controller.data.Ts;

        % Extract state variables
        X = Y(:, 2);
        S = Y(:, 3);

        % Total variation
        TV = sum(trapz(T(1:end-1), abs(diff(U))));

        % Substrate setpoint trajectory
        Ssp = Kp * (Xsp - X);
        Ssp = min(Ssp, 3);

        % Performance indices
        Js = trapz(T, abs(Ssp - S));
        Jx = trapz(T, abs(Xsp - X));

        % Controller parameters
        p = simulation_data.run(i).controller.data.p;
        m = simulation_data.run(i).controller.data.m;
        Q_S = simulation_data.run(i).controller.data.Q_S;

        % Append to results table
        results = [results; {i, p, m, Q_S, Ts, Jx, Js, TV}];
    end
end


%% Normalising
% Reference condition
ref_idx = results.p == 60 & results.m == 3 & results.Q_S == 2;

% Extract reference values
Jx_ref = results.Jx(ref_idx);
Js_ref = results.Js(ref_idx);
TV_ref = results.TV(ref_idx);

% Compute normalised columns
results.Jx_norm = 100*results.Jx / Jx_ref;
results.Js_norm = 100*results.Js / Js_ref;
results.TV_norm = 100*results.TV / TV_ref;

% Display updated table
results = sortrows(results, {'Ts', 'Q_S', 'p', 'm'}, {'ascend', 'descend', 'descend', 'descend'});
disp(results)