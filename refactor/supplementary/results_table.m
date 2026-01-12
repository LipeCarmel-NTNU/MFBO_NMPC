% Ensure we have runs and determine number of experiments
load results\simulation_data.mat
if isfield(simulation_data, 'run')
    n = numel(simulation_data.run);
elseif n == 0
    warning('simulation_data.run is missing or empty. Nothing to process.');
    return
end

% Preallocate results table
results = table('Size', [n, 5], ...
    'VariableTypes', {'double', 'string', 'double', 'double', 'double'}, ...
    'VariableNames', {'Index', 'ControllerType', 'Adj', 'Ideal', 'TV'});

% Common scalars
Jx_ideal = simulation_data.Jx_ideal;
Jx_adj   = simulation_data.Jx_adj;
Xsp      = simulation_data.Xsp;       % setpoint

for i = 1:n
    % Extract trajectories and sampling
    Y  = simulation_data.run(i).results.Y;     % output array, columns correspond to signals
    U  = simulation_data.run(i).results.U;     % input trajectory
    T  = simulation_data.run(i).results.T;     % time vector (monotone)
    Ts = simulation_data.run(i).controller.data.Ts;  % target sample time

    % Ensure uniformity
    t_interp = 0:Ts:T(end);
    X_interp = interp1(T, Y(:,2), t_interp);
    Jx = trapz(t_interp, abs(Xsp - X_interp));
    simulation_data.run(i).results.Jx = Jx;

    % Total variation
    TV = sum(trapz(T(1:end-1), abs(diff(U))));
    simulation_data.run(i).results.TV = TV;

    % Normalised performance indices
    adj   = 100 * Jx / Jx_adj;
    ideal = 100 * Jx / Jx_ideal;

    % Write
    results(i, :) = {i, simulation_data.run(i).controller.type, adj, ideal, TV};
end

% Display
results_disp = results;
cols = {'Adj','Ideal','TV'};
for c = 1:numel(cols)
    results_disp.(cols{c}) = round(results_disp.(cols{c}), 3);
end

disp(results_disp)

save results\simulation_data.mat simulation_data