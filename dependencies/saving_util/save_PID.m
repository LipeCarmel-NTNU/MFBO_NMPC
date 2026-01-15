load results/simulation_data.mat

if isfield(simulation_data, 'run')
    n = length(simulation_data.run);
    n = n + 1;
else
    n = 1;
end

% Run results
simulation_data.run(n).results.Y = Y;
simulation_data.run(n).results.T = T;
simulation_data.run(n).results.U = U;
simulation_data.run(n).results.Jx = Jx;

% Controller information
simulation_data.run(n).controller.type = 'PID';
simulation_data.run(n).controller.data.details  = '';

simulation_data.run(n).controller.data.substrate_pid = substrate_pid;
simulation_data.run(n).controller.data.volume_pid = volume_pid;
simulation_data.run(n).controller.data.growth_pid = growth_pid;
simulation_data.run(n).controller.data.dilution_pid = dilution_pid;

simulation_data.run(n).controller.data.Ts = Ts;
simulation_data.run(n).controller.data.deadband = deadband;

save_simu_data(simulation_data)