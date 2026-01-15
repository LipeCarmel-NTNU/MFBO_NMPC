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

% Performance was moved to the results table script
% simulation_data.run(n).results.Jx = Jx; 


% Controller information
simulation_data.run(n).controller.type = 'MPC';
simulation_data.run(n).controller.data.p  = MPC.p;
simulation_data.run(n).controller.data.m  = MPC.m;
simulation_data.run(n).controller.data.Rdu  = MPC.Rdu;
simulation_data.run(n).controller.data.Q_S  = MPC.Q_S;
simulation_data.run(n).controller.data.Ts = MPC.Ts;

save_simu_data(simulation_data)