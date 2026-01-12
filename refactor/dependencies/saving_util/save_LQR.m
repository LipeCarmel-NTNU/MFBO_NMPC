
load results/simulation_data.mat
n = 3;
% if isfield(simulation_data, 'run')
%     n = length(simulation_data.run);
%     n = n + 1;
% else
%     n = 1;
% end

% Run results
simulation_data.run(n).results.Y = Y;
simulation_data.run(n).results.T = T;
simulation_data.run(n).results.U = U;
simulation_data.run(n).results.Jx = Jx;

% Controller information
simulation_data.run(n).controller.type = 'LQR';
simulation_data.run(n).controller.data.K  = K;
simulation_data.run(n).controller.data.A  = A;
simulation_data.run(n).controller.data.B  = B;
simulation_data.run(n).controller.data.Q  = Q;
simulation_data.run(n).controller.data.R  = R;
simulation_data.run(n).controller.data.Ts = Ts;

save_simu_data(simulation_data)