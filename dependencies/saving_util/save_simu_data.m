function save_simu_data(simulation_data)
% Save main file
save('results/simulation_data.mat', 'simulation_data');

% Generate timestamp
t = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
backup_name = "results/simulation_data_backup_" + string(t) + ".mat";

% Save backup
save(backup_name, 'simulation_data');

fprintf('Data saved to "results/simulation_data.mat" and backup created as "%s".\n', backup_name);
end