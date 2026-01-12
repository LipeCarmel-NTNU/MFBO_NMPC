clear all; clc
 
%% Data
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))
 
 
load simulation_data.mat
nRuns  = numel(simulation_data.run);
 
get_par
xsp = [1; 20; 0];
 
% Indices / constants
Jx_ideal    = simulation_data.Jx_ideal;
Jx_adj = simulation_data.Jx_adj;
 
% Routing indexes for the simulation_data struct
LQR_A = 1;
LQR_B = 2;
LQR_C = 3;
 
PID = 4;
NMPC_A = 5;
NMPC_B = 7;
NMPC_A_star = 12;
 
%% Labels and configuration
 
state_labels = {'V (L)', 'X (g/L)', 'S (g/L)'};
 
 
% Predefine markers for each run
markers = {'o','s','d','^','v','>', 'x'};
marker_size = 7;
line_width  = 1.5;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PID vs NMPC
filename = 'results/PID_NMPC_closedloop.pdf';
include_run     = [PID NMPC_A NMPC_B]; % NMPC A* is 12
YMAX = [1.02; 21; 3.3];
YMIN = [0.98; 0; 0];
controllers  = {'PID','NMPC A','NMPC B'};
legend_labels = {'PID','NMPC A','NMPC B','Limit','Setpoint'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % LQR vs NMPC
% filename = 'results/LQR_NMPC_closedloop.pdf';
% include_run     = [LQR_A, LQR_B, LQR_C, NMPC_A];
% YMAX = [1.5; 21; 50];
% YMIN = [0.98; 0; 0];
% controllers  = {'LQR A','LQR B','LQR C', 'NMPC A'};
% legend_labels = {'LQR A','LQR B','LQR C', 'NMPC A','Limit','Setpoint'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
colors = good_colors(length(include_run) + 1); % +1 for the limit
colors = [
    0.00 0.30 0.50  % 1 dark blue
    0.55 0.25 0.00  % 2 dark orange
    0.20 0.45 0.60  % 3 darker light blue
    0.65 0.40 0.00  % 4 dark yellow-orange
    0.00 0.40 0.35  % 5 dark teal/green
    0.70 0.65 0.10  % 6 dark yellow
    0.00 0.30 0.00  % 7 dark green
    0.55 0.40 0.55  % 8 dark purple
    0.50 0.00 0.05  % 9 dark red-brown
    0.35 0.35 0.35  % 10 dark grey
];
%% Figure
figure(1); clf
ax = gobjects(3,1); % Store axes handles
 
for i = 1:3 % Plotting state i
 
    ax(i) = subplot(3,1,i); 
    hold on
    % loop over runs; each run gets a colour and a marker
    T = simulation_data.run(PID).results.T;
    Y = simulation_data.run(PID).results.Y(:, i);
    plot(T, Y, "-o", 'LineWidth', line_width, ...
        'MarkerIndices', 1:400:length(T), ...
        'MarkerSize', marker_size, 'LineWidth', line_width, ...
        'Color', colors(1,:), 'HandleVisibility', 'on');
    T = simulation_data.run(NMPC_A).results.T;
    Y = simulation_data.run(NMPC_A).results.Y(:, i);
    plot(T, Y, "--s", 'LineWidth', line_width, ...
        'MarkerIndices', 1:400:length(T), ...
        'MarkerSize', marker_size, 'LineWidth', line_width, ...
        'Color', 'red', 'HandleVisibility', 'on');   
   T = simulation_data.run(NMPC_B).results.T;
   Y = simulation_data.run(NMPC_B).results.Y(:, i);
   plot(T, Y, "--d", 'LineWidth', line_width, ...
        'MarkerIndices', 1:400:length(T), ...
        'MarkerSize', marker_size, 'LineWidth', line_width, ...
        'Color', colors(7,:), 'HandleVisibility', 'on');
 
    if i == 2
        % For the biomass plot 
        [Jx, t, y] = best_adjusted(Y(1), xsp(i), par);
        y = y(t< T(end));
        t = t(t< T(end));
        plot(t, y, "-.x", 'LineWidth', line_width, ...
        'MarkerIndices', 1:4:length(t), 'MarkerSize', 10, 'Color', 'k', 'DisplayName', 'Limit');
    end
    % Setpoint line (common across runs)
    t_sp = linspace(T(1), T(end), 5);
    plot(t_sp, xsp(i) * ones(size(t_sp)), "-.", 'LineWidth', 1, ...
         'Color', 'k', 'DisplayName', 'Setpoint');
    % this part is to have the ylim a bit above the setpoint line, so it
    % will be visible 
    lines = findall(ax, 'Type', 'Line');
    ymax = max(arrayfun(@(h) max(h.YData), lines));
    ylim([ax(i).YLim(1), ymax + (ymax>2)*0.1*ymax]);
    % Remove X tick labels for top and middle plots
    if i < 3
        set(gca, 'XTickLabel', [])
    end
    if i == 2
        lgd = legend(legend_labels,'Location','northoutside', 'Interpreter','latex', 'Orientation','horizontal');
    end
        if i == 3
            xlabel('Time (h)', 'Interpreter','latex')
        end
 
    grid off; box off
    ylabel(state_labels{i}, 'Interpreter','latex')
    %ylim([YMIN(i) YMAX(i)])
end
 
 
set_font_size(16)
set_fig_size(650, 500)
 
% Adjust spacing between subplots
subplot_spacing = 0.04;  % Adjust this value to control spacing
set(gcf, 'Units', 'normalized');
pos1 = get(ax(1), 'Position');
pos2 = get(ax(2), 'Position');
pos3 = get(ax(3), 'Position');
 
% Position plots with minimal spacing
pos2(2) = pos3(2) + pos3(4) + subplot_spacing;
pos1(2) = pos2(2) + pos2(4) + subplot_spacing;
 
set(ax(1), 'Position', pos1);
set(ax(2), 'Position', pos2);
set(ax(3), 'Position', pos3);
 
% Reposition legend
lgd.Position(2) = ax(1).Position(2) + ax(1).Position(4) + 0.04; % Place above first subplot
if lgd.Position(1) < 0
    lgd.Position(1) = 0.035;
end
ylim([YMIN(i) YMAX(i)])
 
 
save_figure(filename)
 
%% Data of the shown controllers
legends = controllers;  % unified source of names
k = 0;
for r = 1:nRuns
    if any(include_run == r)
        k = k + 1;
        disp(controllers{k})
        disp(simulation_data.run(r).controller.type)
        fprintf('%.2f %% adjusted\n', 100*simulation_data.run(r).results.Jx / Jx_adj)
        fprintf('%.2f %% ideal\n',     100*simulation_data.run(r).results.Jx / Jx_ideal)
        fprintf('%.3f TV\n',     simulation_data.run(r).results.TV)
        if contains(simulation_data.run(r).controller.type, 'MPC')
            disp(['p: ', num2str(simulation_data.run(r).controller.data.p)]);
            disp(['m: ', num2str(simulation_data.run(r).controller.data.m)]);
            disp(['Ts: ', num2str(simulation_data.run(r).controller.data.Ts)]);
            disp(['Q_S: ', num2str(simulation_data.run(r).controller.data.Q_S)]);
        end
        fprintf('\n\n')
    end
end
 
 
function [Jx, t, y] = best_adjusted(X0, Xsp, par)
    tspan = [0 100];
    [t, y] = ode45(@(t, biomass) ideal_model(biomass, Xsp, par), tspan, X0);
    y(y > Xsp) = Xsp;
 
    Jx = trapz(t, Xsp - y);
 
    function dxdt = ideal_model(X, Xsp, par)
        V = 1;
 
        % Steady-state Fin
        Fin = (1/par.Sin)*par.mu_max*X*par.Y_XSinv*V;
        mu_net = par.mu_max - par.kd - Fin;
        if X >= Xsp
            dxdt = 0;
        else
            dxdt = mu_net*X;
        end
    end
end