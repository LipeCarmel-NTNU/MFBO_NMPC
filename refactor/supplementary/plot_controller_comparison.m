clear all; clc
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir))

%% Data
get_par
load simulation_data.mat
nRuns  = numel(simulation_data.run);
for i = 1 : nRuns
    disp(simulation_data.run(i).controller.type)
end
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
marker_size = 5;
line_width  = 1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PID vs NMPC
filename = 'results/PID_NMPC_closedloop.pdf';
include_run     = [PID NMPC_A NMPC_B]; % NMPC A* is 12
YMAX = [1.02; 21; 3.3];
YMIN = [0.98; 0; 0];
controllers  = {'PID','NMPC A','NMPC B'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LQR vs NMPC
% filename = 'results/LQR_NMPC_closedloop.pdf';
% include_run     = [LQR_A, LQR_B, LQR_C, NMPC_A];
% YMAX = [1.5; 21; 50];
% YMIN = [0.98; 0; 0];
% controllers  = {'LQR A','LQR B','LQR C', 'NMPC A'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors = good_colors(length(include_run) + 1); % +1 for the limit


%% Figure
figure(1); clf
ax = gobjects(3,1); % Store axes handles

for i = 1:3 % Plotting state i

    k = 0; % controller/color idx

    ax(i) = subplot(3,1,i); 
    hold on
    
    % loop over runs; each run gets a colour and a marker
    for r = 1:nRuns
        if any(include_run == r)
            k = k + 1;
            % Extract series
            T = simulation_data.run(r).results.T;
            Y = simulation_data.run(r).results.Y(:, i);
            
            % Continuous line (hidden from legend)
            if contains(controllers{k}, 'NMPC')
                line_type = '--';
            else
                line_type = '-';
            end
            plot(T, Y, line_type, 'LineWidth', line_width, ...
                'Color', colors(k,:), 'HandleVisibility', 'off');
            
            %% Markers with random jitter
            n_markers = 5;
            jitter_percent = 50;
            [marker_t, marker_y] = getMarkerPoints(T, Y, n_markers, jitter_percent);
            
            plot(marker_t, marker_y, markers{k}, ...
                'MarkerSize', marker_size, 'LineWidth', line_width, ...
                'Color', colors(k,:), 'DisplayName', controllers{k});

        end
    end
        
    % Setpoint line (common across runs)
    t_sp = linspace(T(1), T(end), 5);
    plot(t_sp, xsp(i) * ones(size(t_sp)), '.--', 'LineWidth', line_width, ...
        'MarkerSize', 10, 'Color', 'k', 'DisplayName', 'Setpoint');
    
    if i == 2
        % For the biomass plot 
        [Jx, t, y] = best_adjusted(Y(1), xsp(i), par);
        y = y(t< T(end));
        t = t(t< T(end));
        plot(t, y, '-', 'LineWidth', line_width, ...
        'MarkerSize', 10, 'Color', 'g', 'DisplayName', 'Limit');
    end

    % Remove X tick labels for top and middle plots
    if i < 3
        set(gca, 'XTickLabel', [])
    end
    
    if i == 2
        lgd = legend('Location','northoutside', 'Interpreter','latex', 'Orientation','horizontal');
    end
    
        if i == 3
            xlabel('Time (h)', 'Interpreter','latex')
        end

    grid off; box off
    ylabel(state_labels{i}, 'Interpreter','latex')
    ylim([YMIN(i) YMAX(i)])
end


set_font_size(13)
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