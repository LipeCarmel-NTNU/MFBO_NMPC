load results/simulation_data.mat

% RUN_IDX definitions
PID = 4;
NMPC_A = 5;

% Constants
filename = 'results/F_compared.pdf';
X_INSET = [13 13.6];
LINEWIDTH = 1;
YLIM = [0 0.4];

% Load data for both controllers
T_PID = simulation_data.run(PID).results.T;       % time [h]
U_PID = simulation_data.run(PID).results.U;       % [F_in, F_m, F_out] (L/h)

T_NMPC = simulation_data.run(NMPC_A).results.T;   % time [h]
U_NMPC = simulation_data.run(NMPC_A).results.U;   % [F_in, F_m, F_out] (L/h)

% Plot
cols = [
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

col_PID = cols(1,:);
col_NMPC = cols(9,:);

figure(1);
hold on;
stairs(T_PID, U_PID(:,1), 'Color', col_PID, 'LineWidth', LINEWIDTH, 'DisplayName', 'PID');
stairs(T_NMPC, U_NMPC(:,1), 'Color', col_NMPC, 'LineWidth', LINEWIDTH, 'DisplayName', 'NMPC');
hold off;
grid on;
ylabel('$F_{\mathrm{in}}$ (L/h)', 'Interpreter','latex');
xlabel('Time (h)','Interpreter','latex');
ylim(YLIM);
legend('Location','best');

% Inset
ax_inset = axes('Parent', gcf, 'Position', [0.20 0.4 0.30 0.30]); % [left bottom width height]
hold on;
stairs(T_PID, U_PID(:,1), 'Color', col_PID, 'LineWidth', LINEWIDTH);
stairs(T_NMPC, U_NMPC(:,1), 'Color', col_NMPC, 'LineWidth', LINEWIDTH);
hold off;
grid on;
xlim(ax_inset, X_INSET);
ylim(ax_inset, [0 0.14]);
box(ax_inset, 'on');

set_fig_size
set_font_size
save_figure(filename)