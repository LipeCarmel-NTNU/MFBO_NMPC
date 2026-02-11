clear; close all; clc

%% --- Paths / options ---
folder = fullfile("results", "run1");
csvPath         = fullfile(folder, "results.csv");
outScatterPath  = fullfile(folder, "sse_vs_ssdu.png");
outRuntimePath  = fullfile(folder, "runtime_vs_iteration.png");

if ~isfile(csvPath)
    error("CSV not found: %s", csvPath);
end
if ~isfolder("results")
    mkdir("results");
end

%% --- Load CSV (expects columns: timestamp,SSE,SSdU,J,runtime_s,theta_1..theta_12) ---
T = readtable(csvPath, "TextType", "string");

requiredCols = ["timestamp","SSE","SSdU","J","runtime_s"];
for c = requiredCols
    if ~ismember(c, string(T.Properties.VariableNames))
        error("Missing required column '%s' in %s", c, csvPath);
    end
end
for i = 1:12
    nm = "theta_" + i;
    if ~ismember(nm, string(T.Properties.VariableNames))
        error("Missing required column '%s' in %s", nm, csvPath);
    end
end

%% --- Parse / derive fields ---
T.timestamp_dt = datetime(T.timestamp, "InputFormat","yyyyMMdd_HHmmss");

T.iteration   = (1:height(T)).';
T.runtime_min = T.runtime_s / 60;

T.f        = double(T.theta_1);
T.theta_p  = round(double(T.theta_2));
T.theta_m  = round(double(T.theta_3));

T.q1_log10 = double(T.theta_4);
T.q2_log10 = double(T.theta_5);
T.q3_log10 = double(T.theta_6);

T.r_u1_log10  = double(T.theta_7);
T.r_u2_log10  = double(T.theta_8);
T.r_u3_log10  = double(T.theta_9);

T.r_du1_log10 = double(T.theta_10);
T.r_du2_log10 = double(T.theta_11);
T.r_du3_log10 = double(T.theta_12);

T.Q_x1    = 10.^T.q1_log10;
T.Q_x2    = 10.^T.q2_log10;
T.Q_x3    = 10.^T.q3_log10;

T.R_u_x1  = 10.^T.r_u1_log10;
T.R_u_x2  = 10.^T.r_u2_log10;
T.R_u_x3  = 10.^T.r_u3_log10;

T.R_du_x1 = 10.^T.r_du1_log10;
T.R_du_x2 = 10.^T.r_du2_log10;
T.R_du_x3 = 10.^T.r_du3_log10;

% Horizon mapping
T.m = T.theta_m + 1;
T.p = T.theta_p + T.m;

%% --- Pareto mask (minimise J_track=SSE and J_TV=SSdU) ---
J_track = double(T.SSE);
J_TV    = double(T.SSdU);

n = numel(J_track);
isPareto = true(n,1);

for i = 1:n
    dominated = (J_track <= J_track(i)) & (J_TV <= J_TV(i)) & ...
                ( (J_track < J_track(i)) | (J_TV < J_TV(i)) );
    dominated(i) = false;
    if any(dominated)
        isPareto(i) = false;
    end
end

%% --- Global text interpreter (LaTeX) ---
set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");

%% --- Plot: J_track vs J_TV (no size encoding) ---
fig1 = figure("Color","w");
ax1 = axes(fig1); hold(ax1,"on"); grid(ax1,"on"); box(ax1,"on");

c = double(T.J);
markerSize = 80;

scatter(ax1, J_TV, J_track, markerSize, c, ...
    "filled", ...
    "MarkerEdgeColor","k", ...
    "LineWidth",0.7);

set(ax1, "XScale","log", "YScale","log");

title(ax1, "Pareto frontier");
xlabel(ax1, "$J_{\mathrm{TV}}$");
ylabel(ax1, "$J_{\mathrm{track}}$");

ax1.GridLineStyle = "--";
ax1.GridAlpha = 0.4;

cb = colorbar(ax1);
cb.Label.String = "$J$";

set(ax1, "FontSize", 12);

% Highlight Pareto points
paretoMarkerSize = 150;
scatter(ax1, J_TV(isPareto), J_track(isPareto), paretoMarkerSize, ...
    "MarkerEdgeColor","r", ...
    "MarkerFaceColor","none", ...
    "LineWidth",1.2);

exportgraphics(fig1, outScatterPath, "Resolution", 300);

%% --- Plot: runtime vs iteration (minutes), with iter=20 line ---
fig2 = figure("Color","w");
ax2 = axes(fig2); hold(ax2,"on"); grid(ax2,"on"); box(ax2,"on");

plot(ax2, T.iteration, T.runtime_min, "-o", ...
    "LineWidth", 2.0, "MarkerSize", 6);

xline(ax2, 20, "--", "Optimisation start ($k=20$)", ...
    "LineWidth", 1.3, ...
    "LabelVerticalAlignment","middle", ...
    "LabelHorizontalAlignment","left", ...
    "FontSize", 14);

xlabel(ax2, "$k$ (iteration)", "FontSize", 14);
ylabel(ax2, "$t_{\mathrm{run}}$ (min)", "FontSize", 14);
title(ax2, "Iteration runtime");

ax2.GridLineStyle = "--";
ax2.GridAlpha = 0.4;

set(ax2, "FontSize", 12);

exportgraphics(fig2, outRuntimePath, "Resolution", 300);

%% --- Correlation check: Ru_i vs Rdu_i (i = 1..3) ---
Ru_log = [T.r_u1_log10,  T.r_u2_log10,  T.r_u3_log10];
Rdu_log = [T.r_du1_log10, T.r_du2_log10, T.r_du3_log10];
Ru_lin = [T.R_u_x1, T.R_u_x2, T.R_u_x3];
Rdu_lin = [T.R_du_x1, T.R_du_x2, T.R_du_x3];

pairLabel = ["(1)","(2)","(3)"];
pearson_r = NaN(3,1);
pearson_p = NaN(3,1);
spearman_rho = NaN(3,1);
spearman_p = NaN(3,1);
pearson_r_lin = NaN(3,1);
pearson_p_lin = NaN(3,1);
spearman_rho_lin = NaN(3,1);
spearman_p_lin = NaN(3,1);

for i = 1:3
    [pearson_r(i), pearson_p(i)] = corr(Ru_log(:,i), Rdu_log(:,i), ...
        "Type","Pearson", "Rows","complete");
    [spearman_rho(i), spearman_p(i)] = corr(Ru_log(:,i), Rdu_log(:,i), ...
        "Type","Spearman", "Rows","complete");
    [pearson_r_lin(i), pearson_p_lin(i)] = corr(Ru_lin(:,i), Rdu_lin(:,i), ...
        "Type","Pearson", "Rows","complete");
    [spearman_rho_lin(i), spearman_p_lin(i)] = corr(Ru_lin(:,i), Rdu_lin(:,i), ...
        "Type","Spearman", "Rows","complete");
end

RuRduCorrTable = table(pairLabel(:), pearson_r, pearson_p, spearman_rho, spearman_p, ...
    VariableNames=["pair","pearson_r","pearson_p","spearman_rho","spearman_p"]);
RuRduCorrTableLinear = table(pairLabel(:), pearson_r_lin, pearson_p_lin, spearman_rho_lin, spearman_p_lin, ...
    VariableNames=["pair","pearson_r","pearson_p","spearman_rho","spearman_p"]);

disp("Correlation check between Ru_i and Rdu_i (log10-space):");
disp(RuRduCorrTable);
disp("Correlation check between Ru_i and Rdu_i (linear-space):");
disp(RuRduCorrTableLinear);

fig3 = figure("Color","w");
tiledlayout(fig3, 1, 3, "Padding","compact", "TileSpacing","compact");
for i = 1:3
    ax = nexttile; hold(ax, "on"); grid(ax, "on"); box(ax, "on");
    scatter(ax, Ru_log(:,i), Rdu_log(:,i), 40, double(T.J), "filled", ...
        "MarkerEdgeColor","k", "LineWidth",0.5);
    xlabel(ax, sprintf("$\\log_{10}(R_{u,%d})$", i));
    ylabel(ax, sprintf("$\\log_{10}(R_{\\Delta u,%d})$", i));
    title(ax, sprintf("$r=%.3f,\\ \\rho=%.3f$", pearson_r(i), spearman_rho(i)));
    ax.GridLineStyle = "--";
    ax.GridAlpha = 0.35;
end

%% --- Display Pareto table including tuning weights (Q, R, Rdu) ---
Tp = T(isPareto, :);
[~, ord_tune] = sort(double(Tp.SSE), "descend");
Tp = Tp(ord_tune, :);

Tp.runtime_per_f = NaN(height(Tp),1);
if ismember("runtime_s", string(Tp.Properties.VariableNames)) && ismember("f", string(Tp.Properties.VariableNames))
    Tp.runtime_per_f = double(Tp.runtime_s) ./ max(double(Tp.f), eps);
elseif ismember("runtime_min", string(Tp.Properties.VariableNames)) && ismember("f", string(Tp.Properties.VariableNames))
    Tp.runtime_per_f = double(Tp.runtime_min) ./ max(double(Tp.f), eps);
end

tuningCols = ["timestamp","SSE","SSdU","J","p","m","f", ...
    "runtime_min","runtime_per_f", ...
    "Q_x1","Q_x2","Q_x3", ...
    "R_u_x1","R_u_x2","R_u_x3", ...
    "R_du_x1","R_du_x2","R_du_x3"];

tuningCols = tuningCols(ismember(tuningCols, string(Tp.Properties.VariableNames)));
ParetoTuningTable = Tp(:, tuningCols);

disp("Pareto frontier with tuning weights (Q, R, Rdu):");
disp(ParetoTuningTable);

% Optional: write table
% writetable(ParetoTuningTable, fullfile("results","pareto_points.csv"));
