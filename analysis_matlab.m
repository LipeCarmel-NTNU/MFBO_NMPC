% analyze_mfbo_results.m
% Load MFBO optimisation results from results/results.csv (same schema as matlab_interface.read_results)
% Decode theta columns and reproduce:
%   (1) SSE vs SSdU scatter, coloured by J, sized by theta_p
%   (2) Runtime (minutes) vs iteration with vertical line at iter=20
%
% Outputs (PNG):
%   results/sse_vs_ssdu.png
%   results/runtime_vs_iteration.png

clear; close all; clc

%% --- Paths / options ---
csvPath         = fullfile("results", "results.csv");
outScatterPath  = fullfile("results", "sse_vs_ssdu.png");
outRuntimePath  = fullfile("results", "runtime_vs_iteration.png");
showFigures     = false;     % set true to keep figures open

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

%% --- Parse / derive fields (equivalent to Python build_results_dataframe) ---
% timestamp format: YYYYmmdd_HHMMSS
T.timestamp_dt = datetime(T.timestamp, "InputFormat","yyyyMMdd_HHmmss");

T.iteration   = (1:height(T)).';
T.runtime_min = T.runtime_s / 60;

% Rename theta columns (keep originals too)
T.f        = double(T.theta_1);
T.theta_p  = double(T.theta_2);
T.theta_m  = double(T.theta_3);
T.q1_log10 = double(T.theta_4);
T.q2_log10 = double(T.theta_5);
T.q3_log10 = double(T.theta_6);

T.r_u1_log10  = double(T.theta_7);
T.r_u2_log10  = double(T.theta_8);
T.r_u3_log10  = double(T.theta_9);

T.r_du1_log10 = double(T.theta_10);
T.r_du2_log10 = double(T.theta_11);
T.r_du3_log10 = double(T.theta_12);

% Integral parameters
T.theta_p = round(T.theta_p);
T.theta_m = round(T.theta_m);

% Derived diagonal weights (SERIES_LABELS = {"x1","x2","x3"})
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

%% --- Plot: SSE vs SSdU (no size encoding) ---
fig1 = figure("Color","w");
ax1 = axes(fig1); hold(ax1,"on"); grid(ax1,"on"); box(ax1,"on");

x = double(T.SSE);
y = double(T.SSdU);
c = double(T.J);

markerSize = 80;   % fixed marker area (points^2)

sc = scatter(ax1, x, y, markerSize, c, ...
    "filled", ...
    "MarkerEdgeColor","k", ...
    "LineWidth",0.7);

set(ax1, "XScale","log", "YScale","log");

title(ax1, "SSE vs SSdU");
xlabel(ax1, "SSE (State Tracking Error)");
ylabel(ax1, "SSdU (Control Effort)");

ax1.GridLineStyle = "--";
ax1.GridAlpha = 0.4;

cb = colorbar(ax1);
cb.Label.String = "J (Cost)";

set(ax1, "FontSize", 12);


% %--- Pareto frontier extraction (minimise SSE and SSdU) and table display ---
% Assumes table T already exists in workspace from the earlier script, with:
%   T.SSE, T.SSdU, T.p, T.m, T.J, T.iteration, T.timestamp_dt, T.runtime_min

% Ensure numeric vectors
SSE  = double(T.SSE);
SSdU = double(T.SSdU);

n = numel(SSE);
isPareto = true(n,1);

% Nondominated set for 2 objectives (minimise both):
% Point i is dominated if ∃ j with (SSE_j <= SSE_i) & (SSdU_j <= SSdU_i) and at least one strict.
for i = 1:n
    dom = (SSE <= SSE(i)) & (SSdU <= SSdU(i)) & ( (SSE < SSE(i)) | (SSdU < SSdU(i)) );
    dom(i) = false; % ignore self
    if any(dom)
        isPareto(i) = false;
    end
end

% Extract and sort frontier for readability (decreasing SSE)
Tp = T(isPareto, :);
[~, ord] = sort(double(Tp.SSE), "descend");
Tp = Tp(ord, :);

% Select columns to display (edit as needed)
cols = ["iteration","timestamp_dt","SSE","SSdU","J","p","m","theta_p","theta_m","runtime_min","f"];
cols = cols(ismember(cols, string(Tp.Properties.VariableNames))); % keep only those that exist

ParetoTable = Tp(:, cols);

% Display in Command Window
disp("Pareto frontier points (minimise SSE and SSdU):");
disp(ParetoTable);

% If you want indices in the original table:
paretoIdx = find(isPareto);
disp("Row indices in T for Pareto points:");
disp(paretoIdx);

% Optionally write to CSV:
% writetable(ParetoTable, fullfile("results","pareto_points.csv"));

% Highlight Pareto frontier on the log-log scatter plot with open circles
paretoMarkerSize = 150;
scatter(ax1, x(isPareto), y(isPareto), paretoMarkerSize, ...
    "MarkerEdgeColor","r", ...
    "MarkerFaceColor","none", ...
    "LineWidth",1.2);

exportgraphics(fig1, fullfile("results","sse_vs_ssdu.png"), "Resolution", 300);
