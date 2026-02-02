% --- Global text interpreter (LaTeX) ---
set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");

%% --- Plot: J_track vs J_TV (no size encoding) ---
fig1 = figure("Color","w");
ax1 = axes(fig1); hold(ax1,"on"); grid(ax1,"on"); box(ax1,"on");

x = double(T.SSE);   % J_track
y = double(T.SSdU);  % J_TV
c = double(T.J);

markerSize = 80;     % fixed marker area (points^2)

scatter(ax1, x, y, markerSize, c, ...
    "filled", ...
    "MarkerEdgeColor","k", ...
    "LineWidth",0.7);

set(ax1, "XScale","log", "YScale","log");

title(ax1, "Pareto Frontier");
xlabel(ax1, "$J_{\mathrm{track}}$");
ylabel(ax1, "$J_{\mathrm{TV}}$");

ax1.GridLineStyle = "--";
ax1.GridAlpha = 0.4;

cb = colorbar(ax1);
cb.Label.String = "J";

set(ax1, "FontSize", 12);

% Highlight Pareto frontier with open red circles
paretoMarkerSize = 150;
scatter(ax1, x(isPareto), y(isPareto), paretoMarkerSize, ...
    "MarkerEdgeColor","r", ...
    "MarkerFaceColor","none", ...
    "LineWidth",1.2);

exportgraphics(fig1, outScatterPath, "Resolution", 300);

%% --- Plot: runtime vs iteration (minutes), with iter=20 line ---
fig2 = figure("Color","w");
ax2 = axes(fig2); hold(ax2,"on"); grid(ax2,"on"); box(ax2,"on");

plot(ax2, T.iteration, T.runtime_min, "-o", ...
    "LineWidth", 2.0, "MarkerSize", 6);

xline(ax2, 20, "--", "Optimisation start (k=20)", ...
    "LineWidth", 1.3, ...
    "LabelVerticalAlignment","middle", ...
    "LabelHorizontalAlignment","left", ...
    "FontSize",14);

xlabel(ax2, "$k$ (iteration)", "FontSize",14);
ylabel(ax2, "$t_{\mathrm{run}}$ (min)", "FontSize",14);
title(ax2, "Iteration runtime");

ax2.GridLineStyle = "--";
ax2.GridAlpha = 0.4;


set(ax2, "FontSize", 12);

exportgraphics(fig2, outRuntimePath, "Resolution", 300);
