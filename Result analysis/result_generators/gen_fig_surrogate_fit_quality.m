
function gen_fig_surrogate_fit_quality(ctx)
%GEN_FIG_SURROGATE_FIT_QUALITY Fit loss per objective and refit cost against vintage.
%   Result: results/graphical_results/surrogate_fit_quality.png/.pdf
%   Diagnostic, not a paper figure. Reads only the published vintages, so a
%   case that never refits its surrogate has nothing to contribute here.
[S, caseNames, caseLabels] = ra_surrogate_table(ctx);
nCases = numel(caseNames);
targets      = ["SSE", "SSdU"];
targetLabels = ["J_{\mathrm{track}}", "J_{\mathrm{TV}}"];
plotColors = ctx.plotColors; caseMarkers = ctx.caseMarkers;
fontSize = ctx.fontSize; seqMap = ctx.seqMap; graphics_dir = ctx.graphicsDir;

%% Figure 2: fit quality (loss) and refit cost (wall time)
fig2 = figure('Color', 'w', 'Name', 'Surrogate Fit Quality by Vintage');
tiledlayout(fig2, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
panelLabels2 = ["a", "b", "c"];
lossCols   = ["SSE_loss", "SSdU_loss", "fit_wall_s"];
lossLabels = ["$\mathcal{L}$ ($J_{\mathrm{track}}$)", ...
              "$\mathcal{L}$ ($J_{\mathrm{TV}}$)", ...
              "$t_{\mathrm{fit}}$ (s)"];

for ip = 1:numel(lossCols)
    ax = nexttile; hold(ax, 'on');
    for k = 1:nCases
        T = S(S.case == caseNames{k}, :);
        plot(ax, double(T.vintage), double(T.(lossCols(ip))), '-', ...
            'Marker', caseMarkers(k), 'MarkerSize', 7, ...
            'MarkerFaceColor', 'w', 'LineWidth', 2.0, ...
            'Color', plotColors(k, :), 'DisplayName', caseLabels(k));
    end
    xlabel(ax, '$v$ (vintage)');
    ylabel(ax, lossLabels(ip));
    title(ax, "$\mathbf{" + panelLabels2(ip) + "}$", 'Interpreter', 'latex');
    ax.TitleHorizontalAlignment = 'left';
    set(ax, 'FontSize', fontSize);
    grid(ax, 'off'); box(ax, 'off');
    if ip == 1
        legend(ax, 'Location', 'best');
    end
end
save_plot_outputs(fig2, fullfile(graphics_dir, 'surrogate_fit_quality.png'), fontSize, 1400, 420);
end
