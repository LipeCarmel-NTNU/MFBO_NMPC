
function gen_fig_objective_space_z(ctx)
%GEN_FIG_OBJECTIVE_SPACE_Z Objective space per case, shaded by fidelity z.
%   Result: results/graphical_results/sse_vs_ssdu_side_by_side_z.png/.pdf
%   Paper:  fig:obj_z. One panel per case, log-log, shared axis limits, the
%   case frontier as open accent rings over a continuum guide.
%   Both arms are drawn: SF comes out uniformly at the top of the colorbar,
%   which is the visual statement that it never varied its horizon.
F = ra_require(ctx, "frontier");
E = F.E; isPareto = F.isPareto; nCases = F.nCases;
xLimAll = F.xLimAll; yLimAll = F.yLimAll; zLo = F.zLo;
fontSize = ctx.fontSize; seqMap = ctx.seqMap; accentColor = ctx.accentColor;
graphicsDir = ctx.graphicsDir;

%% Figure 1: SSE vs SSdU per case, color mapped by fidelity z
fig1 = figure('Color', 'w', 'Name', 'Pareto SSE vs SSdU by Case');
tiledlayout(fig1, 1, nCases, 'Padding', 'compact', 'TileSpacing', 'compact');
for k = 1:nCases
    T  = E{k};
    ax = nexttile; hold(ax, 'on');
    plot_pareto_continuum(ax, double(T.SSdU(isPareto{k})), double(T.SSE(isPareto{k})), ...
        accentColor, xLimAll, yLimAll);
    scatter(ax, double(T.SSdU), double(T.SSE), 80, double(T.z), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
    scatter(ax, double(T.SSdU(isPareto{k})), double(T.SSE(isPareto{k})), 170, accentColor, ...
        'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', accentColor, 'LineWidth', 1.2);

    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', fontSize);
    xlim(ax, xLimAll); ylim(ax, yLimAll);
    colormap(ax, seqMap);
    caxis(ax, [zLo, 1]);
    xlabel(ax, '$J_{\mathrm{TV}}$');
    ylabel(ax, '$J_{\mathrm{track}}$');
    title(ax, "$\mathbf{" + char('a' + k - 1) + "}$", 'Interpreter', 'latex');
    ax.TitleHorizontalAlignment = 'left';
    grid(ax, 'off'); box(ax, 'off');
    cb = colorbar(ax);
    cb.Label.String = '$z$ (dimensionless)';
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = fontSize;
end
save_plot_outputs(fig1, fullfile(graphicsDir, 'sse_vs_ssdu_side_by_side_z.png'), ...
    fontSize, 600 * nCases, 460);
end
