
function gen_fig_runtime_iteration(ctx)
%GEN_FIG_RUNTIME_ITERATION Per-evaluation runtime and fidelity over iterations.
%   Result: results/graphical_results/runtime_vs_iteration_side_by_side.png/.pdf
%   Paper:  fig:runtime_iter. Left axis t_iter in hours, right axis z, dashed
%   line at the DOE cutoff. DOE and BO run on one iteration axis per case.
F = ra_require(ctx, "frontier");
F = ra_select_cases(F, "mf");   % baseline is shown only in the combined frontier and the cumulative runtime
E = F.E; D = F.D; nCases = F.nCases;
fontSize = ctx.fontSize; plotColors = ctx.plotColors; accentColor = ctx.accentColor;
graphicsDir = ctx.graphicsDir;

%% Figure 2: iteration runtime + fidelity per case (DOE + BO timeline)
tMaxH = 0;
for k = 1:nCases
    tMaxH = max(tMaxH, max([double(D{k}.t_total); double(E{k}.t_total)]) / 3600);
end
fig2 = figure('Color', 'w', 'Name', 'Iteration Runtime and Fidelity by Case');
tiledlayout(fig2, 1, nCases, 'Padding', 'compact', 'TileSpacing', 'compact');
for k = 1:nCases
    nDoe    = height(D{k});
    iterAll = [double(D{k}.iter); nDoe + double(E{k}.iter)];
    tAllH   = [double(D{k}.t_total); double(E{k}.t_total)] / 3600;
    zAll    = [double(D{k}.z); double(E{k}.z)];

    ax = nexttile; hold(ax, 'on');
    yyaxis(ax, 'left');
    plot(ax, iterAll, tAllH, '-', 'LineWidth', 2.0, 'Color', plotColors(k, :));
    ax.YColor = 'k';
    ylim(ax, [0, 1.05 * max(tMaxH, eps)]);
    xline(ax, nDoe, '--', 'LineWidth', 2.0, 'Color', 'k', 'Alpha', 1);
    ylabel(ax, '$t_{\mathrm{iter}}$ (h)');
    yyaxis(ax, 'right');
    plot(ax, iterAll, zAll, 'o', 'LineWidth', 2.0, 'MarkerSize', 4, 'Color', accentColor);
    ax.YColor = 'k';
    ylim(ax, [0, 1]);
    ylabel(ax, '$z$ (dimensionless)');
    xlabel(ax, '$k$ (iteration)');
    xlim(ax, [1, max(1, max(iterAll))]);
    title(ax, "$\mathbf{" + char('a' + k - 1) + "}$", 'Interpreter', 'latex');
    ax.TitleHorizontalAlignment = 'left';
    set(ax, 'FontSize', fontSize);
    grid(ax, 'off'); box(ax, 'off');
    axes(ax); %#ok<LAXES>
    yyaxis(ax, 'left');  format_tick(0, 1);
    yyaxis(ax, 'right'); yticks(ax, 0:0.2:1); format_tick(0, 1);
end
save_plot_outputs(fig2, fullfile(graphicsDir, 'runtime_vs_iteration_side_by_side.png'), ...
    fontSize, 600 * nCases, 460);
end
