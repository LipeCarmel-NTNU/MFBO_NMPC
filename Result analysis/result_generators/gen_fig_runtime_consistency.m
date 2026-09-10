
function gen_fig_runtime_consistency(ctx)
%GEN_FIG_RUNTIME_CONSISTENCY Is a runtime jump solver time or clock artifact.
%   Result: results/graphical_results/runtime_consistency.png/.pdf
%   Diagnostic, not a paper figure. Panel a plots t_total against t_nmpc with
%   the y = x reference; b the gap over execution order; c the gap against
%   failed solves. A resumed evaluation reports wall time for its last segment
%   only, so it falls below the reference in a and negative in b.
T = ra_require(ctx, "timeline");
T = ra_select_cases(T, "mf");   % baseline is shown only in the combined frontier and the cumulative runtime
A = T.A; doeCount = T.doeCount; nCases = T.nCases; caseLabels = T.caseLabels;
fontSize = ctx.fontSize; plotColors = ctx.plotColors; accentColor = ctx.accentColor;
caseMarkers = ctx.caseMarkers; graphicsDir = ctx.graphicsDir;

%% Figure: solve time vs wall time, gap timeline, gap vs crashes
fig = figure('Color', 'w');
tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- (a) t_nmpc vs t_total, y = x reference -----------------------------
ax1 = nexttile; hold(ax1, 'on'); grid(ax1, 'off'); box(ax1, 'off');
lim = [0, 1.05 * max(cellfun(@(T) max([T.totalH; T.nmpcH]), A))];
plot(ax1, lim, lim, '-', 'Color', accentColor, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
for k = 1:nCases
    scatter(ax1, A{k}.nmpcH, A{k}.totalH, 46, plotColors(k, :), ...
        caseMarkers(k), 'LineWidth', 1.1, 'DisplayName', caseLabels(k));
end
xlim(ax1, lim); ylim(ax1, lim); axis(ax1, 'square');
xlabel(ax1, 'NMPC solve time $t_{\mathrm{nmpc}}$ (h)');
ylabel(ax1, 'Wall time $t_{\mathrm{total}}$ (h)');
title(ax1, '$\mathbf{a}$');
ax1.TitleHorizontalAlignment = 'left';
legend(ax1, 'Location', 'northwest', 'Box', 'off');

% --- (b) gap vs execution order ------------------------------------------
ax2 = nexttile; hold(ax2, 'on'); grid(ax2, 'off'); box(ax2, 'off');
for k = 1:nCases
    plot(ax2, A{k}.gidx, A{k}.gapMin, '-', 'Color', plotColors(k, :), ...
        'LineWidth', 1.1, 'Marker', caseMarkers(k), 'MarkerSize', 4.5, ...
        'HandleVisibility', 'off');
end
for c = unique(doeCount(doeCount > 0))'
    xline(ax2, c + 0.5, 'k--', 'LineWidth', 2, 'HandleVisibility', 'off');
end
yline(ax2, 0, 'k-', 'LineWidth', 0.75, 'HandleVisibility', 'off');
xlabel(ax2, 'Evaluation (execution order, DOE then BO)');
ylabel(ax2, '$t_{\mathrm{total}} - t_{\mathrm{nmpc}}$ (min)');
title(ax2, '$\mathbf{b}$');
ax2.TitleHorizontalAlignment = 'left';

% --- (c) gap vs number of failed solves ----------------------------------
ax3 = nexttile; hold(ax3, 'on'); grid(ax3, 'off'); box(ax3, 'off');
for k = 1:nCases
    scatter(ax3, A{k}.nBad, A{k}.gapMin, 46, plotColors(k, :), ...
        caseMarkers(k), 'LineWidth', 1.1, 'HandleVisibility', 'off');
end
yline(ax3, 0, 'k-', 'LineWidth', 0.75, 'HandleVisibility', 'off');
xlabel(ax3, 'Failed solves per eval. (flag $\le 0$)');
ylabel(ax3, '$t_{\mathrm{total}} - t_{\mathrm{nmpc}}$ (min)');
title(ax3, '$\mathbf{c}$');
ax3.TitleHorizontalAlignment = 'left';

save_plot_outputs(fig, fullfile(graphicsDir, 'runtime_consistency.png'), ...
    fontSize, 1400, 460);
end
