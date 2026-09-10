
function gen_fig_runtime_cumulative(ctx)
%GEN_FIG_RUNTIME_CUMULATIVE Cumulative wall-clock runtime per case.
%   Result: results/graphical_results/runtime_cumulative.png/.pdf
%   Paper:  fig:runtime_cum, which still includes the older stem
%           runtime_cumulative_run1_run2.pdf.
%   A campaign stopped on a time budget ends its curve early, which is the
%   comparison this figure exists to make.
F = ra_require(ctx, "frontier");
E = F.E; D = F.D; nCases = F.nCases; caseLabels = F.caseLabels;
fontSize = ctx.fontSize; plotColors = ctx.plotColors; caseLines = ctx.caseLines;
graphicsDir = ctx.graphicsDir;

%% Figure 4: cumulative runtime per case
fig4 = figure('Color', 'w', 'Name', 'Cumulative Runtime by Case');
ax = axes(fig4); hold(ax, 'on');
xMax = 1;
for k = 1:nCases
    nDoe    = height(D{k});
    iterAll = [double(D{k}.iter); nDoe + double(E{k}.iter)];
    cumH    = cumsum([double(D{k}.t_total); double(E{k}.t_total)], 'omitnan') / 3600;
    plot(ax, iterAll, cumH, caseLines(k), 'LineWidth', 2.0, ...
        'Color', plotColors(k, :), 'DisplayName', caseLabels(k));
    xMax = max(xMax, max(iterAll));
end
xline(ax, height(D{1}), '--', 'LineWidth', 2.0, 'Color', 'k', 'Alpha', 1, ...
    'HandleVisibility', 'off');
xlabel(ax, '$k$ (iteration)');
ylabel(ax, '$t_{\mathrm{run}}$ (h)');
xlim(ax, [1, xMax]);
legend(ax, 'Location', 'northwest');
set(ax, 'FontSize', fontSize);
grid(ax, 'off'); box(ax, 'off');
save_plot_outputs(fig4, fullfile(graphicsDir, 'runtime_cumulative.png'), fontSize, 920, 520);
end
