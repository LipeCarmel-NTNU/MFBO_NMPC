
function gen_fig_surrogate_coefficients(ctx)
%GEN_FIG_SURROGATE_COEFFICIENTS Fitted a, b and lambda against vintage.
%   Result: results/graphical_results/surrogate_coefficients.png/.pdf
%   Diagnostic, not a paper figure. Reads only the published vintages, so a
%   case that never refits its surrogate has nothing to contribute here.
[S, caseNames, caseLabels] = ra_surrogate_table(ctx);
nCases = numel(caseNames);
targets      = ["SSE", "SSdU"];
targetLabels = ["J_{\mathrm{track}}", "J_{\mathrm{TV}}"];
plotColors = ctx.plotColors; caseMarkers = ctx.caseMarkers;
fontSize = ctx.fontSize; seqMap = ctx.seqMap; graphics_dir = ctx.graphicsDir;

%% Figure 1: coefficient evolution (a, b, lambda) per target
fig1 = figure('Color', 'w', 'Name', 'Surrogate Coefficients by Vintage');
tiledlayout(fig1, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
coeffs      = ["a", "b", "lambda"];
coeffLabels = ["a", "b", "\lambda"];
panelLabels = ["a", "b", "c"; "d", "e", "f"];

for it = 1:numel(targets)
    for ic = 1:numel(coeffs)
        ax = nexttile; hold(ax, 'on');
        col = targets(it) + "_" + coeffs(ic);
        for k = 1:nCases
            T = S(S.case == caseNames{k}, :);
            plot(ax, double(T.vintage), double(T.(col)), '-', ...
                'Marker', caseMarkers(k), 'MarkerSize', 7, ...
                'MarkerFaceColor', 'w', 'LineWidth', 2.0, ...
                'Color', plotColors(k, :), 'DisplayName', caseLabels(k));
        end
        xlabel(ax, '$v$ (vintage)');
        ylabel(ax, "$" + coeffLabels(ic) + "$ ($" + targetLabels(it) + "$)");
        title(ax, "$\mathbf{" + panelLabels(it, ic) + "}$", 'Interpreter', 'latex');
        ax.TitleHorizontalAlignment = 'left';
        set(ax, 'FontSize', fontSize);
        grid(ax, 'off'); box(ax, 'off');
        if it == 1 && ic == 1
            legend(ax, 'Location', 'best');
        end
    end
end
save_plot_outputs(fig1, fullfile(graphics_dir, 'surrogate_coefficients.png'), fontSize, 1400, 700);
end
