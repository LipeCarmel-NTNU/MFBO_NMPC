
function gen_fig_surrogate_phi_curves(ctx)
%GEN_FIG_SURROGATE_PHI_CURVES The cost fraction phi(z) = I_z(a,b) per published vintage.
%   Result: results/graphical_results/surrogate_phi_curves.png/.pdf
%   Diagnostic, not a paper figure. Reads only the published vintages, so a
%   case that never refits its surrogate has nothing to contribute here.
[S, caseNames, caseLabels] = ra_surrogate_table(ctx);
nCases = numel(caseNames);
targets      = ["SSE", "SSdU"];
targetLabels = ["J_{\mathrm{track}}", "J_{\mathrm{TV}}"];
plotColors = ctx.plotColors; caseMarkers = ctx.caseMarkers;
fontSize = ctx.fontSize; seqMap = ctx.seqMap; graphics_dir = ctx.graphicsDir;

%% Figure 3: phi(z) curves per vintage, colored by vintage (navia)
fig3 = figure('Color', 'w', 'Name', 'phi(z) Curves by Vintage');
tiledlayout(fig3, numel(targets), nCases, 'Padding', 'compact', 'TileSpacing', 'compact');
zq = linspace(0, 1, 201);
vMax = max(double(S.vintage));   % shared scale so all colorbars agree
panelIdx = 0;

for it = 1:numel(targets)
    for k = 1:nCases
        panelIdx = panelIdx + 1;
        ax = nexttile; hold(ax, 'on');
        T    = S(S.case == caseNames{k}, :);
        vAll = double(T.vintage);
        for r = 1:height(T)
            a = double(T.(targets(it) + "_a")(r));
            b = double(T.(targets(it) + "_b")(r));
            if vMax > 0
                cIdx = 1 + round((vAll(r) / vMax) * (size(seqMap, 1) - 1));
            else
                cIdx = 1;
            end
            plot(ax, zq, phi_eval(zq, a, b), '-', 'LineWidth', 1.8, ...
                'Color', seqMap(cIdx, :));
        end
        plot(ax, [0, 1], [0, 1], '--', 'LineWidth', 1.2, 'Color', 'k'); % phi(z) = z reference
        xlabel(ax, '$z$ (dimensionless)');
        ylabel(ax, "$\varphi(z)$ ($" + targetLabels(it) + "$)");
        title(ax, "$\mathbf{" + char('a' + panelIdx - 1) + "}$ " + caseLabels(k), ...
            'Interpreter', 'latex');
        ax.TitleHorizontalAlignment = 'left';
        xlim(ax, [0, 1]); ylim(ax, [0, 1]);
        set(ax, 'FontSize', fontSize);
        grid(ax, 'off'); box(ax, 'off');
        colormap(ax, seqMap);
        caxis(ax, [0, max(vMax, 1)]);
        if k == nCases
            cb = colorbar(ax);
            cb.Label.String = '$v$ (vintage)';
            cb.Label.Interpreter = 'latex';
            cb.TickLabelInterpreter = 'latex';
            cb.FontSize = fontSize;
        end
    end
end
save_plot_outputs(fig3, fullfile(graphics_dir, 'surrogate_phi_curves.png'), fontSize, 1200, 850);
end
