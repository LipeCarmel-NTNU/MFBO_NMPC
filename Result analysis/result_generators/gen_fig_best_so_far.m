function gen_fig_best_so_far(ctx)
%GEN_FIG_BEST_SO_FAR Running best of each objective over the BO iterations.
%   Result: results/graphical_results/best_so_far_side_by_side.png/.pdf
%
%   Panel a: the lowest J_track seen up to iteration k. Panel b: the lowest
%   J_TV. DOE evaluations are excluded, so k = 1 is the first point the
%   optimiser proposed, and a campaign that starts flat started from a design
%   that already held its best value.
%
%   The two panels track the objectives independently, so the evaluation
%   holding the record in a is usually not the one holding it in b. Neither
%   curve is a frontier: this is a convergence diagnostic, not a Pareto
%   statement. Read it next to pareto_samples_combined, which is.
%
%   Log y on both panels: J_TV spans about two decades across the campaigns
%   and the steps down would be invisible on a linear axis.
%
%   Styled after runtime_cumulative — line per campaign, the baseline dotted
%   in Wong Vermillion, ending where its budget ran out.

    F = ra_require(ctx, "frontier");
    isBase = contains(lower(string(F.caseNames)), "baseline");

    cols   = ["SSE", "SSdU"];
    labels = ["$\min J_{\mathrm{track}}$ so far", "$\min J_{\mathrm{TV}}$ so far"];
    panels = ["a", "b"];

    fig = figure('Color', 'w', 'Name', 'Best Objective So Far');
    tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    for ip = 1:numel(cols)
        ax = nexttile; hold(ax, 'on');
        kMax = 1;
        for k = 1:F.nCases
            T    = F.E{k};
            it   = double(T.iter);
            best = cummin(double(T.(char(cols(ip)))));
            if isBase(k)
                style = ':';
                col   = ctx.baselineColor;
            else
                style = ctx.caseLines(k);
                col   = ctx.plotColors(k, :);
            end
            plot(ax, it, best, style, 'LineWidth', 2.0, 'Color', col, ...
                'DisplayName', F.caseLabels(k));
            kMax = max(kMax, max(it));
        end
        set(ax, 'YScale', 'log', 'FontSize', ctx.fontSize);
        xlim(ax, [1, kMax]);
        xlabel(ax, '$k$ (BO iteration)');
        ylabel(ax, labels(ip));
        title(ax, "$\mathbf{" + panels(ip) + "}$", 'Interpreter', 'latex');
        ax.TitleHorizontalAlignment = 'left';
        grid(ax, 'off'); box(ax, 'off');
        if ip == 1
            legend(ax, 'Location', 'northeast');
        end
    end

    save_plot_outputs(fig, fullfile(ctx.graphicsDir, 'best_so_far_side_by_side.png'), ...
        ctx.fontSize, 1200, 460);
end
