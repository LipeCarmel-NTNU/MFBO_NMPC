function gen_fig_pareto_combined(ctx)
%GEN_FIG_PARETO_COMBINED Pooled samples and per-case frontiers in one axes.
%   Result: results/graphical_results/pareto_samples_combined.png/.pdf
%   Paper:  fig:pareto_samples, which still includes the older stem
%           pareto_samples_run1_run2.pdf.
%
%   The multi-fidelity cases are pooled and their merged frontier is drawn in
%   the accent colour, exactly as before. A single-fidelity baseline case is
%   NOT pooled: it is overlaid as orange crosses, its frontier points larger.
%   Keeping it out of the pooling means the merged frontier stays the
%   published two-case result and the figure gains the comparison without
%   restating it, and it is also the honest call while the baseline has fewer
%   BO evaluations than the cases it is being shown against.

    F = ra_require(ctx, "frontier");

    isBase = contains(lower(string(F.caseNames)), "baseline");
    mf     = find(~isBase);
    base   = find(isBase);

    fontSize = ctx.fontSize;
    xLimAll  = F.xLimAll;
    yLimAll  = F.yLimAll;

    fig3 = figure('Color', 'w', 'Name', 'Combined Pareto Samples');
    ax = axes(fig3); hold(ax, 'on');

    if ~isempty(mf)
        plot_combined_pareto_base(ax, F.E(mf), F.Tp(mf), ctx.plotColors(mf, :), ...
            ctx.caseMarkers(mf), ctx.accentColor, xLimAll, yLimAll);
    end

    baselineColor = ctx.baselineColor;   % Wong Vermillion, see ra_context
    for k = base(:)'
        Tb = F.E{k};
        scatter(ax, double(Tb.SSdU), double(Tb.SSE), 70, baselineColor, 'x', ...
            'LineWidth', 1.3, 'DisplayName', F.caseLabels(k) + " samples");
        scatter(ax, double(F.Tp{k}.SSdU), double(F.Tp{k}.SSE), 200, baselineColor, 'x', ...
            'LineWidth', 2.4, 'DisplayName', F.caseLabels(k) + " frontier");
        fprintf(['%s overlaid as orange crosses: %d samples, %d on its own frontier. ' ...
                 'It is excluded from the merged frontier.\n'], ...
            F.caseLabels(k), height(Tb), height(F.Tp{k}));
    end

    apply_combined_axes_style(ax, fontSize, xLimAll, yLimAll);
    save_plot_outputs(fig3, fullfile(ctx.graphicsDir, 'pareto_samples_combined.png'), ...
        fontSize, 920, 520);
end
