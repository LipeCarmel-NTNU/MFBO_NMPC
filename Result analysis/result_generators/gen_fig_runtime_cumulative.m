function gen_fig_runtime_cumulative(ctx)
%GEN_FIG_RUNTIME_CUMULATIVE Cumulative wall-clock runtime per campaign.
%   Result: results/graphical_results/runtime_cumulative.png/.pdf
%   Paper:  fig:runtime_cum, which still includes the older stem
%           runtime_cumulative_run1_run2.pdf.
%
%   One of the two figures the single-fidelity baseline appears in, drawn in
%   Wong Vermillion and dotted. Its curve ends where its wall-clock budget ran
%   out, above both multi-fidelity curves at their full iteration count, which
%   is the comparison this figure exists to make.
%
%   The quantity is t_total, the published definition of this figure. t_nmpc
%   is the quantity that compares campaigns cleanly (see
%   gen_tbl_campaign_runtime); swapping the two cumsum terms below would
%   change what the paper's caption claims, so it is left alone.

    F = ra_require(ctx, "frontier");
    E = F.E; D = F.D; nCases = F.nCases; caseLabels = F.caseLabels;
    isBase = contains(lower(string(F.caseNames)), "baseline");

    fig4 = figure('Color', 'w', 'Name', 'Cumulative Runtime by Case');
    ax = axes(fig4); hold(ax, 'on');
    xMax = 1;
    for k = 1:nCases
        nDoe    = height(D{k});
        iterAll = [double(D{k}.iter); nDoe + double(E{k}.iter)];
        cumH    = cumsum([double(D{k}.t_total); double(E{k}.t_total)], 'omitnan') / 3600;
        if isBase(k)
            style = ':';
            col   = ctx.baselineColor;
        else
            style = ctx.caseLines(k);
            col   = ctx.plotColors(k, :);
        end
        plot(ax, iterAll, cumH, style, 'LineWidth', 2.0, ...
            'Color', col, 'DisplayName', caseLabels(k));
        xMax = max(xMax, max(iterAll));
    end
    xline(ax, height(D{1}), '--', 'LineWidth', 2.0, 'Color', 'k', 'Alpha', 1, ...
        'HandleVisibility', 'off');
    xlabel(ax, '$k$ (iteration)');
    ylabel(ax, '$t_{\mathrm{run}}$ (h)');
    xlim(ax, [1, xMax]);
    legend(ax, 'Location', 'northwest');
    set(ax, 'FontSize', ctx.fontSize);
    grid(ax, 'off'); box(ax, 'off');
    save_plot_outputs(fig4, fullfile(ctx.graphicsDir, 'runtime_cumulative.png'), ...
        ctx.fontSize, 920, 520);
end
