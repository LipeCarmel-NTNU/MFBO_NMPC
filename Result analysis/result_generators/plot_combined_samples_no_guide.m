function plot_combined_samples_no_guide(ax, E, Tp, plotColors, caseMarkers)
%PLOT_COMBINED_SAMPLES_NO_GUIDE Sample cloud + per-case frontier markers, no guideline/ring.
Tall = vertcat(E{:});
scatter(ax, double(Tall.SSdU), double(Tall.SSE), 18, ...
    'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
for k = 1:numel(Tp)
    scatter(ax, double(Tp{k}.SSdU), double(Tp{k}.SSE), 80 + 10 * (k - 1), plotColors(k, :), ...
        caseMarkers(k), 'MarkerFaceColor', plotColors(k, :), ...
        'MarkerEdgeColor', plotColors(k, :), 'LineWidth', 1.4);
end
end
