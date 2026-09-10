function [finalMask, Tall] = plot_combined_pareto_base(ax, E, Tp, plotColors, caseMarkers, accentColor, xBounds, yBounds)
%PLOT_COMBINED_PARETO_BASE All-case sample cloud + per-case frontiers + combined frontier.
Tall = vertcat(E{:});   % BO evaluations only; DOE excluded structurally
scatter(ax, double(Tall.SSdU), double(Tall.SSE), 18, ...
    'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'Optimization samples');

finalMask = compute_pareto_mask(double(Tall.SSE), double(Tall.SSdU));
Tf = Tall(finalMask, :);
plot_pareto_continuum(ax, double(Tf.SSdU), double(Tf.SSE), accentColor, xBounds, yBounds);

for k = 1:numel(Tp)
    scatter(ax, double(Tp{k}.SSdU), double(Tp{k}.SSE), 80 + 10 * (k - 1), plotColors(k, :), ...
        caseMarkers(k), 'MarkerFaceColor', plotColors(k, :), ...
        'MarkerEdgeColor', plotColors(k, :), 'LineWidth', 1.4);
end

scatter(ax, double(Tf.SSdU), double(Tf.SSE), 300, accentColor, ...
    'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', accentColor, 'LineWidth', 2);
end
