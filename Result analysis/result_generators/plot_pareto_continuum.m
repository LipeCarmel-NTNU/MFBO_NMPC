function h = plot_pareto_continuum(ax, x, y, curveColor, xBounds, yBounds)
%PLOT_PARETO_CONTINUUM Smooth monotone visual guide of the frontier trend.
[xSort, ord] = sort(x(:), 'ascend');
ySort = y(ord);

if isempty(xSort)
    h = plot(ax, xSort, ySort, '-', 'Color', curveColor, 'LineWidth', 2.0);
    return
end

if numel(xSort) >= 3
    % Interpolate in log-log domain for a smooth non-segmented visual guide.
    lx = log10(xSort);
    ly = log10(ySort);
    lxqIn = linspace(min(lx), max(lx), 220);
    lyqIn = pchip(lx, ly, lxqIn);
    xq = (10.^lxqIn).';
    yq = (10.^lyqIn).';
else
    % Fewer than three points is not enough to interpolate, so the guide is
    % the raw segment. It used to return here, which meant a two-point
    % frontier -- the single-fidelity arm has exactly two -- was drawn without
    % the extensions below and so showed no dominated region, while a
    % multi-point frontier in the neighbouring panel did. The extensions are
    % what mark the dominated region, not the interpolation, so they have to
    % apply at every frontier size.
    xq = xSort(:);
    yq = ySort(:);
end

% Left extension: vertical at the leftmost Pareto point.
if nargin >= 6 && ~isempty(yBounds)
    yTop = max(yBounds);
    if yTop > yq(1)
        xq = [xq(1); xq];
        yq = [yTop; yq];
    end
end

% Right extension: horizontal from the rightmost Pareto point.
if nargin >= 5 && ~isempty(xBounds)
    xRight = max(xBounds);
    if xRight > xq(end)
        xq = [xq; xRight];
        yq = [yq; yq(end)];
    end
end

h = plot(ax, xq, yq, '-', 'Color', curveColor, 'LineWidth', 2.0);
plot(ax, xSort, ySort, 'o', ...
    'Color', curveColor, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', curveColor, 'LineWidth', 1.4);
end
