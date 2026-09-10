function h = plot_pareto_continuum_line_only(ax, x, y, curveColor, xBounds, yBounds)
%PLOT_PARETO_CONTINUUM_LINE_ONLY Smooth pchip Pareto guide without point markers.
[xSort, ord] = sort(x(:), 'ascend');
ySort = y(ord);

if numel(xSort) < 3
    h = plot(ax, xSort, ySort, '-', 'Color', curveColor, 'LineWidth', 2.0);
    return
end

lx = log10(xSort);
ly = log10(ySort);
lxqIn = linspace(min(lx), max(lx), 220);
lyqIn = pchip(lx, ly, lxqIn);
xq = (10.^lxqIn).';
yq = (10.^lyqIn).';

if nargin >= 6 && ~isempty(yBounds)
    yTop = max(yBounds);
    if yTop > yq(1)
        xq = [xq(1); xq];
        yq = [yTop; yq];
    end
end
if nargin >= 5 && ~isempty(xBounds)
    xRight = max(xBounds);
    if xRight > xq(end)
        xq = [xq; xRight];
        yq = [yq; yq(end)];
    end
end

h = plot(ax, xq, yq, '-', 'Color', curveColor, 'LineWidth', 2.0);
end
