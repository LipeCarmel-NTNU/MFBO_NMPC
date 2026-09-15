%% PICK_SEQ_COLORMAP  Render sequential-colormap candidates side by side.
%
% A chooser, not a result. It writes nothing and is not registered in
% main_results. Run it, look at the figure, pick a row, and tell me which.
%
% Each row is one candidate, drawn four ways:
%   1. the ramp itself, left = low value
%   2. the objective-space scatter with the colour axis floored at 0.5, which
%      is what gen_fig_objective_space_z does today
%   3. the same scatter with the colour axis spanning the data
%   4. the phi(z) vintage curves - thin lines on white, the harder case
%
% Columns 2 and 3 are the same points twice. z runs 0.833 to 1.0 in case2_v3
% and is exactly 1 for every baseline row, so a floor at 0.5 squeezes every
% marker into the top third of the ramp and they all come out one colour. That
% is a separate defect from the choice of map, and the pair of columns is there
% to keep the two apart while judging.
%
% Maps come from dependencies/plot_utils/colormaps/seq_maps.mat - see the
% README beside it for provenance and the citations each one needs.
%
% Every candidate is truncated to its own usable span: the script finds where
% the map's contrast against white falls through crFloor and cuts there, so the
% rows differ in hue path rather than in how much was cut. Orientation is
% normalised first, so position 0 is the dark end of every map, and the result
% is reversed so the HIGH value takes the dark end.

clear; close all; clc;

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(fullfile(root, 'dependencies')));
addpath(fullfile(here, 'result_generators'));

%% ------------------------------------------------------------------ knobs
%  map name     nLevels   crFloor      nLevels 0 = continuous
%                                      crFloor 0 = no truncation
CANDIDATES = {
    'navia'         0       0
    'navia'         0     3.0
    'batlow'        0     3.0
    'bamako'        0     3.0
    'tokyo'         0     3.0
    'acton'         0     3.0
    'cividis'       0     3.0
    'oslo'          0     3.0
    'grayC'         0     3.0
    'navia'         7     3.0
};

MAPFILE = fullfile(root, 'dependencies', 'plot_utils', 'colormaps', 'seq_maps.mat');
CLIM_LO = 0.5;                  % the floor gen_fig_objective_space_z uses today

%% -------------------------------------------------------------- the data
assert(isfile(MAPFILE), 'Colormap file not found: %s', MAPFILE);
M = load(MAPFILE);

[zs, xs, ys] = local_scatter_data(root);
[zc, curves] = local_curve_data(here);

tFloor = (zs - min(CLIM_LO, min(zs))) / max(1 - min(CLIM_LO, min(zs)), eps);
tData  = (zs - min(zs))              / max(max(zs) - min(zs),          eps);

%% ------------------------------------------------------------ the figure
nC = size(CANDIDATES, 1);
fig = figure('Color', 'w', 'Name', 'Sequential colormap candidates');
tiledlayout(fig, nC, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

fprintf('\n%-22s %7s %7s %9s %9s %8s\n', 'candidate', 'L* lo', 'L* hi', ...
    'cr worst', 'cr best', 'min dL*');
fprintf('%s\n', repmat('-', 1, 70));

for k = 1:nC
    mapName = CANDIDATES{k, 1};
    nLevels = CANDIDATES{k, 2};
    crFloor = CANDIDATES{k, 3};
    isDisc  = nLevels > 0;

    assert(isfield(M, mapName), '%s is not in %s', mapName, MAPFILE);
    [cmap, hi] = local_build(M.(mapName), crFloor, nLevels);

    label = mapName;
    if crFloor > 0; label = sprintf('%s [0, %.2f]', label, hi); end
    if isDisc;      label = sprintf('%s, %d lvl', label, nLevels); end

    [Llo, Lhi, crWorst, crBest, dLmin] = local_metrics(cmap);
    fprintf('%-22s %7.1f %7.1f %8.2f:1 %7.2f:1 %8.1f\n', ...
        label, Llo, Lhi, crWorst, crBest, dLmin);

    ax = nexttile((k-1)*4 + 1);
    image(ax, reshape(cmap, [1 size(cmap,1) 3]));
    set(ax, 'YTick', [], 'XTick', []);
    ylabel(ax, label, 'Interpreter', 'none', 'FontSize', 8, ...
        'Rotation', 0, 'HorizontalAlignment', 'right');
    if k == 1; title(ax, 'ramp', 'FontSize', 9, 'Interpreter', 'none'); end

    local_scatter_panel(nexttile((k-1)*4 + 2), xs, ys, ...
        local_colors(cmap, isDisc, tFloor));
    if k == 1; title(gca, 'markers, clim floor 0.5', 'FontSize', 9, 'Interpreter', 'none'); end

    local_scatter_panel(nexttile((k-1)*4 + 3), xs, ys, ...
        local_colors(cmap, isDisc, tData));
    if k == 1; title(gca, 'markers, clim = data', 'FontSize', 9, 'Interpreter', 'none'); end

    ax = nexttile((k-1)*4 + 4); hold(ax, 'on');
    nv = size(curves, 2);
    Cc = local_colors(cmap, isDisc, (0:nv-1).' / max(nv - 1, 1));
    for j = 1:nv
        plot(ax, zc, curves(:, j), '-', 'LineWidth', 1.6, 'Color', Cc(j, :));
    end
    set(ax, 'FontSize', 8); box(ax, 'off'); xlim(ax, [0 1]);
    if k == 1; title(ax, 'lines', 'FontSize', 9, 'Interpreter', 'none'); end
end

fprintf(['\ncr = contrast ratio against white paper; 3:1 is the usual floor for a\n' ...
         'graphical element. min dL* is the smallest lightness step between\n' ...
         'adjacent levels, and about 6 is where adjacent levels stop separating.\n' ...
         'A continuous ramp is measured over 9 evenly spaced samples so it is\n' ...
         'comparable with the quantised rows.\n\n']);

set(fig, 'Position', [60 40 1250 140 * nC]);


%% ===================== local functions =====================

function [cmap, hiPos] = local_build(base, crFloor, nLevels)
%LOCAL_BUILD Orient dark-first, truncate at the contrast floor, quantise, reverse.
    base = min(max(double(base), 0), 1);
    L = local_lstar(base);
    if L(1) > L(end)
        base = flipud(base);            % some maps ship light-first
    end

    cr = local_contrast(base);
    if crFloor > 0
        bad = find(cr < crFloor, 1, 'first');
        if isempty(bad)
            hiIdx = size(base, 1);
        else
            hiIdx = max(2, bad - 1);
        end
    else
        hiIdx = size(base, 1);
    end
    hiPos = (hiIdx - 1) / (size(base, 1) - 1);

    x = linspace(0, 1, size(base, 1));
    if nLevels > 0
        tq = linspace(0, hiPos, nLevels);
    else
        tq = linspace(0, hiPos, 256);
    end
    cmap = interp1(x, base, tq(:), 'linear');
    cmap = flipud(min(max(cmap, 0), 1));   % row 1 palest, last row darkest
end

function local_scatter_panel(ax, x, y, C)
%LOCAL_SCATTER_PANEL One objective-space panel with explicit per-point colour.
    hold(ax, 'on');
    scatter(ax, x, y, 70, C, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.7);
    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 8);
    box(ax, 'off');
end

function C = local_colors(cmap, isDiscrete, t)
%LOCAL_COLORS Colour for each normalised value in t, one row each.
    t = min(max(t(:), 0), 1);
    n = size(cmap, 1);
    if isDiscrete
        idx = min(n, max(1, ceil(t * n)));
        idx(t == 0) = 1;
        C = cmap(idx, :);
    else
        C = interp1(linspace(0, 1, n), cmap, t, 'linear');
    end
end

function [Llo, Lhi, crWorst, crBest, dLmin] = local_metrics(cmap)
%LOCAL_METRICS Lightness span, contrast against white, smallest adjacent step.
    L  = local_lstar(cmap);
    cr = local_contrast(cmap);
    Llo = min(L); Lhi = max(L);
    crWorst = min(cr); crBest = max(cr);
    m   = min(9, size(cmap, 1));
    idx = round(linspace(1, size(cmap, 1), m));
    dLmin = min(abs(diff(L(idx))));
end

function L = local_lstar(rgb)
%LOCAL_LSTAR CIE L* of each row of an sRGB matrix.
    Y = local_luminance(rgb);
    L = 903.3 * Y;
    big = Y > 0.008856;
    L(big) = 116 * Y(big).^(1/3) - 16;
end

function cr = local_contrast(rgb)
%LOCAL_CONTRAST WCAG contrast ratio of each row against white.
    cr = 1.05 ./ (local_luminance(rgb) + 0.05);
end

function Y = local_luminance(rgb)
%LOCAL_LUMINANCE Relative luminance of each row of an sRGB matrix.
    s = min(max(rgb, 0), 1);
    lin = s / 12.92;
    big = s > 0.04045;
    lin(big) = ((s(big) + 0.055) / 1.055).^2.4;
    Y = lin * [0.2126; 0.7152; 0.0722];
end

function [z, x, y] = local_scatter_data(root)
%LOCAL_SCATTER_DATA BO rows of every campaign on disk, else a stand-in.
    z = []; x = []; y = [];
    d = dir(fullfile(root, 'results'));
    for i = 1:numel(d)
        if ~d(i).isdir || startsWith(d(i).name, '.'); continue; end
        f = fullfile(root, 'results', d(i).name, 'results.csv');
        if ~isfile(f); continue; end
        try
            T = readtable(f);
        catch
            continue
        end
        if ~all(ismember({'z', 'SSE', 'SSdU'}, T.Properties.VariableNames)); continue; end
        z = [z; double(T.z)];        %#ok<AGROW>
        y = [y; double(T.SSE)];      %#ok<AGROW>
        x = [x; double(T.SSdU)];     %#ok<AGROW>
    end
    ok = isfinite(z) & isfinite(x) & isfinite(y) & x > 0 & y > 0;
    z = z(ok); x = x(ok); y = y(ok);
    if isempty(z)
        warning('No results.csv found; using synthetic scatter data.');
        rng(0);
        z = [0.833 + 0.167 * rand(40, 1); ones(9, 1)];
        x = 10.^(-1.5 + 2.2 * rand(numel(z), 1));
        y = 10.^(4.0 + 0.9 * rand(numel(z), 1));
    end
end

function [z, curves] = local_curve_data(here)
%LOCAL_CURVE_DATA phi(z) = I_z(a, b) per published vintage, else a stand-in.
    z = linspace(0, 1, 200).';
    a = []; b = [];
    f = fullfile(here, 'storage', 'surrogate_vintages.csv');
    if isfile(f)
        try
            T = readtable(f);
            if all(ismember({'SSE_a', 'SSE_b'}, T.Properties.VariableNames))
                a = double(T.SSE_a); b = double(T.SSE_b);
            end
        catch
        end
    end
    if isempty(a)
        warning('No surrogate_vintages.csv; using synthetic phi curves.');
        a = linspace(0.84, 0.92, 9).'; b = linspace(3.20, 2.90, 9).';
    end
    curves = zeros(numel(z), numel(a));
    for j = 1:numel(a)
        curves(:, j) = betainc(z, a(j), b(j));
    end
end
