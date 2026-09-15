function [cmap, hiPos] = ra_seq_colormap(name, crFloor)
%RA_SEQ_COLORMAP One sequential map, truncated for marks on white and reversed.
%
%   cmap = ra_seq_colormap("navia") returns a 256-by-3 sRGB map in which row 1
%   is the PALEST colour and the last row the darkest, so a value axis running
%   low to high puts the high value on the dark end.
%
%   Two things happen to the stored map before it is returned.
%
%   Orientation. Some maps ship dark-first and some light-first, so the map is
%   flipped if needed to put the dark end at position 0 before anything else.
%
%   Truncation. A perceptually uniform sequential map is uniform in LIGHTNESS,
%   which is exactly what drives one end of it to paper-white: navia ends at
%   #FCF4D9, a contrast ratio of 1.10:1 against white. That is correct for a
%   filled field, which is what these maps are designed for, and useless for
%   discrete markers and thin lines over white paper. The pale end is therefore
%   cut at the point where contrast against white falls through crFloor
%   (default 3.0, the usual floor for a graphical element). Every map cut this
%   way ends near L* = 61, so maps differ in hue path rather than in how much
%   was trimmed.
%
%   Maps come from dependencies/plot_utils/colormaps/seq_maps.mat. The README
%   beside that file carries the provenance and the citation each map needs;
%   whatever reaches a paper figure has to be cited in the caption or methods.
%
%   Result analysis/pick_seq_colormap.m renders the candidates side by side.

arguments
    name (1,1) string
    crFloor (1,1) double {mustBeNonnegative} = 3.0
end

persistent MAPS
if isempty(MAPS)
    genDir  = fileparts(mfilename('fullpath'));
    mapFile = fullfile(fileparts(fileparts(genDir)), 'dependencies', ...
        'plot_utils', 'colormaps', 'seq_maps.mat');
    if ~isfile(mapFile)
        error('ra_seq_colormap:missing', 'Colormap file not found: %s', mapFile);
    end
    MAPS = load(mapFile);
end
if ~isfield(MAPS, name)
    error('ra_seq_colormap:unknown', '%s is not in seq_maps.mat. Available: %s', ...
        name, strjoin(string(fieldnames(MAPS)).', ', '));
end

base = min(max(double(MAPS.(name)), 0), 1);
if local_lstar(base(1, :)) > local_lstar(base(end, :))
    base = flipud(base);
end

n = size(base, 1);
if crFloor > 0
    bad = find(local_contrast(base) < crFloor, 1, 'first');
    if isempty(bad)
        hiIdx = n;
    else
        hiIdx = max(2, bad - 1);
    end
else
    hiIdx = n;
end
hiPos = (hiIdx - 1) / (n - 1);

cmap = interp1(linspace(0, 1, n), base, linspace(0, hiPos, 256).', 'linear');
cmap = flipud(min(max(cmap, 0), 1));
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
