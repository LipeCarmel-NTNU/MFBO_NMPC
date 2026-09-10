%% surrogate_analysis
% Evolution of the phi(z) fidelity-cost surrogate across vintages.
%
% Inputs
%   storage/preprocessed.mat -> surrogate : one row per phi(z) vintage with
%   fitted shape parameters (a, b), regularization (lambda), fit loss, and
%   refit wall time, per target (SSE = J_track, SSdU = J_TV) and per case.
%
% Figures (saved to results/graphical_results)
%   1. surrogate_coefficients.png/.pdf : a, b, lambda vs vintage, per target
%   2. surrogate_fit_quality.png/.pdf  : fit loss and refit wall time vs vintage
%   3. surrogate_phi_curves.png/.pdf   : phi(z) = I_z(a, b) per vintage,
%                                        colored by vintage with the navia map
%
% Color scheme (see COLOR_SCHEME.md)
%   Case 1 -> Wong Blue        [0 114 178], marker "o"
%   Case 2 -> Wong BluishGreen [0 158 115], marker "^"
%   Accent -> Wong ReddishPurple [204 121 167]
%   Sequential (vintage progression) -> navia (Scientific Colour Maps)

%% Paths and data
here      = fileparts(mfilename('fullpath'));
repo_root = fileparts(here);
addpath(genpath(fullfile(repo_root, 'dependencies')));

load(fullfile(here, 'storage', 'preprocessed.mat'), 'surrogate');
% Older preprocessed.mat builds carry 'xCase' (readtable's default naming rule
% renames the reserved word 'case'). Normalize so the script can rely on 'case'.
if ismember('xCase', surrogate.Properties.VariableNames)
    surrogate = renamevars(surrogate, 'xCase', 'case');
end
if isempty(surrogate)
    error('surrogate_analysis:noData', ...
        'surrogate table is empty. Run parse_registry.py, then initial_preprocessing.m.');
end

graphics_dir = fullfile(repo_root, 'results', 'graphical_results');
if ~exist(graphics_dir, 'dir'); mkdir(graphics_dir); end

%% Style (COLOR_SCHEME.md)
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
fontSize = 20;

plotColors  = nature_methods_colors(3);   % Blue (Case 1), BluishGreen (Case 2), ReddishPurple (accent)
caseMarkers = ["o", "^", "d"];
seqMap      = load_navia_colormap(256);   % sequential map for vintage progression

%% Split per case, ordered by vintage
S = sortrows(surrogate, {'case', 'vintage'});
S.case    = removecats(categorical(S.case));
caseNames = categories(S.case);
nCases    = numel(caseNames);
if nCases > numel(caseMarkers)
    error('surrogate_analysis:tooManyCases', ...
        'Only %d case styles defined; found %d cases.', numel(caseMarkers), nCases);
end
caseLabels = arrayfun(@pretty_case, string(caseNames));

targets      = ["SSE", "SSdU"];
targetLabels = ["J_{\mathrm{track}}", "J_{\mathrm{TV}}"];

%% Figure 1: coefficient evolution (a, b, lambda) per target
fig1 = figure('Color', 'w', 'Name', 'Surrogate Coefficients by Vintage');
tiledlayout(fig1, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
coeffs      = ["a", "b", "lambda"];
coeffLabels = ["a", "b", "\lambda"];
panelLabels = ["a", "b", "c"; "d", "e", "f"];

for it = 1:numel(targets)
    for ic = 1:numel(coeffs)
        ax = nexttile; hold(ax, 'on');
        col = targets(it) + "_" + coeffs(ic);
        for k = 1:nCases
            T = S(S.case == caseNames{k}, :);
            plot(ax, double(T.vintage), double(T.(col)), '-', ...
                'Marker', caseMarkers(k), 'MarkerSize', 7, ...
                'MarkerFaceColor', 'w', 'LineWidth', 2.0, ...
                'Color', plotColors(k, :), 'DisplayName', caseLabels(k));
        end
        xlabel(ax, '$v$ (vintage)');
        ylabel(ax, "$" + coeffLabels(ic) + "$ ($" + targetLabels(it) + "$)");
        title(ax, "$\mathbf{" + panelLabels(it, ic) + "}$", 'Interpreter', 'latex');
        ax.TitleHorizontalAlignment = 'left';
        set(ax, 'FontSize', fontSize);
        grid(ax, 'off'); box(ax, 'off');
        if it == 1 && ic == 1
            legend(ax, 'Location', 'best');
        end
    end
end
save_plot_outputs(fig1, fullfile(graphics_dir, 'surrogate_coefficients.png'), fontSize, 1400, 700);

%% Figure 2: fit quality (loss) and refit cost (wall time)
fig2 = figure('Color', 'w', 'Name', 'Surrogate Fit Quality by Vintage');
tiledlayout(fig2, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
panelLabels2 = ["a", "b", "c"];
lossCols   = ["SSE_loss", "SSdU_loss", "fit_wall_s"];
lossLabels = ["$\mathcal{L}$ ($J_{\mathrm{track}}$)", ...
              "$\mathcal{L}$ ($J_{\mathrm{TV}}$)", ...
              "$t_{\mathrm{fit}}$ (s)"];

for ip = 1:numel(lossCols)
    ax = nexttile; hold(ax, 'on');
    for k = 1:nCases
        T = S(S.case == caseNames{k}, :);
        plot(ax, double(T.vintage), double(T.(lossCols(ip))), '-', ...
            'Marker', caseMarkers(k), 'MarkerSize', 7, ...
            'MarkerFaceColor', 'w', 'LineWidth', 2.0, ...
            'Color', plotColors(k, :), 'DisplayName', caseLabels(k));
    end
    xlabel(ax, '$v$ (vintage)');
    ylabel(ax, lossLabels(ip));
    title(ax, "$\mathbf{" + panelLabels2(ip) + "}$", 'Interpreter', 'latex');
    ax.TitleHorizontalAlignment = 'left';
    set(ax, 'FontSize', fontSize);
    grid(ax, 'off'); box(ax, 'off');
    if ip == 1
        legend(ax, 'Location', 'best');
    end
end
save_plot_outputs(fig2, fullfile(graphics_dir, 'surrogate_fit_quality.png'), fontSize, 1400, 420);

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

fprintf('Saved surrogate figures for %d case(s), %d vintage row(s) -> %s\n', ...
    nCases, height(S), graphics_dir);

%% Local functions
function label = pretty_case(name)
%PRETTY_CASE Turn a folder name like results_case1 into "Case 1".
tok = regexp(name, '(\d+)\s*$', 'tokens', 'once');
if isempty(tok)
    label = string(name);
else
    label = "Case " + string(tok{1});
end
end


function cmap = load_navia_colormap(n)
%LOAD_NAVIA_COLORMAP Load the sequential map used for ordered-color plots.
matPath = which('navia.mat');
if strlength(string(matPath)) == 0
    error('Unable to locate navia.mat on MATLAB path.');
end
S = load(matPath, 'navia');
if ~isfield(S, 'navia')
    error('File %s does not contain variable ''navia''.', matPath);
end
cmap = double(S.navia);
if size(cmap, 2) ~= 3
    error('Variable ''navia'' in %s must have 3 columns (RGB).', matPath);
end
if size(cmap, 1) ~= n
    x  = linspace(0, 1, size(cmap, 1));
    xq = linspace(0, 1, n);
    cmap = interp1(x, cmap, xq, 'linear');
end
cmap = min(max(cmap, 0), 1);
end


function save_plot_outputs(figHandle, pngPath, fontSize, figWidthPx, figHeightPx)
%SAVE_PLOT_OUTPUTS Export standardized figure files (PNG 300 dpi + vector PDF).
arguments
    figHandle
    pngPath
    fontSize (1,1) double = 14
    figWidthPx (1,1) double = 900
    figHeightPx (1,1) double = 500
end
figure(figHandle);
set_fig_size(figWidthPx, figHeightPx);
set_font_size(fontSize);
exportgraphics(figHandle, pngPath, 'Resolution', 300);
[folderPath, fileStem] = fileparts(pngPath);
pdfPath = fullfile(folderPath, strcat(fileStem, '.pdf'));
save_figure(pdfPath, NaN, false);
end
