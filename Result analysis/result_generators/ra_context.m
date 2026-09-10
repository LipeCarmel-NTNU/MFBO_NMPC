function ctx = ra_context()
%RA_CONTEXT Shared paths, palette and figure style for every result generator.
%
%   ctx = ra_context() resolves the project layout, puts the dependencies and
%   the generator folder on the path, applies the LaTeX/typography defaults
%   from COLOR_SCHEME.md, and returns the handful of values every generator
%   needs. Generators take ctx as their only configuration argument, so the
%   look of the figure set is defined in one place.
%
%   Fields:
%     repo_root, resultsRoot, graphicsDir, numericalDir, storageDir
%     fontSize, plotColors, accentColor, caseMarkers, caseLines, seqMap
%
%   Palette note: plotColors carries Blue, BluishGreen, Vermillion and the
%   accent is ReddishPurple. The previous scripts took the accent from
%   plotColors(3,:), which is the same colour a third case would have been
%   drawn in, so with three cases a series was indistinguishable from the
%   Pareto continuum and from the y = x reference in runtime_consistency.
%   Cases one and two are unchanged.

    genDir  = fileparts(mfilename('fullpath'));
    here    = fileparts(genDir);              % Result analysis
    ctx.here      = here;
    ctx.repo_root = fileparts(here);

    addpath(genpath(fullfile(ctx.repo_root, 'dependencies')));
    addpath(genDir);

    ctx.resultsRoot  = fullfile(ctx.repo_root, 'results');
    ctx.graphicsDir  = fullfile(ctx.resultsRoot, 'graphical_results');
    ctx.numericalDir = fullfile(ctx.resultsRoot, 'numerical results');
    ctx.storageDir   = fullfile(here, 'storage');
    for d = [string(ctx.graphicsDir), string(ctx.numericalDir), string(ctx.storageDir)]
        if ~isfolder(d); mkdir(d); end
    end

    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    ctx.fontSize = 20;

    NC = nature_methods_colors();
    ctx.plotColors  = [NC.Blue; NC.BluishGreen; NC.Vermillion];
    ctx.accentColor = NC.ReddishPurple;
    % The single-fidelity baseline is drawn in Wong Vermillion wherever it
    % appears, deliberately outside the case sequence: it is a reference run,
    % not a third case. Wong Orange (#E69F00) is the lighter amber and sits at
    % 2.19:1 on a white ground, too faint for unlabelled markers.
    ctx.baselineColor = NC.Vermillion;
    ctx.caseMarkers = ["o", "^", "d"];
    ctx.caseLines   = ["-", "-.", ":"];
    ctx.seqMap      = load_navia_colormap(256);
end
