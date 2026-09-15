function col = ra_selected_color(ctx, arm)
%RA_SELECTED_COLOR One colour per arm for the selected_* diagnostic figures.
%   The single-fidelity arm keeps the Vermillion it has everywhere else; the
%   multi-fidelity arm keeps its Blue; the damped-Rdu blend takes the accent
%   ReddishPurple, which no case or continuum uses. Columns are labelled, so
%   the colour only has to separate the three sources.
switch string(arm)
    case "BO"
        col = ctx.baselineColor;
    case "RDU"
        col = ctx.accentColor;
    otherwise
        col = ctx.plotColors(1, :);
end
end
