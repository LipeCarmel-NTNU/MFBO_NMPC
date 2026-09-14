function col = ra_selected_color(ctx, arm)
%RA_SELECTED_COLOR One colour per arm for the selected_* diagnostic figures.
%   The single-fidelity arm keeps the Vermillion it has everywhere else;
%   the multi-fidelity arm keeps its Blue. Columns are labelled, so the
%   colour only has to separate the two arms.
if arm == "BO"
    col = ctx.baselineColor;
else
    col = ctx.plotColors(1, :);
end
end
