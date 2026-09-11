function F = ra_select_cases(F, mode)
%RA_SELECT_CASES Keep the runtime-aware arm, or only the single-fidelity one.
%
%   F = ra_select_cases(F, "mf")       drops every campaign whose folder name
%                                      contains "baseline"
%   F = ra_select_cases(F, "baseline") keeps only SF
%
%   Works on either storage artifact: the frontier struct (E, D, Tp,
%   isPareto) or the timeline struct (A, doeCount). Shared axis limits are
%   recomputed from the retained evaluations, so a figure never pads its axes
%   for points it does not draw.
%
%   Use this only where the split is structural. With the cost-aware arm and
%   SF as the two arms of one experiment, SF belongs in every figure that
%   compares them;
%   the remaining callers are the z=1 refinement figure and table, where a run
%   that never left z = 1 has nothing to refine and would plot as an empty
%   panel. Styling, not filtering, is what separates the arms elsewhere --
%   see ra_case_style.

    names  = string(F.caseNames);
    isBase = contains(lower(names), "baseline");
    switch lower(string(mode))
        case "mf",       keep = ~isBase;
        case "baseline", keep = isBase;
        otherwise
            error('ra_select_cases:mode', 'mode must be "mf" or "baseline", got "%s".', mode);
    end

    if ~any(keep)
        error('ra_select_cases:empty', ...
            'No case left after selecting "%s" from: %s', mode, strjoin(names, ', '));
    end

    F.caseNames  = F.caseNames(keep);
    F.caseLabels = F.caseLabels(keep);
    for f = ["E", "D", "Tp", "isPareto", "A"]
        if isfield(F, f)
            F.(f) = F.(f)(keep);
        end
    end
    if isfield(F, 'doeCount')
        F.doeCount = F.doeCount(keep);
    end
    F.nCases = numel(F.caseNames);

    if isfield(F, 'E') && isfield(F, 'xLimAll')
        kept = vertcat(F.E{:});
        F.xLimAll = padded_log_limits(double(kept.SSdU), 1.25);
        F.yLimAll = padded_log_limits(double(kept.SSE),  1.25);
        F.zLo     = min(0.5, min(double(kept.z)));
    end
end
