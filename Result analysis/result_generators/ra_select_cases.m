function F = ra_select_cases(F, mode)
%RA_SELECT_CASES Keep the multi-fidelity cases, or only the baseline.
%
%   F = ra_select_cases(F, "mf")       drops every campaign whose folder name
%                                      contains "baseline"
%   F = ra_select_cases(F, "baseline") keeps only those
%
%   Works on either storage artifact: the frontier struct (E, D, Tp,
%   isPareto) or the timeline struct (A, doeCount). Shared axis limits are
%   recomputed from the retained evaluations, so a figure never pads its axes
%   for points it does not draw.
%
%   The baseline is a single-fidelity reference run, not a third case: it
%   belongs in the two figures that compare cost and frontier position, and
%   nowhere else, or it would imply a like-for-like third arm.

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
