function S = ra_case_style(ctx, caseNames)
%RA_CASE_STYLE Colour, marker and line style per campaign, keyed by name.
%
%   S = ra_case_style(ctx, F.caseNames) returns
%       S.color   nCases-by-3
%       S.marker  1-by-nCases string
%       S.line    1-by-nCases string
%       S.isBase  1-by-nCases logical, true for the single-fidelity arm
%
%   The single-fidelity arm (SF) is always Wong Vermillion, cross marker,
%   dotted line, wherever it is drawn. Every other campaign consumes the case
%   sequence -- Blue, BluishGreen, ... -- in order, so the cost-aware arm keeps
%   the blue it has in the submitted figures whether or not SF shares the axes.
%
%   Indexing ctx.plotColors by position instead would give SF whatever colour
%   its slot fell on and shift the others the moment the list changes length,
%   which is exactly the failure ra_order_cases was written to avoid. Every
%   generator that draws a per-case series asks here.

    names  = string(caseNames);
    n      = numel(names);
    isBase = contains(lower(names), "baseline");

    S.color  = zeros(n, 3);
    S.marker = strings(1, n);
    S.line   = strings(1, n);
    S.isBase = isBase(:)';

    nSeq = size(ctx.plotColors, 1);
    j    = 0;
    for k = 1:n
        if isBase(k)
            S.color(k, :) = ctx.baselineColor;
            S.marker(k)   = "x";
            S.line(k)     = ":";
        else
            j = min(j + 1, nSeq);
            S.color(k, :) = ctx.plotColors(j, :);
            S.marker(k)   = ctx.caseMarkers(j);
            S.line(k)     = ctx.caseLines(j);
        end
    end
end
