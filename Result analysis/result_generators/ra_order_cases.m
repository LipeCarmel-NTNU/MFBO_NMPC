function ord = ra_order_cases(names)
%RA_ORDER_CASES Canonical case order: multi-fidelity cases first, baseline last.
%
%   categories() returns its levels alphabetically, so "results_baseline"
%   would sort ahead of "results_case1" and every index-based style lookup
%   (plotColors(k,:), caseMarkers(k), caseLines(k)) would shift: Case 1 would
%   be drawn in the colour Case 2 had in the submitted figures. Ordering here
%   keeps colour attached to the case and leaves the two-case figures
%   byte-identical to what the paper already shows.
%
%   Within the multi-fidelity group the order is natural: a trailing number
%   sorts numerically, everything else alphabetically after it.

    names = string(names);
    isBase = contains(lower(names), "baseline");

    num = nan(numel(names), 1);
    for i = 1:numel(names)
        tok = regexp(names(i), '(\d+)\s*$', 'tokens', 'once');
        if ~isempty(tok)
            num(i) = str2double(tok{1});
        end
    end

    key = [double(isBase(:)), isnan(num), num, double(1:numel(names))'];
    key(isnan(key)) = 0;
    [~, ord] = sortrows(key, [1 2 3 4]);
    ord = ord(:)';
end
