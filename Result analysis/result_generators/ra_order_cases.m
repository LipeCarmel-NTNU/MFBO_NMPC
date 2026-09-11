function ord = ra_order_cases(names)
%RA_ORDER_CASES Canonical arm order: the studied campaigns first, SF last.
%
%   categories() returns its levels alphabetically, so "results_baseline"
%   would sort ahead of "results_cost_aware" and the cost-aware arm would lose
%   the first slot in every ordered list -- legends, table rows, tiled panels.
%   Ordering here puts it first and SF last wherever both are drawn.
%
%   Colour no longer rides on this order (ra_case_style keys on the folder
%   name), but reading order still does, so it stays.
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
