function label = pretty_case(name)
%PRETTY_CASE Turn a folder name like results_case1 into "Case 1".
tok = regexp(name, 'case(\d+)\s*$', 'tokens', 'once');
if isempty(tok)
    label = string(strrep(strrep(name, 'results_', ''), '_', ' '));
    if strlength(label) > 0
        % Keep the fallback in the same register as "Case 1": a folder
        % called results_baseline reads "Baseline" in a legend.
        label = upper(extractBefore(label, 2)) + extractAfter(label, 1);
    end
else
    label = "Case " + string(tok{1});
end
end
