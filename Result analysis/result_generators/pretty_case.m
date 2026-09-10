function label = pretty_case(name)
%PRETTY_CASE Turn a folder name like results_case1 into "Case 1".
tok = regexp(name, 'case(\d+)\s*$', 'tokens', 'once');
if isempty(tok)
    label = string(strrep(strrep(name, 'results_', ''), '_', ' '));
else
    label = "Case " + string(tok{1});
end
end
