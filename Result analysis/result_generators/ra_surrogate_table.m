function [S, caseNames, caseLabels] = ra_surrogate_table(ctx)
%RA_SURROGATE_TABLE The published surrogate vintages, one row per fit.
%   Shared by the three surrogate figures. Normalises the 'case' column name,
%   which readtable renames to 'xCase' in older preprocessed.mat builds
%   because 'case' is a MATLAB keyword.
    P = ra_require(ctx, "preprocessed");
    S = P.surrogate;
    if ismember('xCase', S.Properties.VariableNames)
        S = renamevars(S, 'xCase', 'case');
    end
    if isempty(S)
        error('ra_surrogate_table:noData', ...
            'surrogate table is empty. Run parse_registry.py, then initial_preprocessing.');
    end
    S = sortrows(S, {'case', 'vintage'});
    S.case    = removecats(categorical(S.case));
    caseNames = categories(S.case);
    caseLabels = arrayfun(@pretty_case, string(caseNames));
end
